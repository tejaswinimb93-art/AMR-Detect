
import os
import re
import json
import ast
import shutil
import subprocess
import tempfile
from pathlib import Path

import pandas as pd
import streamlit as st

st.set_page_config(
    page_title="AMRIVA",
    page_icon="🧬",
    layout="wide",
)

st.markdown("""
<style>
.block-container {max-width: 1500px; padding-top: 2rem;}
.hero {
    padding: 28px 32px; border-radius: 20px; margin-bottom: 22px;
    background: linear-gradient(135deg,#123c69,#0f766e); color: white;
}
.hero h1 {font-size: 48px; margin: 0 0 5px 0;}
.hero p {font-size: 17px; margin: 4px 0;}
.card {padding: 16px; border: 1px solid #ddd; border-radius: 14px;}
</style>
""", unsafe_allow_html=True)

st.markdown("""
<div class="hero">
<h1>🧬 AMRIVA</h1>
<p><b>AMRFinderPlus Results Organizer & Explorer</b></p>
<p>Organize genomic AMR findings across samples, organisms and antimicrobial-associated determinants.</p>
</div>
""", unsafe_allow_html=True)

# -----------------------------
# Helpers
# -----------------------------
DISPLAY_COLUMNS = [
    "Sample ID", "Organism", "Antibiotic / antimicrobial association",
    "AMR determinant", "Class", "Subclass", "Mechanism / Product",
    "Method", "Identity", "Coverage", "Sequence / Contig",
    "Nucleotide Sequence", "Start", "Stop", "Reference",
    "Interpretation", "Raw Result"
]

def clean(v):
    if v is None:
        return ""
    s = str(v).strip()
    if s.lower() in {"nan", "none", "null"}:
        return ""
    return s

def norm_col(c):
    return re.sub(r"[^a-z0-9]+", "", str(c).strip().lower())

def find_col(df, *names):
    lookup = {norm_col(c): c for c in df.columns}
    for name in names:
        n = norm_col(name)
        if n in lookup:
            return lookup[n]
    for c in df.columns:
        nc = norm_col(c)
        for name in names:
            n = norm_col(name)
            if n and (n in nc or nc in n):
                return c
    return None

def read_table(uploaded):
    """Read real AMRFinderPlus TSV/CSV robustly, even when the extension is wrong."""
    name = str(getattr(uploaded, "name", "")).lower()
    try:
        uploaded.seek(0)
    except Exception:
        pass
    raw = uploaded.read()
    if isinstance(raw, str):
        raw = raw.encode("utf-8")
    text = raw.decode("utf-8-sig", errors="replace")
    # Detect delimiter from the actual content instead of trusting the file extension.
    first_nonempty = next((line for line in text.splitlines() if line.strip()), "")
    if "\t" in first_nonempty:
        sep = "\t"
    elif "," in first_nonempty:
        sep = ","
    elif ";" in first_nonempty:
        sep = ";"
    else:
        sep = "\t"
    from io import StringIO
    return pd.read_csv(StringIO(text), sep=sep, dtype=str, keep_default_na=False)

def fasta_records(text):
    records = []
    header = None
    seq = []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(seq)))
            header = line[1:].strip()
            seq = []
        else:
            seq.append(line)
    if header is not None:
        records.append((header, "".join(seq)))
    return records

def sample_id_from_header(header, fallback):
    # Prefer explicit sample identifiers in common header forms.
    patterns = [
        r"(?:sample[_\s-]*id|sample)[=:]([A-Za-z0-9_.-]+)",
        r"^([A-Za-z0-9_.-]+)"
    ]
    for p in patterns:
        m = re.search(p, header, flags=re.I)
        if m:
            return m.group(1)
    return fallback

def organism_from_header(header):
    # Use organism metadata only when it is explicitly present.
    patterns = [
        r"(?:organism|species|taxon)[=:]([^|;]+)",
        r"\[([A-Z][A-Za-z]+(?:\s+[a-z][A-Za-z0-9_.-]+){0,2})\]"
    ]
    for p in patterns:
        m = re.search(p, header, flags=re.I)
        if m:
            return m.group(1).strip()
    return ""

def parse_fasta_files(files):
    rows = []
    for f in files or []:
        text = f.getvalue().decode("utf-8", errors="replace")
        records = fasta_records(text)
        for i, (header, seq) in enumerate(records, start=1):
            sid = sample_id_from_header(header, Path(f.name).stem)
            # Prefer an explicit contig token from common NCBI/FASTA headers.
            contig = ""
            for pat in [
                r"(?:^|[|;\s])contig(?:[_\s-]*id)?[=:|]([A-Za-z0-9_.-]+)",
                r"\b(contig_[A-Za-z0-9_.-]+)\b",
            ]:
                m = re.search(pat, header, flags=re.I)
                if m:
                    contig = m.group(1)
                    break
            if not contig:
                token = header.split()[0] if header else f"record_{i}"
                # Handle headers such as SAMPLE|contig_001|organism=E. coli
                parts = token.split("|")
                contig = next((x for x in parts if x.lower().startswith("contig")), token)

            rows.append({
                "Sample ID": sid,
                "Organism": organism_from_header(header),
                "Sequence / Contig": contig,
                "_full_sequence": re.sub(r"\s+", "", seq).upper(),
                "_header": header,
            })
    return pd.DataFrame(rows)

def extract_nucleotide(row, fasta_map):
    contig = clean(row.get("Sequence / Contig", ""))
    start = clean(row.get("Start", ""))
    stop = clean(row.get("Stop", ""))
    sid = clean(row.get("Sample ID", ""))

    if not contig or not start or not stop:
        return ""

    try:
        a, b = int(float(start)), int(float(stop))
    except Exception:
        return ""

    # AMRFinderPlus coordinates are 1-based inclusive.
    lo, hi = min(a, b), max(a, b)
    key_options = [
        (sid, contig),
        ("", contig),
        (sid, ""),
        ("", ""),
    ]
    for key in key_options:
        seq = fasta_map.get(key)
        if seq:
            if lo >= 1 and hi <= len(seq):
                fragment = seq[lo - 1:hi]
                if a > b:
                    # Reverse-complement when coordinates are on the reverse strand.
                    fragment = reverse_complement(fragment)
                return fragment
    return ""

def reverse_complement(seq):
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(table)[::-1]

def derive_antimicrobial(row):
    """Create a searchable antimicrobial association from supplied fields.

    AMRFinderPlus generally reports resistance class/subclass rather than a
    literal drug name. The association below is therefore a class-level
    search aid, not a phenotypic susceptibility call.
    """
    for key in ["Antibiotic", "Drug", "Antimicrobial", "Antimicrobial(s)",
                "Antibiotic / antimicrobial association"]:
        if key in row and clean(row[key]):
            return clean(row[key])

    sub = clean(row.get("Subclass", "")).upper()
    cls = clean(row.get("Class", "")).upper()
    mapping = {
        "BETA-LACTAM": "Beta-lactams (penicillins, cephalosporins, carbapenems)",
        "PENICILLIN": "Penicillins",
        "CEPHALOSPORIN": "Cephalosporins",
        "CARBAPENEM": "Carbapenems",
        "GLYCOPEPTIDE": "Glycopeptides",
        "AMINOGLYCOSIDE": "Aminoglycosides",
        "QUINOLONE": "Quinolones / fluoroquinolones",
        "SULFONAMIDE": "Sulfonamides",
        "TETRACYCLINE": "Tetracyclines",
        "MACROLIDE": "Macrolides",
        "PHENICOL": "Phenicols",
        "LINCOSAMIDE": "Lincosamides",
        "RIFAMYCIN": "Rifamycins",
        "TRIMETHOPRIM": "Trimethoprim",
        "FOSFOMYCIN": "Fosfomycin",
        "POLYMYXIN": "Polymyxins",
    }
    if sub in mapping:
        return mapping[sub]
    if cls in mapping:
        return mapping[cls]
    if sub:
        return f"Associated class: {sub}"
    if cls:
        return f"Associated class: {cls}"
    return "Not specified in AMRFinderPlus output"

def interpretation(row):
    det = clean(row.get("AMR determinant", ""))
    cls = clean(row.get("Class", ""))
    sub = clean(row.get("Subclass", ""))
    ident = clean(row.get("Identity", ""))
    cov = clean(row.get("Coverage", ""))
    if not det or det == "Not reported":
        return "The uploaded record did not contain an AMR determinant field. No determinant was invented."
    parts = [f"Genomic AMR determinant {det} was detected"]
    if cls:
        parts.append(f"in the {cls} class")
    if sub:
        parts.append(f"({sub})")
    if ident or cov:
        qc = []
        if ident:
            qc.append(f"identity {ident}%")
        if cov:
            qc.append(f"coverage {cov}%")
        parts.append("with " + " and ".join(qc))
    return " ".join(parts) + ". Genomic detection does not by itself establish phenotypic susceptibility or resistance."


def parse_raw_result(value):
    """Parse a JSON/dict-like raw_result cell when a pre-organized file contains one."""
    s = clean(value)
    if not s:
        return {}
    # Some exports wrap the JSON in whitespace or other text.
    candidates = [s]
    if "{" in s and "}" in s:
        candidates.append(s[s.find("{"):s.rfind("}") + 1])

    for candidate in candidates:
        try:
            obj = json.loads(candidate)
            if isinstance(obj, dict):
                return obj
        except Exception:
            pass
        try:
            obj = ast.literal_eval(candidate)
            if isinstance(obj, dict):
                return obj
        except Exception:
            pass
    return {}


def flatten_raw_results(df):
    """
    Support both:
      1) genuine AMRFinderPlus TSV/CSV columns, and
      2) AMRIVA/demo files where technical AMRFinderPlus fields are stored
         inside a raw_result JSON column.

    Explicit columns always win. Parsed raw_result values only fill missing
    fields.
    """
    df = df.copy()

    raw_col = find_col(df, "raw_result", "raw result", "raw AMRFinderPlus result")
    if not raw_col:
        return df

    aliases = {
        "Protein identifier": ["Protein identifier", "Protein ID", "protein_id"],
        "Contig id": ["Contig id", "Contig ID", "contig", "Sequence name", "Sequence ID"],
        "Start": ["Start", "start", "Target start"],
        "Stop": ["Stop", "stop", "Target stop"],
        "Strand": ["Strand", "strand", "Target strand"],
        "Element symbol": ["Element symbol", "Gene symbol", "AMR determinant", "determinant", "Gene", "Element"],
        "Element name": ["Element name", "Gene name", "Product", "Mechanism / Product", "Mechanism"],
        "Scope": ["Scope"],
        "Element type": ["Element type", "Type"],
        "Element subtype": ["Element subtype", "Subtype"],
        "Class": ["Class"],
        "Subclass": ["Subclass", "Sub-class"],
        "Method": ["Method"],
        "Target length": ["Target length", "Target length (bp)", "Element length"],
        "Reference sequence length": ["Reference sequence length", "Reference length"],
        "% Coverage of reference sequence": [
            "% Coverage of reference sequence", "% coverage", "Coverage",
            "Percent coverage", "Identity coverage"
        ],
        "% Identity to reference sequence": [
            "% Identity to reference sequence", "% identity", "Identity",
            "Percent identity"
        ],
        "Alignment length": ["Alignment length", "Alignment"],
        "Accession of closest sequence": [
            "Accession of closest sequence", "Closest reference accession",
            "Reference accession", "Accession", "Reference"
        ],
        "Name of closest sequence": [
            "Name of closest sequence", "Closest reference name",
            "Reference name"
        ],
        "HMM id": ["HMM id", "HMM ID"],
        "HMM description": ["HMM description", "HMM Description"],
    }

    parsed_rows = []
    for value in df[raw_col]:
        obj = parse_raw_result(value)
        parsed_rows.append(obj)

    # Build columns only where the raw JSON actually contains useful values.
    for target, possible_keys in aliases.items():
        if target in df.columns:
            continue

        values = []
        found_any = False
        for obj in parsed_rows:
            value = ""
            if obj:
                lookup = {norm_col(k): k for k in obj.keys()}
                for alias in possible_keys:
                    key = lookup.get(norm_col(alias))
                    if key is not None and clean(obj.get(key)):
                        value = obj.get(key)
                        found_any = True
                        break
            values.append(value)

        if found_any:
            df[target] = values

    return df


def standardize_amrfinder(df, fasta_df=None):
    df = df.copy()
    df.columns = [str(c).strip() for c in df.columns]

    # Important: some earlier AMRIVA demo files stored the actual technical
    # AMRFinderPlus fields in a raw_result JSON cell. Flatten that first.
    df = flatten_raw_results(df)

    # These aliases cover current AMRFinderPlus field names as well as
    # common AMRIVA/demo exports.
    c_sample = find_col(df, "Sample ID", "Sample", "SampleID", "Name", "name", "Input name", "assembly")
    c_org = find_col(df, "Organism", "Species", "Taxon", "Organism name", "Organism_Name")
    c_det = find_col(
        df, "Element symbol", "Gene symbol", "AMR determinant",
        "determinant", "AMR gene", "Gene", "Element", "Resistance gene",
        "Resistance determinant"
    )
    c_class = find_col(df, "Class")
    c_sub = find_col(df, "Subclass", "Sub-class")
    c_mech = find_col(
        df, "Mechanism", "Mechanism / Product", "Product", "Element name",
        "Sequence name", "Gene name", "Protein name", "Mechanism/Product"
    )
    c_method = find_col(df, "Method")
    c_ident = find_col(
        df, "% Identity to reference sequence", "% Identity to reference",
        "% identity", "Identity", "Percent identity"
    )
    c_cov = find_col(
        df, "% Coverage of reference sequence", "% Coverage of reference",
        "% coverage", "Coverage", "Percent coverage"
    )
    c_contig = find_col(
        df, "Contig id", "Contig ID", "Contig", "Sequence name",
        "Sequence ID", "Sequence", "Sequence/Contig"
    )
    c_start = find_col(df, "Start", "Target start")
    c_stop = find_col(df, "Stop", "Target stop")
    c_ref_acc = find_col(
        df, "Accession of closest sequence",
        "Closest reference accession", "Reference accession", "Accession"
    )
    c_ref_name = find_col(
        df, "Name of closest sequence",
        "Closest reference name", "Reference name"
    )
    c_ref = find_col(df, "Reference")
    c_ant = find_col(
        df, "Antibiotic", "Drug", "Antimicrobial", "Antimicrobial(s)",
        "Antibiotic / antimicrobial association"
    )

    output = []
    fasta_map = {}
    if fasta_df is not None and not fasta_df.empty:
        for _, r in fasta_df.iterrows():
            sid = clean(r.get("Sample ID"))
            contig = clean(r.get("Sequence / Contig"))
            seq = clean(r.get("_full_sequence"))
            if contig:
                fasta_map[(sid, contig)] = seq
                fasta_map.setdefault(("", contig), seq)

    for idx, r in df.iterrows():
        raw = " | ".join(
            f"{c}={clean(r[c])}" for c in df.columns if clean(r[c]) != ""
        )

        sid = clean(r[c_sample]) if c_sample else ""
        org = clean(r[c_org]) if c_org else ""

        if not sid:
            sid = f"Sample_{idx + 1}"

        # Combine closest-reference accession + closest-reference name
        # when both are present, rather than losing one of them.
        ref_acc = clean(r[c_ref_acc]) if c_ref_acc else ""
        ref_name = clean(r[c_ref_name]) if c_ref_name else ""
        ref_generic = clean(r[c_ref]) if c_ref else ""
        if ref_acc and ref_name:
            reference = f"{ref_acc} — {ref_name}"
        else:
            reference = ref_acc or ref_name or ref_generic

        rec = {
            "Sample ID": sid,
            "Organism": org,
            "Antibiotic / antimicrobial association":
                clean(r[c_ant]) if c_ant else "",
            "AMR determinant":
                clean(r[c_det]) if c_det else "",
            "Class":
                clean(r[c_class]) if c_class else "",
            "Subclass":
                clean(r[c_sub]) if c_sub else "",
            "Mechanism / Product":
                clean(r[c_mech]) if c_mech else "",
            "Method":
                clean(r[c_method]) if c_method else "",
            "Identity":
                clean(r[c_ident]) if c_ident else "",
            "Coverage":
                clean(r[c_cov]) if c_cov else "",
            "Sequence / Contig":
                clean(r[c_contig]) if c_contig else "",
            "Nucleotide Sequence": "",
            "Start":
                clean(r[c_start]) if c_start else "",
            "Stop":
                clean(r[c_stop]) if c_stop else "",
            "Reference": reference,
            "Interpretation": "",
            "Raw Result": raw,
        }

        # Do not invent an organism. If FASTA metadata contains it, use that.
        if not org and fasta_df is not None and not fasta_df.empty:
            matches = fasta_df[
                fasta_df["Sample ID"].astype(str) == sid
            ]
            if not matches.empty:
                orgs = [
                    clean(x) for x in matches["Organism"].tolist()
                    if clean(x)
                ]
                if orgs:
                    rec["Organism"] = orgs[0]

        if not rec["Antibiotic / antimicrobial association"]:
            rec["Antibiotic / antimicrobial association"] = derive_antimicrobial(rec)

        rec["Nucleotide Sequence"] = extract_nucleotide(rec, fasta_map)

        # Never silently leave display cells looking broken. These labels
        # distinguish "not supplied/reported" from a real empty result.
        fallback = {
            "Organism": "Not provided in input metadata",
            "AMR determinant": "Not reported",
            "Class": "Not reported",
            "Subclass": "Not reported",
            "Mechanism / Product": "Not reported",
            "Method": "Not reported",
            "Identity": "Not reported",
            "Coverage": "Not reported",
            "Sequence / Contig": "Not provided",
            "Nucleotide Sequence":
                "Upload corresponding FASTA to extract sequence"
                if not rec["Nucleotide Sequence"] else rec["Nucleotide Sequence"],
            "Start": "Not reported",
            "Stop": "Not reported",
            "Reference": "Not reported",
        }

        for key, replacement in fallback.items():
            if not clean(rec[key]):
                rec[key] = replacement

        rec["Interpretation"] = interpretation(rec)
        output.append(rec)

    return pd.DataFrame(output, columns=DISPLAY_COLUMNS)


def run_amrfinder_on_fasta(files):
    exe = shutil.which("amrfinder")
    if not exe:
        return None, "AMRFinderPlus executable was not found in the deployed environment."

    combined = []
    temp_root = Path(tempfile.mkdtemp(prefix="amriva_"))
    try:
        for f in files:
            path = temp_root / Path(f.name).name
            path.write_bytes(f.getvalue())
            cmd = [exe, "-n", str(path), "--plus"]
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=180)
            if result.returncode != 0:
                return None, result.stderr.strip() or "AMRFinderPlus returned an error."
            lines = result.stdout.splitlines()
            if lines:
                if not combined:
                    combined.extend(lines)
                else:
                    # skip repeated header
                    header = combined[0]
                    combined.extend(x for x in lines[1:] if x.strip() and x != header)
        if not combined:
            return None, "AMRFinderPlus returned no records."
        return "\n".join(combined), None
    except subprocess.TimeoutExpired:
        return None, "AMRFinderPlus timed out while analysing the uploaded FASTA."
    finally:
        shutil.rmtree(temp_root, ignore_errors=True)

# -----------------------------
# State
# -----------------------------
if "results" not in st.session_state:
    st.session_state.results = pd.DataFrame(columns=DISPLAY_COLUMNS)
if "fasta_df" not in st.session_state:
    st.session_state.fasta_df = pd.DataFrame()

# -----------------------------
# Sidebar
# -----------------------------
st.sidebar.title("AMRIVA")
st.sidebar.caption("Focused AMRFinderPlus results organizer")

mode = st.sidebar.radio(
    "Input mode",
    ["AMRFinderPlus result file", "FASTA + AMRFinderPlus", "FASTA only"],
)

if mode == "AMRFinderPlus result file":
    uploaded_result = st.file_uploader(
        "Upload AMRFinderPlus TSV/CSV",
        type=["tsv", "txt", "csv"],
        accept_multiple_files=False,
        key="result_upload",
    )
    fasta_files = st.file_uploader(
        "Optional FASTA files (for nucleotide sequence extraction and organism metadata)",
        type=["fa", "fasta", "fna"],
        accept_multiple_files=True,
        key="optional_fasta",
    )
    if st.button("Organize AMRFinderPlus results", type="primary"):
        if not uploaded_result:
            st.error("Upload an AMRFinderPlus TSV/CSV first.")
        else:
            try:
                raw_df = read_table(uploaded_result)
                fdf = parse_fasta_files(fasta_files) if fasta_files else pd.DataFrame()
                st.session_state.fasta_df = fdf
                st.session_state.results = standardize_amrfinder(raw_df, fdf)
                st.success(f"Organized {len(st.session_state.results)} AMRFinderPlus record(s).")
            except Exception as e:
                st.exception(e)

elif mode == "FASTA + AMRFinderPlus":
    fasta_files = st.file_uploader(
        "Upload one or more FASTA files",
        type=["fa", "fasta", "fna"],
        accept_multiple_files=True,
    )
    result_file = st.file_uploader(
        "Upload the corresponding AMRFinderPlus TSV/CSV",
        type=["tsv", "txt", "csv"],
    )
    if st.button("Process genomic results", type="primary"):
        if not fasta_files or not result_file:
            st.error("Upload both FASTA file(s) and the corresponding AMRFinderPlus result.")
        else:
            try:
                fdf = parse_fasta_files(fasta_files)
                rdf = read_table(result_file)
                st.session_state.fasta_df = fdf
                st.session_state.results = standardize_amrfinder(rdf, fdf)
                st.success(f"Processed {len(st.session_state.results)} genomic finding(s).")
            except Exception as e:
                st.exception(e)

else:
    fasta_files = st.file_uploader(
        "Upload one or more FASTA files",
        type=["fa", "fasta", "fna"],
        accept_multiple_files=True,
    )
    st.info("FASTA-only mode is for sequence metadata/organization. AMR findings require an AMRFinderPlus result or an installed AMRFinderPlus executable.")
    if st.button("Run AMRFinderPlus", type="primary"):
        if not fasta_files:
            st.error("Upload at least one FASTA file.")
        else:
            fdf = parse_fasta_files(fasta_files)
            st.session_state.fasta_df = fdf
            tsv, err = run_amrfinder_on_fasta(fasta_files)
            if err:
                st.error(err)
                st.warning("No AMR findings were invented. Upload the real AMRFinderPlus TSV/CSV if the executable is unavailable.")
            else:
                try:
                    rdf = pd.read_csv(
                        __import__("io").StringIO(tsv),
                        sep="\t",
                        dtype=str,
                        keep_default_na=False,
                    )
                    st.session_state.results = standardize_amrfinder(rdf, fdf)
                    st.success(f"AMRFinderPlus returned {len(st.session_state.results)} finding(s).")
                except Exception as e:
                    st.exception(e)

# -----------------------------
# Main explorer
# -----------------------------
results = st.session_state.results

tabs = st.tabs([
    "🏠 Home",
    "🧬 Genomic Analysis",
    "🔎 Antibiotic Explorer",
    "🔬 Cross-Sample Comparison",
    "🧫 Sample / Organism Profile",
])

with tabs[0]:
    st.markdown("## What problem does AMRIVA solve?")
    st.write(
        "AMRFinderPlus can produce detailed genomic AMR findings. When many isolates are analysed, "
        "researchers may need to inspect many records to understand which samples and organisms contain "
        "which resistance-associated determinants. AMRIVA organizes those records into a searchable explorer."
    )
    st.markdown("### Core relationship")
    st.code("Sample ↔ Organism ↔ Antimicrobial association ↔ AMR determinant ↔ Genomic location")
    st.markdown("### What AMRIVA adds")
    st.write(
        "AMRIVA is a results-organization and exploration layer. It does not replace AMRFinderPlus "
        "and does not claim that a genomic determinant alone proves phenotypic resistance."
    )

with tabs[1]:
    st.header("Genomic Analysis")
    if results.empty:
        st.info("Upload an AMRFinderPlus result file to populate this table.")
    else:
        st.metric("AMRFinderPlus findings", len(results))
        st.caption(
            "Fields are mapped from AMRFinderPlus output. "
            "Organism requires metadata, and nucleotide sequence extraction "
            "requires the corresponding FASTA file. 'Not reported' means the "
            "source file did not contain that field; AMRIVA does not invent values."
        )
        st.dataframe(results, use_container_width=True, hide_index=True)
        csv = results.to_csv(index=False).encode("utf-8")
        st.download_button(
            "⬇️ Download organized genomic results",
            csv,
            "AMRIVA_organized_genomic_results.csv",
            "text/csv",
        )

with tabs[2]:
    st.header("Antibiotic Explorer")
    if results.empty:
        st.info("Upload and organize results first.")
    else:
        st.caption("Search the antimicrobial association, resistance class/subclass, determinant, sample or organism. A match here means a genomic AMR-associated finding was reported; it is not a phenotypic AST interpretation.")
        search_text = st.text_input(
            "Search antibiotic / antimicrobial group / AMR determinant",
            placeholder="e.g. penicillin, beta-lactam, blaTEM, E. coli",
        )
        col1, col2 = st.columns(2)
        with col1:
            organism_filter = st.selectbox("Organism", ["All"] + sorted([x for x in results["Organism"].astype(str).unique() if x.strip()]))
        with col2:
            sample_filter = st.selectbox("Sample ID", ["All"] + sorted(results["Sample ID"].astype(str).unique()))

        filtered = results.copy()
        if search_text.strip():
            q = search_text.strip().lower()
            searchable = ["Antibiotic / antimicrobial association", "Class", "Subclass", "AMR determinant", "Mechanism / Product", "Organism", "Sample ID"]
            mask = filtered[searchable].astype(str).apply(lambda c: c.str.lower().str.contains(q, regex=False)).any(axis=1)
            filtered = filtered[mask]
        if organism_filter != "All":
            filtered = filtered[filtered["Organism"].astype(str) == organism_filter]
        if sample_filter != "All":
            filtered = filtered[filtered["Sample ID"].astype(str) == sample_filter]

        st.write(f"**{len(filtered)} matching genomic finding(s)**")
        st.dataframe(filtered[["Sample ID", "Organism", "Antibiotic / antimicrobial association", "AMR determinant", "Class", "Subclass", "Identity", "Coverage", "Interpretation"]], use_container_width=True, hide_index=True)

with tabs[3]:
    st.header("Cross-Sample Comparison")
    if results.empty:
        st.info("Upload and organize results first.")
    else:
        st.caption("Compare genomic AMR-associated findings across samples. Filters below are based on the uploaded AMRFinderPlus records; they do not convert genomic detection into a phenotypic AST result.")
        c1, c2, c3, c4 = st.columns(4)
        class_choices = sorted([x for x in results["Class"].astype(str).unique() if x.strip()])
        sub_choices = sorted([x for x in results["Subclass"].astype(str).unique() if x.strip()])
        drug_choices = sorted([x for x in results["Antibiotic / antimicrobial association"].astype(str).unique() if x.strip()])
        org_choices = sorted([x for x in results["Organism"].astype(str).unique() if x.strip()])
        with c1:
            selected_class = st.selectbox("AMR class", ["All"] + class_choices)
        with c2:
            selected_sub = st.selectbox("Subclass / antibiotic group", ["All"] + sub_choices + [x for x in drug_choices if x not in sub_choices])
        with c3:
            selected_org = st.selectbox("Organism", ["All"] + org_choices)
        with c4:
            determinant_query = st.text_input("AMR determinant", placeholder="e.g. blaTEM")

        view = results.copy()
        if selected_class != "All":
            view = view[view["Class"].astype(str) == selected_class]
        if selected_sub != "All":
            view = view[(view["Subclass"].astype(str) == selected_sub) | (view["Antibiotic / antimicrobial association"].astype(str) == selected_sub)]
        if selected_org != "All":
            view = view[view["Organism"].astype(str) == selected_org]
        if determinant_query.strip():
            q = determinant_query.strip().lower()
            view = view[view["AMR determinant"].astype(str).str.lower().str.contains(q, regex=False)]

        st.write(f"**{len(view)} matching genomic finding(s) across {view['Sample ID'].nunique()} sample(s)**")
        summary = (
            view.groupby(["AMR determinant", "Class", "Subclass", "Antibiotic / antimicrobial association"], dropna=False)
            .agg(
                Samples=("Sample ID", lambda x: ", ".join(sorted(set(map(str, x))))),
                Organisms=("Organism", lambda x: ", ".join(sorted(set([str(v) for v in x if str(v).strip()])))),
                Count=("Sample ID", "count"),
            )
            .reset_index()
            .sort_values("Count", ascending=False)
        )
        st.dataframe(summary, use_container_width=True, hide_index=True)
        st.markdown("### Matching sample-level records")
        st.dataframe(view, use_container_width=True, hide_index=True)

with tabs[4]:
    st.header("Sample / Organism Profile")
    if results.empty:
        st.info("Upload and organize results first.")
    else:
        sample_options = sorted(results["Sample ID"].astype(str).unique())
        selected_sample = st.selectbox("Select Sample ID", sample_options)
        profile = results[results["Sample ID"].astype(str) == selected_sample]

        organisms = [x for x in profile["Organism"].astype(str).unique() if x.strip()]
        st.write("**Organism:**", organisms[0] if organisms else "Not provided in the uploaded metadata/result")

        st.dataframe(profile, use_container_width=True, hide_index=True)

st.markdown("---")
st.caption(
    "AMRIVA is a research/educational results-organizer prototype. "
    "AMRFinderPlus remains the genomic AMR detection engine. "
    "Genomic detection alone does not establish phenotypic antimicrobial susceptibility."
)
