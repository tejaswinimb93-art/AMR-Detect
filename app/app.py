
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
    name = uploaded.name.lower()
    if name.endswith(".csv"):
        return pd.read_csv(uploaded)
    # AMRFinderPlus standard output is TSV.
    return pd.read_csv(uploaded, sep="\t", dtype=str, keep_default_na=False)

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
        r"(?:target[_\s-]*acc|target)[=:]([A-Za-z0-9_.-]+)",
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
            contig_match = re.search(r"(?:contig|sequence|seq)[=:]([^|;\s]+)", header, flags=re.I)
            # MicroBIGG-E/NCBI-style FASTA headers may use the first token as
            # the contig accession and carry sample/organism metadata after |.
            # Use that first token when no explicit contig= field exists.
            contig = contig_match.group(1) if contig_match else (header.split("|")[0].split()[0] if header else f"record_{i}")
            rows.append({
                "Sample ID": sid,
                "Organism": organism_from_header(header),
                "Sequence / Contig": contig,
                "_full_sequence": re.sub(r"\s+", "", seq).upper(),
                "_header": header,
            })
    return pd.DataFrame(rows)


def parse_metadata_files(files):
    """Read optional NCBI Isolate Browser / metadata exports.
    Only explicit metadata is used; AMRIVA never guesses an organism.
    """
    frames = []
    for f in files or []:
        try:
            name = f.name.lower()
            if name.endswith(".csv"):
                frames.append(pd.read_csv(f, dtype=str, keep_default_na=False))
            else:
                frames.append(pd.read_csv(f, sep="\t", dtype=str, keep_default_na=False))
        except Exception:
            continue
    if not frames:
        return pd.DataFrame()
    df = pd.concat(frames, ignore_index=True)
    df.columns = [str(c).strip() for c in df.columns]
    return df

def metadata_value(df, row, *aliases):
    c = find_col(df, *aliases) if df is not None and not df.empty else None
    return clean(row[c]) if c else ""

def enrich_fasta_with_metadata(fasta_df, metadata_df):
    if fasta_df is None or fasta_df.empty or metadata_df is None or metadata_df.empty:
        return fasta_df

    out = fasta_df.copy()
    sample_col = find_col(metadata_df, "Sample ID", "sample_id", "sample", "isolate",
                           "isolate_id", "sample name", "sample_name", "strain")
    org_col = find_col(metadata_df, "Organism", "organism_name", "taxon_name",
                       "species", "taxname", "taxgroup_name")
    asm_col = find_col(metadata_df, "Assembly", "asm_acc", "assembly_accession",
                       "assembly accession", "accession")
    if not org_col:
        return out

    lookup = {}
    for _, m in metadata_df.iterrows():
        sid = clean(m[sample_col]) if sample_col else ""
        asm = clean(m[asm_col]) if asm_col else ""
        org = clean(m[org_col])
        if org:
            if sid:
                lookup[("sample", sid.lower())] = org
            if asm:
                lookup[("assembly", asm.lower())] = org

    for i, r in out.iterrows():
        if clean(r.get("Organism")):
            continue
        sid = clean(r.get("Sample ID")).lower()
        header = clean(r.get("_header")).lower()
        seqname = clean(r.get("Sequence / Contig")).lower()
        org = lookup.get(("sample", sid), "")
        if not org:
            for key, value in lookup.items():
                if key[0] == "assembly" and key[1] and key[1] in header:
                    org = value
                    break
        if org:
            out.at[i, "Organism"] = org
    return out

def antibiotic_association(determinant, cls, subclass):
    """Human-readable antibiotic groups associated with the detected determinant.
    These are associations for exploration, not phenotype calls.
    """
    d = clean(determinant).lower()
    exact = {
        "blatem": "Penicillins / cephalosporins (beta-lactams)",
        "blaoxa": "Beta-lactams",
        "blapdc": "Cephalosporins (beta-lactams)",
        "vang": "Vancomycin / glycopeptides",
        "qnrs": "Fluoroquinolones",
        "sul1": "Sulfonamides",
        "teta": "Tetracyclines",
        "aac(6')-ib": "Aminoglycosides",
    }
    if d in exact:
        return exact[d]

    c = clean(subclass or cls).lower()
    if "beta-lactam" in c:
        return "Beta-lactams"
    if "glycopeptide" in c:
        return "Glycopeptides"
    if "quinolone" in c:
        return "Fluoroquinolones"
    if "aminoglycoside" in c:
        return "Aminoglycosides"
    if "tetracycline" in c:
        return "Tetracyclines"
    if "sulfonamide" in c:
        return "Sulfonamides"
    if "macrolide" in c:
        return "Macrolides"
    if "phenicol" in c:
        return "Phenicol antibiotics"
    if "trimethoprim" in c:
        return "Trimethoprim"
    return clean(subclass) or clean(cls) or "Not reported by AMRFinderPlus"

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
    # AMRFinderPlus does not normally report a literal drug name for each hit.
    # Use an explicitly supplied drug first, otherwise use a transparent
    # determinant/class association for exploration.
    for key in ["Antibiotic", "Drug", "Antimicrobial", "Antimicrobial(s)",
                "Antibiotic / antimicrobial association"]:
        if key in row and clean(row[key]):
            return clean(row[key])
    return antibiotic_association(
        row.get("AMR determinant", ""),
        row.get("Class", ""),
        row.get("Subclass", ""),
    )


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


def standardize_amrfinder(df, fasta_df=None, metadata_df=None):
    df = df.copy()
    df.columns = [str(c).strip() for c in df.columns]

    # Important: some earlier AMRIVA demo files stored the actual technical
    # AMRFinderPlus fields in a raw_result JSON cell. Flatten that first.
    df = flatten_raw_results(df)

    # These aliases cover current AMRFinderPlus field names as well as
    # common AMRIVA/demo exports.
    c_sample = find_col(df, "Sample ID", "Sample", "SampleID", "sample_name", "sample name", "isolate", "isolate_id", "target_acc", "Target accession", "target accession", "biosample_acc", "BioSample", "biosample accession")
    c_org = find_col(df, "Organism", "Scientific name", "scientific_name", "Taxonomy", "Taxon", "tax_name", "tax-name", "Scientific_name")
    c_det = find_col(
        df, "Element symbol", "element_symbol", "elem_symbol", "Gene symbol", "AMR determinant",
        "determinant", "AMR gene", "Gene", "Element", "element symbol"
    )
    c_class = find_col(df, "Class")
    c_sub = find_col(df, "Subclass", "Sub-class")
    c_mech = find_col(
        df, "Mechanism", "Product", "Element name", "element_name", "Sequence name",
        "Gene name", "Mechanism/Product", "Element description"
    )
    c_method = find_col(df, "Method", "amr_method", "AMR method")
    c_ident = find_col(
        df, "% Identity to reference sequence", "% Identity to reference", "pct_ref_identity", "Percent identity",
        "% identity", "Identity"
    )
    c_cov = find_col(
        df, "% Coverage of reference sequence", "% Coverage of reference", "pct_ref_coverage", "Percent coverage",
        "% coverage", "Coverage"
    )
    c_contig = find_col(
        df, "Contig id", "Contig", "contig_acc", "contig accession", "Sequence name",
        "Sequence ID", "Sequence", "Sequence/Contig", "sequence accession"
    )
    c_start = find_col(df, "Start", "start_on_contig", "target start", "Target start")
    c_stop = find_col(df, "Stop", "end_on_contig", "stop_on_contig", "target stop", "Target stop")
    c_ref_acc = find_col(
        df, "Accession of closest sequence", "Closest reference accession", "closest_reference_acc", "closest-ref-accession", "Reference accession", "Accession"
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

    if metadata_df is not None and not metadata_df.empty:
        fasta_df = enrich_fasta_with_metadata(fasta_df, metadata_df)

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

    # If the AMRFinderPlus table has no sample column, map contigs to a
    # FASTA-derived sample only when that contig belongs to exactly one sample.
    contig_to_samples = {}
    if fasta_df is not None and not fasta_df.empty:
        for _, fr in fasta_df.iterrows():
            ct = clean(fr.get("Sequence / Contig"))
            sidf = clean(fr.get("Sample ID"))
            if ct and sidf:
                contig_to_samples.setdefault(ct, set()).add(sidf)

    for idx, r in df.iterrows():
        raw = " | ".join(
            f"{c}={clean(r[c])}" for c in df.columns if clean(r[c]) != ""
        )

        sid = clean(r[c_sample]) if c_sample else ""
        org = clean(r[c_org]) if c_org else ""

        if not sid and c_contig:
            ct = clean(r[c_contig])
            candidates = sorted(contig_to_samples.get(ct, set()))
            if len(candidates) == 1:
                sid = candidates[0]

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
    metadata_files = st.file_uploader(
        "Optional NCBI isolate metadata CSV/TSV (for Sample ID + Organism)",
        type=["csv", "tsv", "txt"],
        accept_multiple_files=True,
        key="ncbi_metadata",
        help="Use this when the AMRFinderPlus TSV itself does not contain organism/sample metadata.",
    )
    if st.button("Organize AMRFinderPlus results", type="primary"):
        if not uploaded_result:
            st.error("Upload an AMRFinderPlus TSV/CSV first.")
        else:
            try:
                raw_df = read_table(uploaded_result)
                fdf = parse_fasta_files(fasta_files) if fasta_files else pd.DataFrame()
                mdf = parse_metadata_files(metadata_files) if metadata_files else pd.DataFrame()
                st.session_state.fasta_df = fdf
                st.session_state.results = standardize_amrfinder(raw_df, fdf, mdf)
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
    metadata_files = st.file_uploader(
        "Optional NCBI isolate metadata CSV/TSV",
        type=["csv", "tsv", "txt"],
        accept_multiple_files=True,
    )
    if st.button("Process genomic results", type="primary"):
        if not fasta_files or not result_file:
            st.error("Upload both FASTA file(s) and the corresponding AMRFinderPlus result.")
        else:
            try:
                fdf = parse_fasta_files(fasta_files)
                rdf = read_table(result_file)
                mdf = parse_metadata_files(metadata_files) if metadata_files else pd.DataFrame()
                st.session_state.fasta_df = fdf
                st.session_state.results = standardize_amrfinder(rdf, fdf, mdf)
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
            "Organism can come from explicit NCBI metadata or FASTA metadata, and nucleotide sequence extraction "
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
        st.caption("Use an antibiotic/group filter to see which samples contain associated genomic AMR determinants. This is a genomic association view, not an AST phenotype call.")
        col1, col2, col3 = st.columns(3)
        antibiotic_options = sorted([x for x in results["Antibiotic / antimicrobial association"].astype(str).unique() if x.strip()])
        class_options = sorted([x for x in results["Class"].astype(str).unique() if x.strip()])
        organism_options = sorted([x for x in results["Organism"].astype(str).unique() if x.strip() and "Not provided" not in x])

        with col1:
            selected_antibiotic = st.selectbox("Antibiotic / group", ["All"] + antibiotic_options)
        with col2:
            selected_class = st.selectbox("AMR class", ["All"] + class_options)
        with col3:
            selected_organism = st.selectbox("Organism", ["All"] + organism_options)

        search_text = st.text_input(
            "Search determinant, mechanism, sample ID or contig",
            placeholder="e.g. blaTEM, Sample_001",
        )

        filtered = results.copy()
        if selected_antibiotic != "All":
            filtered = filtered[filtered["Antibiotic / antimicrobial association"] == selected_antibiotic]
        if selected_class != "All":
            filtered = filtered[filtered["Class"] == selected_class]
        if selected_organism != "All":
            filtered = filtered[filtered["Organism"] == selected_organism]
        if search_text.strip():
            q = search_text.strip().lower()
            mask = filtered.apply(
                lambda row: row.astype(str).str.lower().str.contains(q, regex=False).any(),
                axis=1,
            )
            filtered = filtered[mask]

        st.write(f"**{len(filtered)} matching genomic finding(s)**")
        if not filtered.empty:
            sample_list = sorted(set(filtered["Sample ID"].astype(str)))
            st.info("Samples in this result: " + ", ".join(sample_list))
        st.dataframe(filtered, use_container_width=True, hide_index=True)


with tabs[3]:
    st.header("Cross-Sample Comparison")
    if results.empty:
        st.info("Upload and organize results first.")
    else:
        st.caption("Filter across samples by antibiotic association, AMR class, organism, determinant and method.")
        c1, c2, c3 = st.columns(3)
        antibiotics = sorted([x for x in results["Antibiotic / antimicrobial association"].astype(str).unique() if x.strip()])
        classes = sorted([x for x in results["Class"].astype(str).unique() if x.strip()])
        organisms = sorted([x for x in results["Organism"].astype(str).unique() if x.strip() and "Not provided" not in x])
        with c1:
            fa = st.selectbox("Filter by antibiotic / group", ["All"] + antibiotics, key="cross_antibiotic")
        with c2:
            fc = st.selectbox("Filter by AMR class", ["All"] + classes, key="cross_class")
        with c3:
            fo = st.selectbox("Filter by organism", ["All"] + organisms, key="cross_organism")

        c4, c5 = st.columns(2)
        determinants = sorted([x for x in results["AMR determinant"].astype(str).unique() if x.strip() and x != "Not reported"])
        methods = sorted([x for x in results["Method"].astype(str).unique() if x.strip() and x != "Not reported"])
        with c4:
            fd = st.selectbox("Filter by determinant", ["All"] + determinants, key="cross_det")
        with c5:
            fm = st.selectbox("Filter by method", ["All"] + methods, key="cross_method")

        view = results.copy()
        if fa != "All":
            view = view[view["Antibiotic / antimicrobial association"] == fa]
        if fc != "All":
            view = view[view["Class"] == fc]
        if fo != "All":
            view = view[view["Organism"] == fo]
        if fd != "All":
            view = view[view["AMR determinant"] == fd]
        if fm != "All":
            view = view[view["Method"] == fm]

        summary = (
            view.groupby(
                ["Antibiotic / antimicrobial association", "AMR determinant", "Class", "Subclass"],
                dropna=False
            )
            .agg(
                Samples=("Sample ID", lambda x: ", ".join(sorted(set(map(str, x))))),
                Organisms=("Organism", lambda x: ", ".join(sorted(set(
                    [str(v) for v in x if str(v).strip()]
                )))),
                Count=("Sample ID", "count"),
            )
            .reset_index()
            .sort_values("Count", ascending=False)
        )
        st.write(f"**{len(view)} genomic finding(s) across {view['Sample ID'].nunique()} sample(s)**")
        st.dataframe(summary, use_container_width=True, hide_index=True)


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
