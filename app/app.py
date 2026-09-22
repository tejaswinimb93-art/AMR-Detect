import io
import os
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

import pandas as pd
import streamlit as st

st.set_page_config(page_title="AMRIVA", page_icon="🧬", layout="wide")

st.markdown("""
<style>
.block-container {max-width: 1500px; padding-top: 2rem;}
.hero {padding:28px 32px;border-radius:20px;margin-bottom:22px;background:linear-gradient(135deg,#123c69,#0f766e);color:white;}
.hero h1 {font-size:48px;margin:0 0 5px 0}.hero p {font-size:17px;margin:4px 0}
</style>
""", unsafe_allow_html=True)

st.markdown("""
<div class="hero"><h1>🧬 AMRIVA</h1>
<p><b>AMRFinderPlus Results Organizer & Explorer</b></p>
<p>Organize genomic AMR findings across samples, organisms, antibiotics and resistance determinants.</p></div>
""", unsafe_allow_html=True)

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
    return "" if s.lower() in {"nan", "none", "null", "na", "n/a"} else s


def norm_col(c):
    return re.sub(r"[^a-z0-9]+", "", str(c).strip().lower())


def find_col(df, *names):
    lookup = {norm_col(c): c for c in df.columns}
    wanted = [norm_col(x) for x in names]
    for n in wanted:
        if n in lookup:
            return lookup[n]
    # Only use substring matching when it is unambiguous.
    for n in wanted:
        hits = [c for c in df.columns if n and (n in norm_col(c) or norm_col(c) in n)]
        if len(hits) == 1:
            return hits[0]
    return None


def read_table(uploaded):
    uploaded.seek(0)
    raw = uploaded.getvalue()
    name = uploaded.name.lower()
    if name.endswith(".csv"):
        return pd.read_csv(io.BytesIO(raw), dtype=str, keep_default_na=False)
    return pd.read_csv(io.BytesIO(raw), sep="\t", dtype=str, keep_default_na=False)


def fasta_records(text):
    records, header, seq = [], None, []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(seq)))
            header, seq = line[1:].strip(), []
        else:
            seq.append(line)
    if header is not None:
        records.append((header, "".join(seq)))
    return records


def sample_id_from_header(header, fallback):
    patterns = [
        r"(?:sample[_\s-]*id|sample|isolate|assembly)[=:]([A-Za-z0-9_.:-]+)",
        r"(?:sample[_\s-]*id|sample|isolate|assembly)\s+([A-Za-z0-9_.:-]+)",
    ]
    for p in patterns:
        m = re.search(p, header, re.I)
        if m:
            return m.group(1)
    return fallback


def organism_from_header(header):
    patterns = [
        r"(?:organism|species|taxon)[=:]([^|;]+)",
        r"\[([A-Z][A-Za-z]+(?:\s+[a-z][A-Za-z0-9_.-]+){0,2})\]",
    ]
    for p in patterns:
        m = re.search(p, header, re.I)
        if m:
            return m.group(1).strip()
    return ""


def parse_fasta_files(files):
    rows = []
    for f in files or []:
        text = f.getvalue().decode("utf-8", errors="replace")
        fallback = Path(f.name).stem
        for i, (header, seq) in enumerate(fasta_records(text), 1):
            sid = sample_id_from_header(header, fallback)
            contig = header.split()[0] if header else f"record_{i}"
            rows.append({
                "Sample ID": sid,
                "Organism": organism_from_header(header),
                "Sequence / Contig": contig,
                "_full_sequence": re.sub(r"\s+", "", seq).upper(),
                "_header": header,
                "_source_file": f.name,
            })
    return pd.DataFrame(rows)


def reverse_complement(seq):
    return seq.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def extract_nucleotide(row, fasta_df):
    if fasta_df is None or fasta_df.empty:
        return ""
    contig = clean(row.get("Sequence / Contig"))
    sid = clean(row.get("Sample ID"))
    try:
        start = int(float(clean(row.get("Start"))))
        stop = int(float(clean(row.get("Stop"))))
    except Exception:
        return ""
    if not contig:
        return ""

    candidates = fasta_df.copy()
    # Exact sample+contig first; then contig-only if unique.
    exact = candidates[(candidates["Sample ID"].astype(str) == sid) & (candidates["Sequence / Contig"].astype(str) == contig)]
    if exact.empty:
        exact = candidates[candidates["Sequence / Contig"].astype(str) == contig]
    if exact.empty:
        # Some files contain a pipe-delimited AMRFinder contig identifier.
        short = contig.split("|")[0]
        exact = candidates[candidates["Sequence / Contig"].astype(str).str.split("|").str[0] == short]
    if exact.empty:
        return ""

    seq = clean(exact.iloc[0]["_full_sequence"])
    lo, hi = min(start, stop), max(start, stop)
    if lo < 1 or hi > len(seq):
        return ""
    frag = seq[lo - 1:hi]
    return reverse_complement(frag) if start > stop else frag


def choose_sample(row, c_sample, fasta_df, fallback):
    if c_sample:
        value = clean(row[c_sample])
        if value:
            return value
    # Match contig to FASTA metadata where possible.
    c_contig = find_col(pd.DataFrame([row]), "Contig id", "Contig", "Sequence id", "Sequence")
    contig = clean(row[c_contig]) if c_contig else ""
    if fasta_df is not None and not fasta_df.empty and contig:
        m = fasta_df[fasta_df["Sequence / Contig"].astype(str) == contig]
        if len(m) == 1 and clean(m.iloc[0]["Sample ID"]):
            return clean(m.iloc[0]["Sample ID"])
    return fallback


def derive_antimicrobial(row):
    # AMRFinderPlus Subclass is the most useful antibiotic-level association.
    # Do not pretend a class is a specific drug.
    for key in ["Antibiotic", "Drug", "Antimicrobial", "Antimicrobial(s)"]:
        if clean(row.get(key, "")):
            return clean(row[key])
    sub = clean(row.get("Subclass", ""))
    cls = clean(row.get("Class", ""))
    if sub:
        return sub.replace("_", " ")
    if cls:
        return cls.replace("_", " ")
    return "Not specified in AMRFinderPlus output"


def interpretation(row):
    det, cls, sub = clean(row.get("AMR determinant")), clean(row.get("Class")), clean(row.get("Subclass"))
    ident, cov = clean(row.get("Identity")), clean(row.get("Coverage"))
    if not det:
        return "No AMR determinant reported for this record."
    s = f"Genomic AMR determinant {det} was detected"
    if cls:
        s += f" in the {cls} class"
    if sub:
        s += f" ({sub})"
    qc = []
    if ident: qc.append(f"identity {ident}%")
    if cov: qc.append(f"coverage {cov}%")
    if qc: s += " with " + " and ".join(qc)
    return s + ". Genomic detection does not by itself establish phenotypic susceptibility or resistance."


def standardize_amrfinder(df, fasta_df=None, source_name="uploaded_result"):
    df = df.copy()
    df.columns = [str(c).strip() for c in df.columns]

    # These aliases match the actual NCBI AMRFinderPlus TSV, including the current
    # test output fields such as Element symbol, Element name, Contig id,
    # % Coverage of reference and % Identity to reference.
    c_sample = find_col(df, "Sample ID", "SampleID", "Sample", "Name", "Isolate", "Assembly")
    c_org = find_col(df, "Organism", "Species", "Scientific name", "Taxon")
    c_det = find_col(df, "Element symbol", "AMR determinant", "AMR gene", "Gene symbol", "Gene", "Element")
    c_name = find_col(df, "Element name", "Mechanism / Product", "Mechanism", "Product")
    c_class = find_col(df, "Class")
    c_sub = find_col(df, "Subclass", "Sub-class")
    c_method = find_col(df, "Method", "AMR method", "amr_method")
    c_ident = find_col(df, "% Identity to reference", "% identity", "Identity", "Percent identity", "pct_ref_identity")
    c_cov = find_col(df, "% Coverage of reference", "% coverage", "Coverage", "Percent coverage", "pct_ref_coverage")
    c_contig = find_col(df, "Contig id", "Contig", "Sequence id", "Sequence/Contig", "Sequence")
    c_start = find_col(df, "Start", "start_on_contig")
    c_stop = find_col(df, "Stop", "End", "end_on_contig")
    c_ref = find_col(df, "Closest reference accession", "Reference accession", "Reference", "Accession", "closest_ref_accession")
    c_type = find_col(df, "Type", "Element type")
    c_subtype = find_col(df, "Subtype", "Element subtype", "Element Subtype")
    c_ant = find_col(df, "Antibiotic", "Drug", "Antimicrobial", "Antimicrobial(s)")

    fallback_sample = Path(source_name).stem
    out = []
    for idx, r in df.iterrows():
        raw = " | ".join(f"{c}={clean(r[c])}" for c in df.columns if clean(r[c]) != "")
        sid = choose_sample(r, c_sample, fasta_df, fallback_sample)
        org = clean(r[c_org]) if c_org else ""

        # If AMRFinderPlus was run with --name, Name becomes Sample ID.
        # Otherwise, a matching FASTA header/file supplies the sample identity.
        contig = clean(r[c_contig]) if c_contig else ""
        if not org and fasta_df is not None and not fasta_df.empty and contig:
            matches = fasta_df[fasta_df["Sequence / Contig"].astype(str) == contig]
            if not matches.empty:
                vals = [clean(x) for x in matches["Organism"].tolist() if clean(x)]
                if vals: org = vals[0]

        rec = {
            "Sample ID": sid,
            "Organism": org or "Not provided in AMRFinderPlus/FASTA metadata",
            "Antibiotic / antimicrobial association": clean(r[c_ant]) if c_ant else "",
            "AMR determinant": clean(r[c_det]) if c_det else "",
            "Class": clean(r[c_class]) if c_class else "",
            "Subclass": clean(r[c_sub]) if c_sub else "",
            "Mechanism / Product": clean(r[c_name]) if c_name else "",
            "Method": clean(r[c_method]) if c_method else "",
            "Identity": clean(r[c_ident]) if c_ident else "",
            "Coverage": clean(r[c_cov]) if c_cov else "",
            "Sequence / Contig": contig,
            "Nucleotide Sequence": "",
            "Start": clean(r[c_start]) if c_start else "",
            "Stop": clean(r[c_stop]) if c_stop else "",
            "Reference": clean(r[c_ref]) if c_ref else "",
            "Interpretation": "",
            "Raw Result": raw,
        }
        if not rec["Antibiotic / antimicrobial association"]:
            rec["Antibiotic / antimicrobial association"] = derive_antimicrobial(rec)
        rec["Nucleotide Sequence"] = extract_nucleotide(rec, fasta_df)
        rec["Interpretation"] = interpretation(rec)
        out.append(rec)

    return pd.DataFrame(out, columns=DISPLAY_COLUMNS)


def run_amrfinder_on_fasta(files):
    exe = shutil.which("amrfinder")
    if not exe:
        return None, "AMRFinderPlus executable was not found in this deployment. Upload the real AMRFinderPlus TSV instead."
    temp_root = Path(tempfile.mkdtemp(prefix="amriva_"))
    try:
        chunks = []
        for f in files:
            path = temp_root / Path(f.name).name
            path.write_bytes(f.getvalue())
            cmd = [exe, "-n", str(path), "--plus"]
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=180)
            if result.returncode != 0:
                return None, result.stderr.strip() or "AMRFinderPlus returned an error."
            lines = [x for x in result.stdout.splitlines() if x.strip()]
            if lines:
                if not chunks:
                    chunks.extend(lines)
                else:
                    chunks.extend(x for x in lines[1:] if x != chunks[0])
        return ("\n".join(chunks), None) if chunks else (None, "AMRFinderPlus returned no records.")
    except subprocess.TimeoutExpired:
        return None, "AMRFinderPlus timed out."
    finally:
        shutil.rmtree(temp_root, ignore_errors=True)


if "results" not in st.session_state:
    st.session_state.results = pd.DataFrame(columns=DISPLAY_COLUMNS)
if "fasta_df" not in st.session_state:
    st.session_state.fasta_df = pd.DataFrame()

st.sidebar.title("AMRIVA")
st.sidebar.caption("AMRFinderPlus results organizer")
mode = st.sidebar.radio("Input mode", ["AMRFinderPlus result file", "FASTA + AMRFinderPlus", "FASTA only"])

if mode == "AMRFinderPlus result file":
    uploaded_result = st.file_uploader("Upload AMRFinderPlus TSV/CSV", type=["tsv", "txt", "csv"], key="result_upload")
    fasta_files = st.file_uploader("Optional corresponding FASTA files (recommended for organism + nucleotide sequence)", type=["fa", "fasta", "fna"], accept_multiple_files=True, key="optional_fasta")
    if st.button("Organize AMRFinderPlus results", type="primary"):
        if not uploaded_result:
            st.error("Upload the AMRFinderPlus TSV/CSV first.")
        else:
            try:
                raw = read_table(uploaded_result)
                fdf = parse_fasta_files(fasta_files) if fasta_files else pd.DataFrame()
                st.session_state.fasta_df = fdf
                st.session_state.results = standardize_amrfinder(raw, fdf, uploaded_result.name)
                st.success(f"Organized {len(st.session_state.results)} AMRFinderPlus record(s).")
            except Exception as e:
                st.error(f"Could not read this result file: {e}")

elif mode == "FASTA + AMRFinderPlus":
    fasta_files = st.file_uploader("Upload corresponding FASTA file(s)", type=["fa", "fasta", "fna"], accept_multiple_files=True)
    result_file = st.file_uploader("Upload the AMRFinderPlus TSV/CSV", type=["tsv", "txt", "csv"])
    if st.button("Process genomic results", type="primary"):
        if not fasta_files or not result_file:
            st.error("Upload both FASTA file(s) and the corresponding AMRFinderPlus result.")
        else:
            try:
                fdf = parse_fasta_files(fasta_files)
                rdf = read_table(result_file)
                st.session_state.fasta_df = fdf
                st.session_state.results = standardize_amrfinder(rdf, fdf, result_file.name)
                st.success(f"Processed {len(st.session_state.results)} genomic finding(s).")
            except Exception as e:
                st.error(f"Could not process the files: {e}")

else:
    fasta_files = st.file_uploader("Upload one or more FASTA files", type=["fa", "fasta", "fna"], accept_multiple_files=True)
    st.info("FASTA-only mode runs the installed AMRFinderPlus executable if available; otherwise use a real AMRFinderPlus TSV.")
    if st.button("Run AMRFinderPlus", type="primary"):
        if not fasta_files:
            st.error("Upload at least one FASTA file.")
        else:
            fdf = parse_fasta_files(fasta_files)
            st.session_state.fasta_df = fdf
            tsv, err = run_amrfinder_on_fasta(fasta_files)
            if err:
                st.error(err)
            else:
                try:
                    rdf = pd.read_csv(io.StringIO(tsv), sep="\t", dtype=str, keep_default_na=False)
                    st.session_state.results = standardize_amrfinder(rdf, fdf, "FASTA analysis")
                    st.success(f"AMRFinderPlus returned {len(st.session_state.results)} finding(s).")
                except Exception as e:
                    st.error(f"Could not parse AMRFinderPlus output: {e}")

results = st.session_state.results

tabs = st.tabs(["🏠 Home", "🧬 Genomic Analysis", "🔎 Antibiotic Explorer", "🔬 Cross-Sample Comparison", "🧫 Sample / Organism Profile"])

with tabs[0]:
    st.markdown("## What problem does AMRIVA solve?")
    st.write("AMRFinderPlus produces detailed genomic findings. AMRIVA organizes those findings by sample, organism, antimicrobial association, determinant, genomic location and quality metrics, then makes the results searchable and comparable across samples.")
    st.markdown("### Core relationship")
    st.code("Sample ↔ Organism ↔ Antibiotic association ↔ AMR determinant ↔ Genomic location")
    st.info("AMRIVA is an organization/exploration layer around AMRFinderPlus. Genomic detection alone does not establish phenotypic resistance.")

with tabs[1]:
    st.header("Genomic Analysis")
    if results.empty:
        st.info("Upload an AMRFinderPlus result file first.")
    else:
        st.metric("AMRFinderPlus findings", len(results))
        st.dataframe(results, use_container_width=True, hide_index=True)
        st.download_button("⬇️ Download organized genomic results", results.to_csv(index=False).encode("utf-8"), "AMRIVA_organized_genomic_results.csv", "text/csv")

with tabs[2]:
    st.header("Antibiotic Explorer")
    if results.empty:
        st.info("Upload and organize results first.")
    else:
        col1, col2 = st.columns(2)
        with col1:
            antibiotics = sorted([x for x in results["Antibiotic / antimicrobial association"].astype(str).unique() if x.strip()])
            antibiotic = st.selectbox("Filter by antibiotic / antimicrobial", ["All"] + antibiotics)
        with col2:
            organism = st.selectbox("Filter by organism", ["All"] + sorted(results["Organism"].astype(str).unique().tolist()))
        search = st.text_input("Search determinant, class, subclass or mechanism", placeholder="e.g. penicillin, BETA-LACTAM, blaTEM")
        filtered = results.copy()
        if antibiotic != "All": filtered = filtered[filtered["Antibiotic / antimicrobial association"] == antibiotic]
        if organism != "All": filtered = filtered[filtered["Organism"] == organism]
        if search.strip():
            q = search.lower().strip()
            filtered = filtered[filtered.apply(lambda r: r.astype(str).str.lower().str.contains(q, regex=False).any(), axis=1)]
        st.write(f"**{len(filtered)} matching finding(s)**")
        st.dataframe(filtered[["Sample ID","Organism","Antibiotic / antimicrobial association","AMR determinant","Class","Subclass","Method","Identity","Coverage"]], use_container_width=True, hide_index=True)
        st.caption("The explorer shows which uploaded samples/organisms contain findings associated with the selected antimicrobial class/subclass. It does not infer AST resistance.")

with tabs[3]:
    st.header("Cross-Sample Comparison")
    if results.empty:
        st.info("Upload and organize results first.")
    else:
        c1, c2, c3 = st.columns(3)
        with c1:
            antibiotic = st.selectbox("Filter by antibiotic / antimicrobial", ["All"] + sorted([x for x in results["Antibiotic / antimicrobial association"].astype(str).unique() if x.strip()]), key="cross_antibiotic")
        with c2:
            amr_class = st.selectbox("Filter by AMR class", ["All"] + sorted([x for x in results["Class"].astype(str).unique() if x.strip()]), key="cross_class")
        with c3:
            organism = st.selectbox("Filter by organism", ["All"] + sorted(results["Organism"].astype(str).unique().tolist()), key="cross_organism")
        c4, c5 = st.columns(2)
        with c4:
            method = st.selectbox("Filter by method", ["All"] + sorted([x for x in results["Method"].astype(str).unique() if x.strip()]), key="cross_method")
        with c5:
            determinant = st.selectbox("Filter by AMR determinant", ["All"] + sorted([x for x in results["AMR determinant"].astype(str).unique() if x.strip()]), key="cross_det")

        view = results.copy()
        if antibiotic != "All": view = view[view["Antibiotic / antimicrobial association"] == antibiotic]
        if amr_class != "All": view = view[view["Class"] == amr_class]
        if organism != "All": view = view[view["Organism"] == organism]
        if method != "All": view = view[view["Method"] == method]
        if determinant != "All": view = view[view["AMR determinant"] == determinant]

        summary = (view.groupby(["AMR determinant","Class","Subclass"], dropna=False)
                   .agg(Samples=("Sample ID", lambda x: ", ".join(sorted(set(map(str,x))))),
                        Organisms=("Organism", lambda x: ", ".join(sorted(set(map(str,x))))),
                        Finding_Count=("Sample ID","count"))
                   .reset_index().sort_values("Finding_Count", ascending=False))
        st.write(f"**{len(view)} finding(s) across {view['Sample ID'].nunique()} sample(s)**")
        st.dataframe(summary, use_container_width=True, hide_index=True)
        st.download_button("⬇️ Download comparison", summary.to_csv(index=False).encode("utf-8"), "AMRIVA_cross_sample_comparison.csv", "text/csv")

with tabs[4]:
    st.header("Sample / Organism Profile")
    if results.empty:
        st.info("Upload and organize results first.")
    else:
        selected_sample = st.selectbox("Select Sample ID", sorted(results["Sample ID"].astype(str).unique()))
        profile = results[results["Sample ID"].astype(str) == selected_sample]
        organisms = [x for x in profile["Organism"].astype(str).unique() if x.strip()]
        st.write("**Organism:**", organisms[0] if organisms else "Not provided")
        st.write("**Finding count:**", len(profile))
        st.dataframe(profile, use_container_width=True, hide_index=True)

st.markdown("---")
st.caption("AMRIVA is a research/educational results-organizer prototype. AMRFinderPlus remains the genomic AMR detection engine. Genomic detection alone does not establish phenotypic antimicrobial susceptibility.")
