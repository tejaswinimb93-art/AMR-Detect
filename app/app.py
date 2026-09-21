"""
AMRIVA — AMRFinderPlus Results Organizer & Explorer
===================================================

Run:
    streamlit run amriva_app.py

Purpose:
    Organize and explore actual AMRFinderPlus output across one or more samples.

Important scientific behavior:
    - AMRIVA does not invent resistance genes, organisms, antibiotics, or sequences.
    - AMRFinderPlus source columns are preserved.
    - "Antibiotic" is populated only from an explicit antibiotic/drug column if one
      exists in the uploaded data. Otherwise the AMRFinderPlus subclass/class remains
      the antimicrobial association supplied by AMRFinderPlus.
    - Nucleotide Sequence is extracted from the uploaded FASTA using the reported
      contig + start/stop coordinates. If the required source sequence is missing,
      it stays unavailable.
"""

from __future__ import annotations

import io
import json
import re
from pathlib import Path
from typing import Iterable

import pandas as pd
import streamlit as st


# ---------------------------------------------------------------------------
# AMRFinderPlus field handling
# ---------------------------------------------------------------------------

ALIASES = {
    "sample": ["sample", "sample id", "sample_id"],
    "organism": ["organism", "scientific name", "species"],
    "antibiotic": ["antibiotic", "drug", "antimicrobial"],
    "determinant": [
        "element symbol", "amr determinant", "amr_determinant",
        "determinant", "gene",
    ],
    "element_name": ["element name"],
    "mechanism": [
        "mechanism", "mechanism/product", "mechanism / product",
        "resistance mechanism",
    ],
    "class": ["class"],
    "subclass": ["subclass"],
    "method": ["method"],
    "identity": [
        "% identity to reference", "identity (%)", "identity",
        "pct_ref_identity",
    ],
    "coverage": [
        "% coverage of reference", "coverage (%)", "coverage",
        "pct_ref_coverage",
    ],
    "contig": ["contig id", "contig", "sequence / contig", "sequence/contig",
               "contig_acc"],
    "start": ["start", "start on contig", "start_on_contig"],
    "stop": ["stop", "end", "end on contig", "end_on_contig"],
    "reference": [
        "closest reference accession", "reference", "reference accession",
        "closest_reference_acc",
    ],
    "raw_result": ["raw result", "raw_result"],
    "strand": ["strand", "orientation"],
}


def clean(value) -> str:
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except Exception:
        pass
    return str(value).strip()


def norm(value: str) -> str:
    return re.sub(r"\s+", " ", str(value).strip().lower())


def normalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    normalized = {norm(c): c for c in out.columns}
    rename = {}
    for target, aliases in ALIASES.items():
        for alias in aliases:
            source = normalized.get(norm(alias))
            if source is not None:
                rename[source] = target
                break
    return out.rename(columns=rename)


def parse_tsv(uploaded_file) -> pd.DataFrame:
    data = uploaded_file.getvalue()
    return pd.read_csv(
        io.BytesIO(data),
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )


# ---------------------------------------------------------------------------
# FASTA handling
# ---------------------------------------------------------------------------

def parse_fasta(uploaded_file) -> dict[str, str]:
    """Return {FASTA record ID: nucleotide sequence}."""
    text = uploaded_file.getvalue().decode("utf-8", errors="ignore")
    records = {}
    current_id = None
    chunks = []

    def flush():
        if current_id is not None:
            records[current_id] = re.sub(r"\s+", "", "".join(chunks)).upper()

    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            flush()
            header = line[1:].strip()
            current_id = header.split()[0]
            chunks = []
        else:
            chunks.append(line)

    flush()
    return records


def fasta_lookup(records: dict[str, str], contig: str) -> str:
    """Match AMRFinderPlus contig to FASTA record ID without guessing."""
    if not contig:
        return ""
    if contig in records:
        return records[contig]

    # Exact accession before whitespace/version differences.
    short = contig.split()[0]
    if short in records:
        return records[short]

    for key, seq in records.items():
        if key.split()[0] == short:
            return seq

    return ""


def reverse_complement(seq: str) -> str:
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(table)[::-1]


def extract_nucleotide_sequence(row, fasta_records: dict[str, str]) -> str:
    """
    Extract the actual determinant region from the analyzed FASTA.

    AMRFinderPlus reports Start/Stop on the contig. Coordinates are treated as
    1-based inclusive. For a '-' strand, the returned sequence is reverse
    complemented so it is oriented with the detected element.
    """
    contig = clean(row.get("contig", ""))
    start = clean(row.get("start", ""))
    stop = clean(row.get("stop", ""))

    if not contig or not start or not stop:
        return ""

    try:
        start_i = int(float(start))
        stop_i = int(float(stop))
    except ValueError:
        return ""

    if start_i < 1 or stop_i < 1:
        return ""

    seq = fasta_lookup(fasta_records, contig)
    if not seq:
        return ""

    lo, hi = min(start_i, stop_i), max(start_i, stop_i)
    if hi > len(seq):
        return ""

    extracted = seq[lo - 1:hi]
    strand = clean(row.get("strand", "")).lower()
    if strand in {"-", "minus", "reverse", "reverse strand"}:
        extracted = reverse_complement(extracted)

    return extracted


# ---------------------------------------------------------------------------
# Organization
# ---------------------------------------------------------------------------

def associated_group(row) -> str:
    """
    Preserve source-supported antimicrobial association only.

    No determinant-to-drug mapping is invented. If an explicit antibiotic
    field exists, it is used. Otherwise AMRFinderPlus subclass/class is used
    as the supplied category.
    """
    antibiotic = clean(row.get("antibiotic", ""))
    subclass = clean(row.get("subclass", ""))
    cls = clean(row.get("class", ""))

    if antibiotic:
        return antibiotic
    if subclass:
        return subclass
    if cls:
        return cls
    return ""


def organize_result(
    source: pd.DataFrame,
    sample_id: str = "",
    organism: str = "",
    fasta_records: dict[str, str] | None = None,
) -> pd.DataFrame:
    source = source.copy()
    source_columns = list(source.columns)
    df = normalize_columns(source)

    if sample_id:
        df["sample"] = sample_id
    if organism:
        df["organism"] = organism

    for col in [
        "sample", "organism", "antibiotic", "determinant", "element_name",
        "mechanism", "class", "subclass", "method", "identity", "coverage",
        "contig", "start", "stop", "reference", "strand",
    ]:
        if col not in df.columns:
            df[col] = ""

    if fasta_records is not None:
        df["Nucleotide Sequence"] = df.apply(
            lambda row: extract_nucleotide_sequence(row, fasta_records),
            axis=1,
        )
    else:
        df["Nucleotide Sequence"] = ""

    df["Associated antibiotic/group"] = df.apply(associated_group, axis=1)

    # Raw Result is a faithful copy of the original uploaded AMRFinderPlus row.
    if "raw_result" not in df.columns:
        df["raw_result"] = [
            json.dumps(
                {str(c): clean(row[c]) for c in source_columns},
                ensure_ascii=False,
            )
            for _, row in source.iterrows()
        ]

    preferred = [
        "sample", "organism", "antibiotic", "determinant", "class", "subclass",
        "mechanism", "element_name", "method", "identity", "coverage",
        "contig", "Nucleotide Sequence", "start", "stop", "reference",
        "raw_result",
    ]
    cols = [c for c in preferred if c in df.columns]
    cols += [c for c in df.columns if c not in cols]
    return df[cols]


def load_batch(uploaded_tsvs, fasta_by_sample=None) -> pd.DataFrame:
    frames = []
    fasta_by_sample = fasta_by_sample or {}

    for uploaded in uploaded_tsvs:
        # Filename stem is only a fallback display/sample ID. The user can
        # override it in the UI; no biological identity is inferred from it.
        sample_id = Path(uploaded.name).stem
        source = parse_tsv(uploaded)
        records = fasta_by_sample.get(sample_id, {})
        frames.append(organize_result(source, sample_id=sample_id,
                                       fasta_records=records))
    return pd.concat(frames, ignore_index=True, sort=False) if frames else pd.DataFrame()


# ---------------------------------------------------------------------------
# Search / profiles / comparison
# ---------------------------------------------------------------------------

SEARCH_COLUMNS = [
    "sample", "organism", "antibiotic", "Associated antibiotic/group",
    "determinant", "element_name", "mechanism", "class", "subclass",
    "method", "contig", "reference",
]


def search_all(df: pd.DataFrame, query: str) -> pd.DataFrame:
    query = clean(query)
    if not query or df.empty:
        return df.iloc[0:0].copy()

    mask = pd.Series(False, index=df.index)
    for col in SEARCH_COLUMNS:
        if col in df.columns:
            mask |= df[col].astype(str).str.contains(
                query, case=False, regex=False, na=False
            )
    return df.loc[mask].copy()


def exact_profile(df: pd.DataFrame, column: str, value: str) -> pd.DataFrame:
    if df.empty or column not in df.columns:
        return df.iloc[0:0].copy()
    value = clean(value).casefold()
    return df[df[column].astype(str).str.strip().str.casefold() == value].copy()


def comparison_table(df: pd.DataFrame, query: str) -> pd.DataFrame:
    matches = search_all(df, query)
    if matches.empty:
        return matches

    group = "determinant" if "determinant" in matches.columns else "element_name"
    association = (
        "antibiotic" if "antibiotic" in matches.columns
        and matches["antibiotic"].astype(str).str.strip().any()
        else "Associated antibiotic/group"
    )

    fields = ["sample", "organism", association, group]
    fields = [x for x in fields if x in matches.columns]

    return (
        matches.groupby(fields, dropna=False)
        .size()
        .reset_index(name="finding_count")
        .sort_values(fields, kind="stable")
        .reset_index(drop=True)
    )


# ---------------------------------------------------------------------------
# UI
# ---------------------------------------------------------------------------

st.set_page_config(
    page_title="AMRIVA — AMRFinderPlus Explorer",
    page_icon="🧬",
    layout="wide",
)

st.title("🧬 AMRIVA")
st.caption("AMRFinderPlus Results Organizer & Explorer")

if "results" not in st.session_state:
    st.session_state.results = pd.DataFrame()

st.sidebar.header("Data input")
tsv_files = st.sidebar.file_uploader(
    "Upload AMRFinderPlus TSV result files",
    type=["tsv", "txt"],
    accept_multiple_files=True,
    help="Upload one or more actual AMRFinderPlus tabular output files.",
)

fasta_files = st.sidebar.file_uploader(
    "Optional FASTA files for nucleotide-sequence extraction",
    type=["fasta", "fa", "fna"],
    accept_multiple_files=True,
    help="Used only to extract the reported Start–Stop region from the analyzed contig.",
)

if st.sidebar.button("Organize / Reload Results", type="primary"):
    if not tsv_files:
        st.sidebar.error("Upload at least one AMRFinderPlus TSV result.")
    else:
        fasta_map = {}
        for f in fasta_files or []:
            fasta_map[Path(f.name).stem] = parse_fasta(f)

        try:
            st.session_state.results = load_batch(tsv_files, fasta_map)
            st.sidebar.success(
                f"Loaded {len(st.session_state.results)} AMRFinderPlus findings."
            )
        except Exception as exc:
            st.sidebar.error(f"Could not organize the uploaded results: {exc}")

df = st.session_state.results

# ---------------- Home ----------------
tab_home, tab_genomic, tab_antibiotic, tab_compare, tab_profile, tab_search = st.tabs(
    [
        "Home",
        "Genomic Analysis",
        "Antibiotic Explorer",
        "Cross-Sample Comparison",
        "Sample / Organism Profile",
        "Search / Filter",
    ]
)

with tab_home:
    st.header("AMRIVA")
    st.write(
        "AMRIVA organizes detailed genomic output from AMRFinderPlus into a "
        "searchable multi-sample workspace."
    )
    st.info(
        "AMRIVA is an organizer/explorer of AMRFinderPlus findings. "
        "It does not independently diagnose phenotypic resistance or invent "
        "biological results."
    )
    st.markdown(
        """
**Workflow**

`AMRFinderPlus TSV → AMRIVA → organized genomic results → search → profiles → comparison`

Upload the AMRFinderPlus TSV results in the sidebar. If the actual analyzed FASTA/contigs
are also supplied, AMRIVA can extract the nucleotide segment corresponding to the reported
contig coordinates.
"""
    )
    if not df.empty:
        c1, c2, c3 = st.columns(3)
        c1.metric("Findings", len(df))
        c2.metric("Samples", df["sample"].nunique() if "sample" in df else 0)
        c3.metric("Organisms", df["organism"].replace("", pd.NA).nunique()
                  if "organism" in df else 0)

# ---------------- Genomic Analysis ----------------
with tab_genomic:
    st.header("Genomic Analysis")
    st.caption("Organized AMRFinderPlus output. The source fields remain traceable.")

    if df.empty:
        st.info("Upload AMRFinderPlus TSV results from the sidebar.")
    else:
        display_cols = [
            "sample", "organism", "antibiotic", "determinant", "class",
            "subclass", "mechanism", "method", "identity", "coverage",
            "contig", "Nucleotide Sequence", "start", "stop", "reference",
            "raw_result",
        ]
        display_cols = [c for c in display_cols if c in df.columns]
        st.dataframe(
            df[display_cols],
            use_container_width=True,
            hide_index=True,
        )

        st.download_button(
            "Download organized AMRIVA results (CSV)",
            df.to_csv(index=False).encode("utf-8"),
            file_name="AMRIVA_organized_results.csv",
            mime="text/csv",
        )

# ---------------- Antibiotic Explorer ----------------
with tab_antibiotic:
    st.header("Antibiotic Explorer")
    query = st.text_input(
        "Search antibiotic / antimicrobial association",
        placeholder="Example: Penicillin",
        key="antibiotic_query",
    )

    if df.empty:
        st.info("Upload results first.")
    elif query:
        matches = search_all(df, query)
        st.write(f"**{len(matches)} matching findings**")
        if matches.empty:
            st.warning("No matching finding was present in the uploaded dataset.")
        else:
            cols = [
                "sample", "organism", "antibiotic", "Associated antibiotic/group",
                "determinant", "class", "subclass", "mechanism", "contig",
                "Nucleotide Sequence", "start", "stop", "reference",
            ]
            cols = [c for c in cols if c in matches.columns]
            st.dataframe(matches[cols], use_container_width=True, hide_index=True)

# ---------------- Cross-Sample Comparison ----------------
with tab_compare:
    st.header("Cross-Sample Comparison")
    query = st.text_input(
        "Compare findings for an antibiotic/group or AMR determinant",
        placeholder="Example: Penicillin or blaTEM",
        key="comparison_query",
    )

    if df.empty:
        st.info("Upload results first.")
    elif query:
        comp = comparison_table(df, query)
        if comp.empty:
            st.warning("No matching findings were present in the uploaded dataset.")
        else:
            st.dataframe(comp, use_container_width=True, hide_index=True)

            matches = search_all(df, query)
            st.subheader("Underlying genomic findings")
            cols = [
                "sample", "organism", "determinant", "class", "subclass",
                "mechanism", "contig", "Nucleotide Sequence", "start",
                "stop", "reference",
            ]
            cols = [c for c in cols if c in matches.columns]
            st.dataframe(matches[cols], use_container_width=True, hide_index=True)

# ---------------- Profiles ----------------
with tab_profile:
    st.header("Sample / Organism Profile")

    if df.empty:
        st.info("Upload results first.")
    else:
        p1, p2 = st.columns(2)
        with p1:
            samples = sorted(
                [x for x in df["sample"].astype(str).unique() if x.strip()]
            ) if "sample" in df.columns else []
            selected_sample = st.selectbox(
                "Sample",
                ["— Select —"] + samples,
                key="profile_sample",
            )
        with p2:
            organisms = sorted(
                [x for x in df["organism"].astype(str).unique() if x.strip()]
            ) if "organism" in df.columns else []
            selected_organism = st.selectbox(
                "Organism",
                ["— Select —"] + organisms,
                key="profile_organism",
            )

        if selected_sample != "— Select —":
            profile = exact_profile(df, "sample", selected_sample)
            st.subheader(f"Sample: {selected_sample}")
            st.dataframe(profile, use_container_width=True, hide_index=True)

        if selected_organism != "— Select —":
            profile = exact_profile(df, "organism", selected_organism)
            st.subheader(f"Organism: {selected_organism}")
            st.dataframe(profile, use_container_width=True, hide_index=True)

# ---------------- Search / Filter ----------------
with tab_search:
    st.header("Search / Filter")

    if df.empty:
        st.info("Upload results first.")
    else:
        f1, f2, f3 = st.columns(3)

        with f1:
            sample_options = ["All"] + sorted(df["sample"].astype(str).unique())
            organism_options = ["All"] + sorted(df["organism"].astype(str).unique())
            sample_filter = st.selectbox("Sample ID", sample_options)
            organism_filter = st.selectbox("Organism", organism_options)

        with f2:
            determinant_options = ["All"] + sorted(
                df["determinant"].astype(str).unique()
            )
            class_options = ["All"] + sorted(df["class"].astype(str).unique())
            determinant_filter = st.selectbox("AMR determinant", determinant_options)
            class_filter = st.selectbox("Class", class_options)

        with f3:
            subclass_options = ["All"] + sorted(
                df["subclass"].astype(str).unique()
            )
            antibiotic_options = ["All"] + sorted(
                df["Associated antibiotic/group"].astype(str).unique()
            )
            subclass_filter = st.selectbox("Subclass", subclass_options)
            antibiotic_filter = st.selectbox(
                "Antibiotic / group", antibiotic_options
            )

        filtered = df.copy()

        filters = [
            ("sample", sample_filter),
            ("organism", organism_filter),
            ("determinant", determinant_filter),
            ("class", class_filter),
            ("subclass", subclass_filter),
            ("Associated antibiotic/group", antibiotic_filter),
        ]
        for column, value in filters:
            if value != "All" and column in filtered.columns:
                filtered = filtered[
                    filtered[column].astype(str) == value
                ]

        st.write(f"**{len(filtered)} findings match the selected filters.**")
        st.dataframe(filtered, use_container_width=True, hide_index=True)

        st.download_button(
            "Download filtered results (CSV)",
            filtered.to_csv(index=False).encode("utf-8"),
            file_name="AMRIVA_filtered_results.csv",
            mime="text/csv",
        )
