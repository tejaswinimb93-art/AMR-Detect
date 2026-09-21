import os
import re
import csv
import io
import json
import zipfile
import shutil
import tempfile
import subprocess
from pathlib import Path

import pandas as pd
import gradio as gr

APP_TITLE = "AMRIVA — Genomic AMR Explorer"

# Render supplies PORT. Binding to 0.0.0.0 is essential.
PORT = int(os.environ.get("PORT", "7860"))

# Optional environment variables for AMRFinderPlus.
AMRFINDER_BIN = os.environ.get("AMRFINDER_BIN", "amrfinder")
AMRFINDER_DB = os.environ.get("AMRFINDER_DB", "")
AMRFINDER_UPDATE_ON_START = os.environ.get("AMRFINDER_UPDATE_ON_START", "0") == "1"

# AMRFinderPlus reports resistance-associated classes/subclasses rather than
# phenotypic AST results. These aliases are used only for dataset exploration.
SEARCH_ALIASES = {
    "penicillin": ["BETA-LACTAM"],
    "penicillins": ["BETA-LACTAM"],
    "cephalosporin": ["CEPHALOSPORIN", "BETA-LACTAM"],
    "cephalosporins": ["CEPHALOSPORIN", "BETA-LACTAM"],
    "carbapenem": ["CARBAPENEM", "BETA-LACTAM"],
    "carbapenems": ["CARBAPENEM", "BETA-LACTAM"],
    "fluoroquinolone": ["FLUOROQUINOLONE"],
    "fluoroquinolones": ["FLUOROQUINOLONE"],
    "ciprofloxacin": ["FLUOROQUINOLONE"],
    "tetracycline": ["TETRACYCLINE"],
    "tetracyclines": ["TETRACYCLINE"],
    "aminoglycoside": ["AMINOGLYCOSIDE"],
    "aminoglycosides": ["AMINOGLYCOSIDE"],
    "macrolide": ["MACROLIDE"],
    "macrolides": ["MACROLIDE"],
    "glycopeptide": ["GLYCOPEPTIDE"],
    "glycopeptides": ["GLYCOPEPTIDE"],
    "vancomycin": ["GLYCOPEPTIDE"],
    "sulfonamide": ["SULFONAMIDE"],
    "sulfonamides": ["SULFONAMIDE"],
    "trimethoprim": ["DIAMINOPYRIMIDINE"],
    "fosfomycin": ["FOSFOMYCIN"],
    "chloramphenicol": ["PHENICOL"],
    "rifampicin": ["RIFAMYCIN"],
    "rifampin": ["RIFAMYCIN"],
    "polymyxin": ["POLYMYXIN"],
    "colistin": ["POLYMYXIN"],
}


def norm(x):
    return re.sub(r"\s+", " ", str(x or "").strip())


def first_present(row, names):
    for n in names:
        if n in row and norm(row[n]):
            return norm(row[n])
    return ""


def infer_sample_id(path):
    return Path(path).stem


def infer_organism_from_header(header):
    h = header.strip().lstrip(">")
    # Common simple forms:
    # >S01|Escherichia coli
    # >S01 Escherichia coli
    if "|" in h:
        parts = [p.strip() for p in h.split("|") if p.strip()]
        if len(parts) >= 2:
            candidate = parts[1]
            if candidate and not re.fullmatch(r"(sample|isolate)?[_ -]?\d+", candidate, re.I):
                return candidate
    # Explicit organism labels
    m = re.search(r"(?:organism|species)\s*[:=]\s*([^|;]+)", h, re.I)
    if m:
        return m.group(1).strip()
    return ""


def fasta_records(path):
    records = []
    name = None
    seq = []
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(seq)))
                name = line[1:].strip()
                seq = []
            else:
                seq.append(re.sub(r"\s+", "", line))
    if name is not None:
        records.append((name, "".join(seq)))
    return records


def sequence_summary(path):
    records = fasta_records(path)
    seq = "".join(s for _, s in records).upper()
    counts = {b: seq.count(b) for b in "ACGTN"}
    total = len(seq)
    gc = ((counts["G"] + counts["C"]) / total * 100) if total else 0
    organism = ""
    for h, _ in records:
        organism = infer_organism_from_header(h)
        if organism:
            break
    return {
        "FASTA records": len(records),
        "Total bases": total,
        "A": counts["A"],
        "C": counts["C"],
        "G": counts["G"],
        "T": counts["T"],
        "N": counts["N"],
        "GC %": round(gc, 2),
        "Organism": organism or "Not provided",
    }


def canonicalize_amrfinder_row(raw):
    # AMRFinderPlus versions have used slightly different column names.
    # Normalize every known spelling into our stable display schema.
    lower = {str(k).strip().lower(): v for k, v in raw.items()}

    def get(*keys):
        for k in keys:
            if k.lower() in lower:
                return norm(lower[k.lower()])
        return ""

    gene = get("Gene symbol", "Element symbol", "AMR determinant", "Gene")
    sequence_name = get("Sequence name", "Element name")
    cls = get("Class", "AMR class")
    subclass = get("Subclass", "AMR subclass")
    method = get("Method")
    identity = get("% Identity to reference", "Identity (%)", "Identity")
    coverage = get("% Coverage of reference", "Coverage (%)", "Coverage")
    contig = get("Contig id", "Contig", "Sequence / Contig")
    start = get("Start")
    stop = get("Stop")
    strand = get("Strand")
    ref_acc = get("Closest reference accession", "Reference", "Reference accession")
    ref_name = get("Closest reference name", "Reference name")
    scope = get("Scope")
    elem_type = get("Element type", "Type")
    elem_subtype = get("Element subtype", "Subtype")
    target_len = get("Target length")
    ref_len = get("Reference sequence length")
    align_len = get("Alignment length")
    hmm_acc = get("HMM accession")
    hmm_desc = get("HMM description")

    # A readable AMR determinant: prefer gene symbol, then sequence name.
    determinant = gene or sequence_name or "Detected AMR-associated element"

    return {
        "AMR determinant": determinant,
        "Class": cls,
        "Subclass": subclass,
        "Mechanism / Product": sequence_name,
        "Method": method,
        "Identity (%)": identity,
        "Coverage (%)": coverage,
        "Sequence / Contig": contig,
        "Start": start,
        "Stop": stop,
        "Reference": ref_acc or ref_name,
        "Scope": scope,
        "Element type": elem_type,
        "Element subtype": elem_subtype,
        "Target length": target_len,
        "Reference length": ref_len,
        "Alignment length": align_len,
        "HMM accession": hmm_acc,
        "HMM description": hmm_desc,
    }


def parse_tsv(text, sample_id, organism=""):
    lines = [x for x in text.splitlines() if x.strip()]
    if not lines:
        return []

    # AMRFinderPlus is TSV. Keep the parser tolerant of BOM/whitespace.
    header_idx = None
    for i, line in enumerate(lines):
        if "\t" in line and (
            "Contig id" in line or "Gene symbol" in line or
            "Element symbol" in line or "Class" in line
        ):
            header_idx = i
            break

    if header_idx is None:
        return []

    reader = csv.DictReader(
        io.StringIO("\n".join(lines[header_idx:])),
        delimiter="\t"
    )
    out = []
    for raw in reader:
        if not any(norm(v) for v in raw.values()):
            continue
        item = canonicalize_amrfinder_row(raw)
        item["Sample ID"] = sample_id
        item["Organism"] = organism or "Not provided"
        item["Raw result"] = " | ".join(
            f"{k}={norm(v)}" for k, v in raw.items() if norm(v)
        )
        out.append(item)
    return out


def find_amrfinder():
    return shutil.which(AMRFINDER_BIN) or shutil.which("amrfinder")


def maybe_update_database():
    if not AMRFINDER_UPDATE_ON_START:
        return
    updater = shutil.which("amrfinder_update")
    if updater:
        cmd = [updater]
        if AMRFINDER_DB:
            cmd += ["-d", AMRFINDER_DB]
        subprocess.run(cmd, capture_output=True, text=True, timeout=900)


def run_amrfinder(fasta_path, workdir):
    binary = find_amrfinder()
    if not binary:
        raise RuntimeError(
            "AMRFinderPlus executable was not found. Install AMRFinderPlus in "
            "the Render environment and make sure 'amrfinder' is on PATH."
        )

    output = Path(workdir) / (Path(fasta_path).stem + ".amrfinder.tsv")
    cmd = [binary, "-n", str(fasta_path), "-o", str(output)]

    if AMRFINDER_DB:
        cmd += ["-d", AMRFINDER_DB]

    proc = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        timeout=900
    )

    if proc.returncode != 0:
        raise RuntimeError(
            "AMRFinderPlus failed.\n\nSTDERR:\n" +
            (proc.stderr[-5000:] if proc.stderr else "No stderr returned.")
        )

    if output.exists():
        return output.read_text(encoding="utf-8", errors="replace"), proc.stdout

    # Some installations may write TSV to stdout despite -o.
    return proc.stdout, proc.stdout


def prepare_inputs(files):
    paths = []
    tempdir = tempfile.mkdtemp(prefix="amriva_inputs_")
    for f in files or []:
        p = Path(f)
        if p.suffix.lower() == ".zip":
            with zipfile.ZipFile(p, "r") as z:
                for name in z.namelist():
                    if name.lower().endswith((".fa", ".fasta", ".fna")) and not name.endswith("/"):
                        safe = Path(name).name
                        target = Path(tempdir) / safe
                        with z.open(name) as src, open(target, "wb") as dst:
                            shutil.copyfileobj(src, dst)
                        paths.append(str(target))
        elif p.suffix.lower() in (".fa", ".fasta", ".fna"):
            target = Path(tempdir) / p.name
            shutil.copy2(p, target)
            paths.append(str(target))
    return tempdir, paths


def analyze_files(files, metadata_file=None):
    if not files:
        return (
            "⚠️ Upload at least one FASTA/FNA file or a ZIP containing FASTA files.",
            pd.DataFrame(), pd.DataFrame(), "No data"
        )

    try:
        maybe_update_database()
        tempdir, paths = prepare_inputs(files)
    except Exception as e:
        return f"❌ Input preparation error: {e}", pd.DataFrame(), pd.DataFrame(), "Error"

    metadata = {}
    if metadata_file:
        try:
            mp = Path(metadata_file)
            if mp.suffix.lower() == ".xlsx":
                mdf = pd.read_excel(mp)
            else:
                mdf = pd.read_csv(mp)
            cols = {str(c).strip().lower(): c for c in mdf.columns}
            sid_col = cols.get("sample id") or cols.get("sample_id") or cols.get("sample")
            org_col = cols.get("organism") or cols.get("species")
            if sid_col and org_col:
                for _, r in mdf.iterrows():
                    metadata[norm(r[sid_col])] = norm(r[org_col])
        except Exception as e:
            shutil.rmtree(tempdir, ignore_errors=True)
            return f"⚠️ Metadata file could not be read: {e}", pd.DataFrame(), pd.DataFrame(), "Error"

    summary_rows = []
    finding_rows = []
    errors = []

    for path in paths:
        sample_id = infer_sample_id(path)
        organism = metadata.get(sample_id, "")
        if not organism:
            recs = fasta_records(path)
            for h, _ in recs:
                organism = infer_organism_from_header(h)
                if organism:
                    break

        s = sequence_summary(path)
        s["Sample ID"] = sample_id
        s["Organism"] = organism or "Not provided"

        try:
            raw, _ = run_amrfinder(path, tempdir)
            findings = parse_tsv(raw, sample_id, organism)
            s["AMR Findings"] = len(findings)
            finding_rows.extend(findings)
        except Exception as e:
            s["AMR Findings"] = 0
            errors.append(f"{sample_id}: {e}")

        summary_rows.append(s)

    shutil.rmtree(tempdir, ignore_errors=True)

    summary_cols = [
        "Sample ID", "Organism", "FASTA records", "Total bases",
        "A", "C", "G", "T", "N", "GC %", "AMR Findings"
    ]
    finding_cols = [
        "Sample ID", "Organism", "AMR determinant", "Class", "Subclass",
        "Mechanism / Product", "Method", "Identity (%)", "Coverage (%)",
        "Sequence / Contig", "Start", "Stop", "Reference", "Raw result"
    ]

    summary_df = pd.DataFrame(summary_rows)
    for c in summary_cols:
        if c not in summary_df.columns:
            summary_df[c] = ""
    summary_df = summary_df[summary_cols]

    findings_df = pd.DataFrame(finding_rows)
    if findings_df.empty:
        findings_df = pd.DataFrame(columns=finding_cols)
    else:
        for c in finding_cols:
            if c not in findings_df.columns:
                findings_df[c] = ""
        findings_df = findings_df[finding_cols]

    status = f"✅ Analysis completed for {len(summary_df)} sample(s). {len(findings_df)} AMRFinderPlus finding(s) returned."
    if errors:
        status += "\n\n⚠️ Some samples could not be analyzed:\n" + "\n".join(errors)

    return status, summary_df, findings_df, "Ready"


def explorer_search(findings_df, query):
    if findings_df is None or len(findings_df) == 0:
        return "No AMRFinderPlus findings are available yet.", pd.DataFrame()

    q = norm(query).lower()
    if not q:
        return "Enter a resistance-associated target, gene, class, subclass, organism or sample ID.", pd.DataFrame()

    aliases = SEARCH_ALIASES.get(q, [])
    terms = [q] + [x.lower() for x in aliases]

    mask = pd.Series(False, index=findings_df.index)
    searchable_cols = [
        "AMR determinant", "Class", "Subclass", "Mechanism / Product",
        "Organism", "Sample ID", "Sequence / Contig", "Reference"
    ]

    for col in searchable_cols:
        if col in findings_df.columns:
            colmask = findings_df[col].fillna("").astype(str).str.lower().apply(
                lambda x: any(t in x for t in terms)
            )
            mask = mask | colmask

    result = findings_df.loc[mask].copy()
    cols = [
        "Sample ID", "Organism", "AMR determinant", "Class", "Subclass",
        "Sequence / Contig", "Start", "Stop", "Identity (%)",
        "Coverage (%)", "Reference"
    ]
    for c in cols:
        if c not in result.columns:
            result[c] = ""
    result = result[cols]

    if result.empty:
        return f"No AMRFinderPlus finding matched '{query}'.", result

    return (
        f"🔎 {len(result)} matching finding(s) for **{query}**. "
        "Results are genomic AMR-associated findings, not phenotypic AST results.",
        result
    )


def export_findings(findings_df):
    if findings_df is None or len(findings_df) == 0:
        return None
    out = Path(tempfile.gettempdir()) / "AMRIVA_AMRFinderPlus_organized_results.csv"
    findings_df.to_csv(out, index=False)
    return str(out)


with gr.Blocks(title=APP_TITLE, theme=gr.themes.Soft()) as demo:
    gr.Markdown(
        """
# 🧬 AMRIVA
### Genomic AMR Explorer for Multi-Isolate Research

**From raw AMRFinderPlus output → organized genomic evidence → searchable relationships across samples.**

AMRIVA does **not** replace AMRFinderPlus and does **not** perform phenotypic susceptibility testing.
It organizes and explores the genomic findings returned by AMRFinderPlus.
"""
    )

    with gr.Row():
        with gr.Column(scale=2):
            gr.Markdown(
                """
### Why AMRIVA?

When a researcher has dozens or hundreds of isolates, the important question is often not
only *“did AMRFinderPlus find an AMR-associated element?”* but:

**Which samples and organisms carry it? Which determinant was detected? What class/subclass is it associated with? Where is it located?**

AMRIVA brings those relationships into one searchable workspace.
"""
            )
        with gr.Column(scale=1):
            gr.Markdown(
                """
**Workflow**

`FASTA / ZIP`
→ `AMRFinderPlus`
→ `Structured findings`
→ `Resistance Explorer`

📌 Organism names are taken from supplied metadata or FASTA headers. AMRIVA never invents an organism name.
"""
            )

    with gr.Tab("🧬 Genomic Analysis"):
        gr.Markdown("### Upload isolates for real AMRFinderPlus analysis")
        fasta_files = gr.File(
            label="FASTA / FNA files or ZIP containing multiple FASTA files",
            file_count="multiple",
            file_types=[".fa", ".fasta", ".fna", ".zip"],
            type="filepath"
        )
        metadata_file = gr.File(
            label="Optional sample metadata CSV/XLSX — columns: Sample ID, Organism",
            file_count="single",
            file_types=[".csv", ".xlsx"],
            type="filepath"
        )
        run_btn = gr.Button("▶ Run AMRFinderPlus Analysis", variant="primary")
        status = gr.Markdown()
        summary = gr.Dataframe(label="Sample Genomic Summary", interactive=False, wrap=True)
        findings = gr.Dataframe(label="AMRFinderPlus Findings — Organized", interactive=False, wrap=True)
        download = gr.File(label="Download organized findings")

        run_btn.click(
            analyze_files,
            inputs=[fasta_files, metadata_file],
            outputs=[status, summary, findings, gr.State()]
        ).then(
            export_findings,
            inputs=[findings],
            outputs=[download]
        )

    with gr.Tab("🔎 Resistance Explorer"):
        gr.Markdown(
            """
### Search across all analyzed isolates

Enter a **gene, AMRFinderPlus class/subclass, organism, sample ID, genomic contig, or common antibiotic name**.

For antibiotic names, AMRIVA searches the corresponding AMRFinderPlus resistance class where an established mapping is available.
It does **not** claim AST susceptibility/resistance for an individual drug.
"""
        )
        explorer_query = gr.Textbox(
            label="Search target",
            placeholder="Examples: beta-lactam, penicillin, ciprofloxacin, blaTEM, Escherichia coli, S01"
        )
        explorer_btn = gr.Button("🔎 Explore", variant="primary")
        explorer_status = gr.Markdown()
        explorer_table = gr.Dataframe(label="Cross-Sample Genomic Relationships", interactive=False, wrap=True)

        explorer_btn.click(
            explorer_search,
            inputs=[findings, explorer_query],
            outputs=[explorer_status, explorer_table]
        )

    with gr.Tab("🧫 Organism Profile"):
        gr.Markdown("Select a sample from the analyzed dataset to inspect all of its genomic AMR-associated findings.")
        sample_dropdown = gr.Dropdown(label="Sample ID", choices=[])
        refresh_btn = gr.Button("Refresh sample list")
        organism_table = gr.Dataframe(label="Sample-specific AMRFinderPlus findings", interactive=False, wrap=True)

        def refresh_samples(df):
            if df is None or len(df) == 0 or "Sample ID" not in df.columns:
                return gr.update(choices=[])
            return gr.update(choices=sorted(df["Sample ID"].dropna().astype(str).unique().tolist()))

        def sample_profile(df, sid):
            if df is None or len(df) == 0 or not sid:
                return pd.DataFrame()
            return df[df["Sample ID"].astype(str) == str(sid)].copy()

        refresh_btn.click(refresh_samples, inputs=[findings], outputs=[sample_dropdown])
        sample_dropdown.change(sample_profile, inputs=[findings, sample_dropdown], outputs=[organism_table])

    with gr.Tab("ℹ️ Interpretation"):
        gr.Markdown(
            """
### How to interpret AMRIVA

- **AMR determinant** — the AMRFinderPlus-detected gene/element symbol when available.
- **Class / Subclass** — resistance-associated classification reported by AMRFinderPlus.
- **Method** — AMRFinderPlus detection method.
- **Identity / Coverage** — sequence-match metrics reported by AMRFinderPlus.
- **Sequence / Contig + Start / Stop** — genomic location reported by AMRFinderPlus.
- **Reference** — closest reference accession/name when supplied.
- **Organism** — supplied by the user through metadata or parsed from the FASTA header.

**Important:** genomic detection is not the same thing as an AST result. A detected resistance-associated determinant should be interpreted with the organism, gene/context, database version, and appropriate laboratory evidence.

### What makes the explorer useful?

For a large isolate collection, the researcher can move from:

**“I have 100 AMRFinderPlus outputs.”**

to:

**“Show me every isolate in this dataset carrying this resistance-associated class/gene, which organism it belongs to, which determinant was detected, and where it occurs.”**

That is the organizational layer AMRIVA adds around the underlying AMRFinderPlus analysis.
"""
        )

    gr.Markdown(
        """
---
**AMRIVA is a research/educational genomic analysis interface. It is not a clinical diagnostic or AST replacement.**
"""
    )


if __name__ == "__main__":
    demo.launch(
        server_name="0.0.0.0",
        server_port=PORT,
        share=False,
        show_error=True
    )
