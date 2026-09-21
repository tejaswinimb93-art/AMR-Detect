
import os, re, json, zipfile, shutil, subprocess, tempfile, uuid
from pathlib import Path
import pandas as pd
import gradio as gr

APP_NAME = "AMRIVA"
WORK = Path(tempfile.gettempdir()) / "amriva"
WORK.mkdir(parents=True, exist_ok=True)
EVIDENCE = {}

FASTA_EXTS = {".fa", ".fasta", ".fna", ".fas"}

def clean_sample_id(name):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", Path(name).stem)

def read_fasta_records(path):
    records = []
    header = None
    seq = []
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq)))
                header = line[1:].strip()
                seq = []
            else:
                seq.append(re.sub(r"\s+", "", line).upper())
        if header is not None:
            records.append((header, "".join(seq)))
    return records

def sequence_stats(path):
    records = read_fasta_records(path)
    seq = "".join(s for _, s in records).upper()
    counts = {b: seq.count(b) for b in "ACGTN"}
    valid = counts["A"] + counts["C"] + counts["G"] + counts["T"]
    gc = ((counts["G"] + counts["C"]) / valid * 100) if valid else 0
    return len(records), len(seq), counts, gc

def safe_extract_zip(zip_path, outdir):
    outdir.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(zip_path, "r") as z:
        for info in z.infolist():
            target = (outdir / info.filename).resolve()
            if not str(target).startswith(str(outdir.resolve()) + os.sep):
                raise ValueError("Unsafe ZIP path detected.")
            if info.is_dir():
                target.mkdir(parents=True, exist_ok=True)
            else:
                target.parent.mkdir(parents=True, exist_ok=True)
                with z.open(info) as src, open(target, "wb") as dst:
                    shutil.copyfileobj(src, dst)

def collect_fastas(files):
    root = WORK / str(uuid.uuid4())
    inputdir = root / "inputs"
    inputdir.mkdir(parents=True)
    paths = []
    for item in files or []:
        p = Path(item)
        if p.suffix.lower() == ".zip":
            zdir = inputdir / clean_sample_id(p.name)
            safe_extract_zip(p, zdir)
            paths.extend([x for x in zdir.rglob("*") if x.is_file() and x.suffix.lower() in FASTA_EXTS])
        elif p.suffix.lower() in FASTA_EXTS:
            dest = inputdir / p.name
            shutil.copy2(p, dest)
            paths.append(dest)
    # Deduplicate by resolved path
    return root, sorted(set(x.resolve() for x in paths))

def amrfinder_available():
    try:
        r = subprocess.run(["amrfinder", "--version"], capture_output=True, text=True, timeout=30)
        return r.returncode == 0, (r.stdout + r.stderr).strip()
    except Exception as e:
        return False, str(e)

def parse_amr_tsv(tsv_path):
    if not tsv_path.exists() or tsv_path.stat().st_size == 0:
        return pd.DataFrame()
    # AMRFinderPlus emits tab-delimited output; comment lines are ignored.
    try:
        df = pd.read_csv(tsv_path, sep="\t", dtype=str, comment="#", keep_default_na=False)
    except Exception:
        df = pd.read_csv(tsv_path, sep="\t", dtype=str, keep_default_na=False)
    df.columns = [str(c).strip() for c in df.columns]
    return df

def first_col(df, names):
    low = {str(c).lower(): c for c in df.columns}
    for n in names:
        if n.lower() in low:
            return low[n.lower()]
    return None

def normalize(df, sample_id):
    if df.empty:
        return pd.DataFrame(columns=[
            "Sample ID","Organism","AMR determinant","Class","Subclass",
            "Mechanism / Product","Method","Identity (%)","Coverage (%)",
            "Sequence / Contig","Start","Stop","Reference","Raw result"
        ])
    def val(names):
        c = first_col(df, names)
        return df[c] if c else pd.Series([""] * len(df), index=df.index)
    out = pd.DataFrame(index=df.index)
    out["Sample ID"] = sample_id
    out["Organism"] = val(["organism","Organism"])
    out["AMR determinant"] = val(["gene_symbol","Gene symbol","gene","Gene"])
    out["Class"] = val(["class","Class"])
    out["Subclass"] = val(["subclass","Subclass"])
    out["Mechanism / Product"] = val(["product","Product","mechanism","Mechanism"])
    out["Method"] = val(["method","Method","amr_method"])
    out["Identity (%)"] = val(["%identity","% Identity","identity"])
    out["Coverage (%)"] = val(["%coverage","% Coverage","coverage"])
    out["Sequence / Contig"] = val(["contig_name","Contig","sequence","Sequence name","name"])
    out["Start"] = val(["start","Start","start_on_contig"])
    out["Stop"] = val(["stop","Stop","end","stop_on_contig"])
    out["Reference"] = val(["closest_reference","Closest reference","accession","Accession","closest_reference_name"])
    out["Raw result"] = df.apply(lambda r: " | ".join(f"{c}={r[c]}" for c in df.columns if str(r[c]).strip()), axis=1)
    return out.reset_index(drop=True)

def find_sequence_evidence(nuc_fasta, sample_id):
    if not nuc_fasta or not nuc_fasta.exists():
        return {}
    evidence = {}
    current = None
    seq = []
    for h, s in read_fasta_records(nuc_fasta):
        evidence[h] = s
    return evidence

def run_one(path, organism):
    sample_id = clean_sample_id(path.name)
    outdir = WORK / str(uuid.uuid4())
    outdir.mkdir(parents=True)
    report = outdir / "amrfinder.tsv"
    matched_nuc = outdir / "matched_nucleotide.fasta"
    cmd = ["amrfinder", "-n", str(path), "-o", str(report),
           "--nucleotide_output", str(matched_nuc)]
    if organism and organism.strip():
        cmd += ["-O", organism.strip()]
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
    if proc.returncode != 0:
        return None, f"AMRFinderPlus failed for {sample_id}:\n{proc.stderr[-4000:]}"
    raw = parse_amr_tsv(report)
    norm = normalize(raw, sample_id)
    nrec, seqlen, counts, gc = sequence_stats(path)
    stats = {
        "Sample ID": sample_id, "FASTA records": nrec, "Total bases": seqlen,
        "A": counts["A"], "C": counts["C"], "G": counts["G"], "T": counts["T"],
        "N": counts["N"], "GC %": round(gc, 2),
        "AMR findings": len(norm)
    }
    return (norm, stats, matched_nuc, raw), ""

def analyze(files, organism):
    ok, ver = amrfinder_available()
    if not ok:
        return (
            "❌ AMRFinderPlus is not available in this deployment. Deploy the included Dockerfile so the server installs AMRFinderPlus and its database.",
            pd.DataFrame(), pd.DataFrame(), "", pd.DataFrame()
        )
    try:
        root, paths = collect_fastas(files)
    except Exception as e:
        return f"❌ Upload error: {e}", pd.DataFrame(), pd.DataFrame(), "", pd.DataFrame()
    if not paths:
        return "❌ No FASTA files found. Upload FASTA/FA/FNA files or a ZIP containing them.", pd.DataFrame(), pd.DataFrame(), "", pd.DataFrame()
    all_rows, stat_rows = [], []
    errors = []
    evidence_map = {}
    for p in paths:
        result, err = run_one(p, organism)
        if err:
            errors.append(err)
            continue
        norm, stats, matched_nuc, _ = result
        evidence_map[stats["Sample ID"]] = str(matched_nuc)
        all_rows.append(norm)
        stat_rows.append(stats)
    findings = pd.concat(all_rows, ignore_index=True) if all_rows else pd.DataFrame()
    stats = pd.DataFrame(stat_rows)
    if findings.empty:
        msg = f"✅ Analysis completed for {len(paths)} sample(s). No AMR-associated findings were reported by AMRFinderPlus."
    else:
        msg = f"✅ Analysis completed for {len(paths)} sample(s). {len(findings)} AMRFinderPlus finding(s) reported."
    if errors:
        msg += "\n\n⚠️ Some files failed:\n" + "\n\n".join(errors)
    global EVIDENCE
    EVIDENCE = evidence_map
    return msg, findings, stats, ver, findings

def explorer(query, findings):
    if findings is None or len(findings) == 0:
        return "Run genomic analysis first.", pd.DataFrame()
    q = (query or "").strip().lower()
    if not q:
        return "Enter an antibiotic, antibiotic class, gene, or AMRFinderPlus term.", pd.DataFrame()
    # This is deliberately a transparent text search over reported AMRFinderPlus annotations.
    # It avoids inventing antibiotic relationships not present in the tool output.
    cols = ["Class","Subclass","AMR determinant","Mechanism / Product","Reference","Raw result"]
    mask = pd.Series(False, index=findings.index)
    for c in cols:
        if c in findings.columns:
            mask |= findings[c].astype(str).str.lower().str.contains(re.escape(q), na=False)
    hits = findings.loc[mask].copy()
    if hits.empty:
        return f"No AMRFinderPlus result contains '{query}'. Try an antibiotic class, gene/determinant, or term shown in the results.", hits
    return f"Found {len(hits)} matching genomic finding(s) for '{query}'. These are AMRFinderPlus-associated findings, not phenotypic AST results.", hits

def build_sequence_evidence(sample, determinant, files):
    sample_id = clean_sample_id(sample or "")
    path = EVIDENCE.get(sample_id)
    if not path or not Path(path).exists():
        return "Run genomic analysis first, then enter the exact Sample ID from the results."
    records = read_fasta_records(Path(path))
    if not records:
        return "AMRFinderPlus did not return a matched nucleotide sequence for this finding."
    text = "\n".join(f">{h}\n{s}" for h, s in records)
    return (
        f"### {sample} — AMRFinderPlus matched nucleotide evidence\n\n"
        f"**Requested determinant:** {determinant or 'not specified'}\n\n"
        "The sequence below comes from AMRFinderPlus `--nucleotide_output`. "
        "It is not a fabricated sequence or a sequence reconstructed from the gene name.\n\n"
        "```text\n" + text[:20000] + "\n```"
    )

def load_demo():
    return None

css = """
body {background: #f7fafc;}
.gradio-container {max-width: 1250px !important;}
.hero {padding: 26px; border-radius: 20px; background: linear-gradient(135deg,#0f172a,#164e63); color:white; margin-bottom:18px;}
.hero h1 {font-size: 42px; margin:0 0 8px;}
.hero p {font-size:17px; opacity:.92;}
"""

with gr.Blocks(css=css, title="AMRIVA — Genomic AMR Explorer") as demo:
    findings_state = gr.State(pd.DataFrame())
    files_state = gr.State([])

    gr.HTML("""
    <div class="hero">
      <h1>🧬 AMRIVA</h1>
      <p>Genomic Antimicrobial Resistance Explorer</p>
      <p>Run NCBI AMRFinderPlus on assembled nucleotide sequences, then explore resistance-associated genomic findings across samples and organisms.</p>
    </div>
    """)

    with gr.Tab("Home"):
        gr.Markdown("""
### Why AMRIVA?

AMRIVA is a research and learning interface built around **NCBI AMRFinderPlus**. It does not replace the underlying detection engine. Instead, it provides a focused workflow for analyzing multiple genomic samples and exploring the resulting resistance-associated determinants.

**Workflow**

`FASTA / ZIP → AMRFinderPlus → genomic findings → antibiotic/class exploration → sequence evidence`

**Important:** a genomic finding is not the same as phenotypic susceptibility. AMRIVA does not provide clinical diagnosis or AST interpretation.
        """)

    with gr.Tab("Genomic Analysis"):
        gr.Markdown("### Upload genomic samples")
        uploads = gr.File(file_count="multiple", file_types=[".fasta",".fa",".fna",".fas",".zip"], label="FASTA files or ZIP of FASTA files")
        organism = gr.Textbox(label="Optional organism for AMRFinderPlus mutation screening", placeholder="e.g. Escherichia coli")
        analyze_btn = gr.Button("🔬 Run AMRFinderPlus Analysis", variant="primary")
        status = gr.Markdown()
        summary = gr.Dataframe(label="Sample statistics", interactive=False)
        results = gr.Dataframe(label="AMRFinderPlus findings", interactive=False, wrap=True)
        version = gr.Textbox(label="AMRFinderPlus runtime information", interactive=False)
        analyze_btn.click(analyze, [uploads, organism], [status, results, summary, version, findings_state])

        gr.Markdown("""
### What the results mean

The findings table is built from the **actual AMRFinderPlus TSV output**. Columns such as identity, coverage, method, coordinates and reference are shown only when AMRFinderPlus provides corresponding fields.

The uploaded nucleotide sequence is the genomic input; AMRFinderPlus is responsible for detecting AMR-associated genes and supported resistance-associated point mutations.
        """)

    with gr.Tab("Antibiotic Explorer"):
        gr.Markdown("""
### ⭐ Explore genomic findings by antibiotic/class

Enter a term such as an antibiotic class, gene/determinant, or another term that appears in the AMRFinderPlus annotation.

**Example research question:**  
*Which samples and organisms in my dataset have genomic findings associated with β-lactam resistance, and which determinants were detected?*

The explorer searches the reported AMRFinderPlus annotations. It does **not** infer phenotypic resistance from the genome.
        """)
        query = gr.Textbox(label="Search antibiotic / class / determinant", placeholder="e.g. beta-lactam, fluoroquinolone, blaTEM")
        explore = gr.Button("Explore")
        explore_status = gr.Markdown()
        explore_table = gr.Dataframe(interactive=False, wrap=True)
        explore.click(explorer, [query, findings_state], [explore_status, explore_table])

        gr.Markdown("### Finding details")
        sample = gr.Textbox(label="Sample ID")
        determinant = gr.Textbox(label="AMR determinant")
        evidence = gr.Markdown()
        evidence_btn = gr.Button("View uploaded sequence context")
        evidence_btn.click(build_sequence_evidence, [sample, determinant, uploads], evidence)


if __name__ == "__main__":
    demo.launch(server_name="0.0.0.0", server_port=int(os.environ.get("PORT","7860")))
