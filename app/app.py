import os
import zipfile
import tempfile
import shutil
import subprocess
from io import StringIO
import pandas as pd
import gradio as gr

APP_NAME = "AMRIVA"

DETAIL_COLUMNS = [
    "Sample ID", "Organism", "AMR determinant", "Sequence / Element name",
    "Class", "Subclass", "Mechanism / Product", "Method",
    "Identity (%)", "Coverage (%)", "Contig", "Start", "Stop",
    "Reference", "Raw result"
]

RAW_FIELDS = [
    "Protein id", "Contig id", "Start", "Stop", "Strand", "Element symbol",
    "Element name", "Scope", "Type", "Subtype", "Class", "Subclass", "Method",
    "Target length", "Reference sequence length", "% Coverage of reference",
    "% Identity to reference", "Alignment length", "Closest reference name",
    "HMM accession", "HMM description", "Reference accession"
]

ORGANISM_MAP = {
    "escherichia coli": "Escherichia coli", "e. coli": "Escherichia coli",
    "klebsiella pneumoniae": "Klebsiella pneumoniae",
    "pseudomonas aeruginosa": "Pseudomonas aeruginosa",
    "acinetobacter baumannii": "Acinetobacter baumannii",
    "staphylococcus aureus": "Staphylococcus aureus",
    "enterococcus faecalis": "Enterococcus faecalis",
    "enterococcus faecium": "Enterococcus faecium",
    "salmonella": "Salmonella spp.", "campylobacter": "Campylobacter spp.",
    "serratia marcescens": "Serratia marcescens",
}


def clean(v):
    if v is None:
        return ""
    s = str(v).strip()
    return "" if s.lower() in {"nan", "none", "<na>", "na"} else s


def read_fasta_records(path):
    records = []
    header = None
    parts = []
    with open(path, "r", encoding="utf-8") as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(parts).upper()))
                header = line[1:].strip()
                parts = []
            else:
                parts.append(line)
    if header is not None:
        records.append((header, "".join(parts).upper()))
    return records


def fasta_stats(seq):
    valid = "ACGTN"
    seq = "".join(x for x in seq.upper() if x in valid)
    if not seq:
        raise ValueError("No valid nucleotide sequence found.")
    n = len(seq)
    return {
        "length": n, "A": seq.count("A"), "C": seq.count("C"),
        "G": seq.count("G"), "T": seq.count("T"), "N": seq.count("N"),
        "GC": (seq.count("G") + seq.count("C")) / n * 100,
    }


def sample_id_from_header(header, filename=""):
    return (header.split()[0] if header else os.path.splitext(os.path.basename(filename))[0]).strip()


def infer_organism(header, filename=""):
    text = f"{header} {filename}".lower()
    for key in sorted(ORGANISM_MAP, key=len, reverse=True):
        if key in text:
            return ORGANISM_MAP[key]
    return "Not provided"


def load_metadata(path):
    if not path:
        return {}
    ext = os.path.splitext(str(path))[1].lower()
    if ext == ".csv":
        df = pd.read_csv(path, dtype=str).fillna("")
    elif ext in {".xlsx", ".xls"}:
        df = pd.read_excel(path, dtype=str).fillna("")
    else:
        raise ValueError("Metadata must be CSV or Excel.")
    cols = {str(c).strip().lower(): c for c in df.columns}
    sid_col = next((cols[k] for k in ["sample id", "sample_id", "sample", "id"] if k in cols), None)
    org_col = next((cols[k] for k in ["organism", "organism name", "species"] if k in cols), None)
    if not sid_col or not org_col:
        raise ValueError("Metadata requires columns: Sample ID and Organism.")
    return {clean(r[sid_col]): clean(r[org_col]) for _, r in df.iterrows() if clean(r[sid_col])}


def run_amrfinder(path, organism=None):
    cmd = ["amrfinder", "-n", path]
    if organism and organism != "Not provided":
        # AMRFinderPlus expects its supported organism name; do not guess it from arbitrary text.
        normalized = organism.strip().replace(" ", "_")
        cmd += ["-O", normalized]
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    except FileNotFoundError:
        return False, "AMRFinderPlus was not found on this server."
    except subprocess.TimeoutExpired:
        return False, "AMRFinderPlus analysis timed out after 300 seconds."
    except Exception as e:
        return False, f"AMRFinderPlus execution error: {e}"
    if p.returncode != 0:
        return False, p.stderr.strip() or "AMRFinderPlus returned a non-zero exit code."
    return True, p.stdout.strip()


def parse_amrfinder_tsv(raw, sample_id, organism):
    """Parse the real AMRFinderPlus tabular output, preserving the complete raw line."""
    if not raw or not raw.strip():
        return pd.DataFrame(columns=DETAIL_COLUMNS)
    lines = [x for x in raw.splitlines() if x.strip()]
    # AMRFinderPlus writes a header line beginning with '#'. Remove the comment marker only for parsing.
    header_idx = next((i for i, x in enumerate(lines) if x.startswith("#") and "Element symbol" in x), None)
    if header_idx is None:
        # Some versions emit the header without #.
        header_idx = next((i for i, x in enumerate(lines) if "Element symbol" in x and "Class" in x), None)
    if header_idx is None:
        return pd.DataFrame(columns=DETAIL_COLUMNS)

    header = lines[header_idx].lstrip("#").strip().split("\t")
    data_lines = lines[header_idx + 1:]
    data_lines = [x for x in data_lines if not x.startswith("#")]
    if not data_lines:
        return pd.DataFrame(columns=DETAIL_COLUMNS)

    rows = []
    for line in data_lines:
        vals = line.split("\t")
        if len(vals) < len(header):
            vals += [""] * (len(header) - len(vals))
        row = {header[i].strip(): clean(vals[i]) for i in range(len(header))}
        row["_raw"] = line
        rows.append(row)

    def first(row, *names):
        for name in names:
            if clean(row.get(name, "")):
                return clean(row.get(name))
        return ""

    out = []
    for r in rows:
        gene = first(r, "Element symbol", "Gene symbol")
        elem_name = first(r, "Element name", "Sequence name")
        cls = first(r, "Class")
        subclass = first(r, "Subclass")
        method = first(r, "Method")
        ident = first(r, "% Identity to reference", "% Identity")
        cov = first(r, "% Coverage of reference", "% Coverage")
        contig = first(r, "Contig id", "Sequence / Contig")
        start = first(r, "Start")
        stop = first(r, "Stop")
        ref_acc = first(r, "Closest reference accession", "Reference accession")
        ref_name = first(r, "Closest reference name", "Name of closest sequence")
        reference = " — ".join([x for x in [ref_acc, ref_name] if x])
        mechanism = first(r, "Element name", "Sequence name", "HMM description")
        out.append({
            "Sample ID": sample_id,
            "Organism": organism,
            "AMR determinant": gene,
            "Sequence / Element name": elem_name,
            "Class": cls,
            "Subclass": subclass,
            "Mechanism / Product": mechanism,
            "Method": method,
            "Identity (%)": ident,
            "Coverage (%)": cov,
            "Contig": contig,
            "Start": start,
            "Stop": stop,
            "Reference": reference,
            "Raw result": r.get("_raw", ""),
        })
    return pd.DataFrame(out, columns=DETAIL_COLUMNS)


def analyze_one(path, metadata=None):
    metadata = metadata or {}
    records = read_fasta_records(path)
    all_details = []
    summaries = []
    if not records:
        raise ValueError("The FASTA file contains no records.")
    for header, seq in records:
        sid = sample_id_from_header(header, path)
        organism = metadata.get(sid) or infer_organism(header, path)
        stats = fasta_stats(seq)
        # Write this record to a temporary FASTA so multi-record files are handled one sample at a time.
        with tempfile.NamedTemporaryFile("w", suffix=".fasta", delete=False, encoding="utf-8") as tmp:
            tmp.write(f">{header}\n{seq}\n")
            tmp_path = tmp.name
        try:
            success, raw = run_amrfinder(tmp_path)
        finally:
            try: os.unlink(tmp_path)
            except OSError: pass
        details = parse_amrfinder_tsv(raw if success else "", sid, organism)
        if not details.empty:
            all_details.append(details)
        determinants = ", ".join(dict.fromkeys(details["AMR determinant"].tolist())) if not details.empty else ""
        targets = ", ".join(dict.fromkeys(details["Subclass"].replace("", pd.NA).dropna().tolist())) if not details.empty else ""
        summaries.append({
            "Sample ID": sid, "Organism": organism,
            "Sequence length": stats["length"], "GC (%)": f'{stats["GC"]:.2f}',
            "AMR status": "AMR determinant(s) detected" if determinants else ("No known AMR determinant reported" if success else "Analysis failed"),
            "AMR determinants": determinants,
            "Resistance targets": targets,
            "Analysis message": "" if success else raw,
        })
    detail_df = pd.concat(all_details, ignore_index=True) if all_details else pd.DataFrame(columns=DETAIL_COLUMNS)
    return pd.DataFrame(summaries), detail_df


def batch_analyze(zip_path, metadata_path=None):
    if not zip_path:
        return pd.DataFrame(), pd.DataFrame(), "Upload a ZIP containing FASTA files."
    metadata = load_metadata(metadata_path) if metadata_path else {}
    work = tempfile.mkdtemp(prefix="amriva_")
    try:
        with zipfile.ZipFile(zip_path, "r") as z:
            members = [m for m in z.infolist() if not m.is_dir() and m.filename.lower().endswith((".fa", ".fasta", ".fna"))]
            if not members:
                return pd.DataFrame(), pd.DataFrame(), "No FASTA files were found in the ZIP."
            summaries, details = [], []
            for i, m in enumerate(members):
                safe_name = os.path.basename(m.filename) or f"sample_{i}.fasta"
                path = os.path.join(work, f"{i}_{safe_name}")
                with z.open(m) as src, open(path, "wb") as dst:
                    shutil.copyfileobj(src, dst)
                try:
                    s, d = analyze_one(path, metadata)
                    summaries.append(s)
                    if not d.empty:
                        details.append(d)
                except Exception as e:
                    summaries.append(pd.DataFrame([{
                        "Sample ID": os.path.splitext(safe_name)[0], "Organism": metadata.get(os.path.splitext(safe_name)[0], "Not provided"),
                        "Sequence length": "", "GC (%)": "", "AMR status": "Analysis failed",
                        "AMR determinants": "", "Resistance targets": "", "Analysis message": str(e)
                    }]))
            summary = pd.concat(summaries, ignore_index=True) if summaries else pd.DataFrame()
            detail = pd.concat(details, ignore_index=True) if details else pd.DataFrame(columns=DETAIL_COLUMNS)
            return summary, detail, f"Analysis completed for {len(summary)} sample record(s)."
    except zipfile.BadZipFile:
        return pd.DataFrame(), pd.DataFrame(), "The uploaded file is not a valid ZIP archive."
    finally:
        shutil.rmtree(work, ignore_errors=True)


def search_explorer(details, query):
    if not isinstance(details, pd.DataFrame) or details.empty:
        return pd.DataFrame(columns=DETAIL_COLUMNS), "No AMRFinderPlus findings are available yet. Upload and analyze samples first."
    q = clean(query).lower()
    if not q:
        return details.copy(), f"Showing all {len(details)} AMRFinderPlus finding(s)."
    searchable = ["Class", "Subclass", "AMR determinant", "Sequence / Element name", "Mechanism / Product", "Reference", "Organism", "Sample ID"]
    mask = pd.Series(False, index=details.index)
    for col in searchable:
        mask = mask | details[col].astype(str).str.lower().str.contains(q, regex=False, na=False)
    result = details.loc[mask].copy()
    if result.empty:
        msg = (f'No exact text match for "{query}". AMRFinderPlus may classify the target at a broader class/subclass level. '
               "Try a class/subclass such as BETA-LACTAM or CEPHALOSPORIN.")
    else:
        msg = f'Explorer found {len(result)} finding(s) matching "{query}".'
    return result, msg


def organism_profile(details, sample_id):
    if not isinstance(details, pd.DataFrame) or details.empty:
        return pd.DataFrame(columns=DETAIL_COLUMNS), "No findings available."
    sid = clean(sample_id)
    if not sid:
        return details.copy(), "Select a sample to see its resistance-associated genomic profile."
    out = details[details["Sample ID"].astype(str) == sid].copy()
    return out, f"{sid}: {len(out)} AMRFinderPlus finding(s)."


def explorer_choices(details):
    if not isinstance(details, pd.DataFrame) or details.empty:
        return gr.update(choices=[], value=None)
    vals = []
    for col in ["Subclass", "Class"]:
        for x in details[col].astype(str):
            x = clean(x)
            if x and x not in vals:
                vals.append(x)
    return gr.update(choices=vals, value=(vals[0] if vals else None))


def sample_choices(details):
    if not isinstance(details, pd.DataFrame) or details.empty:
        return gr.update(choices=[], value=None)
    vals = list(dict.fromkeys([clean(x) for x in details["Sample ID"].astype(str) if clean(x)]))
    return gr.update(choices=vals, value=(vals[0] if vals else None))


def run_single(path, metadata_path):
    if not path:
        return "Please upload a FASTA file.", pd.DataFrame(columns=DETAIL_COLUMNS), pd.DataFrame(columns=DETAIL_COLUMNS), gr.update(choices=[], value=None)
    try:
        metadata = load_metadata(metadata_path) if metadata_path else {}
        summary, details = analyze_one(path, metadata)
        message = f"✅ {len(summary)} FASTA record(s) analyzed. {len(details)} AMRFinderPlus finding(s) reported."
        return message, details, details, sample_choices(details)
    except Exception as e:
        return f"❌ Analysis failed: {e}", pd.DataFrame(columns=DETAIL_COLUMNS), pd.DataFrame(columns=DETAIL_COLUMNS), gr.update(choices=[], value=None)


def run_batch(zip_path, metadata_path):
    summary, details, msg = batch_analyze(zip_path, metadata_path)
    return msg, summary, details, details, sample_choices(details)

HOME = """# 🧬 AMRIVA\n## Genomic AMR Research Workspace\n\n### The problem we solve\nA researcher may have **100 isolates** and AMRFinderPlus findings spread across many results. AMRIVA organizes those findings around two practical questions:\n\n**Organism → What resistance-associated determinants were detected in this organism?**\n\n**Resistance target → Which samples/organisms have a matching AMRFinderPlus finding, which determinant was detected, and where is it located?**\n\n### Workflow\n**FASTA → AMRFinderPlus → structured genomic findings → Organism Profile / Resistance Explorer → sequence location & reference evidence → further laboratory validation**\n\nAMRIVA is a research/educational interface around AMRFinderPlus. It does not replace AMRFinderPlus or phenotypic AST."""

ABOUT = """## What makes the explorer useful?\n\nAMRIVA does not invent resistance calls. It uses the actual AMRFinderPlus output and reorganizes it into a searchable workspace.\n\nFor each finding you can see:\n- sample ID\n- organism (from uploaded metadata or FASTA header; never guessed from a gene)\n- AMR determinant\n- resistance class/subclass reported by AMRFinderPlus\n- method\n- identity and coverage\n- contig and genomic start/stop\n- closest reference\n- original raw AMRFinderPlus result\n\n**Important:** AMRFinderPlus `Subclass` can provide a more specific antibiotic or antibiotic-class association. A class/subclass annotation is not automatically proof of phenotypic resistance to every individual drug in that class. Phenotypic AST remains separate confirmation.\n"""

AMRIVA_CSS = """
.gradio-container{max-width:1450px!important;margin:auto!important}
#hero{padding:34px;border-radius:24px;margin-bottom:20px;background:linear-gradient(135deg,#123c69,#0f766e);color:white}
#hero h1{font-size:48px!important;margin:0!important}
"""

with gr.Blocks(title=APP_NAME) as demo:
    detail_state = gr.State(pd.DataFrame(columns=DETAIL_COLUMNS))
    gr.HTML("<div id='hero'><h1>🧬 AMRIVA</h1><p>Genomic AMR Research Workspace</p><p>AMRFinderPlus-powered analysis • Organism profiles • Resistance Explorer</p></div>")
    with gr.Tabs():
        with gr.Tab("🏠 Home"):
            gr.Markdown(HOME)
        with gr.Tab("🧬 Genomic Analysis"):
            gr.Markdown("## Genomic Analysis\nUpload FASTA. AMRFinderPlus performs the genomic screening; AMRIVA structures the returned fields into a readable evidence table.")
            with gr.Row():
                fasta = gr.File(label="FASTA / FNA", file_types=[".fa", ".fasta", ".fna"], type="filepath")
                metadata = gr.File(label="Optional organism metadata (CSV/XLSX)", file_types=[".csv", ".xlsx", ".xls"], type="filepath")
            gr.Markdown("**Metadata format:** `Sample ID, Organism`. If organism information is not supplied, AMRIVA shows **Not provided** rather than inventing a species.")
            run = gr.Button("🧬 Run AMRFinderPlus", variant="primary")
            status = gr.Markdown()
            results = gr.Dataframe(headers=DETAIL_COLUMNS, datatype=["str"] * len(DETAIL_COLUMNS), interactive=False, wrap=True, label="Structured AMRFinderPlus Findings")
            gr.Markdown("### What the table means")
            gr.Markdown("**AMR determinant** = detected genetic element; **Class/Subclass** = AMRFinderPlus resistance target classification; **Identity/Coverage** = sequence-reference evidence; **Contig + Start/Stop** = genomic location; **Reference** = closest reference evidence; **Raw result** = original AMRFinderPlus row.")
        with gr.Tab("📁 Batch Analysis"):
            gr.Markdown("## Analyze many isolates\nUpload one ZIP containing any number of `.fa`, `.fasta` or `.fna` files. Optional metadata connects Sample IDs to organism names.")
            with gr.Row():
                batch_zip = gr.File(label="ZIP of FASTA files", file_types=[".zip"], type="filepath")
                batch_meta = gr.File(label="Optional Sample ID → Organism CSV/XLSX", file_types=[".csv", ".xlsx", ".xls"], type="filepath")
            batch_run = gr.Button("📊 Analyze All Samples", variant="primary")
            batch_status = gr.Markdown()
            batch_summary = gr.Dataframe(interactive=False, wrap=True, label="Sample overview")
            batch_detail = gr.Dataframe(headers=DETAIL_COLUMNS, datatype=["str"] * len(DETAIL_COLUMNS), interactive=False, wrap=True, label="All AMRFinderPlus Findings")
        with gr.Tab("🔎 Resistance Explorer"):
            gr.Markdown("## Resistance Explorer\nSearch the complete uploaded dataset by **antibiotic/resistance target, class, subclass, gene, organism or sample**.")
            query = gr.Textbox(label="Search resistance target / antibiotic / class / subclass / gene", placeholder="e.g. BETA-LACTAM, CEPHALOSPORIN, blaTEM")
            search_btn = gr.Button("🔎 Explore", variant="primary")
            explorer_msg = gr.Markdown()
            explorer = gr.Dataframe(headers=DETAIL_COLUMNS, datatype=["str"] * len(DETAIL_COLUMNS), interactive=False, wrap=True, label="Matching samples, organisms, genes and genomic locations")
            gr.Markdown("**Note:** If you search an individual antibiotic such as `penicillin` but AMRFinderPlus reports only a broader class such as `BETA-LACTAM`, AMRIVA will not silently convert the broad class into a drug-specific resistance claim. Search the reported class/subclass or use the AMRFinderPlus-supported annotation.")
        with gr.Tab("🧫 Organism Profile"):
            gr.Markdown("## Organism / Sample Resistance Profile\nSee all AMRFinderPlus resistance-associated findings detected in one sample.")
            sample_drop = gr.Dropdown(label="Sample ID", choices=[])
            profile_btn = gr.Button("View Profile", variant="primary")
            profile_msg = gr.Markdown()
            profile = gr.Dataframe(headers=DETAIL_COLUMNS, datatype=["str"] * len(DETAIL_COLUMNS), interactive=False, wrap=True, label="Complete genomic profile")
        with gr.Tab("ℹ️ About"):
            gr.Markdown(ABOUT)

    run.click(run_single, [fasta, metadata], [status, results, detail_state, sample_drop])
    batch_run.click(run_batch, [batch_zip, batch_meta], [batch_status, batch_summary, batch_detail, detail_state, sample_drop])
    search_btn.click(search_explorer, [detail_state, query], [explorer_msg, explorer])
    query.submit(search_explorer, [detail_state, query], [explorer_msg, explorer])
    profile_btn.click(organism_profile, [detail_state, sample_drop], [profile, profile_msg])

if __name__ == "__main__":
    demo.launch(theme=gr.themes.Soft(), css=AMRIVA_CSS)
