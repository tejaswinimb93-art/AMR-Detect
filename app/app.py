import os, re, zipfile, tempfile, shutil, subprocess, hashlib
from concurrent.futures import ThreadPoolExecutor, as_completed
import pandas as pd
import gradio as gr

APP_NAME = "AMRIVA"

# ---------- Constants ----------
AST_COLUMNS = ["Sample ID","AMR Finding","Antibiotic","MIC","MIC unit","AST","pH","Temperature","Other information","Comparison"]
GENOMIC_COLUMNS = ["Sample ID","Detected organism","FASTA file","Sequence length","GC content","AMR status","AMRFinderPlus result","Interpretation"]
COMPARISON_COLUMNS = ["Sample ID","Detected organism","AMR Finding","AMRFinderPlus result","Antibiotic","MIC","MIC unit","AST","pH","Temperature","Other information","Comparison"]
MIC_SERIES_COLUMNS = ["Sample ID","Antibiotic","Concentration","Unit","Growth"]

ORGANISM_MAP = {
    "escherichia coli":"Escherichia_coli", "e. coli":"Escherichia_coli",
    "klebsiella pneumoniae":"Klebsiella_pneumoniae", "pseudomonas aeruginosa":"Pseudomonas_aeruginosa",
    "acinetobacter baumannii":"Acinetobacter_baumannii", "staphylococcus aureus":"Staphylococcus_aureus",
    "staphylococcus epidermidis":"Staphylococcus_epidermidis", "enterococcus faecalis":"Enterococcus_faecalis",
    "enterococcus faecium":"Enterococcus_faecium", "salmonella":"Salmonella", "haemophilus influenzae":"Haemophilus_influenzae",
    "neisseria gonorrhoeae":"Neisseria_gonorrhoeae", "streptococcus agalactiae":"Streptococcus_agalactiae",
    "streptococcus pneumoniae":"Streptococcus_pneumoniae", "streptococcus pyogenes":"Streptococcus_pyogenes",
    "serratia marcescens":"Serratia_marcescens", "vibrio cholerae":"Vibrio_cholerae"
}

# Cache AMRFinder results in memory for the running server. This avoids rerunning the same file.
_AMR_CACHE = {}


def empty_df(columns):
    return pd.DataFrame(columns=columns)


def read_fasta(filepath):
    header, parts = "", []
    with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if not header:
                    header = line[1:].strip()
            else:
                parts.append(line.upper())
    return header, "".join(parts)


def detect_organism(header, filename):
    text = f"{header} {filename}".lower()
    for name in sorted(ORGANISM_MAP, key=len, reverse=True):
        if name in text:
            return name.title(), ORGANISM_MAP[name]
    return "Identification not available from FASTA metadata", None


def sequence_stats(sequence):
    seq = "".join(c for c in sequence.upper() if c in "ACGTN")
    if not seq:
        raise ValueError("No valid nucleotide sequence was found.")
    n = len(seq)
    return {"sequence": seq, "length": n, "A": seq.count("A"), "G": seq.count("G"),
            "C": seq.count("C"), "T": seq.count("T"), "N": seq.count("N"),
            "GC": (seq.count("G") + seq.count("C")) / n * 100}


def file_key(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_amrfinder_table(stdout):
    """Return a clean dataframe from AMRFinderPlus tabular output."""
    lines = [x.rstrip("\n") for x in (stdout or "").splitlines() if x.strip()]
    if not lines:
        return empty_df([]), "No known AMR determinant detected."
    data_lines = [x for x in lines if not x.startswith("#")]
    if not data_lines:
        return empty_df([]), "No known AMR determinant detected."
    try:
        from io import StringIO
        df = pd.read_csv(StringIO("\n".join(data_lines)), sep="\t", dtype=str)
        df = df.fillna("")
        if df.shape[1] <= 1:
            df = pd.read_csv(StringIO("\n".join(data_lines)), sep=r"\s+", engine="python", dtype=str).fillna("")
        return df, ""
    except Exception:
        # Last-resort tab split using first row as header when present.
        rows = [x.split("\t") for x in data_lines]
        if len(rows) >= 2:
            width = len(rows[0])
            rows = [r[:width] + [""] * max(0, width-len(r)) for r in rows[1:]]
            return pd.DataFrame(rows, columns=rows[0] if rows else []), ""
        return empty_df([]), "No known AMR determinant detected."


def run_amrfinder(path, organism=None):
    """Run AMRFinderPlus once per unique sequence file and return parsed data."""
    try:
        key = file_key(path)
        if key in _AMR_CACHE:
            return _AMR_CACHE[key]
        cmd = ["amrfinder", "-n", path]
        if organism:
            cmd += ["-O", organism]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=180)
        if result.returncode != 0:
            out = (False, empty_df([]), "", "AMRFinderPlus failed: " + (result.stderr.strip() or "unknown error"))
        else:
            table, fallback = parse_amrfinder_table(result.stdout)
            if table.empty:
                out = (True, table, "No known AMR determinant detected.", "")
            else:
                out = (True, table, "AMR determinant(s) detected", "")
        _AMR_CACHE[key] = out
        return out
    except subprocess.TimeoutExpired:
        return False, empty_df([]), "", "AMRFinderPlus analysis timed out."
    except FileNotFoundError:
        return False, empty_df([]), "", "AMRFinderPlus was not found on the server."
    except Exception as e:
        return False, empty_df([]), "", f"AMRFinderPlus failed: {e}"


def amr_table_text(df):
    if df is None or df.empty:
        return "No known AMR determinant detected."
    # Keep the actual AMRFinderPlus columns; do not invent values.
    return df.to_json(orient="records")


def genomic_interpretation(success, table):
    if not success:
        return "AMR genomic analysis could not be completed. Phenotypic AST remains necessary for susceptibility confirmation."
    if table is None or table.empty:
        return "No known AMR determinant was detected by this screening. This does not prove susceptibility; phenotypic AST is required."
    return "Detected genomic determinants may be associated with antimicrobial resistance. Genomic detection alone does not establish phenotypic resistance; phenotypic AST is required for confirmation."


def analyse_single(fasta_file):
    if not fasta_file:
        return ["Please upload a FASTA file."] + [""] * 14 + [empty_df([])]
    try:
        header, seq = read_fasta(fasta_file)
        stats = sequence_stats(seq)
        filename = os.path.basename(fasta_file)
        organism, org = detect_organism(header, filename)
        success, table, status, error = run_amrfinder(fasta_file, org)
        if error:
            status = "Analysis failed"
            result_text = error
        else:
            result_text = amr_table_text(table)
        interpretation = genomic_interpretation(success, table)
        sample_id = header.split()[0] if header else os.path.splitext(filename)[0]
        summary = "Sample analysed successfully." if success else "Sample analysed; AMRFinderPlus reported an error."
        return [summary, sample_id, organism, filename, str(stats["length"]), str(stats["A"]), str(stats["G"]),
                str(stats["C"]), str(stats["T"]), str(stats["N"]), f'{stats["GC"]:.2f}%', status,
                result_text, interpretation, header, table]
    except Exception as e:
        return [f"Analysis error: {e}"] + [""] * 14 + [empty_df([])]


def extract_fasta_files(zip_path, folder):
    files = []
    with zipfile.ZipFile(zip_path, "r") as z:
        for member in z.infolist():
            if member.is_dir():
                continue
            name = os.path.basename(member.filename)
            if not name or name.startswith(".") or not name.lower().endswith((".fa", ".fasta", ".fna")):
                continue
            dest = os.path.join(folder, name)
            base, ext = os.path.splitext(name)
            i = 1
            while os.path.exists(dest):
                dest = os.path.join(folder, f"{base}_{i}{ext}")
                i += 1
            with z.open(member) as src, open(dest, "wb") as dst:
                shutil.copyfileobj(src, dst)
            files.append(dest)
    return files


def analyse_one_batch(path, index):
    name = os.path.basename(path)
    try:
        header, seq = read_fasta(path)
        stats = sequence_stats(seq)
        organism, org = detect_organism(header, name)
        success, table, status, error = run_amrfinder(path, org)
        if error:
            status = "Analysis failed"
            result = error
        else:
            result = amr_table_text(table)
        sid = header.split()[0] if header else os.path.splitext(name)[0]
        return {"Sample ID": sid, "Detected organism": organism, "FASTA file": name,
                "Sequence length": stats["length"], "GC content": f'{stats["GC"]:.2f}%',
                "AMR status": status, "AMRFinderPlus result": result,
                "Interpretation": genomic_interpretation(success, table)}
    except Exception as e:
        return {"Sample ID": f"AMR-BATCH-{index:03d}", "Detected organism": "Unknown", "FASTA file": name,
                "Sequence length": "Error", "GC content": "Error", "AMR status": "Analysis failed",
                "AMRFinderPlus result": str(e), "Interpretation": "Analysis could not be completed."}


def analyse_zip(zip_file):
    if not zip_file:
        return empty_df(GENOMIC_COLUMNS), empty_df(GENOMIC_COLUMNS)
    folder = tempfile.mkdtemp(prefix="amriva_batch_")
    try:
        if not zipfile.is_zipfile(zip_file):
            df = pd.DataFrame([{"Sample ID":"ERROR","Detected organism":"N/A","FASTA file":os.path.basename(zip_file),
                "Sequence length":"N/A","GC content":"N/A","AMR status":"Invalid ZIP",
                "AMRFinderPlus result":"The uploaded file is not a valid ZIP archive.","Interpretation":"Upload a ZIP containing FASTA files."}], columns=GENOMIC_COLUMNS)
            return df, df
        files = extract_fasta_files(zip_file, folder)
        if not files:
            df = pd.DataFrame([{"Sample ID":"ERROR","Detected organism":"N/A","FASTA file":"No FASTA files found",
                "Sequence length":"N/A","GC content":"N/A","AMR status":"No FASTA files",
                "AMRFinderPlus result":"No .fa, .fasta or .fna files were found inside the ZIP.","Interpretation":"Upload a ZIP containing FASTA files."}], columns=GENOMIC_COLUMNS)
            return df, df
        # Parallelize independent files. Results remain scientifically identical; this only reduces wall-clock time.
        rows = [None] * len(files)
        workers = min(4, len(files))
        with ThreadPoolExecutor(max_workers=workers) as ex:
            futures = {ex.submit(analyse_one_batch, p, i+1): i for i, p in enumerate(files)}
            for fut in as_completed(futures):
                rows[futures[fut]] = fut.result()
        df = pd.DataFrame(rows, columns=GENOMIC_COLUMNS)
        return df, df
    except Exception as e:
        df = pd.DataFrame([{"Sample ID":"ERROR","Detected organism":"N/A","FASTA file":"ZIP processing",
            "Sequence length":"N/A","GC content":"N/A","AMR status":"Batch analysis failed",
            "AMRFinderPlus result":str(e),"Interpretation":"The ZIP could not be processed."}], columns=GENOMIC_COLUMNS)
        return df, df
    finally:
        shutil.rmtree(folder, ignore_errors=True)


def normalize_columns(df):
    aliases = {
        "sample_id":"Sample ID", "sample":"Sample ID", "Sample_ID":"Sample ID",
        "amr_finding":"AMR Finding", "AMR finding":"AMR Finding",
        "antibiotic":"Antibiotic", "mic":"MIC", "mic value":"MIC", "mic_unit":"MIC unit", "unit":"MIC unit",
        "ast":"AST", "ast result":"AST", "ph":"pH", "temperature":"Temperature", "temp":"Temperature",
        "other information":"Other information", "notes":"Other information", "experimental notes":"Other information"
    }
    return df.rename(columns={c: aliases.get(str(c).strip(), c) for c in df.columns})


def load_ast_file(file):
    if not file:
        return empty_df(AST_COLUMNS), "No file selected."
    try:
        ext = os.path.splitext(file)[1].lower()
        df = pd.read_csv(file) if ext == ".csv" else pd.read_excel(file)
        df = normalize_columns(df)
        for c in AST_COLUMNS:
            if c not in df.columns:
                df[c] = ""
        df = df[AST_COLUMNS].fillna("")
        return df, f"Loaded {len(df)} AST/MIC row(s) from {ext.upper()[1:]} ."
    except Exception as e:
        return empty_df(AST_COLUMNS), f"AST/MIC file could not be loaded: {e}"


def add_ast(data, sample, antibiotic, mic, unit, ast, ph, temp, notes):
    df = data.copy() if isinstance(data, pd.DataFrame) else empty_df(AST_COLUMNS)
    if not str(sample).strip() or not str(antibiotic).strip():
        return df, "Sample ID and antibiotic are required."
    row = {c:"" for c in AST_COLUMNS}
    row.update({"Sample ID":str(sample).strip(), "Antibiotic":str(antibiotic).strip(), "MIC":str(mic).strip(),
                "MIC unit":str(unit), "AST":str(ast), "pH":str(ph).strip(), "Temperature":str(temp).strip(),
                "Other information":str(notes).strip(), "Comparison":""})
    return pd.concat([df, pd.DataFrame([row])], ignore_index=True), f"Added AST/MIC data for {sample}."


def calculate_mic_series(data):
    out = ["Sample ID","Antibiotic","MIC","Unit","Lowest no-growth concentration","Interpretation note"]
    if data is None: return empty_df(out), "No observations entered."
    try:
        df = data.copy() if isinstance(data,pd.DataFrame) else pd.DataFrame(data)
        for c in MIC_SERIES_COLUMNS:
            if c not in df.columns: df[c] = ""
        df = df[MIC_SERIES_COLUMNS].copy()
        df["conc"] = pd.to_numeric(df["Concentration"], errors="coerce")
        df["growth_norm"] = df["Growth"].astype(str).str.strip().str.lower()
        rows=[]
        for (sample, ab), g in df.groupby(["Sample ID","Antibiotic"], dropna=False, sort=False):
            g=g.dropna(subset=["conc"])
            no=g[g["growth_norm"].isin(["no growth","nogrowth","no","absent","none"])]
            unit=next((str(x).strip() for x in g["Unit"] if str(x).strip() and str(x).lower()!="nan"),"")
            if no.empty:
                mic="Not determined"; lowest="No no-growth observation"; note="No tested concentration was recorded as No growth."
            else:
                x=float(no["conc"].min()); mic=f"{x:g} {unit}".strip(); lowest=mic; note="MIC reported as the lowest tested concentration recorded as No growth."
            rows.append({"Sample ID":str(sample),"Antibiotic":str(ab),"MIC":mic,"Unit":unit,"Lowest no-growth concentration":lowest,"Interpretation note":note})
        return pd.DataFrame(rows,columns=out), f"Calculated {len(rows)} sample/antibiotic combination(s)."
    except Exception as e:
        return empty_df(out), f"MIC series error: {e}"


def extract_gene_tokens(text):
    s=str(text).lower()
    return re.findall(r"\b(?:bla[a-z0-9_-]+|gyr[a-z0-9_-]*|par[a-z0-9_-]*|qnr[a-z0-9_-]*|[a-z]*mcr[-_]?\d+[a-z0-9_-]*|van[a-z0-9_-]*|erm[a-z0-9_-]*|tet[a-z0-9_-]*|aac[a-z0-9_-]*|aph[a-z0-9_-]*|aad[a-z0-9_-]*|sul[a-z0-9_-]*|dfr[a-z0-9_-]*|cat[a-z0-9_-]*)\b", s)


def resistance_relevance(amr_text, antibiotic):
    s=str(amr_text).lower(); ab=str(antibiotic).lower()
    genes=extract_gene_tokens(s)
    if not genes: return False
    classes = {
        "ciprofloxacin":["gyr","par","qnr"], "levofloxacin":["gyr","par","qnr"], "moxifloxacin":["gyr","par","qnr"],
        "ampicillin":["bla"], "amoxicillin":["bla"], "piperacillin":["bla"], "ceftriaxone":["bla"], "cefotaxime":["bla"], "ceftazidime":["bla"],
        "meropenem":["bla"], "imipenem":["bla"], "ertapenem":["bla"], "doripenem":["bla"],
        "vancomycin":["van"], "tetracycline":["tet"], "doxycycline":["tet"], "erythromycin":["erm"],
        "clindamycin":["erm"], "gentamicin":["aac","aph","aad"], "amikacin":["aac","aph","aad"], "tobramycin":["aac","aph","aad"],
        "trimethoprim":["dfr"], "sulfamethoxazole":["sul"], "chloramphenicol":["cat"], "colistin":["mcr"]
    }
    keys=[]
    for k,v in classes.items():
        if k in ab: keys.extend(v)
    if not keys:
        # If the antibiotic itself appears in the AMRFinder text, allow a cautious match.
        return ab and ab in s
    return any(any(g.startswith(k) for k in keys) for g in genes)


def run_comparative(genomic_df, ast_df, metadata_df):
    if genomic_df is None: genomic_df=empty_df(GENOMIC_COLUMNS)
    if ast_df is None: ast_df=empty_df(AST_COLUMNS)
    if metadata_df is None: metadata_df=pd.DataFrame()
    g=genomic_df.copy() if isinstance(genomic_df,pd.DataFrame) else pd.DataFrame(genomic_df)
    a=ast_df.copy() if isinstance(ast_df,pd.DataFrame) else pd.DataFrame(ast_df)
    m=metadata_df.copy() if isinstance(metadata_df,pd.DataFrame) else pd.DataFrame(metadata_df)
    if g.empty or a.empty:
        return empty_df(COMPARISON_COLUMNS), "Upload genomic FASTA/ZIP results and AST/MIC data first."
    for c in GENOMIC_COLUMNS:
        if c not in g.columns: g[c]=""
    for c in AST_COLUMNS:
        if c not in a.columns: a[c]=""
    if not m.empty:
        m=normalize_columns(m)
        for c in ["Sample ID","pH","Temperature","Other information"]:
            if c not in m.columns: m[c]=""
    rows=[]
    for _, ar in a.iterrows():
        sid=str(ar.get("Sample ID","")).strip()
        matches=g[g["Sample ID"].astype(str).str.strip().str.lower()==sid.lower()]
        if matches.empty:
            grows=[None]
        else:
            grows=[matches.iloc[0]]
        for gr in grows:
            ph=ar.get("pH",""); temp=ar.get("Temperature",""); other=ar.get("Other information","")
            if not m.empty:
                mm=m[m["Sample ID"].astype(str).str.strip().str.lower()==sid.lower()]
                if not mm.empty:
                    mr=mm.iloc[0]
                    ph=ph or mr.get("pH",""); temp=temp or mr.get("Temperature",""); other=other or mr.get("Other information","")
            organism=gr.get("Detected organism","") if gr is not None else "Genomic data not available"
            amrf=gr.get("AMRFinderPlus result","") if gr is not None else "Genomic data not available"
            amrfinding=gr.get("AMR status","") if gr is not None else "Genomic data not available"
            ast=str(ar.get("AST","")).strip().lower()
            relevant=resistance_relevance(amrf, ar.get("Antibiotic",""))
            if gr is None:
                comparison="Genomic data not available"
            elif not ast or ast=="not provided":
                comparison="Insufficient data"
            elif relevant and ast=="resistant":
                comparison="Potentially concordant"
            elif relevant and ast=="susceptible":
                comparison="Potentially discordant"
            elif not relevant and ast in {"resistant","susceptible"}:
                comparison="Requires interpretation"
            else:
                comparison="Requires interpretation"
            rows.append({"Sample ID":sid,"Detected organism":organism,"AMR Finding":amrfinding,"AMRFinderPlus result":amrf,
                "Antibiotic":ar.get("Antibiotic",""),"MIC":ar.get("MIC",""),"MIC unit":ar.get("MIC unit",""),"AST":ar.get("AST",""),
                "pH":ph,"Temperature":temp,"Other information":other,"Comparison":comparison})
    out=pd.DataFrame(rows,columns=COMPARISON_COLUMNS)
    msg=f"Built {len(out)} integrated comparison row(s). Genomic and AST/MIC data were matched by Sample ID."
    return out,msg


def summary(df):
    if df is None or not isinstance(df,pd.DataFrame) or df.empty: return "No comparative results yet."
    counts=df["Comparison"].value_counts().to_dict()
    return "### Comparison summary\n" + "\n".join([f"- **{k}:** {v}" for k,v in counts.items()])


def make_mic_plot(data):
    try:
        import matplotlib.pyplot as plt
        df=data.copy(); df["MIC_num"]=pd.to_numeric(df["MIC"].astype(str).str.extract(r"([0-9.]+)")[0],errors="coerce"); df=df.dropna(subset=["MIC_num"])
        if df.empty:return None
        fig,ax=plt.subplots(figsize=(9,5)); ax.bar(df["Sample ID"].astype(str),df["MIC_num"]); ax.set_xlabel("Sample ID"); ax.set_ylabel("MIC"); ax.set_title("MIC comparison"); ax.tick_params(axis="x",rotation=45); fig.tight_layout(); return fig
    except Exception:return None


def make_mic_ph_plot(data):
    return scatter_plot(data,"pH","MIC","MIC vs pH")

def make_mic_temp_plot(data):
    return scatter_plot(data,"Temperature","MIC","MIC vs temperature")

def scatter_plot(data,xcol,ycol,title):
    try:
        import matplotlib.pyplot as plt
        df=data.copy(); df[xcol]=pd.to_numeric(df[xcol],errors="coerce"); df[ycol]=pd.to_numeric(df[ycol].astype(str).str.extract(r"([0-9.]+)")[0],errors="coerce"); df=df.dropna(subset=[xcol,ycol])
        if df.empty:return None
        fig,ax=plt.subplots(figsize=(8,5)); ax.scatter(df[xcol],df[ycol]); ax.set_xlabel(xcol); ax.set_ylabel("MIC"); ax.set_title(title); fig.tight_layout(); return fig
    except Exception:return None


def make_genotype_phenotype_plot(data):
    try:
        import matplotlib.pyplot as plt
        if data is None or data.empty:return None
        df=data.copy(); counts=df["Comparison"].value_counts(); fig,ax=plt.subplots(figsize=(8,5)); ax.bar(counts.index.astype(str),counts.values); ax.set_ylabel("Samples"); ax.set_title("Genotype–phenotype comparison"); ax.tick_params(axis="x",rotation=25); fig.tight_layout(); return fig
    except Exception:return None


def dashboard(df,sample_id):
    if df is None or not isinstance(df,pd.DataFrame) or df.empty:return "No integrated data available. Run Comparative Analysis first."
    hit=df[df["Sample ID"].astype(str).str.strip().str.lower()==str(sample_id).strip().lower()]
    if hit.empty:return f"No integrated result found for `{sample_id}`."
    r=hit.iloc[0]
    return "## Sample Dashboard — " + str(r["Sample ID"]) + "\n\n" + "\n".join([f"- **{c}:** {r[c] if str(r[c]).strip() else 'Not provided'}" for c in COMPARISON_COLUMNS[1:]])


HOME="""# 🧬 AMRIVA\n\n**Antimicrobial Resistance Genomic Analysis Platform**\n\nAMRIVA connects genomic AMR screening with laboratory-provided AST/MIC and experimental information for research and education.\n\n**Research prototype — not a clinical diagnostic system.**"""
WOMEN="""# 👩‍🔬 Women's Health & AMR\n\nAMR is relevant to infections affecting women, including urinary tract infections and other areas of reproductive-health research. AMRIVA presents this as a research/education theme and does not diagnose patient samples or prescribe treatment."""
METHOD="""# 📚 Methodology\n\n**Sample → DNA extraction → sequencing → FASTA → AMRIVA → AMRFinderPlus → genomic findings**\n\nSeparately:\n\n**Laboratory AST/MIC + pH + temperature → AMRIVA**\n\nThen AMRIVA matches the datasets using **Sample ID** and performs comparative analysis and visualization. AMRIVA does not perform the wet-lab procedures itself."""
LIMITS="""# ⚠️ Limitations & Privacy\n\n- Research/educational prototype, not a clinical diagnostic system.\n- Does not prescribe antibiotics or replace AST.\n- Genomic detection does not automatically prove phenotypic resistance.\n- Absence of a detected determinant does not prove susceptibility.\n- AST/MIC/pH/temperature are laboratory-provided values.\n- Use authorized, de-identified research data; do not upload unnecessary patient-identifying information."""

CSS=""".gradio-container{max-width:1450px!important;margin:auto!important}#hero{padding:32px;border-radius:24px;margin-bottom:18px;background:linear-gradient(135deg,#123c69,#0f766e);color:white}#hero h1{font-size:48px!important;margin:0}#hero p{font-size:18px!important}"""

with gr.Blocks(title=APP_NAME,css=CSS,theme=gr.themes.Soft()) as demo:
    genomic_state=gr.State(empty_df(GENOMIC_COLUMNS))
    ast_state=gr.State(empty_df(AST_COLUMNS))
    comparison_state=gr.State(empty_df(COMPARISON_COLUMNS))
    gr.HTML("<div id='hero'><h1>🧬 AMRIVA</h1><p>Antimicrobial Resistance Genomic Analysis Platform</p><p>Genomic screening • AST/MIC integration • Comparative analysis</p></div>")
    with gr.Tabs():
        with gr.Tab("🏠 Home"): gr.Markdown(HOME)
        with gr.Tab("🧬 Single Sample"):
            sf=gr.File(label="Upload FASTA",file_types=[".fa",".fasta",".fna"],type="filepath"); sb=gr.Button("🧬 Analyse Sample",variant="primary"); ss=gr.Markdown()
            with gr.Row(): sid=gr.Textbox(label="Sample ID"); sorg=gr.Textbox(label="Detected organism")
            with gr.Row(): sfn=gr.Textbox(label="FASTA file"); sl=gr.Textbox(label="Sequence length"); sgc=gr.Textbox(label="GC content"); sst=gr.Textbox(label="AMR status")
            with gr.Row(): sa=gr.Textbox(label="A"); sg=gr.Textbox(label="G"); sc=gr.Textbox(label="C"); st=gr.Textbox(label="T"); sn=gr.Textbox(label="N")
            samr=gr.Dataframe(label="AMRFinderPlus result",interactive=False,wrap=True); sint=gr.Textbox(label="Interpretation",lines=5); sh=gr.Textbox(label="FASTA header",visible=False)
            sb.click(analyse_single,sf,[ss,sid,sorg,sfn,sl,sa,sg,sc,st,sn,sgc,sst,samr,sint,sh]).then(lambda sid,sorg,fn,sl,gc,status,tab: pd.DataFrame([{"Sample ID":sid,"Detected organism":sorg,"FASTA file":fn,"Sequence length":sl,"GC content":gc,"AMR status":status,"AMRFinderPlus result":amr_table_text(tab),"Interpretation":genomic_interpretation(True,tab)}]) if sid else empty_df(GENOMIC_COLUMNS),[sid,sorg,sfn,sl,sgc,sst,samr],genomic_state)
        with gr.Tab("📁 Multiple Samples"):
            bz=gr.File(label="Upload ZIP containing FASTA files",file_types=[".zip"],type="filepath"); bb=gr.Button("📊 Analyse All Samples",variant="primary"); bt=gr.Dataframe(headers=GENOMIC_COLUMNS,interactive=False,wrap=True,label="Batch genomic results"); bb.click(analyse_zip,bz,[bt,genomic_state])
        with gr.Tab("🧪 AST + MIC"):
            with gr.Row(): ast_sample=gr.Textbox(label="Sample ID"); ast_antibiotic=gr.Textbox(label="Antibiotic"); ast_mic=gr.Textbox(label="MIC value"); ast_unit=gr.Dropdown(["µg/mL","mg/L","other"],value="µg/mL",label="MIC unit")
            with gr.Row(): ast_result=gr.Dropdown(["Susceptible","Intermediate","Resistant","Not provided"],value="Not provided",label="AST result"); ast_ph=gr.Textbox(label="pH"); ast_temp=gr.Textbox(label="Temperature (°C)"); ast_notes=gr.Textbox(label="Other information / notes")
            ast_add=gr.Button("＋ Add AST/MIC Result",variant="primary"); ast_msg=gr.Markdown(); ast_table=gr.Dataframe(headers=AST_COLUMNS,interactive=True,wrap=True,label="AST/MIC raw laboratory data")
            ast_add.click(add_ast,[ast_table,ast_sample,ast_antibiotic,ast_mic,ast_unit,ast_result,ast_ph,ast_temp,ast_notes],[ast_table,ast_msg]).then(lambda x:x,[ast_table],ast_state)
            ast_file=gr.File(label="Upload AST/MIC CSV or Excel",file_types=[".csv",".xlsx",".xls"],type="filepath"); ast_load=gr.Button("📥 Load CSV/Excel"); ast_file_msg=gr.Markdown(); ast_load.click(load_ast_file,ast_file,[ast_table,ast_file_msg]).then(lambda x:x,[ast_table],ast_state)
            gr.Markdown("### MIC concentration series")
            mic_series=gr.Dataframe(headers=MIC_SERIES_COLUMNS,interactive=True,wrap=True); mic_calc=gr.Button("🔬 Calculate MIC"); mic_out=gr.Dataframe(interactive=False,wrap=True); mic_msg=gr.Markdown(); mic_calc.click(calculate_mic_series,mic_series,[mic_out,mic_msg])
        with gr.Tab("🌡️ Experimental Metadata"):
            meta=gr.Dataframe(headers=["Sample ID","pH","Temperature","Other information"],interactive=True,wrap=True,label="Optional additional metadata")
            gr.Markdown("If pH, temperature and notes are already in your AST/MIC spreadsheet, you do not need to enter them again.")
        with gr.Tab("⭐ Comparative Analysis"):
            gr.Markdown("### Match genomic and phenotypic data using the same Sample ID.")
            cb=gr.Button("⭐ Build Comparative Analysis",variant="primary"); cmsg=gr.Markdown(); co=gr.Dataframe(headers=COMPARISON_COLUMNS,interactive=False,wrap=True); cs=gr.Markdown()
            cb.click(run_comparative,[genomic_state,ast_state,meta],[co,cmsg]).then(lambda x:x,[co],comparison_state).then(summary,co,cs)
        with gr.Tab("📋 Sample Dashboard"):
            ds=gr.Textbox(label="Sample ID"); db=gr.Button("📋 View Dashboard",variant="primary"); dout=gr.Markdown(); db.click(dashboard,[comparison_state,ds],dout)
        with gr.Tab("📊 Visualisation"):
            gr.Markdown("Only actual user-provided measurements are plotted.")
            gdata=gr.Dataframe(headers=COMPARISON_COLUMNS,interactive=True,wrap=True); g1b=gr.Button("📊 MIC comparison"); g1=gr.Plot(); g1b.click(make_mic_plot,gdata,g1)
            g2b=gr.Button("📈 MIC vs pH"); g2=gr.Plot(); g2b.click(make_mic_ph_plot,gdata,g2)
            g3b=gr.Button("🌡️ MIC vs temperature"); g3=gr.Plot(); g3b.click(make_mic_temp_plot,gdata,g3)
            g4b=gr.Button("🧬 Genotype vs phenotype"); g4=gr.Plot(); g4b.click(make_genotype_phenotype_plot,gdata,g4)
        with gr.Tab("🤖 AMRIVA Research Assistant"): gr.Markdown("# 🤖 AMRIVA Research Assistant\nEducational explanations only. It does not diagnose or prescribe.")
        with gr.Tab("👩‍🔬 Women's Health & AMR"): gr.Markdown(WOMEN)
        with gr.Tab("📚 Methodology"): gr.Markdown(METHOD)
        with gr.Tab("⚠️ Limitations & Privacy"): gr.Markdown(LIMITS)
    gr.Markdown("---\n**AMRIVA** · Research and educational prototype — not a clinical diagnostic system.")

if __name__ == "__main__":
    demo.launch(server_name="0.0.0.0", server_port=int(os.environ.get("PORT",7860)))
