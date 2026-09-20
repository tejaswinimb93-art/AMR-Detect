import os, zipfile, tempfile, shutil, subprocess
import pandas as pd
import gradio as gr

APP_NAME = "AMRIVA"

ORGANISM_MAP = {
    "escherichia coli":"Escherichia", "e. coli":"Escherichia",
    "klebsiella pneumoniae":"Klebsiella_pneumoniae",
    "pseudomonas aeruginosa":"Pseudomonas_aeruginosa",
    "acinetobacter baumannii":"Acinetobacter_baumannii",
    "staphylococcus aureus":"Staphylococcus_aureus",
    "staphylococcus epidermidis":"Staphylococcus_epidermidis",
    "enterococcus faecalis":"Enterococcus_faecalis",
    "enterococcus faecium":"Enterococcus_faecium",
    "salmonella spp.":"Salmonella", "salmonella":"Salmonella",
    "haemophilus influenzae":"Haemophilus_influenzae",
    "neisseria gonorrhoeae":"Neisseria_gonorrhoeae",
    "streptococcus agalactiae":"Streptococcus_agalactiae",
    "streptococcus pneumoniae":"Streptococcus_pneumoniae",
    "streptococcus pyogenes":"Streptococcus_pyogenes",
    "campylobacter":"Campylobacter", "serratia marcescens":"Serratia_marcescens",
    "vibrio cholerae":"Vibrio_cholerae", "vibrio parahaemolyticus":"Vibrio_parahaemolyticus",
    "vibrio vulnificus":"Vibrio_vulnificus"
}

COMPARISON_COLUMNS = ["Sample ID","AMR Finding","Antibiotic","MIC","AST","pH","Temperature","Comparison"]


def read_fasta(filepath):
    header, parts = "", []
    with open(filepath, "r", encoding="utf-8") as f:
        for line in f:
            line=line.strip()
            if not line: continue
            if line.startswith(">"):
                if not header: header=line[1:].strip()
            else: parts.append(line.upper())
    return header, "".join(parts)


def detect_organism(header, filename):
    text=f"{header} {filename}".lower()
    for name in sorted(ORGANISM_MAP, key=len, reverse=True):
        if name in text: return name.title(), ORGANISM_MAP[name]
    return "Identification not available from FASTA metadata", None


def calculate_statistics(sequence):
    seq="".join(x for x in sequence.upper() if x in set("ACGTN"))
    if not seq: raise ValueError("No valid nucleotide sequence was found.")
    n=len(seq); a=seq.count("A"); g=seq.count("G"); c=seq.count("C"); t=seq.count("T"); unknown=seq.count("N")
    return {"sequence":seq,"length":n,"A":a,"G":g,"C":c,"T":t,"N":unknown,"GC":(g+c)/n*100}


def run_amrfinder(fasta_file, amrfinder_organism=None):
    cmd=["amrfinder","-n",fasta_file]
    if amrfinder_organism: cmd += ["-O",amrfinder_organism]
    try:
        r=subprocess.run(cmd,capture_output=True,text=True,timeout=300)
        if r.returncode != 0:
            return False, "AMRFinderPlus analysis failed.\n\n" + (r.stderr.strip() or "Unknown AMRFinderPlus error.")
        return True, r.stdout.strip() or "No known AMR determinant detected."
    except subprocess.TimeoutExpired: return False,"AMRFinderPlus analysis failed: analysis timed out."
    except FileNotFoundError: return False,"AMRFinderPlus was not found on the server."
    except Exception as e: return False,f"AMRFinderPlus analysis failed: {e}"


def interpret_amr(success, result):
    if not success:
        return "AMR analysis could not be completed", "Susceptibility cannot be inferred because genomic AMR analysis did not complete. Validated phenotypic AST is required for confirmation."
    if result.startswith("No known AMR determinant detected"):
        return "No known AMR determinant detected", "No known genomic AMR determinant was detected by this screening. This does not prove susceptibility. Phenotypic AST is required."
    return "AMR determinant(s) detected", "Detected genomic determinants may be associated with antimicrobial resistance. Absence of a detected resistance determinant does not prove susceptibility. Phenotypic AST is required for clinical confirmation."


def analyse_single(fasta_file):
    empty=("",)*14
    if not fasta_file: return ("Please upload a FASTA file.",)+empty
    try:
        header,seq=read_fasta(fasta_file)
        if not seq: raise ValueError("The FASTA file does not contain a nucleotide sequence.")
        filename=os.path.basename(fasta_file)
        organism,org=detect_organism(header,filename)
        stats=calculate_statistics(seq)
        success,result=run_amrfinder(fasta_file,org)
        status,interpretation=interpret_amr(success,result)
        sample_id=header.split()[0] if header else os.path.splitext(filename)[0]
        return ("Sample analysed successfully.",sample_id,organism,filename,str(stats["length"]),str(stats["A"]),str(stats["G"]),str(stats["C"]),str(stats["T"]),str(stats["N"]),f'{stats["GC"]:.2f}%',status,result,interpretation,header)
    except Exception as e:
        return (f"Analysis error: {e}",)+empty[:10]+("Analysis failed","","")+empty[-1:]


def extract_fasta_files(zip_path,folder):
    exts=(".fa",".fasta",".fna"); files=[]
    with zipfile.ZipFile(zip_path,"r") as z:
        for member in z.infolist():
            if member.is_dir(): continue
            name=os.path.basename(member.filename)
            if not name or name.startswith(".") or not name.lower().endswith(exts): continue
            dest=os.path.join(folder,name); base,ext=os.path.splitext(name); i=1
            while os.path.exists(dest):
                name=f"{base}_{i}{ext}"; dest=os.path.join(folder,name); i+=1
            with z.open(member) as src, open(dest,"wb") as dst: shutil.copyfileobj(src,dst)
            files.append(dest)
    return files


def analyse_zip(zip_file):
    cols=["Sample ID","Detected organism","FASTA file","Sequence length","GC content","AMR status","AMRFinderPlus result","Interpretation"]
    empty=pd.DataFrame(columns=cols)
    if not zip_file: return empty
    folder=tempfile.mkdtemp(prefix="amriva_batch_")
    try:
        if not zipfile.is_zipfile(zip_file):
            return pd.DataFrame([{"Sample ID":"ERROR","Detected organism":"N/A","FASTA file":os.path.basename(zip_file),"Sequence length":"N/A","GC content":"N/A","AMR status":"Invalid ZIP","AMRFinderPlus result":"The uploaded file is not a valid ZIP archive.","Interpretation":"Upload a ZIP containing FASTA files."}],columns=cols)
        files=extract_fasta_files(zip_file,folder)
        if not files:
            return pd.DataFrame([{"Sample ID":"ERROR","Detected organism":"N/A","FASTA file":"No FASTA files found","Sequence length":"N/A","GC content":"N/A","AMR status":"No FASTA files","AMRFinderPlus result":"No .fa, .fasta or .fna files were found inside the ZIP.","Interpretation":"Upload a ZIP containing one or more FASTA files."}],columns=cols)
        rows=[]
        for i,path in enumerate(files,1):
            name=os.path.basename(path)
            try:
                header,seq=read_fasta(path)
                if not seq:
                    rows.append({"Sample ID":f"AMR-BATCH-{i:03d}","Detected organism":"Unknown","FASTA file":name,"Sequence length":"0","GC content":"N/A","AMR status":"Invalid FASTA","AMRFinderPlus result":"No nucleotide sequence found.","Interpretation":"Analysis could not be performed."}); continue
                organism,org=detect_organism(header,name); stats=calculate_statistics(seq); success,result=run_amrfinder(path,org); status,interpretation=interpret_amr(success,result); sid=header.split()[0] if header else os.path.splitext(name)[0]
                rows.append({"Sample ID":sid,"Detected organism":organism,"FASTA file":name,"Sequence length":stats["length"],"GC content":f'{stats["GC"]:.2f}%',"AMR status":status,"AMRFinderPlus result":result,"Interpretation":interpretation})
            except Exception as e:
                rows.append({"Sample ID":f"AMR-BATCH-{i:03d}","Detected organism":"Unknown","FASTA file":name,"Sequence length":"Error","GC content":"Error","AMR status":"Analysis failed","AMRFinderPlus result":str(e),"Interpretation":"Analysis could not be completed."})
        return pd.DataFrame(rows,columns=cols)
    except Exception as e:
        return pd.DataFrame([{"Sample ID":"ERROR","Detected organism":"N/A","FASTA file":"ZIP processing","Sequence length":"N/A","GC content":"N/A","AMR status":"Batch analysis failed","AMRFinderPlus result":str(e),"Interpretation":"The ZIP could not be processed."}],columns=cols)
    finally: shutil.rmtree(folder,ignore_errors=True)



# ============================================================
# MIC CONCENTRATION-SERIES ANALYSIS
# ============================================================

MIC_SERIES_COLUMNS = [
    "Sample ID", "Antibiotic", "Concentration", "Unit", "Growth"
]

def calculate_mic_series(data):
    """Determine the lowest tested concentration marked as No growth.

    This uses only the laboratory observations supplied by the user.
    It does not assign Susceptible/Intermediate/Resistant categories.
    """
    out_columns = [
        "Sample ID", "Antibiotic", "MIC", "Unit",
        "Lowest no-growth concentration", "Interpretation note"
    ]
    empty = pd.DataFrame(columns=out_columns)
    if data is None:
        return empty, "No MIC observations entered."
    try:
        df = data.copy() if isinstance(data, pd.DataFrame) else pd.DataFrame(data)
        if df.empty:
            return empty, "No MIC observations entered."

        for col in MIC_SERIES_COLUMNS:
            if col not in df.columns:
                df[col] = ""
        df = df[MIC_SERIES_COLUMNS].copy()
        df["Concentration_num"] = pd.to_numeric(df["Concentration"], errors="coerce")
        df["Growth_clean"] = df["Growth"].astype(str).str.strip().str.lower()

        results = []
        for (sample, antibiotic), group in df.groupby(
            ["Sample ID", "Antibiotic"], dropna=False, sort=False
        ):
            group = group.dropna(subset=["Concentration_num"])
            group = group[group["Concentration_num"] >= 0]
            if group.empty:
                results.append({
                    "Sample ID": str(sample),
                    "Antibiotic": str(antibiotic),
                    "MIC": "Not determined",
                    "Unit": "",
                    "Lowest no-growth concentration": "",
                    "Interpretation note": "No valid numeric concentrations were supplied."
                })
                continue

            no_growth = group[group["Growth_clean"].isin(
                ["no", "no growth", "nogrowth", "absent", "none"]
            )]
            unit_values = group["Unit"].astype(str).str.strip()
            unit = next((u for u in unit_values if u and u.lower() != "nan"), "")

            if no_growth.empty:
                mic_text = "Not determined"
                lowest_text = "No no-growth observation"
                note = "No tested concentration was entered as No growth."
            else:
                mic = float(no_growth["Concentration_num"].min())
                mic_text = f"{mic:g} {unit}".strip()
                lowest_text = mic_text
                note = "MIC calculated as the lowest tested concentration recorded as No growth."

            results.append({
                "Sample ID": str(sample),
                "Antibiotic": str(antibiotic),
                "MIC": mic_text,
                "Unit": unit,
                "Lowest no-growth concentration": lowest_text,
                "Interpretation note": note
            })

        result_df = pd.DataFrame(results, columns=out_columns)
        return result_df, f"Calculated MIC information for {len(result_df)} sample/antibiotic combination(s)."
    except Exception as e:
        return empty, f"MIC series could not be analysed: {e}"

def add_ast_result(data, sample_id, antibiotic, mic, unit, ast_result, notes):
    columns=["Sample ID","AMR Finding","Antibiotic","MIC","AST","pH","Temperature","Comparison"]
    try:
        df=data.copy() if isinstance(data,pd.DataFrame) else pd.DataFrame(columns=columns)
        if df.empty: df=pd.DataFrame(columns=columns)
        for col in columns:
            if col not in df.columns: df[col]=""
        if not str(sample_id).strip(): return df, "Please enter a Sample ID."
        if not str(antibiotic).strip(): return df, "Please enter an antibiotic."
        if not str(mic).strip(): return df, "Please enter the MIC value."
        try:
            float(mic)
        except ValueError:
            return df, "MIC must be a numeric value."
        row={"Sample ID":str(sample_id).strip(),"AMR Finding":"","Antibiotic":str(antibiotic).strip(),"MIC":f"{str(mic).strip()} {unit}".strip(),"AST":str(ast_result).strip(),"pH":"","Temperature":"","Comparison":"","Notes":str(notes).strip()}
        df=pd.concat([df,pd.DataFrame([row])],ignore_index=True)
        return df[columns], f"Added AST/MIC result for {sample_id}."
    except Exception as e:
        return data, f"Could not add result: {e}"


def load_ast_csv(csv_file):
    columns=["Sample ID","AMR Finding","Antibiotic","MIC","AST","pH","Temperature","Comparison"]
    if not csv_file:
        return pd.DataFrame(columns=columns), "No CSV file selected."
    try:
        df=pd.read_csv(csv_file)
        aliases={"sample_id":"Sample ID","Sample_ID":"Sample ID","sample":"Sample ID","antibiotic":"Antibiotic","MIC":"MIC","mic":"MIC","AST":"AST","ast":"AST","pH":"pH","ph":"pH","temperature":"Temperature","temp":"Temperature","AMR Finding":"AMR Finding","amr_finding":"AMR Finding"}
        df=df.rename(columns={c:aliases.get(c,c) for c in df.columns})
        for col in columns:
            if col not in df.columns: df[col]=""
        df=df[columns].copy()
        return df, f"Loaded {len(df)} AST/MIC row(s) from CSV."
    except Exception as e:
        return pd.DataFrame(columns=columns), f"CSV could not be loaded: {e}"



def add_metadata(data, sample_id, temperature, ph, antibiotic, mic, ast, other):
    columns=["Sample ID","Temperature","pH","Antibiotic","MIC","AST","Other information"]
    try:
        df=data.copy() if isinstance(data,pd.DataFrame) else pd.DataFrame(columns=columns)
        if df.empty: df=pd.DataFrame(columns=columns)
        for col in columns:
            if col not in df.columns: df[col]=""
        sample_id=str(sample_id).strip()
        if not sample_id:
            return df[columns], "Please enter a Sample ID."
        if str(temperature).strip():
            try: float(temperature)
            except ValueError: return df[columns], "Temperature must be numeric."
        if str(ph).strip():
            try: float(ph)
            except ValueError: return df[columns], "pH must be numeric."
        row={"Sample ID":sample_id,"Temperature":str(temperature).strip(),"pH":str(ph).strip(),"Antibiotic":str(antibiotic).strip(),"MIC":str(mic).strip(),"AST":str(ast).strip(),"Other information":str(other).strip()}
        df=pd.concat([df,pd.DataFrame([row])],ignore_index=True)
        return df[columns], f"Added experimental metadata for {sample_id}."
    except Exception as e:
        return data, f"Could not add metadata: {e}"


def load_metadata_csv(csv_file):
    columns=["Sample ID","Temperature","pH","Antibiotic","MIC","AST","Other information"]
    if not csv_file:
        return pd.DataFrame(columns=columns), "No CSV file selected."
    try:
        df=pd.read_csv(csv_file)
        aliases={"sample_id":"Sample ID","Sample_ID":"Sample ID","sample":"Sample ID","temperature":"Temperature","temp":"Temperature","Temperature (C)":"Temperature","pH":"pH","ph":"pH","antibiotic":"Antibiotic","MIC":"MIC","mic":"MIC","AST":"AST","ast":"AST","notes":"Other information","Other Information":"Other information"}
        df=df.rename(columns={c:aliases.get(c,c) for c in df.columns})
        for col in columns:
            if col not in df.columns: df[col]=""
        return df[columns].copy(), f"Loaded {len(df)} metadata row(s) from CSV."
    except Exception as e:
        return pd.DataFrame(columns=columns), f"CSV could not be loaded: {e}"

def compare_samples(data):
    if data is None: return pd.DataFrame(columns=COMPARISON_COLUMNS)
    try:
        df=data.copy() if isinstance(data,pd.DataFrame) else pd.DataFrame(data)
        if df.empty: return pd.DataFrame(columns=COMPARISON_COLUMNS)
        for col in COMPARISON_COLUMNS:
            if col not in df.columns: df[col]=""
        df=df[COMPARISON_COLUMNS].copy(); out=[]
        for _,r in df.iterrows():
            genomic=str(r["AMR Finding"]).lower().strip(); ast=str(r["AST"]).lower().strip()
            if not genomic or not ast: c="Insufficient data"
            elif ("resistance" in genomic or "resistant" in genomic or "detected" in genomic) and ast=="resistant": c="Concordant"
            elif ("resistance" in genomic or "resistant" in genomic or "detected" in genomic) and ast=="susceptible": c="Discordant"
            else: c="Requires interpretation"
            out.append(c)
        df["Comparison"]=out
        return df
    except Exception as e: return pd.DataFrame({"Error":[str(e)]})


def comparison_summary(data):
    df=compare_samples(data)
    if df.empty: return "### No comparison data entered."
    c=(df["Comparison"]=="Concordant").sum(); d=(df["Comparison"]=="Discordant").sum(); i=(df["Comparison"]=="Insufficient data").sum()
    return f"### 📊 Summary\n\n**Samples entered:** {len(df)}\n\n🟢 Concordant: **{c}**  \n🟠 Discordant: **{d}**  \n⚪ Insufficient data: **{i}**\n\nComparison is based only on entered information and does not replace laboratory interpretation."


def dashboard(data,sample_id):
    if data is None: return "No sample data available."
    try:
        df=data.copy() if isinstance(data,pd.DataFrame) else pd.DataFrame(data)
        if df.empty or "Sample ID" not in df.columns: return "No sample data available."
        hit=df[df["Sample ID"].astype(str).str.strip()==str(sample_id).strip()]
        if hit.empty: return f"No information found for sample `{sample_id}`."
        r=hit.iloc[0]
        v=lambda x:str(r[x]) if x in r.index else "Not provided"
        return f"## 📋 Integrated Sample Dashboard\n\n### Sample: {v('Sample ID')}\n\n| Parameter | Result |\n|---|---|\n| AMR genomic finding | {v('AMR Finding')} |\n| Antibiotic | {v('Antibiotic')} |\n| MIC | {v('MIC')} |\n| AST | {v('AST')} |\n| pH | {v('pH')} |\n| Temperature | {v('Temperature')} |\n| Comparison | {v('Comparison')} |\n\n**Note:** Genomic findings do not independently establish phenotypic susceptibility or resistance. Appropriate laboratory AST is required."
    except Exception as e: return f"Dashboard error: {e}"


def make_mic_plot(data):
    try:
        import matplotlib.pyplot as plt
        df=data.copy(); df["MIC"]=pd.to_numeric(df["MIC"],errors="coerce"); df=df.dropna(subset=["MIC"])
        if df.empty: return None
        fig,ax=plt.subplots(figsize=(9,5)); ax.bar(df["Sample ID"].astype(str),df["MIC"]); ax.set_title("MIC Comparison"); ax.set_xlabel("Sample"); ax.set_ylabel("MIC"); ax.tick_params(axis="x",rotation=45); fig.tight_layout(); return fig
    except Exception: return None


def make_mic_ph_plot(data):
    try:
        import matplotlib.pyplot as plt
        df=data.copy(); df["MIC"]=pd.to_numeric(df["MIC"],errors="coerce"); df["pH"]=pd.to_numeric(df["pH"],errors="coerce"); df=df.dropna(subset=["MIC","pH"])
        if df.empty: return None
        fig,ax=plt.subplots(figsize=(8,5)); ax.scatter(df["pH"],df["MIC"]); ax.set_title("MIC vs pH"); ax.set_xlabel("pH"); ax.set_ylabel("MIC"); fig.tight_layout(); return fig
    except Exception: return None


HOME="""# 🧬 AMRIVA\n## Antimicrobial Resistance Genomic Analysis Platform\n\nAMRIVA connects **genomic AMR screening** with laboratory-provided AST/MIC and experimental information for research and education.\n\n### What AMRIVA provides\n- 🧬 FASTA genomic analysis\n- 🔬 AMRFinderPlus-based AMR screening\n- 📁 Batch analysis of any number of FASTA files in a ZIP\n- 🧪 AST and MIC entry\n- 🌡️ pH and temperature metadata\n- ⭐ Comparative analysis\n- 📋 Integrated sample dashboard\n- 📊 Data-based visualisation\n- 👩‍🔬 Women's Health & AMR research/education section\n\n> **AMRIVA is a research and educational prototype, not a clinical diagnostic system.**"""

WOMEN="""# 👩‍🔬 Women's Health & AMR\n\nAMR is relevant to infections affecting women, including urinary tract infections and other areas of reproductive-health research.\n\nAMRIVA addresses the **Women in Science** theme through research and education around genomic AMR surveillance and appropriate public/research datasets.\n\n### Research areas\n- AMR in urinary tract infections\n- AMR surveillance relevant to women's health\n- Genomic surveillance of resistance determinants\n- Genotype–phenotype comparison using research datasets\n- Publicly available research datasets\n\n### Responsible scope\nThis student prototype does **not** diagnose patient vaginal, urine or blood samples, and it does not prescribe antibiotics. Real clinical testing requires appropriate laboratory, institutional, ethical and biosafety approvals."""

METHOD="""# 📚 AMRIVA Methodology\n\n**Sample → DNA preparation → sequencing → FASTA → AMRIVA analysis → AMRFinderPlus screening → genomic findings → AST/MIC → pH/temperature metadata → comparative analysis → visualisation**\n\n### AMRFinderPlus\nAMRIVA uses AMRFinderPlus to screen genomic sequence data for known antimicrobial-resistance determinants supported by the installed software/database.\n\n**Important:** AMRIVA performs genomic screening; it does not experimentally measure antimicrobial susceptibility. Phenotypic AST is required for confirmation."""

LIMITS="""# ⚠️ Limitations & Disclaimer\n\nAMRIVA is a **research and educational prototype**.\n\n- Not a clinical diagnostic system.\n- Does not prescribe antibiotics.\n- Does not replace phenotypic AST.\n- A detected determinant does not automatically establish phenotypic resistance.\n- Absence of a detected determinant does not prove susceptibility.\n- MIC, AST, pH and temperature are user-provided laboratory/experimental values.\n- AMRIVA does not infer pH or temperature from FASTA sequence.\n\nClinical decisions require validated laboratory methods and qualified professionals."""

CSS=""".gradio-container{max-width:1400px!important;margin:auto!important}#hero{padding:36px 30px;border-radius:24px;margin-bottom:22px;background:linear-gradient(135deg,#123c69,#0f766e);color:white}#hero h1{font-size:50px!important;margin-bottom:5px!important}#hero p{font-size:18px!important}"""

with gr.Blocks(title="AMRIVA",css=CSS,theme=gr.themes.Soft()) as demo:
    gr.HTML("<div id='hero'><h1>🧬 AMRIVA</h1><p>Antimicrobial Resistance Genomic Analysis Platform</p><p>Genomic screening • AST/MIC integration • Comparative analysis</p></div>")
    with gr.Tabs():
        with gr.Tab("🏠 Home"): gr.Markdown(HOME)
        with gr.Tab("🧬 Single Sample"):
            gr.Markdown("## Single-Sample Analysis\nUpload one FASTA sequence for genomic AMR screening.")
            sf=gr.File(label="Upload FASTA",file_types=[".fa",".fasta",".fna"],type="filepath")
            sb=gr.Button("🧬 Analyse Sample",variant="primary"); ss=gr.Textbox(label="Analysis summary")
            with gr.Row(): sid=gr.Textbox(label="Sample ID"); sorg=gr.Textbox(label="Organism information")
            sfn=gr.Textbox(label="FASTA file")
            with gr.Row(): sl=gr.Textbox(label="Sequence length"); sgc=gr.Textbox(label="GC content"); sst=gr.Textbox(label="AMR status")
            with gr.Row(): sa=gr.Textbox(label="A"); sg=gr.Textbox(label="G"); sc=gr.Textbox(label="C"); st=gr.Textbox(label="T"); sn=gr.Textbox(label="N")
            samr=gr.Textbox(label="AMRFinderPlus Result",lines=14); sint=gr.Textbox(label="Genomic Interpretation",lines=6); sh=gr.Textbox(label="FASTA Header",visible=False)
            sb.click(analyse_single,sf,[ss,sid,sorg,sfn,sl,sa,sg,sc,st,sn,sgc,sst,samr,sint,sh])
        with gr.Tab("📁 Multiple Samples"):
            gr.Markdown("## Batch Analysis\nUpload **one ZIP containing any number of FASTA files**. AMRIVA analyses every `.fa`, `.fasta` and `.fna` file automatically.")
            bz=gr.File(label="Upload ZIP containing FASTA files",file_types=[".zip"],type="filepath"); bb=gr.Button("📊 Analyse All Samples",variant="primary")
            bt=gr.Dataframe(headers=["Sample ID","Detected organism","FASTA file","Sequence length","GC content","AMR status","AMRFinderPlus result","Interpretation"],interactive=False,wrap=True); bb.click(analyse_zip,bz,bt)
        with gr.Tab("🧪 AST + MIC"):
            gr.Markdown("## 🧪 Phenotypic AST & MIC Data")
            gr.Markdown("Enter **laboratory-generated** susceptibility results. AMRIVA does not perform the wet-lab test or infer AST/MIC from a FASTA sequence.")

            gr.Markdown("### 1️⃣ Quick entry — single or repeated samples")
            with gr.Row():
                ast_sample=gr.Textbox(label="Sample ID",placeholder="e.g. AMR-001")
                ast_antibiotic=gr.Textbox(label="Antibiotic",placeholder="e.g. Ciprofloxacin")
                ast_mic=gr.Textbox(label="MIC value",placeholder="e.g. 0.5")
                ast_unit=gr.Dropdown(["µg/mL","mg/L","other"],value="µg/mL",label="MIC unit")
            with gr.Row():
                ast_result=gr.Dropdown(["Susceptible","Intermediate","Resistant","Not provided"],value="Not provided",label="AST result")
                ast_notes=gr.Textbox(label="Experimental notes",placeholder="Optional")
            ast_add=gr.Button("＋ Add AST/MIC Result",variant="primary")
            ast_message=gr.Markdown()
            ast_table=gr.Dataframe(
                headers=COMPARISON_COLUMNS,
                value=[],
                datatype=["str"]*8,
                interactive=True,
                wrap=True,
                label="AST/MIC results — add as many samples/antibiotics as needed"
            )
            ast_add.click(
                add_ast_result,
                [ast_table,ast_sample,ast_antibiotic,ast_mic,ast_unit,ast_result,ast_notes],
                [ast_table,ast_message]
            )

            gr.Markdown("### 2️⃣ Multiple-sample upload")
            gr.Markdown("For many laboratory results, upload a CSV with **any number of samples and antibiotics**. Recommended columns: `Sample ID, Antibiotic, MIC, AST, pH, Temperature`. AMRIVA does not limit the number of rows to two or three.")
            ast_csv=gr.File(label="Upload AST/MIC CSV",file_types=[".csv"],type="filepath")
            ast_csv_button=gr.Button("📥 Load Multiple-Sample CSV")
            ast_csv_message=gr.Markdown()
            ast_csv_button.click(load_ast_csv,ast_csv,[ast_table,ast_csv_message])

            gr.Markdown("### 3️⃣ MIC concentration series")
            gr.Markdown("Enter the **experimentally observed growth/no-growth result at each tested concentration**. AMRIVA identifies the lowest tested concentration recorded as `No growth`. It does not independently label the sample Susceptible/Intermediate/Resistant.")
            mic_series=gr.Dataframe(
                headers=MIC_SERIES_COLUMNS,
                value=[],
                datatype=["str"]*5,
                interactive=True,
                wrap=True,
                label="MIC observations — supports multiple samples"
            )
            mic_calc=gr.Button("🔬 Calculate MIC from Observations",variant="primary")
            mic_message=gr.Markdown()
            mic_results=gr.Dataframe(
                headers=["Sample ID","Antibiotic","MIC","Unit","Lowest no-growth concentration","Interpretation note"],
                interactive=False,
                wrap=True,
                label="Calculated MIC results"
            )
            mic_calc.click(calculate_mic_series,[mic_series],[mic_results,mic_message])

            gr.Markdown("**Example:** 0.125 µg/mL → Growth; 0.25 → Growth; 0.5 → Growth; 1 → No growth. AMRIVA reports the MIC as 1 µg/mL based on the supplied observations. Clinical AST interpretation still requires the appropriate organism-, drug- and standard-specific breakpoint information.")
        with gr.Tab("🌡️ Experimental Metadata"):
            gr.Markdown("## 🌡️ Experimental Metadata")
            gr.Markdown("Add laboratory/research information such as **pH and temperature**. AMRIVA does not determine these values from FASTA. These values can later be linked to the same Sample ID in Comparative Analysis.")

            gr.Markdown("### 1️⃣ Quick entry — single or repeated samples")
            with gr.Row():
                meta_sample=gr.Textbox(label="Sample ID",placeholder="e.g. AMR-001")
                meta_temp=gr.Textbox(label="Temperature (°C)",placeholder="e.g. 37")
                meta_ph=gr.Textbox(label="pH",placeholder="e.g. 7.0")
            with gr.Row():
                meta_antibiotic=gr.Textbox(label="Antibiotic",placeholder="Optional")
                meta_mic=gr.Textbox(label="MIC",placeholder="Optional")
                meta_ast=gr.Dropdown(["Susceptible","Intermediate","Resistant","Not provided"],value="Not provided",label="AST")
            meta_other=gr.Textbox(label="Other experimental/sample information",placeholder="Optional")
            meta_add=gr.Button("＋ Add Experimental Data",variant="primary")
            meta_message=gr.Markdown()
            meta=gr.Dataframe(headers=["Sample ID","Temperature","pH","Antibiotic","MIC","AST","Other information"],value=[],datatype=["str"]*7,interactive=True,wrap=True,label="Experimental metadata — add as many samples as needed")
            meta_add.click(add_metadata,[meta,meta_sample,meta_temp,meta_ph,meta_antibiotic,meta_mic,meta_ast,meta_other],[meta,meta_message])

            gr.Markdown("### 2️⃣ Multiple-sample upload")
            gr.Markdown("Upload a CSV containing **any number of samples**. Recommended columns: `Sample ID, Temperature, pH, Antibiotic, MIC, AST, Other information`.")
            meta_csv=gr.File(label="Upload experimental metadata CSV",file_types=[".csv"],type="filepath")
            meta_csv_button=gr.Button("📥 Load Multiple-Sample Metadata")
            meta_csv_message=gr.Markdown()
            meta_csv_button.click(load_metadata_csv,meta_csv,[meta,meta_csv_message])
        with gr.Tab("⭐ Comparative Analysis"):
            gr.Markdown("## ⭐ Comparative Analysis\nCombine genomic findings, AST, MIC, pH and temperature to compare samples.")
            ci=gr.Dataframe(headers=COMPARISON_COLUMNS,value=[["","","","","","","",""]],interactive=True,wrap=True); cb=gr.Button("⭐ Run Comparative Analysis",variant="primary"); co=gr.Dataframe(label="Comparative Results",interactive=False,wrap=True); cs=gr.Markdown()
            cb.click(compare_samples,ci,co); cb.click(comparison_summary,ci,cs)
        with gr.Tab("📋 Sample Dashboard"):
            gr.Markdown("## 📋 Integrated Sample Dashboard\nEnter a sample ID to view its combined information.")
            dd=gr.Dataframe(headers=COMPARISON_COLUMNS,interactive=True,wrap=True); ds=gr.Textbox(label="Sample ID"); db=gr.Button("📋 View Dashboard",variant="primary"); dout=gr.Markdown(); db.click(dashboard,[dd,ds],dout)
        with gr.Tab("📊 Visualisation"):
            gr.Markdown("## 📊 Data Visualisation\nOnly actual user-provided measurements are visualised.")
            gd=gr.Dataframe(headers=["Sample ID","AMR Finding","Antibiotic","MIC","AST","pH","Temperature"],interactive=True,wrap=True)
            g1b=gr.Button("📊 Generate MIC Comparison"); g1=gr.Plot(label="MIC Comparison"); g1b.click(make_mic_plot,gd,g1)
            g2b=gr.Button("📈 Generate MIC vs pH"); g2=gr.Plot(label="MIC vs pH"); g2b.click(make_mic_ph_plot,gd,g2)
        with gr.Tab("🤖 AMRIVA Research Assistant"):
            gr.Markdown("## 🤖 AMRIVA Research Assistant")
            gr.Markdown("The planned assistant will explain AMRIVA results in plain language using the entered genomic, AST/MIC and experimental data. It will **not diagnose patients, prescribe antibiotics, or replace laboratory interpretation**.")
            ai_question=gr.Textbox(label="Ask about AMR, MIC, AST or your analysis",placeholder="e.g. What does MIC mean?")
            gr.Markdown("**AI assistant integration is reserved for the next build step so that it can be connected safely without exposing patient-identifying information or hard-coding an API key.**")
        with gr.Tab("👩‍🔬 Women's Health & AMR"): gr.Markdown(WOMEN)
        with gr.Tab("📚 Methodology"): gr.Markdown(METHOD)
        with gr.Tab("⚠️ Limitations"): gr.Markdown(LIMITS)
    gr.Markdown("---\n**AMRIVA** · Antimicrobial Resistance Genomic Analysis Platform\n\n*Research and educational prototype — not a clinical diagnostic system.*")

if __name__ == "__main__":
    demo.launch(server_name="0.0.0.0",server_port=int(os.environ.get("PORT",7860)))
