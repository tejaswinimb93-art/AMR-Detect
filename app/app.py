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

COMPARISON_COLUMNS = ["Sample ID","Detected organism","AMR Finding","AMRFinderPlus result","Antibiotic","MIC","MIC unit","AST","pH","Temperature","Other information","Comparison"]


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


AMRFINDER_DISPLAY_COLUMNS = [
    "Sample ID", "Gene symbol", "Sequence name", "Element type",
    "Element subtype", "Class", "Subclass", "Method",
    "% Coverage", "% Identity", "Accession of closest sequence",
    "Name of closest sequence"
]


def _amrfinder_dataframe(sample_id, raw_result):
    """Parse AMRFinderPlus TSV output into a clean display dataframe.

    AMRFinderPlus writes tab-separated output. We keep the original information
    where available and only select known display columns. If the output is an
    error/no-hit message, return a one-row status table instead of pretending
    that no result was produced.
    """
    text = str(raw_result or "").strip()
    empty = pd.DataFrame(columns=AMRFINDER_DISPLAY_COLUMNS)
    if not text or text.startswith("No known AMR determinant detected"):
        return empty
    if text.startswith("AMRFinderPlus analysis failed") or text.startswith("AMRFinderPlus was not found"):
        return empty

    try:
        from io import StringIO
        df = pd.read_csv(StringIO(text), sep="\t", dtype=str, keep_default_na=False)
        if df.empty:
            return empty
        # Normalise common AMRFinderPlus header spellings.
        aliases = {
            "Gene symbol": "Gene symbol",
            "Sequence name": "Sequence name",
            "Element type": "Element type",
            "Element subtype": "Element subtype",
            "Class": "Class",
            "Subclass": "Subclass",
            "Method": "Method",
            "% Coverage of reference sequence": "% Coverage",
            "% Identity to reference sequence": "% Identity",
            "Accession of closest sequence": "Accession of closest sequence",
            "Name of closest sequence": "Name of closest sequence",
        }
        df = df.rename(columns={c: aliases.get(str(c).strip(), str(c).strip()) for c in df.columns})
        for col in AMRFINDER_DISPLAY_COLUMNS:
            if col not in df.columns:
                df[col] = ""
        df.insert(0, "Sample ID", str(sample_id)) if "Sample ID" not in df.columns else None
        return df[["Sample ID"] + [c for c in AMRFINDER_DISPLAY_COLUMNS if c != "Sample ID"]].copy()
    except Exception:
        # Do not lose a valid result just because a future AMRFinderPlus version
        # changes a header. Try a headerless TSV fallback.
        try:
            from io import StringIO
            raw = pd.read_csv(StringIO(text), sep="\t", header=None, dtype=str, keep_default_na=False)
            if raw.empty:
                return empty
            # Standard AMRFinderPlus output has many columns. Keep the first
            # available fields and expose them with safe labels.
            n = min(raw.shape[1], 12)
            out = pd.DataFrame({
                "Sample ID": [str(sample_id)] * len(raw),
                "Gene symbol": raw.iloc[:, 5] if n > 5 else "",
                "Sequence name": raw.iloc[:, 6] if n > 6 else "",
                "Element type": raw.iloc[:, 8] if n > 8 else "",
                "Element subtype": raw.iloc[:, 9] if n > 9 else "",
                "Class": raw.iloc[:, 10] if n > 10 else "",
                "Subclass": raw.iloc[:, 11] if n > 11 else "",
            })
            for c in AMRFINDER_DISPLAY_COLUMNS:
                if c not in out.columns:
                    out[c] = ""
            return out[AMRFINDER_DISPLAY_COLUMNS]
        except Exception:
            return empty


def _amrfinder_determinants(parsed):
    """Return a concise determinant summary for the AMR status column."""
    if parsed is None or parsed.empty:
        return ""
    vals=[]
    for _, row in parsed.iterrows():
        gene=str(row.get("Gene symbol", "")).strip()
        seq=str(row.get("Sequence name", "")).strip()
        value=gene or seq
        if value and value not in vals:
            vals.append(value)
    return ", ".join(vals)


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


def interpret_amr(success, result, parsed=None):
    if not success:
        return "AMR analysis could not be completed", "AMRFinderPlus did not complete successfully. Susceptibility cannot be inferred from this genomic analysis. Validated phenotypic AST is required for confirmation."
    if parsed is None:
        parsed = pd.DataFrame()
    determinants = _amrfinder_determinants(parsed)
    if parsed.empty:
        return "No known AMR determinant detected", "AMRFinderPlus did not report a known AMR determinant in this sequence. This does not prove susceptibility. Phenotypic AST is required for confirmation."
    return f"AMR determinant(s) detected: {determinants}", "AMRFinderPlus detected the determinant(s) listed in the result table. These genomic findings may be associated with antimicrobial resistance; they do not by themselves establish phenotypic resistance. Phenotypic AST remains important for confirmation."


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
        sample_id=header.split()[0] if header else os.path.splitext(filename)[0]
        parsed=_amrfinder_dataframe(sample_id,result) if success else pd.DataFrame(columns=AMRFINDER_DISPLAY_COLUMNS)
        status,interpretation=interpret_amr(success,result,parsed)
        # The result field is now a concise determinant summary; the complete
        # AMRFinderPlus output is shown separately as a dataframe in the UI.
        result_summary = _amrfinder_determinants(parsed) if not parsed.empty else ("No known AMR determinant detected." if success else result)
        return ("Sample analysed successfully.",sample_id,organism,filename,str(stats["length"]),str(stats["A"]),str(stats["G"]),str(stats["C"]),str(stats["T"]),str(stats["N"]),f'{stats["GC"]:.2f}%',status,result_summary,interpretation,header,parsed)
    except Exception as e:
        return (f"Analysis error: {e}",)+empty[:10]+("Analysis failed","","")+empty[-1:]+(pd.DataFrame(columns=AMRFINDER_DISPLAY_COLUMNS),)


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


_LAST_BATCH_AMRFINDER_DETAILS = pd.DataFrame(columns=AMRFINDER_DISPLAY_COLUMNS)

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
        detail_frames=[]
        for i,path in enumerate(files,1):
            name=os.path.basename(path)
            try:
                header,seq=read_fasta(path)
                sid=header.split()[0] if header else os.path.splitext(name)[0]
                if not seq:
                    rows.append({"Sample ID":sid,"Detected organism":"Unknown","FASTA file":name,"Sequence length":"0","GC content":"N/A","AMR status":"Invalid FASTA","AMRFinderPlus result":"No nucleotide sequence found.","Interpretation":"Analysis could not be performed."}); continue
                organism,org=detect_organism(header,name); stats=calculate_statistics(seq); success,result=run_amrfinder(path,org)
                parsed=_amrfinder_dataframe(sid,result) if success else pd.DataFrame(columns=AMRFINDER_DISPLAY_COLUMNS)
                if isinstance(parsed, pd.DataFrame) and not parsed.empty:
                    detail_frames.append(parsed)
                status,interpretation=interpret_amr(success,result,parsed)
                determinant_summary=_amrfinder_determinants(parsed) if not parsed.empty else ("No known AMR determinant detected." if success else result)
                rows.append({"Sample ID":sid,"Detected organism":organism,"FASTA file":name,"Sequence length":stats["length"],"GC content":f'{stats["GC"]:.2f}%',"AMR status":status,"AMRFinderPlus result":determinant_summary,"Interpretation":interpretation})
            except Exception as e:
                rows.append({"Sample ID":f"AMR-BATCH-{i:03d}","Detected organism":"Unknown","FASTA file":name,"Sequence length":"Error","GC content":"Error","AMR status":"Analysis failed","AMRFinderPlus result":str(e),"Interpretation":"Analysis could not be completed."})
        global _LAST_BATCH_AMRFINDER_DETAILS
        _LAST_BATCH_AMRFINDER_DETAILS = pd.concat(detail_frames, ignore_index=True) if detail_frames else pd.DataFrame(columns=AMRFINDER_DISPLAY_COLUMNS)
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

def add_ast_result(data, genomic_data, sample_id, antibiotic, mic, unit, ast_result, ph, temperature, notes):
    columns=["Sample ID","AMR Finding","Antibiotic","MIC","MIC unit","AST","pH","Temperature","Comparison","Other information"]
    try:
        df=data.copy() if isinstance(data,pd.DataFrame) else pd.DataFrame(columns=columns)
        if df.empty: df=pd.DataFrame(columns=columns)
        for col in columns:
            if col not in df.columns: df[col]=""
        sample_id=str(sample_id).strip()
        if not sample_id: return df[columns], "Please enter a Sample ID."
        if not str(antibiotic).strip(): return df[columns], "Please enter an antibiotic."
        if not str(mic).strip(): return df[columns], "Please enter the MIC value."
        try: float(mic)
        except ValueError: return df[columns], "MIC must be a numeric value."
        for label, value in (("pH", ph), ("Temperature", temperature)):
            if str(value).strip():
                try: float(value)
                except ValueError: return df[columns], f"{label} must be numeric."

        amr_finding = ""
        if isinstance(genomic_data, pd.DataFrame) and not genomic_data.empty and "Sample ID" in genomic_data.columns:
            hits = genomic_data[genomic_data["Sample ID"].astype(str).str.strip().eq(sample_id)]
            if not hits.empty:
                if "AMR status" in hits.columns:
                    amr_finding = str(hits.iloc[0]["AMR status"])
                elif "AMR Finding" in hits.columns:
                    amr_finding = str(hits.iloc[0]["AMR Finding"])

        row={"Sample ID":sample_id,"AMR Finding":amr_finding,"Antibiotic":str(antibiotic).strip(),
             "MIC":str(mic).strip(),"MIC unit":str(unit).strip(),"AST":str(ast_result).strip(),
             "pH":str(ph).strip(),"Temperature":str(temperature).strip(),"Comparison":"",
             "Other information":str(notes).strip()}
        df=pd.concat([df,pd.DataFrame([row])],ignore_index=True)
        note = f"Added AST/MIC result for {sample_id}."
        if not amr_finding:
            note += " AMR Finding will populate after genomic analysis for the same Sample ID."
        return df[columns], note
    except Exception as e:
        return data, f"Could not add result: {e}"


def load_ast_csv(csv_file):
    columns=["Sample ID","AMR Finding","Antibiotic","MIC","MIC unit","AST","pH","Temperature","Comparison","Other information"]
    if not csv_file:
        return pd.DataFrame(columns=columns), "No file selected."
    try:
        ext=os.path.splitext(csv_file)[1].lower()
        if ext in {".xlsx", ".xls"}:
            try:
                df=pd.read_excel(csv_file)
            except ImportError as e:
                return pd.DataFrame(columns=columns), "Excel support is not installed on the server. Add openpyxl to requirements.txt and redeploy AMRIVA."
        else:
            df=pd.read_csv(csv_file)
        aliases={
            "sample_id":"Sample ID","Sample_ID":"Sample ID","sample":"Sample ID","Sample ID":"Sample ID",
            "antibiotic":"Antibiotic","Antibiotic":"Antibiotic",
            "MIC":"MIC","mic":"MIC","mic_value":"MIC","MIC value":"MIC",
            "MIC unit":"MIC unit","mic_unit":"MIC unit","Unit":"MIC unit","unit":"MIC unit",
            "AST":"AST","ast":"AST","AST result":"AST","ast_result":"AST",
            "pH":"pH","ph":"pH","PH":"pH",
            "temperature":"Temperature","temp":"Temperature","Temperature (C)":"Temperature","Temperature (°C)":"Temperature",
            "AMR Finding":"AMR Finding","amr_finding":"AMR Finding","AMR status":"AMR Finding",
            "Comparison":"Comparison","comparison":"Comparison",
            "Other information":"Other information","Other Information":"Other information","Notes":"Other information","notes":"Other information","Experimental notes":"Other information"
        }
        df=df.rename(columns={c:aliases.get(str(c).strip(),str(c).strip()) for c in df.columns})
        # Also support a headerless spreadsheet using the documented column order.
        required_basic={"Sample ID","Antibiotic","MIC"}
        if not required_basic.issubset(set(df.columns)) and len(df.columns) >= 3:
            positional=["Sample ID","Antibiotic","MIC","MIC unit","AST","pH","Temperature","Other information"]
            df.columns=positional[:len(df.columns)]
        for col in columns:
            if col not in df.columns: df[col]=""
        df=df[columns].copy().fillna("")
        return df, f"Loaded {len(df)} AST/MIC row(s) from {ext.upper().replace('.', '')} file."
    except Exception as e:
        return pd.DataFrame(columns=columns), f"AST/MIC spreadsheet could not be loaded: {e}"


def parse_pasted_ast_data(pasted_text):
    columns=["Sample ID","AMR Finding","Antibiotic","MIC","MIC unit","AST","pH","Temperature","Comparison","Other information"]
    if not pasted_text or not str(pasted_text).strip():
        return pd.DataFrame(columns=columns), "Nothing was pasted."
    try:
        # Excel copy/paste normally produces tab-separated rows. Accept tabs, commas,
        # or semicolons so users can paste from common spreadsheet exports.
        text=str(pasted_text).strip()
        from io import StringIO
        try:
            pasted=pd.read_csv(StringIO(text), sep="\t")
            if len(pasted.columns) == 1:
                pasted=pd.read_csv(StringIO(text), sep=None, engine="python")
        except Exception:
            pasted=pd.read_csv(StringIO(text), sep=None, engine="python")
        aliases={"sample_id":"Sample ID","Sample_ID":"Sample ID","sample":"Sample ID",
                 "antibiotic":"Antibiotic","Antibiotic":"Antibiotic","MIC":"MIC","mic":"MIC",
                 "MIC unit":"MIC unit","mic_unit":"MIC unit","Unit":"MIC unit","unit":"MIC unit",
                 "AST":"AST","ast":"AST","AST result":"AST","pH":"pH","ph":"pH",
                 "temperature":"Temperature","temp":"Temperature","AMR Finding":"AMR Finding",
                 "amr_finding":"AMR Finding","Other information":"Other information","Notes":"Other information","notes":"Other information"}
        pasted=pasted.rename(columns={c:aliases.get(str(c).strip(),str(c).strip()) for c in pasted.columns})
        # If the pasted spreadsheet has no header row, use the documented column order.
        if not any(c in pasted.columns for c in ["Sample ID","Antibiotic","MIC"]):
            pasted=pd.read_csv(StringIO(text), sep="\t", header=None)
            pasted.columns=(['Sample ID','Antibiotic','MIC','MIC unit','AST','pH','Temperature','Other information'][:len(pasted.columns)])
        for col in columns:
            if col not in pasted.columns: pasted[col]=""
        pasted=pasted[columns].copy().fillna("")
        return pasted, f"Pasted {len(pasted)} AST/MIC row(s)."
    except Exception as e:
        return pd.DataFrame(columns=columns), f"Pasted data could not be read: {e}"


def merge_ast_sources(file_df, paste_df):
    columns=["Sample ID","AMR Finding","Antibiotic","MIC","MIC unit","AST","pH","Temperature","Comparison","Other information"]
    frames=[]
    for df in (file_df,paste_df):
        if isinstance(df,pd.DataFrame) and not df.empty:
            x=df.copy()
            for c in columns:
                if c not in x.columns: x[c]=""
            frames.append(x[columns])
    if not frames:
        return pd.DataFrame(columns=columns)
    return pd.concat(frames,ignore_index=True).fillna("")


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
        aliases={"sample_id":"Sample ID","Sample_ID":"Sample ID","sample":"Sample ID","temperature":"Temperature","temp":"Temperature","Temperature (C)":"Temperature","pH":"pH","ph":"pH","antibiotic":"Antibiotic","MIC":"MIC","mic":"MIC","MIC unit":"MIC unit","mic_unit":"MIC unit","Unit":"MIC unit","unit":"MIC unit","AST":"AST","ast":"AST","notes":"Other information","Other Information":"Other information"}
        df=df.rename(columns={c:aliases.get(c,c) for c in df.columns})
        for col in columns:
            if col not in df.columns: df[col]=""
        return df[columns].copy(), f"Loaded {len(df)} metadata row(s) from CSV."
    except Exception as e:
        return pd.DataFrame(columns=columns), f"CSV could not be loaded: {e}"

def analyse_single_with_state(fasta_file):
    result = analyse_single(fasta_file)
    empty_genomic = pd.DataFrame(columns=["Sample ID","Detected organism","AMR status","AMR Finding","AMRFinderPlus result"])
    empty_detail = pd.DataFrame(columns=AMRFINDER_DISPLAY_COLUMNS)
    if not fasta_file or len(result) < 16 or result[1] == "":
        base = list(result[:15])
        if len(base) >= 13: base[12] = empty_detail
        return (*base, empty_genomic)
    parsed=result[15]
    determinant_summary=_amrfinder_determinants(parsed) if isinstance(parsed,pd.DataFrame) and not parsed.empty else result[12]
    genomic = pd.DataFrame([{
        "Sample ID": result[1],
        "Detected organism": result[2],
        "AMR status": result[11],
        "AMR Finding": determinant_summary,
        "AMRFinderPlus result": determinant_summary
    }])
    base=list(result[:15])
    base[12]=parsed
    return (*base, genomic)


def analyse_zip_with_state(zip_file):
    global _LAST_BATCH_AMRFINDER_DETAILS
    df = analyse_zip(zip_file)
    return df, df.copy(), _LAST_BATCH_AMRFINDER_DETAILS.copy()


def _normalise_columns(df, aliases):
    if df is None:
        return pd.DataFrame()
    if not isinstance(df, pd.DataFrame):
        df = pd.DataFrame(df)
    rename = {c: aliases.get(str(c).strip(), str(c).strip()) for c in df.columns}
    return df.rename(columns=rename)


def build_comparative_analysis(genomic_data, ast_data, metadata_data):
    columns = ["Sample ID","Detected organism","AMR Finding","AMRFinderPlus result","Antibiotic","MIC","MIC unit","AST","pH","Temperature","Other information","Comparison"]
    empty = pd.DataFrame(columns=columns)
    try:
        g = _normalise_columns(genomic_data, {
            "Sample ID":"Sample ID", "Detected organism":"Detected organism",
            "AMR status":"AMR status", "AMRFinderPlus result":"AMRFinderPlus result",
            "AMR Finding":"AMR Finding"
        })
        a = _normalise_columns(ast_data, {
            "Sample ID":"Sample ID", "Antibiotic":"Antibiotic", "MIC":"MIC",
            "MIC unit":"MIC unit", "AST":"AST", "AST result":"AST",
            "Experimental notes":"Other information", "Other information":"Other information"
        })
        m = _normalise_columns(metadata_data, {
            "Sample ID":"Sample ID", "Temperature":"Temperature", "Temperature (°C)":"Temperature",
            "pH":"pH", "Antibiotic":"Antibiotic", "MIC":"MIC", "AST":"AST",
            "Other information":"Other information", "Other Information":"Other information"
        })

        if a.empty and m.empty:
            return empty, "Add AST/MIC or experimental metadata first."

        # AST is the primary phenotypic table. If metadata has antibiotic-level rows,
        # retain those values as additional information rather than duplicating blindly.
        if a.empty:
            base = m.copy()
        else:
            base = a.copy()

        for c in ["Sample ID","Antibiotic","MIC","MIC unit","AST","Other information"]:
            if c not in base.columns: base[c] = ""

        # Add sample-level metadata by Sample ID. Metadata can be sample-level
        # (blank antibiotic) or antibiotic-specific. Blank values in the AST table
        # are filled from metadata without creating duplicate pH/temperature columns.
        if not m.empty and "Sample ID" in m.columns:
            for c in ["pH", "Temperature"]:
                if c not in m.columns: m[c] = ""
            base["Sample ID"] = base["Sample ID"].astype(str).str.strip()
            m["Sample ID"] = m["Sample ID"].astype(str).str.strip()

            if "Antibiotic" in m.columns and "Antibiotic" in base.columns:
                m["Antibiotic"] = m["Antibiotic"].astype(str).str.strip()
                base["Antibiotic"] = base["Antibiotic"].astype(str).str.strip()

                exact = m[m["Antibiotic"].ne("")][["Sample ID","Antibiotic","pH","Temperature"]].copy()
                exact = exact.drop_duplicates(["Sample ID","Antibiotic"], keep="last")
                exact = exact.rename(columns={"pH":"pH_meta_exact","Temperature":"Temperature_meta_exact"})
                base = base.merge(exact, on=["Sample ID","Antibiotic"], how="left")

                sample_meta = m[m["Antibiotic"].eq("")][["Sample ID","pH","Temperature"]].copy()
                sample_meta = sample_meta.drop_duplicates("Sample ID", keep="last")
                sample_meta = sample_meta.rename(columns={"pH":"pH_meta_sample","Temperature":"Temperature_meta_sample"})
                if not sample_meta.empty:
                    base = base.merge(sample_meta, on="Sample ID", how="left")

                for target, candidates in {
                    "pH":["pH_meta_exact","pH_meta_sample"],
                    "Temperature":["Temperature_meta_exact","Temperature_meta_sample"]
                }.items():
                    if target not in base.columns: base[target] = ""
                    base[target] = base[target].replace({"": pd.NA, "nan": pd.NA, "None": pd.NA})
                    for candidate in candidates:
                        if candidate in base.columns:
                            base[target] = base[target].fillna(base[candidate])
                    base[target] = base[target].fillna("")

                drop_cols=[c for c in ["pH_meta_exact","Temperature_meta_exact","pH_meta_sample","Temperature_meta_sample"] if c in base.columns]
                if drop_cols: base.drop(columns=drop_cols, inplace=True)
            else:
                mm=m[["Sample ID","pH","Temperature"]].drop_duplicates("Sample ID", keep="last").rename(columns={"pH":"pH_meta","Temperature":"Temperature_meta"})
                base=base.merge(mm,on="Sample ID",how="left")
                for target,meta_col in [("pH","pH_meta"),("Temperature","Temperature_meta")]:
                    if target not in base.columns: base[target]=""
                    if meta_col in base.columns:
                        base[target]=base[target].replace({"":pd.NA}).fillna(base[meta_col]).fillna("")
                        base.drop(columns=[meta_col],inplace=True)

        if "pH" not in base.columns: base["pH"] = ""
        if "Temperature" not in base.columns: base["Temperature"] = ""

        # Attach genomic findings by Sample ID.
        if not g.empty and "Sample ID" in g.columns:
            g["Sample ID"] = g["Sample ID"].astype(str).str.strip()
            g["AMR Finding"] = g.get("AMR Finding", "")
            g["AMRFinderPlus result"] = g.get("AMRFinderPlus result", "")
            if g["AMR Finding"].astype(str).str.strip().eq("").all():
                status = g.get("AMR status", "").astype(str)
                result = g.get("AMRFinderPlus result", "").astype(str)
                g["AMR Finding"] = result.where(~result.str.contains("No known AMR determinant detected", case=False, na=False), "No known AMR determinant detected")
            keep = [c for c in ["Sample ID","Detected organism","AMR Finding","AMRFinderPlus result"] if c in g.columns]
            for old_col in ["Detected organism","AMR Finding","AMRFinderPlus result"]:
                if old_col in base.columns:
                    base.drop(columns=[old_col], inplace=True)
            base = base.merge(g[keep].drop_duplicates("Sample ID"), on="Sample ID", how="left")
        else:
            base["Detected organism"] = ""
            base["AMR Finding"] = ""

        if "Detected organism" not in base.columns: base["Detected organism"] = ""
        if "AMR Finding" not in base.columns: base["AMR Finding"] = ""
        if "AMRFinderPlus result" not in base.columns: base["AMRFinderPlus result"] = ""
        if "Other information" not in base.columns: base["Other information"] = ""
        if "MIC unit" not in base.columns: base["MIC unit"] = ""

        # Conservative comparison: only flag a possible genotype/phenotype relationship
        # when the genomic result contains a recognizable antibiotic/class keyword.
        # This is a screening aid, not a clinical susceptibility call.
        CLASS_KEYWORDS = {
            "ciprofloxacin": ["fluoroquinolone", "quinolone", "ciprofloxacin"],
            "levofloxacin": ["fluoroquinolone", "quinolone", "levofloxacin"],
            "moxifloxacin": ["fluoroquinolone", "quinolone", "moxifloxacin"],
            "ampicillin": ["beta-lactam", "beta lactam", "ampicillin", "penicillin"],
            "amoxicillin": ["beta-lactam", "beta lactam", "amoxicillin", "penicillin"],
            "ceftriaxone": ["beta-lactam", "beta lactam", "cephalosporin", "ceftriaxone"],
            "cefotaxime": ["beta-lactam", "beta lactam", "cephalosporin", "cefotaxime"],
            "meropenem": ["carbapenem", "beta-lactam", "beta lactam", "meropenem"],
            "imipenem": ["carbapenem", "beta-lactam", "beta lactam", "imipenem"],
            "gentamicin": ["aminoglycoside", "gentamicin"],
            "amikacin": ["aminoglycoside", "amikacin"],
            "streptomycin": ["aminoglycoside", "streptomycin"],
            "tetracycline": ["tetracycline", "tet"],
            "doxycycline": ["tetracycline", "doxycycline"],
            "trimethoprim": ["trimethoprim", "sulfonamide"],
            "sulfamethoxazole": ["sulfonamide", "sulfamethoxazole"],
            "vancomycin": ["vancomycin", "glycopeptide"],
            "linezolid": ["linezolid", "oxazolidinone"],
            "erythromycin": ["macrolide", "erythromycin"],
            "azithromycin": ["macrolide", "azithromycin"],
            "chloramphenicol": ["chloramphenicol"],
            "rifampicin": ["rifampin", "rifamycin", "rifampicin"],
        }

        def classify(row):
            finding = str(row.get("AMR Finding", "")).strip().lower()
            raw_result = str(row.get("AMRFinderPlus result", "")).strip().lower()
            evidence = f"{finding} {raw_result}"
            ast = str(row.get("AST", "")).strip().lower()
            antibiotic = str(row.get("Antibiotic", "")).strip().lower()
            if not ast or ast in {"not provided", "nan", "none", ""}:
                return "Insufficient data"
            if not evidence or "no known amr determinant detected" in evidence:
                return "Insufficient data"
            keywords = CLASS_KEYWORDS.get(antibiotic, [antibiotic] if antibiotic else [])
            relevant = any(k and k in evidence for k in keywords)

            # Also recognize common AMRFinderPlus gene/family names that imply
            # resistance to the same antibiotic class. This is a screening-level
            # genotype/phenotype comparison, not a clinical susceptibility call.
            gene_class_keywords = {
                "ciprofloxacin": ["gyrA", "parC", "qnr", "aac(6')-Ib-cr", "quinolone"],
                "levofloxacin": ["gyrA", "parC", "qnr", "quinolone"],
                "moxifloxacin": ["gyrA", "parC", "qnr", "quinolone"],
                "ampicillin": ["blaTEM", "blaSHV", "blaCTX-M", "blaOXA", "penicillinase", "beta-lactam"],
                "amoxicillin": ["blaTEM", "blaSHV", "blaCTX-M", "blaOXA", "beta-lactam"],
                "ceftriaxone": ["blaCTX-M", "blaSHV", "blaTEM", "blaCMY", "cephalosporin", "beta-lactam"],
                "cefotaxime": ["blaCTX-M", "blaSHV", "blaTEM", "blaCMY", "cephalosporin", "beta-lactam"],
                "meropenem": ["blaNDM", "blaKPC", "blaVIM", "blaIMP", "blaOXA-48", "carbapenem"],
                "imipenem": ["blaNDM", "blaKPC", "blaVIM", "blaIMP", "blaOXA-48", "carbapenem"],
                "gentamicin": ["aac", "aph", "ant", "aminoglycoside"],
                "amikacin": ["aac", "aph", "ant", "aminoglycoside"],
                "tetracycline": ["tet", "tetracycline"],
                "doxycycline": ["tet", "tetracycline"],
                "trimethoprim": ["dfr", "trimethoprim"],
                "sulfamethoxazole": ["sul", "sulfamethoxazole", "sulfonamide"],
                "vancomycin": ["vanA", "vanB", "vancomycin"],
                "erythromycin": ["erm", "macrolide"],
                "azithromycin": ["erm", "macrolide"],
            }
            gene_keywords = gene_class_keywords.get(antibiotic, [])
            relevant = relevant or any(k.lower() in evidence for k in gene_keywords)

            if not relevant:
                return "Requires interpretation"
            if ast == "resistant":
                return "Potentially concordant"
            if ast == "susceptible":
                return "Potentially discordant"
            if ast == "intermediate":
                return "Requires interpretation"
            return "Insufficient data"

        base["Comparison"] = base.apply(classify, axis=1)
        out = base[columns].fillna("")
        return out, f"Built comparative analysis for {len(out)} antibiotic test row(s)."
    except Exception as e:
        return empty, f"Comparative analysis could not be built: {e}"


def compare_samples(data):
    # Kept for compatibility with older calls.
    if data is None: return pd.DataFrame(columns=COMPARISON_COLUMNS)
    df = data.copy() if isinstance(data,pd.DataFrame) else pd.DataFrame(data)
    for col in COMPARISON_COLUMNS:
        if col not in df.columns: df[col] = ""
    return df[COMPARISON_COLUMNS]


def run_comparative(genomic_data, ast_data, metadata_data):
    result, message = build_comparative_analysis(genomic_data, ast_data, metadata_data)
    summary = comparison_summary(result)
    return result, message, summary

def comparison_summary(data):
    if data is None: return "### No comparison data entered."
    df=data.copy() if isinstance(data,pd.DataFrame) else pd.DataFrame(data)
    if df.empty: return "### No comparison data entered."
    c=(df.get("Comparison",pd.Series(dtype=str)).astype(str)=="Potentially concordant").sum()
    d=(df.get("Comparison",pd.Series(dtype=str)).astype(str)=="Potentially discordant").sum()
    i=(df.get("Comparison",pd.Series(dtype=str)).astype(str)=="Insufficient data").sum()
    r=(df.get("Comparison",pd.Series(dtype=str)).astype(str)=="Requires interpretation").sum()
    return f"### 📊 Summary\n\n**Antibiotic test rows:** {len(df)}\n\n🟢 Potentially concordant: **{c}**  \n🟠 Potentially discordant: **{d}**  \n🟡 Requires interpretation: **{r}**  \n⚪ Insufficient data: **{i}**\n\nThese are screening-level comparisons only. AMRIVA does not infer clinical susceptibility from a genomic determinant alone; organism-, drug- and standard-specific interpretation remains necessary."


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



# ============================================================
# SIMPLIFIED COMPETITION WORKFLOW
# ============================================================

def teacher_antibiotic_comparison(genomic_data, ast_data, selected_antibiotic):
    """Compare samples that are phenotypically Resistant to one selected antibiotic.

    This is the core competition feature: same antibiotic phenotype -> compare
    genomic resistance-associated determinants across different organisms.
    """
    cols=["Sample ID","Organism","AST","AMR determinant","Mechanism / class","Comparison note"]
    empty=pd.DataFrame(columns=cols)
    try:
        g=genomic_data.copy() if isinstance(genomic_data,pd.DataFrame) else pd.DataFrame(genomic_data)
        a=ast_data.copy() if isinstance(ast_data,pd.DataFrame) else pd.DataFrame(ast_data)
        if a.empty:
            return empty, "Add simple AST results first."
        for c in ["Sample ID","Antibiotic","AST"]:
            if c not in a.columns: a[c]=""
        a["Sample ID"]=a["Sample ID"].astype(str).str.strip()
        a["Antibiotic"]=a["Antibiotic"].astype(str).str.strip()
        a["AST"]=a["AST"].astype(str).str.strip()
        selected=str(selected_antibiotic or "").strip()
        if not selected:
            return empty, "Select an antibiotic first."
        resistant=a[(a["Antibiotic"].str.casefold()==selected.casefold()) & (a["AST"].str.casefold()=="resistant")].copy()
        if resistant.empty:
            return empty, f"No samples are currently marked Resistant to {selected}."

        if g.empty:
            g=pd.DataFrame(columns=["Sample ID","Detected organism","AMR Finding","AMRFinderPlus result"])
        for c in ["Sample ID","Detected organism","AMR Finding","AMRFinderPlus result"]:
            if c not in g.columns: g[c]=""
        g["Sample ID"]=g["Sample ID"].astype(str).str.strip()
        g= g.drop_duplicates("Sample ID", keep="last")
        merged=resistant.merge(g[["Sample ID","Detected organism","AMR Finding","AMRFinderPlus result"]],on="Sample ID",how="left")
        rows=[]
        for _,r in merged.iterrows():
            finding=str(r.get("AMR Finding","")).strip()
            raw=str(r.get("AMRFinderPlus result","")).strip()
            if not finding: finding=raw or "No genomic finding linked yet"
            mechanism=""
            text=f"{finding} {raw}".lower()
            if any(x in text for x in ["bla", "beta-lactam", "beta lactam"]): mechanism="β-lactam / β-lactamase-associated finding"
            elif any(x in text for x in ["gyrA","parC","qnr","quinolone"]): mechanism="Quinolone-resistance-associated finding"
            elif any(x in text for x in ["tet","tetracycline"]): mechanism="Tetracycline-resistance-associated finding"
            elif any(x in text for x in ["vanA","vanB","vancomycin"]): mechanism="Vancomycin-resistance-associated finding"
            elif any(x in text for x in ["erm","macrolide"]): mechanism="Macrolide-resistance-associated finding"
            else: mechanism="Genomic resistance-associated finding"
            rows.append({"Sample ID":r["Sample ID"],"Organism":r.get("Detected organism","") or "Not available","AST":r["AST"],"AMR determinant":finding,"Mechanism / class":mechanism,"Comparison note":"Same observed phenotype; compare the genomic determinant with the other resistant samples."})
        out=pd.DataFrame(rows,columns=cols)
        determinants=out["AMR determinant"].astype(str).tolist()
        unique=[]
        for x in determinants:
            if x not in unique and x not in ["", "No genomic finding linked yet"]: unique.append(x)
        if len(unique)>1:
            note=f"### Same-antibiotic comparison\n**{selected}**: {len(out)} resistant sample(s). Different genomic determinant(s) are represented in the uploaded results, so the researcher can investigate whether the observed phenotype has different genomic bases across organisms."
        elif len(unique)==1:
            note=f"### Same-antibiotic comparison\n**{selected}**: {len(out)} resistant sample(s). The uploaded results show the same reported determinant across the linked samples."
        else:
            note=f"### Same-antibiotic comparison\n**{selected}**: {len(out)} resistant sample(s), but genomic findings have not been linked to all of them yet."
        note += "\n\n**Important:** AST provides the observed laboratory phenotype; genomic detection provides resistance-associated evidence. AMRIVA does not claim that a genomic determinant alone proves phenotypic resistance."
        return out,note
    except Exception as e:
        return empty,f"Comparison could not be built: {e}"


def list_antibiotics(ast_data):
    try:
        a=ast_data.copy() if isinstance(ast_data,pd.DataFrame) else pd.DataFrame(ast_data)
        if a.empty or "Antibiotic" not in a.columns: return gr.update(choices=[],value=None)
        vals=[]
        for x in a["Antibiotic"].astype(str):
            x=x.strip()
            if x and x.lower() not in {"nan","none"} and x not in vals: vals.append(x)
        return gr.update(choices=vals,value=(vals[0] if vals else None))
    except Exception:
        return gr.update(choices=[],value=None)


def merge_simple_ast(current, new):
    cols=["Sample ID","AMR Finding","Antibiotic","MIC","MIC unit","AST","pH","Temperature","Comparison","Other information"]
    frames=[]
    for d in [current,new]:
        if isinstance(d,pd.DataFrame) and not d.empty:
            x=d.copy()
            for c in cols:
                if c not in x.columns: x[c]=""
            frames.append(x[cols])
    return pd.concat(frames,ignore_index=True).fillna("") if frames else pd.DataFrame(columns=cols)


def simple_ast_add(data,sample_id,organism,antibiotic,ast_value):
    cols=["Sample ID","Organism","Antibiotic","AST"]
    try:
        df=data.copy() if isinstance(data,pd.DataFrame) else pd.DataFrame(columns=cols)
        if df.empty: df=pd.DataFrame(columns=cols)
        for c in cols:
            if c not in df.columns: df[c]=""
        sid=str(sample_id).strip(); org=str(organism).strip(); drug=str(antibiotic).strip(); ast=str(ast_value).strip()
        if not sid or not drug: return df[cols],"Enter Sample ID and antibiotic."
        row=pd.DataFrame([{"Sample ID":sid,"Organism":org,"Antibiotic":drug,"AST":ast}])
        return pd.concat([df[cols],row],ignore_index=True),f"Added AST result for {sid}."
    except Exception as e: return data,f"Could not add AST result: {e}"


HOME_SIMPLE="""# 🧬 AMRIVA\n## Genomic AMR Research Workspace\n\n**AMRIVA uses AMRFinderPlus for genomic screening and provides a researcher-focused workflow to organize and compare the findings.**\n\n### The central research question\n> When different organisms show resistance to the **same antibiotic**, do they show the **same or different resistance-associated genomic determinants**?\n\n### Workflow\n**FASTA → AMRFinderPlus → genomic finding → simple AST → select one antibiotic → compare resistant samples → inspect the sequence/AMRFinderPlus evidence → further laboratory validation**\n\nAMRIVA is a research/educational prototype. It does not replace AMRFinderPlus, AST, or clinical decision-making."""

METHOD_SIMPLE="""# 📚 Methodology\n\n1. The researcher supplies bacterial FASTA sequence data.\n2. AMRIVA sends the sequence through the existing AMRFinderPlus workflow.\n3. AMRFinderPlus reports known resistance-associated determinants supported by its database/methods.\n4. The researcher can enter a simple laboratory AST result: **Susceptible / Intermediate / Resistant**.\n5. The researcher selects an antibiotic and AMRIVA identifies samples recorded as **Resistant** to that same antibiotic.\n6. AMRIVA places the linked genomic determinants side by side so the researcher can investigate whether the same phenotype has the same or different genomic basis.\n7. Individual AMRFinderPlus findings can be inspected in detail, including the available sequence-analysis evidence.\n\n**Important:** genomic detection is evidence, not a standalone clinical susceptibility result. Appropriate laboratory validation remains necessary."""

LIMITS_SIMPLE="""# ⚠️ Limitations\n\n- AMRIVA is a research/educational prototype, not a clinical diagnostic system.\n- AMRFinderPlus performs the underlying genomic detection; AMRIVA does not claim to replace it.\n- A detected resistance-associated determinant does not by itself prove phenotypic resistance.\n- AST values are laboratory observations supplied by the user.\n- AMRIVA does not prescribe antibiotics.\n- Results depend on the quality of the input sequence and the capabilities/reference database of the underlying AMRFinderPlus installation.\n"""

with gr.Blocks(title="AMRIVA — Genomic AMR Research Workspace", theme=gr.themes.Soft()) as demo:
    genomic_state=gr.State(pd.DataFrame(columns=["Sample ID","Detected organism","AMR status","AMR Finding","AMRFinderPlus result"]))
    ast_state=gr.State(pd.DataFrame(columns=["Sample ID","Organism","Antibiotic","AST"]))
    gr.HTML("<div style='padding:30px;border-radius:22px;background:linear-gradient(135deg,#123c69,#0f766e);color:white'><h1 style='font-size:48px;margin:0'>🧬 AMRIVA</h1><p style='font-size:18px'>Genomic AMR Research Workspace</p><p>Sequence screening • simple AST • same-antibiotic genomic comparison</p></div>")
    with gr.Tabs():
        with gr.Tab("🏠 Home"):
            gr.Markdown(HOME_SIMPLE)
        with gr.Tab("🧬 Genomic Analysis"):
            gr.Markdown("## FASTA → AMRFinderPlus")
            gr.Markdown("Upload one FASTA. Your existing AMRFinderPlus connection is used on the server; the result is displayed as a detailed table and a concise researcher-friendly summary.")
            sf=gr.File(label="Upload FASTA",file_types=[".fa",".fasta",".fna"],type="filepath")
            sb=gr.Button("🧬 Analyse with AMRFinderPlus",variant="primary")
            msg=gr.Markdown()
            with gr.Row(): sid=gr.Textbox(label="Sample ID"); organism=gr.Textbox(label="Detected organism")
            with gr.Row(): status=gr.Textbox(label="AMR status"); summary=gr.Textbox(label="Detected AMR determinant(s)")
            detail=gr.Dataframe(headers=AMRFINDER_DISPLAY_COLUMNS,interactive=False,wrap=True,label="AMRFinderPlus result")
            interpretation=gr.Markdown()
            sb.click(analyse_single_with_state,sf,[msg,sid,organism,gr.Textbox(visible=False),gr.Textbox(visible=False),gr.Textbox(visible=False),gr.Textbox(visible=False),gr.Textbox(visible=False),gr.Textbox(visible=False),gr.Textbox(visible=False),gr.Textbox(visible=False),status,summary,interpretation,gr.Textbox(visible=False),detail,genomic_state])
        with gr.Tab("📁 Batch Samples"):
            gr.Markdown("## Multiple FASTA Samples")
            gr.Markdown("Upload one ZIP containing any number of FASTA files. Every FASTA is analysed with the existing AMRFinderPlus workflow.")
            zf=gr.File(label="ZIP containing FASTA files",file_types=[".zip"],type="filepath")
            zb=gr.Button("📊 Analyse all samples",variant="primary")
            batch_msg=gr.Markdown(); batch=gr.Dataframe(headers=["Sample ID","Detected organism","FASTA file","Sequence length","GC content","AMR status","AMRFinderPlus result","Interpretation"],interactive=False,wrap=True,label="Batch summary")
            batch_detail=gr.Dataframe(headers=AMRFINDER_DISPLAY_COLUMNS,interactive=False,wrap=True,label="Detailed AMRFinderPlus findings")
            zb.click(analyse_zip_with_state,zf,[batch,genomic_state,batch_detail])
        with gr.Tab("🧪 Simple AST"):
            gr.Markdown("## Simple Phenotypic AST Input")
            gr.Markdown("Enter the **observed laboratory result**. We deliberately keep this simple: no breakpoint calculations and no automated clinical interpretation.")
            with gr.Row():
                asid=gr.Textbox(label="Sample ID",placeholder="S01")
                aorg=gr.Textbox(label="Organism",placeholder="E. coli")
                adrug=gr.Textbox(label="Antibiotic",placeholder="Penicillin")
                aast=gr.Dropdown(["Resistant","Intermediate","Susceptible"],value="Resistant",label="AST")
            aadd=gr.Button("＋ Add AST result",variant="primary"); amsg=gr.Markdown()
            ast_table=gr.Dataframe(headers=["Sample ID","Organism","Antibiotic","AST"],interactive=True,wrap=True,label="AST results")
            aadd.click(simple_ast_add,[ast_table,asid,aorg,adrug,aast],[ast_table,amsg])
            ast_table.change(lambda x: x, ast_table, ast_state)
        with gr.Tab("⭐ Same-Antibiotic Comparison"):
            gr.Markdown("## ⭐ Teacher's Core Comparison")
            gr.Markdown("**Main idea:** first use AST to identify different organisms that are resistant to the **same antibiotic**. Then compare the genomic resistance-associated determinants detected in those samples.")
            antibiotic_choice=gr.Dropdown(label="Select antibiotic",choices=[])
            refresh=gr.Button("🔄 Refresh antibiotic list")
            compare_button=gr.Button("⭐ Compare resistant samples",variant="primary")
            compare_msg=gr.Markdown()
            compare_table=gr.Dataframe(headers=["Sample ID","Organism","AST","AMR determinant","Mechanism / class","Comparison note"],interactive=False,wrap=True,label="Same-antibiotic genomic comparison")
            compare_summary=gr.Markdown()
            refresh.click(list_antibiotics,ast_table,antibiotic_choice)
            compare_button.click(teacher_antibiotic_comparison,[genomic_state,ast_table,antibiotic_choice],[compare_table,compare_summary])
            ast_table.change(list_antibiotics,ast_table,antibiotic_choice)
        with gr.Tab("🔎 Finding Details"):
            gr.Markdown("## 🔎 Inspect an individual genomic finding")
            gr.Markdown("The detailed AMRFinderPlus result from the genomic analysis is shown here. Use the underlying result/reference as the traceable evidence for the finding.")
            finding_detail=gr.Dataframe(headers=AMRFINDER_DISPLAY_COLUMNS,interactive=False,wrap=True,label="AMRFinderPlus evidence")
            gr.Markdown("**Sequence-level note:** the exact nucleotide sequence/coordinates available from your AMRFinderPlus output should be displayed here if those columns are enabled in your existing output. AMRIVA does not invent sequence evidence.")
        with gr.Tab("📋 Sample Workspace"):
            gr.Markdown("## 📋 Sample Workspace")
            gr.Markdown("A compact view of the analysed samples and their linked AST observations.")
            workspace=gr.Dataframe(headers=["Sample ID","Organism","Antibiotic","AST","Genomic finding"],interactive=False,wrap=True)
            def make_workspace(g,a):
                g=g.copy() if isinstance(g,pd.DataFrame) else pd.DataFrame(g); a=a.copy() if isinstance(a,pd.DataFrame) else pd.DataFrame(a)
                if a.empty: return pd.DataFrame(columns=["Sample ID","Organism","Antibiotic","AST","Genomic finding"])
                if g.empty: g=pd.DataFrame(columns=["Sample ID","AMR Finding"])
                for c in ["Sample ID","AMR Finding"]:
                    if c not in g.columns:g[c]=""
                return a.merge(g[["Sample ID","AMR Finding"]].drop_duplicates("Sample ID"),on="Sample ID",how="left").rename(columns={"AMR Finding":"Genomic finding"})[["Sample ID","Organism","Antibiotic","AST","Genomic finding"]]
            wb=gr.Button("📋 Build workspace"); wb.click(make_workspace,[genomic_state,ast_table],workspace)
        with gr.Tab("📚 Methodology"):
            gr.Markdown(METHOD_SIMPLE)
        with gr.Tab("⚠️ Limitations"):
            gr.Markdown(LIMITS_SIMPLE)
    gr.Markdown("---\n**AMRIVA** · Research/educational prototype · Not a clinical diagnostic or treatment recommendation system")

if __name__ == "__main__":
    demo.launch(server_name="0.0.0.0",server_port=int(os.environ.get("PORT",7860)))
