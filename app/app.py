import os, zipfile, tempfile, shutil, subprocess
import pandas as pd
import gradio as gr

APP_NAME = "AMRIVA"

ORGANISM_MAP = {
    "escherichia coli":"Escherichia coli", "e. coli":"Escherichia coli", "escherichia_coli":"Escherichia coli", "ecoli":"Escherichia coli", "e_coli":"Escherichia coli",

    "klebsiella pneumoniae":"Klebsiella pneumoniae", "klebsiella_pneumoniae":"Klebsiella pneumoniae", "klebsiella":"Klebsiella pneumoniae",
    "pseudomonas aeruginosa":"Pseudomonas aeruginosa", "pseudomonas_aeruginosa":"Pseudomonas aeruginosa", "pseudomonas":"Pseudomonas aeruginosa",
    "acinetobacter baumannii":"Acinetobacter baumannii", "acinetobacter_baumannii":"Acinetobacter baumannii",
    "staphylococcus aureus":"Staphylococcus aureus", "staphylococcus_aureus":"Staphylococcus aureus",
    "staphylococcus epidermidis":"Staphylococcus_epidermidis",
    "enterococcus faecalis":"Enterococcus faecalis", "enterococcus_faecalis":"Enterococcus faecalis",
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
    """Resolve organism from explicit FASTA metadata/header or filename.

    AMRFinderPlus is an AMR-determinant detector, not a species-identification
    engine. Therefore this field is populated only from supplied FASTA metadata
    (header/filename); it is never invented from the AMR gene result.
    """
    text=f"{header} {filename}".lower().replace("-", "_")
    for name in sorted(ORGANISM_MAP, key=len, reverse=True):
        if name in text:
            return ORGANISM_MAP[name], ORGANISM_MAP[name]
    return "Identification not available from FASTA metadata", None


def calculate_statistics(sequence):
    seq="".join(x for x in sequence.upper() if x in set("ACGTN"))
    if not seq: raise ValueError("No valid nucleotide sequence was found.")
    n=len(seq); a=seq.count("A"); g=seq.count("G"); c=seq.count("C"); t=seq.count("T"); unknown=seq.count("N")
    return {"sequence":seq,"length":n,"A":a,"G":g,"C":c,"T":t,"N":unknown,"GC":(g+c)/n*100}


AMRFINDER_DISPLAY_COLUMNS = [
    "Sample ID", "AMR Gene", "Protein / Function", "Resistance Class",
    "Identity (%)", "Coverage (%)"
]


def _amrfinder_dataframe(sample_id, raw_result):
    """Convert AMRFinderPlus TSV output into the compact AMRIVA table.

    AMRIVA keeps the detailed AMRFinderPlus evidence separate from the batch
    summary. When AMRFinderPlus reports no hit, a visible placeholder row is
    returned instead of leaving the table completely blank.
    """
    sid = str(sample_id or "").strip()
    text = str(raw_result or "").strip()
    empty_cols = AMRFINDER_DISPLAY_COLUMNS

    def no_hit():
        return pd.DataFrame([{
            "Sample ID": sid,
            "AMR Gene": "—",
            "Protein / Function": "—",
            "Resistance Class": "—",
            "Identity (%)": "—",
            "Coverage (%)": "—",
        }], columns=empty_cols)

    if not text:
        return no_hit()
    if text.startswith("No known AMR determinant detected"):
        return no_hit()
    if text.startswith("AMRFinderPlus analysis failed") or text.startswith("AMRFinderPlus was not found"):
        return no_hit()

    try:
        from io import StringIO
        df = pd.read_csv(StringIO(text), sep="\t", dtype=str, keep_default_na=False)
        if df.empty:
            return no_hit()

        # AMRFinderPlus standard TSV names.
        aliases = {
            "Gene symbol": "AMR Gene",
            "Element symbol": "AMR Gene",
            "Sequence name": "Protein / Function",
            "Element name": "Protein / Function",
            "HMM description": "Protein / Function",
            "HMM description": "Protein / Function",
            "Class": "Resistance Class",
            "% Identity to reference sequence": "Identity (%)",
            "% Coverage of reference sequence": "Coverage (%)",
            "% Identity": "Identity (%)",
            "% Coverage": "Coverage (%)",
        }
        df = df.rename(columns={str(c).strip(): aliases.get(str(c).strip(), str(c).strip()) for c in df.columns})

        # If AMRFinderPlus has a different version's naming, use sensible
        # fallbacks rather than dropping the valid hit.
        if "AMR Gene" not in df.columns:
            for c in df.columns:
                lc = str(c).lower()
                if "gene" in lc and "symbol" in lc:
                    df["AMR Gene"] = df[c]
                    break
        if "Protein / Function" not in df.columns:
            for c in df.columns:
                lc = str(c).lower()
                if "element name" in lc or "sequence name" in lc or "description" in lc or "protein" in lc or "function" in lc:
                    df["Protein / Function"] = df[c]
                    break
        if "Resistance Class" not in df.columns and "Class" in df.columns:
            df["Resistance Class"] = df["Class"]

        for c in empty_cols:
            if c not in df.columns:
                df[c] = ""

        out = df[empty_cols].copy()
        out.insert(0, "Sample ID", sid) if "Sample ID" not in out.columns else None
        # Replace genuinely empty cells with an em dash for readable results.
        for c in empty_cols[1:]:
            out[c] = out[c].astype(str).replace({"": "—", "nan": "—", "None": "—"})
        return out[empty_cols]

    except Exception:
        # Headerless fallback for older/variant AMRFinderPlus output.
        try:
            from io import StringIO
            raw = pd.read_csv(StringIO(text), sep="\t", header=None, dtype=str, keep_default_na=False)
            if raw.empty:
                return no_hit()
            # Standard AMRFinderPlus positional layout has Gene symbol around
            # column 5 and Sequence name around column 6. Locate Class and
            # identity/coverage by their standard positions when available.
            def val(i):
                return raw.iloc[:, i].astype(str) if i < raw.shape[1] else pd.Series(["—"] * len(raw))
            out = pd.DataFrame({
                "Sample ID": [sid] * len(raw),
                "AMR Gene": val(5),
                "Protein / Function": val(6),
                "Resistance Class": val(10),
                "Identity (%)": val(16),
                "Coverage (%)": val(15),
            })
            for c in empty_cols[1:]:
                out[c] = out[c].replace({"": "—", "nan": "—", "None": "—"})
            return out[empty_cols]
        except Exception:
            return no_hit()

def _amrfinder_determinants(parsed):
    """Return a concise determinant summary for the AMR status column."""
    if parsed is None or parsed.empty:
        return ""
    vals=[]
    for _, row in parsed.iterrows():
        gene=str(row.get("AMR Gene", row.get("Gene symbol", ""))).strip()
        seq=str(row.get("Sequence name", "")).strip()
        value=gene or seq
        if value and value not in vals:
            vals.append(value)
    return ", ".join(vals)


def run_amrfinder(fasta_file, amrfinder_organism=None):
    # AMRFinderPlus nucleotide mode is authoritative for the AMR result.
    # Do not pass our metadata-derived organism string to -O because that
    # option has its own controlled organism vocabulary and is not needed
    # for ordinary nucleotide AMR screening.
    cmd=["amrfinder","-n",fasta_file]
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
                else:
                    detail_frames.append(pd.DataFrame([{
                        "Sample ID": sid,
                        "AMR Gene": "—",
                        "Protein / Function": "—",
                        "Resistance Class": "—",
                        "Identity (%)": "—",
                        "Coverage (%)": "—",
                    }], columns=AMRFINDER_DISPLAY_COLUMNS))
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

with gr.Blocks(title="AMRIVA",css=CSS,theme=gr.themes.Soft()) as demo:
    genomic_state = gr.State(pd.DataFrame(columns=["Sample ID","Detected organism","AMR status","AMRFinderPlus result"]))
    gr.HTML("<div id='hero'><h1>🧬 AMRIVA</h1><p>Antimicrobial Resistance Genomic Analysis Platform</p><p>Genomic screening • AST/MIC integration • Comparative analysis</p></div>")
    with gr.Tabs():
        with gr.Tab("🏠 Home"): gr.Markdown(HOME)
        with gr.Tab("🧬 Single Sample"):
            gr.Markdown("## Single-Sample Analysis\nUpload one FASTA sequence for genomic AMR screening.")
            sf=gr.File(label="Upload FASTA",file_types=[".fa",".fasta",".fna"],type="filepath")
            sb=gr.Button("🧬 Analyse Sample",variant="primary"); ss=gr.Textbox(label="Analysis summary")
            with gr.Row(): sid=gr.Textbox(label="Sample ID"); sorg=gr.Textbox(label="Detected organism (from FASTA metadata)")
            sfn=gr.Textbox(label="FASTA file")
            with gr.Row(): sl=gr.Textbox(label="Sequence length"); sgc=gr.Textbox(label="GC content"); sst=gr.Textbox(label="AMR status")
            with gr.Row(): sa=gr.Textbox(label="A"); sg=gr.Textbox(label="G"); sc=gr.Textbox(label="C"); st=gr.Textbox(label="T"); sn=gr.Textbox(label="N")
            samr=gr.Dataframe(headers=AMRFINDER_DISPLAY_COLUMNS,interactive=False,wrap=True,label="AMRFinderPlus Result — detailed table (actual AMRFinderPlus output)"); sint=gr.Textbox(label="Genomic Interpretation",lines=6); sh=gr.Textbox(label="FASTA Header",visible=False)
            sb.click(analyse_single_with_state,sf,[ss,sid,sorg,sfn,sl,sa,sg,sc,st,sn,sgc,sst,samr,sint,sh,genomic_state])
        with gr.Tab("📁 Multiple Samples"):
            gr.Markdown("## Batch Analysis\nUpload **one ZIP containing any number of FASTA files**. AMRIVA analyses every `.fa`, `.fasta` and `.fna` file automatically.")
            bz=gr.File(label="Upload ZIP containing FASTA files",file_types=[".zip"],type="filepath"); bb=gr.Button("📊 Analyse All Samples",variant="primary")
            bt=gr.Dataframe(headers=["Sample ID","Detected organism","FASTA file","Sequence length","GC content","AMR status","AMRFinderPlus result","Interpretation"],interactive=False,wrap=True, label="Batch genomic summary")
            bdetail=gr.Dataframe(headers=AMRFINDER_DISPLAY_COLUMNS,interactive=False,wrap=True,label="AMRFinderPlus Result — detailed table (actual AMRFinderPlus output, separate from AMR status)")
            bb.click(analyse_zip_with_state,bz,[bt,genomic_state,bdetail])
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
                ast_ph=gr.Textbox(label="pH",placeholder="e.g. 7.0")
                ast_temperature=gr.Textbox(label="Temperature (°C)",placeholder="e.g. 37")
                ast_notes=gr.Textbox(label="Experimental notes",placeholder="Optional")
            gr.Markdown("**AMR Finding:** AMRIVA will automatically pull the genomic AMR status for this Sample ID when the same Sample ID has already been analysed in the genomic section. **pH and temperature** are entered from laboratory/research observations.")
            ast_add=gr.Button("＋ Add AST/MIC Result",variant="primary")
            ast_message=gr.Markdown()
            ast_table=gr.Dataframe(
                headers=["Sample ID","AMR Finding","Antibiotic","MIC","MIC unit","AST","pH","Temperature","Comparison","Other information"],
                value=[],
                datatype=["str"]*10,
                interactive=True,
                wrap=True,
                column_widths=[120,170,150,90,90,120,70,100,130,180],
                label="AST/MIC results — add as many samples/antibiotics as needed"
            )
            ast_add.click(
                add_ast_result,
                [ast_table,genomic_state,ast_sample,ast_antibiotic,ast_mic,ast_unit,ast_result,ast_ph,ast_temperature,ast_notes],
                [ast_table,ast_message]
            )

            gr.Markdown("### 2️⃣ Multiple-sample upload or spreadsheet paste")
            gr.Markdown("For many laboratory results, you can **upload CSV/Excel** or **copy-paste rows directly from Excel/Google Sheets**. One Sample ID can appear in unlimited antibiotic rows.")
            with gr.Row():
                with gr.Column(scale=1):
                    ast_csv=gr.File(label="📁 Upload AST/MIC CSV or Excel",file_types=[".csv",".xlsx",".xls"],type="filepath")
                    ast_csv_button=gr.Button("📥 Load Uploaded File",variant="primary")
                    ast_paste=gr.Textbox(label="📋 Paste AST/MIC rows from Excel",placeholder="Copy cells from Excel and paste here. Header row is recommended.",lines=7)
                    ast_paste_button=gr.Button("📋 Load Pasted Rows")
                with gr.Column(scale=2):
                    ast_upload_preview=gr.Dataframe(headers=["Sample ID","AMR Finding","Antibiotic","MIC","MIC unit","AST","pH","Temperature","Comparison","Other information"],interactive=False,wrap=True,label="Uploaded / pasted result preview")
            ast_csv_message=gr.Markdown()
            ast_paste_message=gr.Markdown()
            ast_csv_button.click(load_ast_csv,ast_csv,[ast_upload_preview,ast_csv_message])
            ast_csv_button.click(load_ast_csv,ast_csv,[ast_table,ast_csv_message])
            ast_paste_button.click(parse_pasted_ast_data,ast_paste,[ast_upload_preview,ast_paste_message])
            ast_paste_button.click(parse_pasted_ast_data,ast_paste,[ast_table,ast_paste_message])

            gr.Markdown("### 3️⃣ MIC concentration series")
            gr.Markdown("Enter experimentally observed **Growth / No growth** at each tested concentration. You can type into the table, paste spreadsheet rows, or upload CSV/Excel. AMRIVA reports the lowest tested concentration recorded as No growth and does not independently label Susceptible/Intermediate/Resistant.")
            with gr.Row():
                with gr.Column(scale=1):
                    mic_file=gr.File(label="📁 Upload MIC series CSV or Excel",file_types=[".csv",".xlsx",".xls"],type="filepath")
                    mic_file_button=gr.Button("📥 Load MIC File",variant="primary")
                    mic_paste=gr.Textbox(label="📋 Paste MIC series from Excel",placeholder="Sample ID<TAB>Antibiotic<TAB>Concentration<TAB>Unit<TAB>Growth",lines=7)
                    mic_paste_button=gr.Button("📋 Load Pasted MIC Rows")
                with gr.Column(scale=2):
                    mic_series=gr.Dataframe(headers=MIC_SERIES_COLUMNS,value=[],datatype=["str"]*5,interactive=True,wrap=True,column_widths=[130,170,120,90,130],label="MIC observations — supports multiple samples")

            mic_file_message=gr.Markdown()
            mic_paste_message=gr.Markdown()
            mic_calc=gr.Button("🔬 Calculate MIC from Observations",variant="primary")
            mic_message=gr.Markdown()
            mic_results=gr.Dataframe(headers=["Sample ID","Antibiotic","MIC","Unit","Lowest no-growth concentration","Interpretation note"],interactive=False,wrap=True,label="Calculated MIC results")

            def _normalise_mic_dataframe(df):
                """Normalise an uploaded/pasted MIC-series table to the five required columns."""
                if df is None:
                    raise ValueError("No MIC data was supplied.")

                df = df.copy()
                df.columns = [str(c).strip() for c in df.columns]

                aliases = {
                    "sample id":"Sample ID", "sample_id":"Sample ID", "sample_id ":"Sample ID",
                    "sample":"Sample ID", "sample name":"Sample ID",
                    "antibiotic":"Antibiotic", "drug":"Antibiotic", "antimicrobial":"Antibiotic",
                    "concentration":"Concentration", "conc":"Concentration",
                    "concentration value":"Concentration", "mic concentration":"Concentration",
                    "mic":"Concentration", "mic value":"Concentration",
                    "unit":"Unit", "mic unit":"Unit", "mic_unit":"Unit",
                    "growth":"Growth", "growth/no growth":"Growth", "growth / no growth":"Growth",
                    "growth status":"Growth", "observation":"Growth", "result":"Growth"
                }
                renamed = {}
                for c in df.columns:
                    key = str(c).strip().lower()
                    renamed[c] = aliases.get(key, c)
                df = df.rename(columns=renamed)

                required = ["Sample ID", "Antibiotic", "Concentration", "Unit", "Growth"]

                # Headerless five-column data is accepted in the documented order.
                if not set(required).issubset(df.columns):
                    if len(df.columns) == 5:
                        df.columns = required
                    else:
                        missing = [c for c in required if c not in df.columns]
                        raise ValueError(
                            "This is not a MIC concentration-series file. "
                            "Required columns are: Sample ID, Antibiotic, Concentration, Unit, Growth. "
                            f"Missing: {', '.join(missing)}."
                        )

                df = df[required].copy().fillna("")
                # Remove completely empty rows.
                mask = df.astype(str).apply(lambda row: any(v.strip() for v in row), axis=1)
                df = df.loc[mask].reset_index(drop=True)
                if df.empty:
                    raise ValueError("The MIC file contains no data rows.")
                return df

            def load_mic_file(path):
                if not path:
                    return pd.DataFrame(columns=MIC_SERIES_COLUMNS), "No MIC file selected."
                try:
                    # Gradio may return a string path or a one-item list in some versions.
                    if isinstance(path, (list, tuple)):
                        path = path[0] if path else None
                    if not path:
                        return pd.DataFrame(columns=MIC_SERIES_COLUMNS), "No MIC file selected."

                    ext = os.path.splitext(str(path))[1].lower()
                    if ext == ".xlsx":
                        df = pd.read_excel(path, engine="openpyxl")
                    elif ext == ".xls":
                        try:
                            df = pd.read_excel(path, engine="xlrd")
                        except ImportError:
                            return pd.DataFrame(columns=MIC_SERIES_COLUMNS), "Old .xls files need xlrd. Please use .xlsx or CSV."
                    elif ext == ".csv":
                        df = pd.read_csv(path)
                    else:
                        return pd.DataFrame(columns=MIC_SERIES_COLUMNS), "Unsupported file type. Please upload CSV, XLSX or XLS."

                    df = _normalise_mic_dataframe(df)
                    return df, f"✅ Loaded {len(df)} MIC observation row(s) from {ext.upper().replace('.', '')}."
                except Exception as e:
                    return pd.DataFrame(columns=MIC_SERIES_COLUMNS), f"❌ MIC file could not be loaded: {e}"

            def parse_pasted_mic(text):
                if not text or not str(text).strip():
                    return pd.DataFrame(columns=MIC_SERIES_COLUMNS), "Nothing was pasted."
                try:
                    from io import StringIO
                    raw = str(text).strip()

                    # Excel/Google Sheets clipboard data is normally tab-separated.
                    try:
                        df = pd.read_csv(StringIO(raw), sep="\t", header=0)
                    except Exception:
                        df = pd.read_csv(StringIO(raw), sep=None, engine="python", header=0)

                    try:
                        df = _normalise_mic_dataframe(df)
                    except ValueError:
                        # Also support headerless pasted rows in the exact five-column order.
                        df = pd.read_csv(StringIO(raw), sep="\t", header=None)
                        if len(df.columns) == 1:
                            df = pd.read_csv(StringIO(raw), sep=None, engine="python", header=None)
                        if len(df.columns) != 5:
                            raise
                        df.columns = MIC_SERIES_COLUMNS
                        df = _normalise_mic_dataframe(df)

                    return df, f"✅ Pasted {len(df)} MIC observation row(s)."
                except Exception as e:
                    return pd.DataFrame(columns=MIC_SERIES_COLUMNS), f"❌ Pasted MIC data could not be read: {e}"

            mic_file_button.click(load_mic_file,mic_file,[mic_series,mic_file_message])
            mic_file.change(load_mic_file,mic_file,[mic_series,mic_file_message])
            mic_paste_button.click(parse_pasted_mic,mic_paste,[mic_series,mic_paste_message])
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
            meta_csv=gr.File(label="Upload experimental metadata CSV or Excel",file_types=[".csv",".xlsx",".xls"],type="filepath")
            meta_csv_button=gr.Button("📥 Load Multiple-Sample Metadata")
            meta_csv_message=gr.Markdown()
            meta_csv_button.click(load_metadata_csv,meta_csv,[meta,meta_csv_message])
        with gr.Tab("⭐ Comparative Analysis"):
            gr.Markdown("## ⭐ Comparative Analysis")
            gr.Markdown("AMRIVA combines **genomic AMR screening + phenotypic AST/MIC + pH + temperature** using the common **Sample ID**. Add AST/MIC and metadata first, then run the comparison.")
            gr.Markdown("### Data sources")
            gr.Markdown("The genomic source comes from the FASTA/ZIP analysis in this browser session. AST/MIC and experimental metadata come from the tables in their respective tabs.")
            cb=gr.Button("⭐ Build Comparative Analysis",variant="primary")
            cmsg=gr.Markdown()
            co=gr.Dataframe(headers=["Sample ID","Detected organism","AMR Finding","AMRFinderPlus result","Antibiotic","MIC","MIC unit","AST","pH","Temperature","Other information","Comparison"],interactive=False,wrap=True,label="Integrated Genotype–Phenotype Comparison")
            cs=gr.Markdown()
            cb.click(run_comparative,[genomic_state,ast_table,meta],[co,cmsg,cs])
            gr.Markdown("**Comparison status is conservative:** 🟢 Potentially concordant means the supplied genomic evidence contains a recognizable antibiotic/class relationship and the supplied AST says Resistant. 🟠 Potentially discordant means the genomic evidence is relevant but AST says Susceptible. 🟡 Requires interpretation means the relationship cannot be established safely from the supplied data. These are not clinical diagnostic calls.")
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
            gr.Markdown("Ask educational questions about AMR, AST, MIC, genomic screening or the AMRIVA workflow. The assistant is educational and does not diagnose patients, prescribe antibiotics or replace laboratory interpretation.")
            ai_question=gr.Textbox(label="Ask a question",placeholder="e.g. What does MIC mean? Why can genotype and AST disagree?",lines=3)
            ai_answer=gr.Markdown()
            def answer_assistant(q):
                q=str(q or "").strip().lower()
                if not q: return "Please enter a question."
                if "mic" in q and "ast" in q:
                    return "### MIC vs AST\n**MIC** is the lowest tested antimicrobial concentration that prevents visible growth under the test conditions. **AST** reports the susceptibility category or observed susceptibility result. MIC values are interpreted using organism-, drug- and standard-specific breakpoints. AMRIVA records laboratory-provided results; it does not create clinical breakpoints itself."
                if "mic" in q:
                    return "### MIC\nMIC means **minimum inhibitory concentration**. In AMRIVA's concentration-series tool, you enter the experimentally observed growth/no-growth result at each tested concentration. AMRIVA reports the lowest tested concentration recorded as no growth. That value is not automatically labelled Susceptible/Intermediate/Resistant."
                if "ast" in q or "susceptible" in q or "resistant" in q:
                    return "### AST\nAntimicrobial susceptibility testing is a laboratory process used to determine how an organism responds to antimicrobial agents. AMRIVA accepts the laboratory-generated AST result and links it to the sample's genomic findings."
                if "genotype" in q or "genomic" in q or "gene" in q:
                    return "### Genomic AMR screening\nAMRIVA uses AMRFinderPlus to screen sequence data for known resistance determinants. A detected determinant can provide genomic evidence, but it does not automatically prove phenotypic resistance to every antibiotic. Phenotypic AST remains important."
                if "ph" in q or "temperature" in q:
                    return "### Experimental metadata\npH and temperature are laboratory/research measurements supplied by the user. AMRIVA does not calculate them from a FASTA sequence. They can be linked to the same Sample ID for comparative analysis."
                if "compare" in q or "concord" in q or "discord" in q:
                    return "### Comparative analysis\nAMRIVA brings together genomic findings, AST/MIC and experimental metadata using Sample ID. A potential concordance/discordance flag is only a screening-level indication; proper interpretation requires the organism, antimicrobial, breakpoint standard and laboratory context."
                if "women" in q or "uti" in q:
                    return "### Women's Health & AMR\nAMRIVA's Women's Health section is educational/research-focused, including AMR surveillance relevant to urinary and reproductive health. It does not diagnose patient vaginal, urine or blood samples."
                return "### AMRIVA Research Assistant\nI can explain AMR, genomic AMR screening, AST, MIC, pH/temperature metadata, comparative analysis and the AMRIVA workflow. Try asking: **What is MIC?**, **Why can genotype and AST disagree?**, or **What does AMRFinderPlus do?**"
            ai_question.submit(answer_assistant, ai_question, ai_answer)
            gr.Button("💡 Explain").click(answer_assistant, ai_question, ai_answer)
        with gr.Tab("👩‍🔬 Women's Health & AMR"): gr.Markdown(WOMEN)
        with gr.Tab("📚 Methodology"): gr.Markdown(METHOD)
        with gr.Tab("⚠️ Limitations"): gr.Markdown(LIMITS)
    gr.Markdown("---\n**AMRIVA** · Antimicrobial Resistance Genomic Analysis Platform\n\n*Research and educational prototype — not a clinical diagnostic system.*")

if __name__ == "__main__":
    demo.launch(server_name="0.0.0.0",server_port=int(os.environ.get("PORT",7860)))
