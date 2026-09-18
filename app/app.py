import os
import zipfile
import tempfile
import shutil
import subprocess

import pandas as pd
import gradio as gr


# ============================================================
# ORGANISM / AMRFINDERPLUS MAPPING
# ============================================================

ORGANISM_MAP = {
    "escherichia coli": "Escherichia",
    "e. coli": "Escherichia",
    "klebsiella pneumoniae": "Klebsiella_pneumoniae",
    "pseudomonas aeruginosa": "Pseudomonas_aeruginosa",
    "acinetobacter baumannii": "Acinetobacter_baumannii",
    "staphylococcus aureus": "Staphylococcus_aureus",
    "staphylococcus epidermidis": "Staphylococcus_epidermidis",
    "enterococcus faecalis": "Enterococcus_faecalis",
    "enterococcus faecium": "Enterococcus_faecium",
    "salmonella": "Salmonella",
    "salmonella spp.": "Salmonella",
    "haemophilus influenzae": "Haemophilus_influenzae",
    "neisseria gonorrhoeae": "Neisseria_gonorrhoeae",
    "streptococcus agalactiae": "Streptococcus_agalactiae",
    "streptococcus pneumoniae": "Streptococcus_pneumoniae",
    "streptococcus pyogenes": "Streptococcus_pyogenes",
    "campylobacter": "Campylobacter",
    "serratia marcescens": "Serratia_marcescens",
    "vibrio cholerae": "Vibrio_cholerae",
    "vibrio parahaemolyticus": "Vibrio_parahaemolyticus",
    "vibrio vulnificus": "Vibrio_vulnificus",
}


# ============================================================
# FASTA READER
# ============================================================

def read_fasta(filepath):
    """
    Read the first FASTA header and all nucleotide sequence lines.
    """

    header = ""
    sequence_parts = []

    with open(filepath, "r", encoding="utf-8") as file:

        for line in file:

            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):

                if not header:
                    header = line[1:].strip()

            else:
                sequence_parts.append(line.upper())

    sequence = "".join(sequence_parts)

    return header, sequence


# ============================================================
# ORGANISM INFORMATION
# ============================================================

def detect_organism(header, filename):

    text = f"{header} {filename}".lower()

    # Longest names first to avoid partial matching problems.
    sorted_names = sorted(
        ORGANISM_MAP.keys(),
        key=len,
        reverse=True
    )

    for organism_name in sorted_names:

        if organism_name in text:

            return (
                organism_name.title(),
                ORGANISM_MAP[organism_name]
            )

    return (
        "Organism identification pending",
        None
    )


# ============================================================
# SEQUENCE STATISTICS
# ============================================================

def calculate_statistics(sequence):

    allowed = set("ACGTN")

    # Remove whitespace and keep only valid nucleotide symbols.
    sequence = "".join(
        base
        for base in sequence.upper()
        if base in allowed
    )

    if not sequence:

        raise ValueError(
            "No valid nucleotide sequence was found."
        )

    length = len(sequence)

    a_count = sequence.count("A")
    g_count = sequence.count("G")
    c_count = sequence.count("C")
    t_count = sequence.count("T")
    n_count = sequence.count("N")

    gc_content = (
        ((g_count + c_count) / length) * 100
    )

    return {
        "sequence": sequence,
        "length": length,
        "A": a_count,
        "G": g_count,
        "C": c_count,
        "T": t_count,
        "N": n_count,
        "GC": gc_content
    }


# ============================================================
# AMRFINDERPLUS
# ============================================================

def run_amrfinder(
    fasta_file,
    amrfinder_organism=None
):

    command = [
        "amrfinder",
        "-n",
        fasta_file
    ]

    # Only provide -O when we have a supported organism.
    if amrfinder_organism:

        command.extend([
            "-O",
            amrfinder_organism
        ])

    try:

        result = subprocess.run(
            command,
            capture_output=True,
            text=True,
            timeout=300
        )

        if result.returncode != 0:

            error_message = result.stderr.strip()

            if not error_message:
                error_message = (
                    "Unknown AMRFinderPlus error."
                )

            return (
                False,
                "AMRFinderPlus analysis failed.\n\n"
                + error_message
            )

        output = result.stdout.strip()

        if not output:

            return (
                True,
                "No known AMR determinant detected."
            )

        return (
            True,
            output
        )

    except subprocess.TimeoutExpired:

        return (
            False,
            "AMRFinderPlus analysis failed: "
            "analysis timed out."
        )

    except FileNotFoundError:

        return (
            False,
            "AMRFinderPlus was not found on the server."
        )

    except Exception as error:

        return (
            False,
            f"AMRFinderPlus analysis failed: {error}"
        )


# ============================================================
# AMR INTERPRETATION
# ============================================================

def interpret_amr(
    success,
    amr_result
):

    if not success:

        return (
            "AMR analysis could not be completed",

            "Susceptibility cannot be inferred because "
            "the genomic AMR analysis did not complete. "
            "Please review the analysis error and use "
            "validated phenotypic AST for clinical confirmation."
        )

    if amr_result.startswith(
        "No known AMR determinant detected"
    ):

        return (
            "No known AMR determinant detected",

            "No known genomic AMR determinant was detected "
            "by this screening. This does not prove "
            "susceptibility. Phenotypic AST is required."
        )

    return (
        "AMR determinant(s) detected",

        "Detected genomic determinants may be associated "
        "with antimicrobial resistance. Absence of a "
        "detected resistance determinant does not prove "
        "susceptibility. Phenotypic AST is required for "
        "clinical confirmation."
    )


# ============================================================
# SINGLE SAMPLE ANALYSIS
# ============================================================

def analyse_single(fasta_file):

    # Exactly 15 outputs are returned by this function.
    empty_outputs = (
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        ""
    )

    if not fasta_file:

        return (
            "Please upload a FASTA file.",
            *empty_outputs
        )

    try:

        header, sequence = read_fasta(
            fasta_file
        )

        if not sequence:

            raise ValueError(
                "The FASTA file does not contain "
                "a nucleotide sequence."
            )

        filename = os.path.basename(
            fasta_file
        )

        organism, amrfinder_organism = detect_organism(
            header,
            filename
        )

        stats = calculate_statistics(
            sequence
        )

        success, amr_result = run_amrfinder(
            fasta_file,
            amrfinder_organism
        )

        status, interpretation = interpret_amr(
            success,
            amr_result
        )

        if header:

            sample_id = header.split()[0]

        else:

            sample_id = os.path.splitext(
                filename
            )[0]

        return (
            "Sample analysed successfully.",
            sample_id,
            organism,
            filename,
            str(stats["length"]),
            str(stats["A"]),
            str(stats["G"]),
            str(stats["C"]),
            str(stats["T"]),
            str(stats["N"]),
            f'{stats["GC"]:.2f}%',
            status,
            amr_result,
            interpretation,
            header
        )

    except Exception as error:

        return (
            f"Analysis error: {error}",
            "",
            "",
            "",
            "",
            "",
            "",
            "",
            "",
            "",
            "",
            "Analysis failed",
            "",
            "",
            ""
        )


# ============================================================
# ZIP EXTRACTION
# ============================================================

def extract_fasta_files(
    zip_path,
    extraction_folder
):

    fasta_extensions = (
        ".fa",
        ".fasta",
        ".fna"
    )

    extracted_files = []

    with zipfile.ZipFile(
        zip_path,
        "r"
    ) as archive:

        for member in archive.infolist():

            if member.is_dir():
                continue

            original_name = member.filename

            filename = os.path.basename(
                original_name
            )

            if not filename:
                continue

            if filename.startswith("."):
                continue

            if not filename.lower().endswith(
                fasta_extensions
            ):
                continue

            # Use basename only.
            # This prevents ZIP path traversal.
            safe_name = filename

            destination = os.path.join(
                extraction_folder,
                safe_name
            )

            # Handle duplicate filenames.
            base, extension = os.path.splitext(
                safe_name
            )

            counter = 1

            while os.path.exists(destination):

                safe_name = (
                    f"{base}_{counter}{extension}"
                )

                destination = os.path.join(
                    extraction_folder,
                    safe_name
                )

                counter += 1

            with archive.open(
                member
            ) as source:

                with open(
                    destination,
                    "wb"
                ) as target:

                    shutil.copyfileobj(
                        source,
                        target
                    )

            extracted_files.append(
                destination
            )

    return extracted_files


# ============================================================
# MULTIPLE SAMPLE ZIP ANALYSIS
# ============================================================

def analyse_zip(zip_file):

    columns = [
        "Sample ID",
        "Detected organism",
        "FASTA file",
        "Sequence length",
        "GC content",
        "AMR status",
        "AMRFinderPlus result",
        "Interpretation"
    ]

    empty_result = pd.DataFrame(
        columns=columns
    )

    if not zip_file:

        return empty_result

    temporary_folder = tempfile.mkdtemp(
        prefix="amr_detect_"
    )

    try:

        zip_path = zip_file

        # ----------------------------------------------------
        # CHECK ZIP
        # ----------------------------------------------------

        if not zipfile.is_zipfile(
            zip_path
        ):

            return pd.DataFrame(
                [{
                    "Sample ID": "ERROR",
                    "Detected organism": "N/A",
                    "FASTA file":
                        os.path.basename(zip_path),
                    "Sequence length": "N/A",
                    "GC content": "N/A",
                    "AMR status": "Invalid ZIP",
                    "AMRFinderPlus result":
                        "The uploaded file is not a valid ZIP archive.",
                    "Interpretation":
                        "Please upload a ZIP containing FASTA files."
                }],
                columns=columns
            )

        # ----------------------------------------------------
        # EXTRACT EVERY FASTA FILE
        # ----------------------------------------------------

        fasta_files = extract_fasta_files(
            zip_path,
            temporary_folder
        )

        # ----------------------------------------------------
        # CHECK WHETHER FASTA FILES EXIST
        # ----------------------------------------------------

        if not fasta_files:

            return pd.DataFrame(
                [{
                    "Sample ID": "ERROR",
                    "Detected organism": "N/A",
                    "FASTA file":
                        "No FASTA files found",
                    "Sequence length": "N/A",
                    "GC content": "N/A",
                    "AMR status": "No FASTA files",
                    "AMRFinderPlus result":
                        "No .fa, .fasta or .fna files "
                        "were found inside the ZIP.",
                    "Interpretation":
                        "Upload a ZIP containing one or "
                        "more FASTA files."
                }],
                columns=columns
            )

        results = []

        # ----------------------------------------------------
        # ANALYSE EVERY FASTA FILE
        # ----------------------------------------------------

        for index, fasta_path in enumerate(
            fasta_files,
            start=1
        ):

            filename = os.path.basename(
                fasta_path
            )

            try:

                header, sequence = read_fasta(
                    fasta_path
                )

                if not sequence:

                    results.append({
                        "Sample ID":
                            f"AMR-BATCH-{index:03d}",

                        "Detected organism":
                            "Unknown",

                        "FASTA file":
                            filename,

                        "Sequence length":
                            "0",

                        "GC content":
                            "N/A",

                        "AMR status":
                            "Invalid FASTA",

                        "AMRFinderPlus result":
                            "No nucleotide sequence found.",

                        "Interpretation":
                            "Analysis could not be performed."
                    })

                    continue

                organism, amrfinder_organism = detect_organism(
                    header,
                    filename
                )

                stats = calculate_statistics(
                    sequence
                )

                success, amr_result = run_amrfinder(
                    fasta_path,
                    amrfinder_organism
                )

                status, interpretation = interpret_amr(
                    success,
                    amr_result
                )

                if header:

                    sample_id = header.split()[0]

                else:

                    sample_id = os.path.splitext(
                        filename
                    )[0]

                results.append({

                    "Sample ID":
                        sample_id,

                    "Detected organism":
                        organism,

                    "FASTA file":
                        filename,

                    "Sequence length":
                        stats["length"],

                    "GC content":
                        f'{stats["GC"]:.2f}%',

                    "AMR status":
                        status,

                    "AMRFinderPlus result":
                        amr_result,

                    "Interpretation":
                        interpretation
                })

            except Exception as error:

                results.append({

                    "Sample ID":
                        f"AMR-BATCH-{index:03d}",

                    "Detected organism":
                        "Unknown",

                    "FASTA file":
                        filename,

                    "Sequence length":
                        "Error",

                    "GC content":
                        "Error",

                    "AMR status":
                        "Analysis failed",

                    "AMRFinderPlus result":
                        str(error),

                    "Interpretation":
                        "Analysis could not be completed."
                })

        return pd.DataFrame(
            results,
            columns=columns
        )

    except Exception as error:

        return pd.DataFrame(
            [{
                "Sample ID": "ERROR",
                "Detected organism": "N/A",
                "FASTA file": "ZIP processing",
                "Sequence length": "N/A",
                "GC content": "N/A",
                "AMR status":
                    "Batch analysis failed",
                "AMRFinderPlus result":
                    str(error),
                "Interpretation":
                    "The ZIP could not be processed."
            }],
            columns=columns
        )

    finally:

        shutil.rmtree(
            temporary_folder,
            ignore_errors=True
        )


# ============================================================
# WEBSITE
# ============================================================

with gr.Blocks(
    title="AMR-Detect"
) as demo:

    gr.Markdown(
        """
        # 🧬 AMR-Detect

        ### Antimicrobial Resistance Detection and Analysis

        Genomic sequence screening and AMR analysis support.
        """
    )

    # ========================================================
    # SINGLE SAMPLE TAB
    # ========================================================

    with gr.Tab("Single Sample"):

        gr.Markdown(
            """
            ### 🧬 Single Sample Analysis

            Upload one FASTA sequence for analysis.
            """
        )

        single_file = gr.File(
            label="Upload FASTA sequence",
            file_types=[
                ".fa",
                ".fasta",
                ".fna"
            ],
            type="filepath"
        )

        single_button = gr.Button(
            "🧬 Analyse Sample",
            variant="primary"
        )

        single_summary = gr.Textbox(
            label="Analysis summary"
        )

        with gr.Row():

            single_id = gr.Textbox(
                label="Sample ID"
            )

            single_organism = gr.Textbox(
                label="Detected organism"
            )

        single_filename = gr.Textbox(
            label="FASTA file"
        )

        with gr.Row():

            single_length = gr.Textbox(
                label="Sequence length"
            )

            single_gc = gr.Textbox(
                label="GC content"
            )

            single_status = gr.Textbox(
                label="AMR status"
            )

        with gr.Row():

            single_a = gr.Textbox(
                label="A"
            )

            single_g = gr.Textbox(
                label="G"
            )

            single_c = gr.Textbox(
                label="C"
            )

            single_t = gr.Textbox(
                label="T"
            )

            single_n = gr.Textbox(
                label="N"
            )

        single_amr = gr.Textbox(
            label="AMRFinderPlus Result",
            lines=12
        )

        single_interpretation = gr.Textbox(
            label="Susceptibility interpretation",
            lines=5
        )

        single_header = gr.Textbox(
            label="FASTA header",
            visible=False
        )

        single_button.click(
            fn=analyse_single,
            inputs=single_file,
            outputs=[
                single_summary,
                single_id,
                single_organism,
                single_filename,
                single_length,
                single_a,
                single_g,
                single_c,
                single_t,
                single_n,
                single_gc,
                single_status,
                single_amr,
                single_interpretation,
                single_header
            ]
        )

    # ========================================================
    # MULTIPLE SAMPLE TAB
    # ========================================================

    with gr.Tab("Multiple Samples"):

        gr.Markdown(
            """
            ### 🔬 Multiple Sample Analysis

            Upload *one ZIP file* containing any number
            of FASTA files.

            Supported FASTA formats:

            .fa • .fasta • .fna

            You do not need to select the FASTA files
            individually. AMR-Detect automatically finds
            every supported FASTA file inside the ZIP.
            """
        )

        multiple_zip = gr.File(
            label="📦 Upload ZIP containing FASTA files",
            file_types=[
                ".zip"
            ],
            type="filepath"
        )

        multiple_button = gr.Button(
            "🔬 Analyse All Samples",
            variant="primary"
        )

        multiple_results = gr.Dataframe(
            headers=[
                "Sample ID",
                "Detected organism",
                "FASTA file",
                "Sequence length",
                "GC content",
                "AMR status",
                "AMRFinderPlus result",
                "Interpretation"
            ],
            datatype=[
                "str",
                "str",
                "str",
                "number",
                "str",
                "str",
                "str",
                "str"
            ],
            label="Multiple Sample AMR Results",
            interactive=False,
            wrap=True
        )

        multiple_button.click(
            fn=analyse_zip,
            inputs=multiple_zip,
            outputs=multiple_results
        )

    # ========================================================
    # PROFESSIONAL REVIEW
    # ========================================================

    gr.Markdown(
        """
        ---

        ### 🩺 Professional Review

        AMR-Detect provides genomic screening and
        research-support information.

        Genomic AMR findings do not independently confirm
        clinical antimicrobial susceptibility or determine
        treatment.

        Phenotypic antimicrobial susceptibility testing (AST),
        clinical context and qualified healthcare-professional
        interpretation are required for clinical decisions.
        """
    )


# ============================================================
# SERVER
# ============================================================

port = int(
    os.environ.get(
        "PORT",
        "7860"
    )
)

demo.launch(
    server_name="0.0.0.0",
    server_port=port
)
               
