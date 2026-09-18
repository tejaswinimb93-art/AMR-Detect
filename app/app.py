import os
import re
import subprocess
import gradio as gr
import pandas as pd


ORGANISM_MAP = {
    "escherichia coli": "Escherichia",
    "e. coli": "Escherichia",
    "klebsiella pneumoniae": "Klebsiella_pneumoniae",
    "pseudomonas aeruginosa": "Pseudomonas_aeruginosa",
    "acinetobacter baumannii": "Acinetobacter_baumannii",
    "staphylococcus aureus": "Staphylococcus_aureus",
    "enterococcus faecalis": "Enterococcus_faecalis",
    "salmonella": "Salmonella",
    "salmonella spp.": "Salmonella",
}


def read_fasta(filepath):
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

    return header, "".join(sequence_parts)


def detect_organism_from_header(header):
    header_lower = header.lower()

    for name in ORGANISM_MAP:
        if name in header_lower:
            return name.title()

    return "Organism identification pending"


def get_amrfinder_organism(organism_name):
    organism_lower = organism_name.lower()

    for name, amrfinder_name in ORGANISM_MAP.items():
        if name in organism_lower:
            return amrfinder_name

    return None


def run_amrfinder(fasta_file, organism_name):
    command = [
        "amrfinder",
        "-n",
        fasta_file
    ]

    amrfinder_organism = get_amrfinder_organism(organism_name)

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
            error = result.stderr.strip()

            if not error:
                error = "Unknown AMRFinderPlus error."

            return (
                "AMRFinderPlus analysis failed.\n\n"
                + error
            )

        output = result.stdout.strip()

        if not output:
            return "No AMR determinants were detected."

        return output

    except subprocess.TimeoutExpired:
        return "AMRFinderPlus analysis failed: analysis timed out."

    except FileNotFoundError:
        return (
            "AMRFinderPlus analysis failed: "
            "AMRFinderPlus was not found."
        )

    except Exception as error:
        return f"AMRFinderPlus analysis failed: {error}"


def calculate_sequence_statistics(sequence):
    allowed = set("ACGTN")

    if any(base not in allowed for base in sequence):
        raise ValueError(
            "Sequence contains characters other than A, C, G, T or N."
        )

    length = len(sequence)

    a = sequence.count("A")
    c = sequence.count("C")
    g = sequence.count("G")
    t = sequence.count("T")
    n = sequence.count("N")

    gc = (
        ((g + c) / length) * 100
        if length > 0
        else 0
    )

    return length, a, g, c, t, n, gc


def interpret_amr(amr_result):

    if amr_result.startswith("No AMR determinants"):
        status = "No known AMR determinant detected"

        interpretation = (
            "No known genomic AMR determinant was detected "
            "by this screening. This does not confirm clinical "
            "susceptibility. Phenotypic AST is required."
        )

    elif amr_result.startswith("AMRFinderPlus analysis failed"):
        status = "AMR analysis could not be completed"

        interpretation = (
            "Susceptibility cannot be inferred because the "
            "genomic AMR analysis did not complete."
        )

    else:
        status = "AMR determinant(s) detected"

        interpretation = (
            "Detected genomic determinants may be associated "
            "with antimicrobial resistance. Absence of a "
            "detected determinant does not prove susceptibility. "
            "Phenotypic AST is required for clinical confirmation."
        )

    return status, interpretation


def analyse_single(fasta_file):

    if not fasta_file:
        return (
            "Please upload a FASTA file.",
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

    try:
        header, sequence = read_fasta(fasta_file)

        if not sequence:
            raise ValueError(
                "The FASTA file does not contain a nucleotide sequence."
            )

        organism = detect_organism_from_header(header)

        (
            length,
            a,
            g,
            c,
            t,
            n,
            gc
        ) = calculate_sequence_statistics(sequence)

        filename = os.path.basename(fasta_file)

        amr_result = run_amrfinder(
            fasta_file,
            organism
        )

        status, interpretation = interpret_amr(
            amr_result
        )

        sample_id = (
            header.split("_")[0]
            if header
            else "AMR-SAMPLE-001"
        )

        return (
            "Sample analysed successfully.",
            sample_id,
            organism,
            filename,
            str(length),
            str(a),
            str(g),
            str(c),
            str(t),
            str(n),
            f"{gc:.2f}%",
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


def analyse_multiple(files):

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

    if not files:
        return pd.DataFrame(columns=columns)

    results = []

    for index, fasta_file in enumerate(files, start=1):

        filename = os.path.basename(fasta_file)

        try:
            header, sequence = read_fasta(fasta_file)

            if not sequence:
                results.append({
                    "Sample ID": f"AMR-BATCH-{index:03d}",
                    "Detected organism": "Unknown",
                    "FASTA file": filename,
                    "Sequence length": 0,
                    "GC content": "N/A",
                    "AMR status": "Invalid FASTA",
                    "AMRFinderPlus result": "No nucleotide sequence found.",
                    "Interpretation": "Analysis could not be performed."
                })
                continue

            organism = detect_organism_from_header(header)

            (
                length,
                a,
                g,
                c,
                t,
                n,
                gc
            ) = calculate_sequence_statistics(sequence)

            amr_result = run_amrfinder(
                fasta_file,
                organism
            )

            status, interpretation = interpret_amr(
                amr_result
            )

            sample_id = (
                header.split("_")[0]
                if header
                else f"AMR-BATCH-{index:03d}"
            )

            results.append({
                "Sample ID": sample_id,
                "Detected organism": organism,
                "FASTA file": filename,
                "Sequence length": length,
                "GC content": f"{gc:.2f}%",
                "AMR status": status,
                "AMRFinderPlus result": amr_result,
                "Interpretation": interpretation
            })

        except Exception as error:

            results.append({
                "Sample ID": f"AMR-BATCH-{index:03d}",
                "Detected organism": "Unknown",
                "FASTA file": filename,
                "Sequence length": "Error",
                "GC content": "Error",
                "AMR status": "Analysis failed",
                "AMRFinderPlus result": str(error),
                "Interpretation": "Analysis could not be completed."
            })

    return pd.DataFrame(results, columns=columns)


with gr.Blocks(
    title="AMR-Detect"
) as demo:

    gr.Markdown(
        """
        # 🧬 AMR-Detect

        ### Antimicrobial Resistance Detection and Analysis

        Upload genomic FASTA sequences for sequence analysis and
        genomic AMR screening.
        """
    )

    with gr.Tab("Single Sample"):

        single_file = gr.File(
            label="Upload FASTA sequence",
            file_types=[".fa", ".fasta", ".fna"],
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
            single_a = gr.Textbox(label="A")
            single_g = gr.Textbox(label="G")
            single_c = gr.Textbox(label="C")
            single_t = gr.Textbox(label="T")
            single_n = gr.Textbox(label="N")

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

    with gr.Tab("Multiple Samples"):

        gr.Markdown(
            """
            ### 🔬 Batch genomic analysis

            Upload multiple FASTA files. Each sample is analysed
            independently and receives its own result.
            """
        )

        multiple_files = gr.File(
            label="Upload multiple FASTA files",
            file_types=[".fa", ".fasta", ".fna"],
            type="filepath",
            file_count="multiple"
        )

        multiple_button = gr.Button(
            "🔬 Analyse Multiple Samples",
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
            label="Multiple Sample AMR Comparison",
            interactive=False,
            wrap=True
        )

        multiple_button.click(
            fn=analyse_multiple,
            inputs=multiple_files,
            outputs=multiple_results
        )

    gr.Markdown(
        """
        ---

        ### 🩺 Professional review

        AMR-Detect provides genomic screening information and
        research-support analysis. Genomic findings do not independently
        confirm clinical antimicrobial susceptibility. Final clinical
        interpretation and treatment decisions require qualified
        healthcare-professional review together with validated
        antimicrobial susceptibility testing (AST) and clinical context.
        """
    )


port = int(os.environ.get("PORT", "7860"))

demo.launch(
    server_name="0.0.0.0",
    server_port=port
)
