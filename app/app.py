import os
import subprocess
import gradio as gr
import pandas as pd


organisms = [
    "Escherichia coli",
    "Klebsiella pneumoniae",
    "Staphylococcus aureus",
    "Pseudomonas aeruginosa",
    "Acinetobacter baumannii",
    "Enterococcus faecalis",
    "Salmonella spp.",
    "Other / Not listed"
]


organism_groups = {
    "Escherichia coli": "Escherichia",
    "Klebsiella pneumoniae": "Klebsiella_pneumoniae",
    "Staphylococcus aureus": "Staphylococcus_aureus",
    "Pseudomonas aeruginosa": "Pseudomonas_aeruginosa",
    "Acinetobacter baumannii": "Acinetobacter_baumannii",
    "Enterococcus faecalis": "Enterococcus_faecalis",
    "Salmonella spp.": "Salmonella"
}


def read_fasta(filepath):
    sequence = ""

    with open(filepath, "r", encoding="utf-8") as file:
        for line in file:
            line = line.strip()

            if not line or line.startswith(">"):
                continue

            sequence += line.upper()

    return sequence


def run_amrfinder(fasta_file, organism):
    command = ["amrfinder", "-n", fasta_file]

    organism_group = organism_groups.get(organism)

    if organism_group:
        command.extend(["-O", organism_group])

    try:
        result = subprocess.run(
            command,
            capture_output=True,
            text=True,
            timeout=300
        )

        if result.returncode != 0:
            return (
                "AMRFinderPlus analysis failed.\n\n"
                + (result.stderr.strip() or "Unknown error.")
            )

        output = result.stdout.strip()

        if not output:
            return "No AMR determinants were detected."

        return output

    except Exception as error:
        return f"AMRFinderPlus error: {error}"


def analyse_single(organism, sample_id, fasta_file):

    if not fasta_file:
        return (
            "Please upload a FASTA file.",
            "", "", "", "", "", "", "", "",
            "", "", "", ""
        )

    if not organism:
        organism = "Not specified"

    if not sample_id:
        sample_id = "Not provided"

    try:
        sequence = read_fasta(fasta_file)

        if not sequence:
            return (
                "The FASTA file does not contain a nucleotide sequence.",
                "", "", "", "", "", "", "", "",
                "", "", "", ""
            )

        allowed = set("ACGTN")

        if any(base not in allowed for base in sequence):
            return (
                "The FASTA sequence contains invalid nucleotide characters.",
                "", "", "", "", "", "", "", "",
                "", "", "", ""
            )

        length = len(sequence)

        a_count = sequence.count("A")
        c_count = sequence.count("C")
        g_count = sequence.count("G")
        t_count = sequence.count("T")
        n_count = sequence.count("N")

        gc_content = (
            ((g_count + c_count) / length) * 100
            if length > 0
            else 0
        )

        filename = os.path.basename(fasta_file)

        amr_result = run_amrfinder(
            fasta_file,
            organism
        )

        if amr_result.startswith("No AMR determinants"):
            amr_status = "No known AMR determinant detected"

            susceptibility = (
                "No known genomic AMR determinant was detected "
                "by this screening. This does not confirm clinical "
                "susceptibility. Phenotypic AST is required."
            )

        elif amr_result.startswith("AMRFinderPlus analysis failed"):
            amr_status = "AMR analysis could not be completed"

            susceptibility = (
                "Susceptibility cannot be inferred because the "
                "genomic AMR analysis did not complete."
            )

        else:
            amr_status = "AMR determinant(s) detected"

            susceptibility = (
                "Detected genomic determinants may be associated "
                "with antimicrobial resistance. Absence of a "
                "detected determinant does not prove susceptibility. "
                "Phenotypic AST is required for clinical confirmation."
            )

        return (
            "Sample received and analysed successfully.",
            sample_id,
            organism,
            filename,
            str(length),
            str(a_count),
            str(g_count),
            str(c_count),
            str(t_count),
            str(n_count),
            f"{gc_content:.2f}%",
            amr_status,
            amr_result,
            susceptibility
        )

    except Exception as error:
        return (
            f"Error: {error}",
            "", "", "", "", "", "", "", "",
            "", "", "", ""
        )


def analyse_multiple(organism, files):

    if not files:
        return pd.DataFrame(
            columns=[
                "Sample ID",
                "FASTA file",
                "Organism",
                "Sequence length",
                "GC content",
                "AMR status",
                "AMRFinderPlus result"
            ]
        )

    results = []

    for index, fasta_file in enumerate(files, start=1):

        sample_id = f"AMR-BATCH-{index:03d}"
        filename = os.path.basename(fasta_file)

        try:
            sequence = read_fasta(fasta_file)

            if not sequence:
                results.append({
                    "Sample ID": sample_id,
                    "FASTA file": filename,
                    "Organism": organism or "Not specified",
                    "Sequence length": 0,
                    "GC content": "N/A",
                    "AMR status": "Invalid FASTA",
                    "AMRFinderPlus result": "No nucleotide sequence found."
                })
                continue

            length = len(sequence)

            gc_content = (
                ((sequence.count("G") + sequence.count("C")) / length) * 100
                if length > 0
                else 0
            )

            amr_result = run_amrfinder(
                fasta_file,
                organism
            )

            if amr_result.startswith("No AMR determinants"):
                status = "No known AMR determinant detected"

            elif amr_result.startswith("AMRFinderPlus analysis failed"):
                status = "AMR analysis failed"

            else:
                status = "AMR determinant(s) detected"

            results.append({
                "Sample ID": sample_id,
                "FASTA file": filename,
                "Organism": organism or "Not specified",
                "Sequence length": length,
                "GC content": f"{gc_content:.2f}%",
                "AMR status": status,
                "AMRFinderPlus result": amr_result
            })

        except Exception as error:

            results.append({
                "Sample ID": sample_id,
                "FASTA file": filename,
                "Organism": organism or "Not specified",
                "Sequence length": "Error",
                "GC content": "Error",
                "AMR status": "Analysis failed",
                "AMRFinderPlus result": str(error)
            })

    return pd.DataFrame(results)


with gr.Blocks(title="AMR-Detect") as demo:

    gr.Markdown(
        """
        # 🧬 AMR-Detect

        ### Antimicrobial Resistance Detection and Analysis

        Genomic screening using AMRFinderPlus with support for
        single-sample and multiple-sample analysis.
        """
    )

    with gr.Tab("Single Sample"):

        single_organism = gr.Dropdown(
            choices=organisms,
            label="Select organism"
        )

        single_sample_id = gr.Textbox(
            label="Sample ID",
            placeholder="Example: AMR001"
        )

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
            single_id = gr.Textbox(label="Sample ID")
            single_org = gr.Textbox(label="Organism")

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

        single_susceptibility = gr.Textbox(
            label="Susceptibility interpretation",
            lines=5
        )

        single_button.click(
            fn=analyse_single,
            inputs=[
                single_organism,
                single_sample_id,
                single_file
            ],
            outputs=[
                single_summary,
                single_id,
                single_org,
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
                single_susceptibility
            ]
        )

    with gr.Tab("Multiple Samples"):

        batch_organism = gr.Dropdown(
            choices=organisms,
            label="Select organism",
            info="The selected organism will be used for all uploaded samples."
        )

        batch_files = gr.File(
            label="Upload multiple FASTA files",
            file_types=[".fa", ".fasta", ".fna"],
            type="filepath",
            file_count="multiple"
        )

        batch_button = gr.Button(
            "🔬 Analyse Multiple Samples",
            variant="primary"
        )

        batch_results = gr.Dataframe(
            headers=[
                "Sample ID",
                "FASTA file",
                "Organism",
                "Sequence length",
                "GC content",
                "AMR status",
                "AMRFinderPlus result"
            ],
            label="Multiple Sample AMR Comparison",
            interactive=False,
            wrap=True
        )

        batch_button.click(
            fn=analyse_multiple,
            inputs=[
                batch_organism,
                batch_files
            ],
            outputs=batch_results
        )

    gr.Markdown(
        """
        ---
        
        ### 🩺 Professional review

        AMR-Detect provides genomic screening information and
        experimental analysis support. Genomic findings do not
        independently confirm clinical antimicrobial susceptibility.
        Final interpretation and treatment decisions require qualified
        healthcare-professional review together with validated
        antimicrobial susceptibility testing (AST) and clinical context.
        """
    )


port = int(os.environ.get("PORT", "7860"))

demo.launch(
    server_name="0.0.0.0",
    server_port=port
)
