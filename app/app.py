import os
import subprocess
import gradio as gr


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
    "Klebsiella pneumoniae": "Klebsiella",
    "Staphylococcus aureus": "Staphylococcus",
    "Pseudomonas aeruginosa": "Pseudomonas",
    "Acinetobacter baumannii": "Acinetobacter",
    "Enterococcus faecalis": "Enterococcus",
    "Salmonella spp.": "Salmonella"
}


def empty_result(message):
    return (
        message, "", "", "", "", "", "", "", "",
        "", "", "", "", ""
    )


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
                + (result.stderr.strip() or "Unknown AMRFinderPlus error.")
            )

        output = result.stdout.strip()

        if not output:
            return "No AMR determinants were detected."

        return output

    except subprocess.TimeoutExpired:
        return "AMRFinderPlus analysis timed out."

    except FileNotFoundError:
        return "AMRFinderPlus was not found in the deployment environment."

    except Exception as error:
        return f"AMRFinderPlus error: {error}"


def analyse_sample(organism, sample_id, fasta_file):

    if not fasta_file:
        return empty_result("Please upload a FASTA file.")

    if not organism:
        organism = "Not specified"

    if not sample_id:
        sample_id = "Not provided"

    try:
        sequence = read_fasta(fasta_file)

        if not sequence:
            return empty_result(
                "Error: The FASTA file does not contain a DNA sequence."
            )

        allowed = set("ACGTN")

        if any(base not in allowed for base in sequence):
            return empty_result(
                "Error: The FASTA sequence contains invalid DNA characters."
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
            susceptibility_note = (
                "No known genomic AMR determinant was detected "
                "by this screening. This does not confirm clinical "
                "susceptibility; phenotypic antimicrobial susceptibility "
                "testing (AST) is required."
            )
        elif amr_result.startswith("AMRFinderPlus analysis failed"):
            amr_status = "AMR analysis could not be completed"
            susceptibility_note = (
                "Susceptibility cannot be inferred because the "
                "genomic AMR analysis did not complete."
            )
        else:
            amr_status = "AMR determinant(s) detected"
            susceptibility_note = (
                "Antibiotics without a detected known genomic resistance "
                "determinant should not automatically be considered "
                "clinically susceptible. Phenotypic AST is required "
                "for confirmation."
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
            susceptibility_note
        )

    except Exception as error:
        return empty_result(
            f"Error while processing the FASTA file: {error}"
        )


demo = gr.Interface(
    fn=analyse_sample,

    inputs=[
        gr.Dropdown(
            choices=organisms,
            label="Select organism",
            info="Choose the organism associated with the sample."
        ),

        gr.Textbox(
            label="Sample ID",
            placeholder="Example: AMR001"
        ),

        gr.File(
            label="Upload FASTA sequence",
            file_types=[".fa", ".fasta", ".fna"],
            type="filepath"
        )
    ],

    outputs=[
        gr.Textbox(label="Analysis summary"),
        gr.Textbox(label="Sample ID"),
        gr.Textbox(label="Organism"),
        gr.Textbox(label="FASTA file"),
        gr.Textbox(label="Sequence length"),
        gr.Textbox(label="A"),
        gr.Textbox(label="G"),
        gr.Textbox(label="C"),
        gr.Textbox(label="T"),
        gr.Textbox(label="N"),
        gr.Textbox(label="GC content"),
        gr.Textbox(label="AMR status"),
        gr.Textbox(label="AMRFinderPlus Result"),
        gr.Textbox(label="Susceptibility interpretation")
    ],

    title="AMR-Detect",

    description=(
        "Antimicrobial Resistance Detection and Analysis"
    ),

    flagging_mode="never"
)


port = int(os.environ.get("PORT", "7860"))

demo.launch(
    server_name="0.0.0.0",
    server_port=port
)
