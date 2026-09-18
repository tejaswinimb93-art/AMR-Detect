import os
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


def empty_result(message):
    return (
        message,
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


def read_fasta(filepath):
    with open(filepath, "r", encoding="utf-8") as file:
        lines = file.readlines()

    sequence = ""

    for line in lines:
        line = line.strip()

        if not line:
            continue

        if line.startswith(">"):
            continue

        sequence += line.upper()

    return sequence


def analyze_sample(organism, sample_id, fasta_file):

    if not organism:
        return empty_result("Please select an organism.")

    if not fasta_file:
        return empty_result("Please upload a FASTA file.")

    if not sample_id:
        sample_id = "Not provided"

    try:
        sequence = read_fasta(fasta_file)

        if not sequence:
            return empty_result(
                "Error: The FASTA file does not contain a DNA sequence."
            )

        allowed_bases = set("ACGTN")

        invalid_bases = [
            base for base in sequence
            if base not in allowed_bases
        ]

        if invalid_bases:
            return empty_result(
                "Error: The FASTA sequence contains invalid characters."
            )

        length = len(sequence)

        a_count = sequence.count("A")
        c_count = sequence.count("C")
        g_count = sequence.count("G")
        t_count = sequence.count("T")
        n_count = sequence.count("N")

        if length > 0:
            gc_content = ((g_count + c_count) / length) * 100
        else:
            gc_content = 0

        filename = os.path.basename(fasta_file)

        status = "Ready for AMR screening"
        next_stage = "AMR database analysis"

        summary = (
            "Sample received successfully.\n\n"
            "The FASTA sequence has been successfully "
            "validated and analysed."
        )

        return (
            summary,
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
            status,
            next_stage
        )

    except Exception as error:
        return empty_result(
            f"Error while processing the FASTA file: {error}"
        )


demo = gr.Interface(
    fn=analyze_sample,

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
        gr.Textbox(label="Status"),
        gr.Textbox(label="Next stage")
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
