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
    sequence = ""

    with open(filepath, "r", encoding="utf-8") as file:
        for line in file:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                continue

            sequence += line.upper()

    return sequence


def format_amr_result(output):
    if not output.strip():
        return "No AMR determinants were detected by AMRFinderPlus."

    lines = output.strip().splitlines()

    if len(lines) <= 1:
        return output.strip()

    header = lines[0].split("\t")

    results = []

    for line in lines[1:]:
        values = line.split("\t")

        row = dict(zip(header, values))

        element_symbol = row.get(
            "Element symbol",
            row.get("Gene symbol", "Not reported")
        )

        element_name = row.get(
            "Element name",
            "Not reported"
        )

        resistance = row.get(
            "Resistance",
            "Not reported"
        )

        method = row.get(
            "Method",
            "Not reported"
        )

        target = row.get(
            "Target",
            "Not reported"
        )

        identity = row.get(
            "% Identity",
            "Not reported"
        )

        coverage = row.get(
            "% Coverage",
            "Not reported"
        )

        results.append(
            f"AMR determinant: {element_symbol}\n"
            f"Name: {element_name}\n"
            f"Resistance: {resistance}\n"
            f"Method: {method}\n"
            f"Target: {target}\n"
            f"Identity: {identity}\n"
            f"Coverage: {coverage}\n"
        )

    return "\n------------------------------\n\n".join(results)


def run_amrfinder(fasta_file, organism):
    command = [
        "amrfinder",
        "-n",
        fasta_file
    ]

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
            error_message = result.stderr.strip()

            if not error_message:
                error_message = "AMRFinderPlus returned an unknown error."

            return (
                "AMRFinderPlus analysis failed.\n\n"
                + error_message
            )

        return format_amr_result(result.stdout)

    except subprocess.TimeoutExpired:
        return (
            "AMRFinderPlus analysis timed out. "
            "Please try a smaller FASTA file."
        )

    except FileNotFoundError:
        return (
            "AMRFinderPlus is not available in the deployment environment."
        )

    except Exception as error:
        return f"AMRFinderPlus error: {error}"


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

        if any(base not in allowed_bases for base in sequence):
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

        summary = (
            "Sample received and analysed successfully."
        )

        status = "AMR analysis completed"

        next_stage = (
            "Review AMRFinderPlus resistance determinants "
            "and supporting sequence information."
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
            amr_result
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
        gr.Textbox(label="AMRFinderPlus Result")
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
