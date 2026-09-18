import os
import gradio as gr
from Bio import SeqIO


def analyze_sample(organism, sample_id, fasta_file):
    if not fasta_file:
        return "Please upload a FASTA file."

    if not sample_id:
        sample_id = "Not provided"

    if not organism:
        organism = "Not specified"

    try:
        records = list(SeqIO.parse(fasta_file, "fasta"))

        if not records:
            return "Error: The FASTA file does not contain a valid sequence."

        sequence = "".join(str(record.seq).upper() for record in records)

        valid_bases = set("ACGTN")

        if any(base not in valid_bases for base in sequence):
            return (
                "Error: The FASTA sequence contains invalid characters. "
                "Please upload a DNA sequence containing A, C, G, T or N."
            )

        length = len(sequence)

        a_count = sequence.count("A")
        c_count = sequence.count("C")
        g_count = sequence.count("G")
        t_count = sequence.count("T")
        n_count = sequence.count("N")

        gc_count = g_count + c_count
        gc_content = (gc_count / length * 100) if length > 0 else 0

        filename = os.path.basename(fasta_file)

        return f"""AMR-Detect Sequence Analysis

Sample ID: {sample_id}
Organism: {organism}
FASTA file: {filename}

Sequence length: {length:,} bp

Base composition:
A: {a_count:,}
C: {c_count:,}
G: {g_count:,}
T: {t_count:,}
N: {n_count:,}

GC content: {gc_content:.2f}%

Status:
FASTA sequence successfully received and analysed.

Next stage:
AMR detection engine will analyse the sequence for antimicrobial
resistance determinants.
"""

    except Exception as e:
        return f"Error while reading FASTA file: {str(e)}"


demo = gr.Interface(
    fn=analyze_sample,
    inputs=[
        gr.Textbox(
            label="Organism",
            placeholder="Example: Escherichia coli (leave blank if unknown)"
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
    outputs=gr.Textbox(
        label="AMR-Detect Result"
    ),
    title="AMR-Detect",
    description=(
        "Antimicrobial Resistance Detection and Analysis"
    ),
    flagging_mode="never"
)


port = int(os.environ.get("PORT", "10000"))

demo.launch(
    server_name="0.0.0.0",
    server_port=port
)
