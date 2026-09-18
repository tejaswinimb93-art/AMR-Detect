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


def read_fasta(filepath):
    """Read a FASTA file and return the nucleotide sequence."""
    with open(filepath, "r", encoding="utf-8") as f:
        lines = f.readlines()

    sequence = "".join(
        line.strip()
        for line in lines
        if not line.startswith(">")
    ).upper()

    return sequence


def analyze_sample(organism, sample_id, fasta_file):
    if not organism:
        return (
            "Please select an organism.",
            "", "", "", "", "", "", "", "", ""
        )

    if not fasta_file:
        return (
            "Please upload a FASTA file.",
            "", "", "", "", "", "", "", "", ""
        )

    if not sample_id:
        sample_id = "Not provided"

    try:
        sequence = read_fasta(fasta_file)

        if not sequence:
            return (
                "The FASTA file does not contain a nucleotide sequence.",
                "", "", "", "", "", "", "", "", ""
            )

        allowed = set("ACGTN")

        # Keep only nucleotide characters
        sequence = "".join(base for base in sequence if base in allowed)

        length = len(sequence)

        a_count = sequence.count("A")
        c_count = sequence.count
