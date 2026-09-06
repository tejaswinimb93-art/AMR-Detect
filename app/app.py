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


def analyze_sample(organism, sample_id, fasta_file):
    if not organism:
        return "Please select an organism."

    if not fasta_file:
        return "Please upload a FASTA file."

    if not sample_id:
        sample_id = "Not provided"

    filename = os.path.basename(fasta_file)

    return f"""Sample received successfully.

Sample ID: {sample_id}
Organism: {organism}
FASTA file: {filename}

AMR analysis engine will be connected in the next step.
"""


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
    outputs=gr.Textbox(
        label="AMR-Detect Result"
    ),
    title="AMR-Detect",
    description="Antimicrobial Resistance Detection and Analysis",
    flagging_mode="never"
)


port = int(os.environ.get("PORT", "10000"))

demo.launch(
    server_name="0.0.0.0",
    server_port=port
)
