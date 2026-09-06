import os
import gradio as gr

def analyze_sample(organism, sample_id):
    if not organism:
        return "Please select an organism."
    
    if not sample_id:
        sample_id = "Not provided"
    
    return f"""Sample received successfully.

Sample ID: {sample_id}
Organism: {organism}

AMR analysis will be performed here.
"""

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
