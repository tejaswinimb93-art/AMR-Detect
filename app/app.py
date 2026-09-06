import os
import gradio as gr

def check_status():
    return "AMR-Detect is working! ✅"

demo = gr.Interface(
    fn=check_status,
    inputs=[],
    outputs="text",
    title="AMR-Detect",
    description="Antimicrobial Resistance Detection and Analysis"
)

port = int(os.environ.get("PORT", "10000"))

demo.launch(
    server_name="0.0.0.0",
    server_port=port
)
