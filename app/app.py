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

demo.launch()
