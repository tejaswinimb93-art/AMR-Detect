FROM mambaorg/micromamba:latest

WORKDIR /app

COPY --chown=$MAMBA_USER:$MAMBA_USER requirements.txt .

RUN micromamba install -y -n base \
    -c conda-forge \
    -c bioconda \
    python=3.11 \
    ncbi-amrfinderplus \
    && micromamba clean --all --yes

RUN pip install --no-cache-dir -r requirements.txt

RUN amrfinder -u

COPY --chown=$MAMBA_USER:$MAMBA_USER app/ .

EXPOSE 7860

CMD ["python", "app.py"]
