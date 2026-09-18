FROM mambaorg/micromamba:latest

WORKDIR /app

COPY --chown=$MAMBA_USER:$MAMBA_USER requirements.txt .

RUN micromamba install -y -n base \
    -c conda-forge \
    -c bioconda \
    python=3.11 \
    ncbi-amrfinderplus=4.2.7 \
    && micromamba clean --all --yes

RUN /opt/conda/bin/python -m pip install --no-cache-dir -r requirements.txt

RUN /opt/conda/bin/amrfinder -u

COPY --chown=$MAMBA_USER:$MAMBA_USER app/ .

EXPOSE 7860

CMD ["/opt/conda/bin/python", "app.py"]
