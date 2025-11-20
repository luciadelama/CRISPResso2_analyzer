# syntax=docker/dockerfile:1.7
FROM continuumio/miniconda3:24.7.1-0

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    CONDA_HTTP_TIMEOUT=120 \
    CONDA_SOLVER=libmamba

# 1) mamba + curl
RUN --mount=type=cache,target=/opt/conda/pkgs \
    conda install -y -c conda-forge mamba curl && conda clean -afy

# 2) Bio-dependencias pesadas (ajusta según necesitéis)
#    CRISPResso2 está en bioconda; a veces requiere samtools/bowtie2/etc.
RUN --mount=type=cache,target=/opt/conda/pkgs \
    mamba install -y -c conda-forge -c bioconda \
      crispresso2 \
      bowtie2 \
      samtools \
      fastqc \
      pigz \
    && conda clean -afy

# 3) Crear usuario no-root
RUN useradd -m -u 1000 appuser

# 4) Instalar deps Python ligeras con pip
WORKDIR /crispr-shiny
COPY requirements.txt /crispr-shiny/requirements.txt
RUN python -m pip install --no-cache-dir -r /crispr-shiny/requirements.txt

# 5) Copiar código
COPY . /crispr-shiny
RUN chown -R appuser:appuser /crispr-shiny
RUN mkdir -p /outputs && chown -R appuser:appuser /outputs
USER appuser

EXPOSE 8000
HEALTHCHECK --interval=30s --timeout=5s --retries=5 \
  CMD curl -fsS http://localhost:8000/ || exit 1

CMD ["sh", "-c", "shiny run --host 0.0.0.0 --port 8000 /crispr-shiny/app.py"]
