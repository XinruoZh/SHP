# SHP Docker image
# Contains all bioinformatics tools needed for the pipeline.
# Nextflow runs on the HOST and invokes this container per-process.
#
# Tools:
#   Python 3.10, BioPython 1.85, Matplotlib 3.9.4, NumPy 2.0.2, Pandas 2.3.1
#   BLAST+ 2.16.0
#   Clustal Omega (multiple sequence alignment)
#   WebLogo 3.7.9
#   Ghostscript (required by WebLogo for EPS output)
#   CD-HIT (protein clustering)
#   HMMER 3.4 (phmmer, jackhmmer, hmmsearch)
#
# Build:
#   docker build -t shp:latest ./docker/

FROM continuumio/miniconda3:latest

LABEL maintainer="SHP pipeline"
LABEL description="Short Peptide Hunting pipeline dependencies"

# Install libmamba solver for faster dependency resolution
RUN conda install -y -n base conda-libmamba-solver \
    && conda config --set solver libmamba

# Install Clustal Omega via apt (not available in conda-forge/bioconda for this platform)
RUN apt-get update && apt-get install -y --no-install-recommends clustalo gawk \
    && rm -rf /var/lib/apt/lists/*

# Increase conda download retries to handle flaky connections
RUN conda config --set remote_max_retries 10 \
    && conda config --set remote_backoff_factor 2

# Bioinformatics tools (split into layers for Docker cache resilience)
RUN conda install -y -c conda-forge -c bioconda \
    python=3.10 \
    blast=2.16.0 \
    cd-hit \
    && conda clean -afy

RUN conda install -y -c conda-forge -c bioconda \
    biopython=1.85 \
    numpy=2.0.2 \
    pandas=2.3.1 \
    matplotlib=3.9.4 \
    && conda clean -afy

RUN conda install -y -c conda-forge -c bioconda \
    weblogo=3.7.9 \
    ghostscript \
    && conda clean -afy

RUN conda install -y -c conda-forge -c bioconda \
    hmmer=3.4 \
    && conda clean -afy

# Matplotlib non-interactive backend (safeguard; alignment_figure.py also sets it explicitly)
ENV MPLBACKEND=Agg

WORKDIR /work
