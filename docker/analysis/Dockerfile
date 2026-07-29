# Analysis/plotting container for the low-input-nanopore pipeline.
# Deliberately minimal and fully pinned: every number and figure in the
# manuscript is produced inside this image.
FROM python:3.12-slim-bookworm@sha256:d50fb7611f86d04a3b0471b46d7557818d88983fc3136726336b2a4c657aa30b

LABEL org.opencontainers.image.source="https://github.com/CarrCE/low-input-nanopore"
LABEL org.opencontainers.image.description="Metrics + figure generation for depletion-mode adaptive sampling analysis"

# procps is not optional: Nextflow's task wrapper shells out to `ps` to collect
# metrics, and without it the wrapper aborts before ever running the command --
# producing an exit-1 with completely empty stdout/stderr.
RUN apt-get update \
 && apt-get install -y --no-install-recommends \
      zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev libssl-dev gcc \
      procps \
 && rm -rf /var/lib/apt/lists/*

COPY requirements.txt /tmp/requirements.txt
RUN pip install --no-cache-dir -r /tmp/requirements.txt \
 && rm /tmp/requirements.txt

# Headless plotting
ENV MPLBACKEND=Agg

WORKDIR /work
