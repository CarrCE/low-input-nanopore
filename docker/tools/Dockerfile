# Alignment/utility container for the low-input-nanopore pipeline.
#
# Built from source rather than pulled from bioconda because the bioconda
# images are linux/amd64 only, while this study's reference platform is Apple
# Silicon (arm64). Measured on an M3 Max, minimap2 map-ont over a 400 MB ONT
# read slice against lambda + E. coli K-12, 8 threads:
#
#   native arm64 (this image)                        6.64 s
#   amd64 under Docker Desktop Rosetta emulation     9.35 s   (1.41x slower)
#
# So the gain is ~40%, not the order of magnitude one might expect -- Rosetta
# handles this workload well. The native build is kept because 40% is material
# across a ~45 Gbp run, and because building from pinned source makes the image
# reproducible on both architectures without depending on Rosetta being
# available. See docs/benchmarks.md.
FROM debian:bookworm-20241111-slim AS build

ARG MINIMAP2_VERSION=2.28
ARG HTSLIB_VERSION=1.21
ARG SAMTOOLS_VERSION=1.21

RUN apt-get update \
 && apt-get install -y --no-install-recommends \
      build-essential ca-certificates curl \
      zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev libssl-dev \
      libncurses-dev autoconf automake \
 && rm -rf /var/lib/apt/lists/*

WORKDIR /src

# ---- minimap2 -------------------------------------------------------------
# On aarch64 minimap2 needs the NEON/SIMDe path selected explicitly.
RUN curl -fsSL "https://github.com/lh3/minimap2/archive/refs/tags/v${MINIMAP2_VERSION}.tar.gz" \
      -o minimap2.tar.gz \
 && tar xf minimap2.tar.gz \
 && cd "minimap2-${MINIMAP2_VERSION}" \
 && if [ "$(uname -m)" = "aarch64" ]; then \
        make arm_neon=1 aarch64=1 ; \
    else \
        make ; \
    fi \
 && install -m 0755 minimap2 /usr/local/bin/minimap2 \
 && install -m 0755 misc/paftools.js /usr/local/bin/paftools.js

# ---- htslib + samtools ----------------------------------------------------
RUN curl -fsSL "https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2" \
      -o htslib.tar.bz2 \
 && tar xf htslib.tar.bz2 \
 && cd "htslib-${HTSLIB_VERSION}" \
 && ./configure --prefix=/usr/local --enable-libcurl \
 && make -j"$(nproc)" && make install

RUN curl -fsSL "https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2" \
      -o samtools.tar.bz2 \
 && tar xf samtools.tar.bz2 \
 && cd "samtools-${SAMTOOLS_VERSION}" \
 && ./configure --prefix=/usr/local --with-htslib=/usr/local \
 && make -j"$(nproc)" && make install

# ---- seqkit (prebuilt Go binary) -----------------------------------------
ARG SEQKIT_VERSION=2.8.2
RUN case "$(uname -m)" in \
      aarch64) SK_ARCH=arm64 ;; \
      x86_64)  SK_ARCH=amd64 ;; \
      *) echo "unsupported arch $(uname -m)" >&2; exit 1 ;; \
    esac \
 && curl -fsSL "https://github.com/shenwei356/seqkit/releases/download/v${SEQKIT_VERSION}/seqkit_linux_${SK_ARCH}.tar.gz" \
      -o seqkit.tar.gz \
 && tar xf seqkit.tar.gz \
 && install -m 0755 seqkit /usr/local/bin/seqkit

# ---- NCBI datasets CLI ----------------------------------------------------
ARG DATASETS_VERSION=v2
RUN case "$(uname -m)" in \
      aarch64) DS_ARCH=linux-arm64 ;; \
      x86_64)  DS_ARCH=linux-amd64 ;; \
    esac \
 && curl -fsSL "https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/${DATASETS_VERSION}/${DS_ARCH}/datasets" \
      -o /usr/local/bin/datasets \
 && curl -fsSL "https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/${DATASETS_VERSION}/${DS_ARCH}/dataformat" \
      -o /usr/local/bin/dataformat \
 && chmod 0755 /usr/local/bin/datasets /usr/local/bin/dataformat

# ---------------------------------------------------------------------------
FROM debian:bookworm-20241111-slim

LABEL org.opencontainers.image.source="https://github.com/CarrCE/low-input-nanopore"
LABEL org.opencontainers.image.description="minimap2 + samtools + seqkit + NCBI datasets, built natively for amd64 and arm64"

RUN apt-get update \
 && apt-get install -y --no-install-recommends \
      zlib1g libbz2-1.0 liblzma5 libcurl4 libssl3 libncurses6 \
      ca-certificates procps unzip curl \
 && rm -rf /var/lib/apt/lists/*

COPY --from=build /usr/local/bin/ /usr/local/bin/
COPY --from=build /usr/local/lib/ /usr/local/lib/
RUN ldconfig

WORKDIR /work

# Fail fast at build time if anything is broken or mis-linked.
RUN minimap2 --version && samtools --version | head -1 \
 && seqkit version && datasets --version
