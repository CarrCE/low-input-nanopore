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
FROM debian:bookworm-20241111-slim@sha256:ca3372ce30b03a591ec573ea975ad8b0ecaf0eb17a354416741f8001bbcae33d AS build

# Every download below is checked against a recorded SHA-256. A date-tagged
# Debian is still a mutable tag, GitHub tag archives are not byte-stable, and a
# moved tag or a substituted artifact would otherwise be built into the image
# that produces every reference in the paper. A mismatch here is intended to
# stop the build: update the hash deliberately, having looked at what changed.
ARG MINIMAP2_VERSION=2.28
ARG MINIMAP2_SHA256=5ea6683b4184b5c49f6dbaef2bc5b66155e405888a0790d1b21fd3c93e474278
ARG HTSLIB_VERSION=1.21
ARG HTSLIB_SHA256=84b510e735f4963641f26fd88c8abdee81ff4cb62168310ae716636aac0f1823
ARG SAMTOOLS_VERSION=1.21
ARG SAMTOOLS_SHA256=05724b083a6b6f0305fcae5243a056cc36cf826309c3cb9347a6b89ee3fc5ada

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
 && echo "${MINIMAP2_SHA256}  minimap2.tar.gz" | sha256sum -c - \
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
 && echo "${HTSLIB_SHA256}  htslib.tar.bz2" | sha256sum -c - \
 && tar xf htslib.tar.bz2 \
 && cd "htslib-${HTSLIB_VERSION}" \
 && ./configure --prefix=/usr/local --enable-libcurl \
 && make -j"$(nproc)" && make install

RUN curl -fsSL "https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2" \
      -o samtools.tar.bz2 \
 && echo "${SAMTOOLS_SHA256}  samtools.tar.bz2" | sha256sum -c - \
 && tar xf samtools.tar.bz2 \
 && cd "samtools-${SAMTOOLS_VERSION}" \
 && ./configure --prefix=/usr/local --with-htslib=/usr/local \
 && make -j"$(nproc)" && make install

# ---- seqkit (prebuilt Go binary) -----------------------------------------
ARG SEQKIT_VERSION=2.8.2
ARG SEQKIT_SHA256_AMD64=a6515a64bf28ceba6922a3854d8ec898ca2ca637b4b00fc883c83076dfaf1ffa
ARG SEQKIT_SHA256_ARM64=8c12816c0bf6c01cfcda329fdb5f6c1479cf2485ce9413ce36b96237524da9dc
RUN case "$(uname -m)" in \
      aarch64) SK_ARCH=arm64; SK_SHA="${SEQKIT_SHA256_ARM64}" ;; \
      x86_64)  SK_ARCH=amd64; SK_SHA="${SEQKIT_SHA256_AMD64}" ;; \
      *) echo "unsupported arch $(uname -m)" >&2; exit 1 ;; \
    esac \
 && curl -fsSL "https://github.com/shenwei356/seqkit/releases/download/v${SEQKIT_VERSION}/seqkit_linux_${SK_ARCH}.tar.gz" \
      -o seqkit.tar.gz \
 && echo "${SK_SHA}  seqkit.tar.gz" | sha256sum -c - \
 && tar xf seqkit.tar.gz \
 && install -m 0755 seqkit /usr/local/bin/seqkit

# ---- NCBI datasets CLI ----------------------------------------------------
# NCBI serves this path ROLLING: /v2/ is always their current build, so the
# same Dockerfile silently produced different binaries over time. There is no
# versioned URL to pin to, so the version is pinned by hash instead. These are
# the binaries that built every reference in the paper (datasets 18.33.1). When
# NCBI ships a new build the checksum fails and the build stops; that is the
# intent -- record the new hash and the new version together in
# results/summary/software_versions.tsv, deliberately.
ARG DATASETS_VERSION=v2
ARG DATASETS_SHA256_AMD64=8459ef1e87433f7b1198f5703c8cc10b55f1904cd448cf7e996c0892b141cd1f
ARG DATAFORMAT_SHA256_AMD64=8450cf7cbdb0ed7fece567405732cd1ff838b5352faa58abbccdeff56e1ff0e8
ARG DATASETS_SHA256_ARM64=54816c2e8586c667ed8107e1fbf94efa6d924845d46b8c4c12f1be21c8353260
ARG DATAFORMAT_SHA256_ARM64=e6f7e9abcde5f231ab92f8a2748938fa4ddd6f3e8de710681a53c591cb35ddf1
RUN case "$(uname -m)" in \
      aarch64) DS_ARCH=linux-arm64; DS_SHA="${DATASETS_SHA256_ARM64}"; DF_SHA="${DATAFORMAT_SHA256_ARM64}" ;; \
      x86_64)  DS_ARCH=linux-amd64; DS_SHA="${DATASETS_SHA256_AMD64}"; DF_SHA="${DATAFORMAT_SHA256_AMD64}" ;; \
      *) echo "unsupported arch $(uname -m)" >&2; exit 1 ;; \
    esac \
 && curl -fsSL "https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/${DATASETS_VERSION}/${DS_ARCH}/datasets" \
      -o /usr/local/bin/datasets \
 && echo "${DS_SHA}  /usr/local/bin/datasets" | sha256sum -c - \
 && curl -fsSL "https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/${DATASETS_VERSION}/${DS_ARCH}/dataformat" \
      -o /usr/local/bin/dataformat \
 && echo "${DF_SHA}  /usr/local/bin/dataformat" | sha256sum -c - \
 && chmod 0755 /usr/local/bin/datasets /usr/local/bin/dataformat

# ---------------------------------------------------------------------------
FROM debian:bookworm-20241111-slim@sha256:ca3372ce30b03a591ec573ea975ad8b0ecaf0eb17a354416741f8001bbcae33d

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
