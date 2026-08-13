# NCBI's Human Read Removal Tool (HRRT / sra-human-scrubber), plus the one
# thing Nextflow needs that the upstream image does not ship.
#
# Why layer instead of rebuilding from source: the upstream image carries the
# k-mer database (human_filter.db.20231218v2, 1.02 GB) that is not separately
# versioned or downloadable at a stable URL. Pinning the base by DIGEST
# therefore pins the tool and its database together, which a from-source build
# could not do. Everything that decides a masking outcome comes from that
# digest; the layer below adds no analysis capability whatsoever.
#
# HRRT is amd64-only -- there is no arm64 build -- so this image is amd64 and
# runs under Rosetta on Apple Silicon, as breseq already does. Measured cost is
# acceptable: HRRT was 760 s of a 954 s pass over a 17 GB read file, and it is
# the dominant phase on either architecture.
#
# Build with:  docker build --platform linux/amd64 -t low-input-nanopore/scrubber:0.1.0 docker/scrubber
FROM ncbi/sra-human-scrubber@sha256:3e87c77dd1c00a4eec59051aaf819032d9b33c1d78f810b8bf89fadeb7ff9ad7

# procps supplies `ps`. Nextflow's task wrapper collects metrics via
#     done < <(ps -e -o pid= -o ppid=)
# inside .command.run, which runs under `bash -euo pipefail`. With no `ps` that
# process substitution fails, -e aborts the wrapper BEFORE .command.sh runs, and
# the task dies with exit 1, empty stdout and empty stderr -- no indication that
# the missing binary is a Nextflow dependency rather than a scrubber problem.
# docker/tools/Dockerfile installs procps for exactly this reason.
RUN apt-get update \
 && apt-get install -y --no-install-recommends procps \
 && rm -rf /var/lib/apt/lists/*

# Fail the build rather than a run if either half is missing.
RUN ps --version && test -s /opt/scrubber/data/human_filter.db \
 && test -x /opt/scrubber/scripts/scrub.sh

WORKDIR /work
