# Reference genomes for validation

Ground-truth reference sets for the Tier-1 (reference-based) validation harness. Reads
are mapped to these known genomes to score per-bin purity and per-organism completeness.

Provenance: accessions traced to Zymo Research Corporation's own NCBI BioProjects
(PRJNA933688–PRJNA1003957) and cross-checked against the strain table in Nicholls et al.
2019 (GigaScience, PMC6520541). The genomes are **not** tracked in git — fetch them with
the `datasets` command below (they land in the git-ignored `refs/`).

## Set A — `ur_r3`: D6311 Microbial Community DNA Standard II (Log Distribution)

| Species | Strain | Accession | Notes |
|---|---|---|---|
| Listeria monocytogenes | B-33116 = ATCC 19117 | GCF_028743575.1 | complete |
| Pseudomonas aeruginosa | B-3509 = ATCC 15442 | GCF_028743595.1 | complete |
| Bacillus subtilis (→ B. spizizenii) | B-354 = ATCC 6633 | GCF_028743795.1 | complete |
| Escherichia coli | B-1109 | GCF_028743555.1 | complete |
| Salmonella enterica | B-4212 | GCF_028743635.1 | complete |
| Limosilactobacillus fermentum | B-1840 = ATCC 14931 | GCF_030770375.1 | complete |
| Enterococcus faecalis | B-537 = ATCC 7080 | GCF_028743535.1 | complete |
| Staphylococcus aureus | B-41012 | GCF_028743615.1 | complete |
| Saccharomyces cerevisiae | Y-567 = ATCC 9763 | GCA_030867715.1 | **contig-level** (exact strain) |
| Cryptococcus neoformans | Y-2534 = ATCC 32045 | GCA_028975465.1 | **contig-level** (exact strain) |

Complete-genome substitutes for the two yeasts (different strain, chromosome-level), if a
fragmented reference proves problematic for mapping:
`S. cerevisiae` S288C `GCF_000146045.2`; `C. neoformans` H99 `GCF_000149245.1`.

## Set B — `lowinput_s2_r4`: D6321 Spike-in Control II (Low Microbial Load, cells)

| Species | Strain | Accession | Notes |
|---|---|---|---|
| Truepera radiovictrix | CIP 108686 (DSM 17093) | GCF_031201145.1 | complete |
| Imtechella halotolerans | JCM 17677 (type K1) | GCF_028743515.2 | complete |
| Allobacillus halotolerans | BCRC 17939 (type B3A) | GCF_028743495.1 | complete |

## Carrier

Enterobacteria phage lambda — `NC_001416.1` (48,502 bp), assembly `GCF_000840245.1`.

## Download

```bash
mkdir -p refs
docker run --rm -v "$PWD/refs":/work -w /work ensemblorg/datasets-cli:v18.2.0 \
  datasets download genome accession \
    GCF_028743575.1 GCF_028743595.1 GCF_028743795.1 GCF_028743555.1 \
    GCF_028743635.1 GCF_030770375.1 GCF_028743535.1 GCF_028743615.1 \
    GCA_030867715.1 GCA_028975465.1 \
    GCF_031201145.1 GCF_028743515.2 GCF_028743495.1 \
    GCF_000840245.1 \
    GCF_000146045.2 GCF_000149245.1 \
    --include genome --filename zymo_refs.zip
unzip -o refs/zymo_refs.zip -d refs/extracted
```
