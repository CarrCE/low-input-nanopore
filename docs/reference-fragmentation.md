# Reference fragmentation: does a contig-level reference cost reads?

Two `lowinput_s1` references are contig-level assemblies of the exact strains in
the ZymoBIOMICS D6311 community: *Saccharomyces cerevisiae* `GCA_030867715.1`
(51 contigs) and *Cryptococcus neoformans* `GCA_028975465.1` (45 contigs). Every
other reference in the panel is complete. A fragmented reference could lose
alignment competitions at contig boundaries, and the same two yeasts account for
almost all the bases the 10% attribution floor removes. This test measures both
effects. It closes the first part of [`TODO.md`](TODO.md) item 5.

**Result: fragmentation changes nothing, and the floor's yeast losses come from
genome content, not assembly quality.** Cutting complete yeast genomes where the
assemblies break changes the assignment of no yeast read and the floor's removals
by not one base. The exact strains, contig-level as published, assign 7.8% more
*S. cerevisiae* reads than complete genomes of other strains, so the published
choice is the better one.

## Design

Three arms on the three D6311 replicates (`lowinput_s1_r1`–`r3`), competitive
mode. Only the two yeast references differ between arms.

| Arm | Yeast references | Pieces | Reference set |
|---|---|---|---|
| A | the published run: exact strains, contig-level | 51 and 45 contigs | `assets/references/lowinput_s1.tsv` |
| B | complete genomes of other strains: S288C `GCF_000146045.2`, H99 `GCF_000149245.1` | 17 and 15 sequences | `assets/fragmentation/lowinput_s1_complete.tsv` |
| C | B's genomes cut where A's assemblies break | 64 and 54 pieces | `assets/fragmentation/lowinput_s1_fragmented.tsv` |

- **B against C isolates fragmentation.** The sequence is the same; only the
  breaks differ.
- **A against B is the strain swap.** It mixes fragmentation with strain
  divergence. The exact strains are about 95% identical to S288C (*S. cerevisiae*)
  and about 90% identical to H99 (*C. neoformans*, whose assembly is 26.3 Mb
  against H99's 18.9 Mb).
- **Each arm runs with the floor on (`--min_aln_frac 0.10`) and off (`0`).** The
  difference is what the floor removes.

**How the cuts are placed** (`bin/fragment_reference.py`). Each assembly is
aligned to the complete genome of its species with minimap2 (`-x asm20`).
Alignments shorter than 10 kb or with MAPQ 0 are ignored. Walking along each
complete chromosome, a cut is placed wherever the aligned contig changes, midway
between the two alignments. A cut within 5 kb of a chromosome end or of the
previous cut is dropped. This gives 47 cuts in S288C and 39 in H99. Pieces are
named `<chromosome>_frag<N>`, and the cut positions are written to
`cuts_<accession>.tsv`.

**C is the conservative test.** Its genomes are more fragmented than the real
assemblies: N50 529 kb against 595 kb for *S. cerevisiae*, and 556 kb against
1.45 Mb for *C. neoformans*.

## Running it

```sh
make all nofloor        # arm A, floor on and off
make fragmentation      # arms B and C, the cuts, and the comparison
```

`make fragmentation` runs `bin/fragmentation_test.sh`. It takes `NF_ARGS`,
`REFDIR` and `PROFILE` as the other targets do, for example
`make fragmentation NF_ARGS="--fetch_from_sra --sradir DIR"` to read the
archive's copies. It takes about 75 minutes: about 35 for each of B and C, and
about 4 for each floor-off run, which reuses the mapping.

The pipeline downloads genomes by accession and cannot download a cut one. The
script therefore copies B's genomes, cuts the two yeasts, and places the result
in the reference cache (`<refdir>/lowinput_s1_fragmented`) before arm C runs.
The genome download is then skipped.

Everything is written under `results/fragmentation/`:

| Path | Contents |
|---|---|
| `complete/`, `complete_nofloor/` | arm B, floor on and off |
| `fragmented/`, `fragmented_nofloor/` | arm C, floor on and off |
| `build/` | the assembly-to-genome alignments (`*.paf`) and the cut positions (`cuts_*.tsv`) |
| `summary/fragmentation.json` | the comparison: yeast reads per arm, every read that changes between B and C, the reads whose alignment in B crosses a cut, floor removals, enrichment |
| `summary/fragmentation_by_organism.tsv` | reads and bases per organism in each arm |

## Results

Run end to end on 5 Nov 2026 from a clean checkout, on ENA's copies of the reads.
The run reproduced a first, hand-run version of 4 Nov exactly, and the cut
genomes were byte-identical.

### Floor on

| Arm | *S. cerevisiae* reads, r1 / r2 / r3 | *C. neoformans* reads | D6311 enrichment (bases, mean) |
|---|---|---|---:|
| A, as published | 1,467 / 996 / 991 | 4 / 0 / 9 | 53.6311 |
| B, complete genomes | 1,351 / 937 / 898 | 0 / 0 / 6 | 53.6262 |
| C, B cut at A's breaks | 1,351 / 937 / 898 | 0 / 0 / 6 | 53.6265 |

**B against C:**
- Every yeast count is identical, read for read and base for base.
- 39 *S. cerevisiae* reads have their primary alignment in B across one of C's
  cuts. All 39 are still assigned to *S. cerevisiae* in C. No *C. neoformans*
  read crosses a cut.
- 44 of 25,039,598 reads (0.00018%) change in total, none of them involving a
  yeast. All are near-ties among *E. coli*, *E. coli* K-12 and λ. Changing the
  reference set changes minimap2's index, and this is the same near-tie noise
  seen when the analysis was rerun from the archive (`TODO.md` item 1).

**A against B:** the complete genomes of other strains assign 268 fewer
*S. cerevisiae* reads (−7.8%) and 7 fewer *C. neoformans* reads. Those reads become
unassigned. The exact strains are the better reference, contig-level or not.

### Floor off: what the floor removes

| Arm | Community bases the floor removes | Yeast share | *S. cerevisiae* | *C. neoformans* | D6311 enrichment, floor off → on |
|---|---:|---:|---:|---:|---|
| A, as published | 23,356,872 | 96.5% | 21,545,310 | 982,423 | 54.58 → 53.63 |
| B, complete genomes | 9,912,829 | 90.5% | 7,758,673 | 1,215,467 | 54.03 → 53.63 |
| C, B cut at A's breaks | 9,912,829 | 90.5% | 7,758,673 | 1,215,467 | 54.03 → 53.63 |

B and C are identical to the base. Even against complete chromosome-level
genomes, the yeasts account for 90% of what the floor removes. The cause is the
yeasts' genome content, not the assemblies' fragmentation.

## What this does not test

- **Two genomes, one community.** The test covers the panel's only contig-level
  references. It says nothing about references fragmented in other ways, such as
  short-read assemblies with thousands of contigs.
- **The cuts are where A's contigs change along B.** Breaks inside sequence that B
  lacks, or that the alignment filter drops, are not reproduced.
