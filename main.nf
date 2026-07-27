#!/usr/bin/env nextflow
/*
 * low-input-nanopore
 *
 * Quantify how much depletion-mode adaptive sampling enriches sample DNA that
 * was carried into library prep on a genomic carrier (default: lambda, which
 * brings E. coli K-12 with it from its production host).
 *
 * Two assignment modes:
 *   competitive  one minimap2 pass against a combined index; best organism wins
 *                by alignment-score margin, ties reported as ambiguous.
 *   sequential   the classic subtraction chain: remove carrier, then remove
 *                contaminant, then map what survives to the community.
 *   both         run both and emit a per-organism delta table.
 */

nextflow.enable.dsl = 2

def helpMessage() {
    log.info """
    low-input-nanopore ${workflow.manifest.version}

    Usage:
      nextflow run . -profile docker --samplesheet assets/samplesheets/lowinput_s2.csv

    Required:
      --samplesheet     CSV describing samples (see assets/samplesheets/)

    Options:
      --mode            competitive | sequential | both   [${params.mode}]
      --outdir          results directory                 [${params.outdir}]
      --refdir          reference cache directory         [${params.refdir}]
      --breseq_consensus  build a reference-guided consensus of the contaminant
                          actually present in the carrier prep [${params.breseq_consensus}]
      --min_mapq        MAPQ floor before scoring         [${params.min_mapq}]
      --coverage_window coverage bin size in bp           [${params.coverage_window}]
      --keep_bams       retain alignment BAMs             [${params.keep_bams}]
      --max_cpus        cap on CPUs per process           [${params.max_cpus}]
    """.stripIndent()
}

if (params.help) { helpMessage(); exit 0 }
if (!params.samplesheet) { helpMessage(); exit 1, "error: --samplesheet is required" }
if (!(params.mode in ['competitive', 'sequential', 'both'])) {
    exit 1, "error: --mode must be one of competitive, sequential, both"
}
if (params.breseq_consensus && params.mode == 'competitive') {
    exit 1, "error: --breseq_consensus has no effect with --mode competitive. " +
            "Competitive assignment never subtracts, so there is no subtraction " +
            "step for a consensus to improve. Use --mode sequential or --mode both."
}

// ---------------------------------------------------------------------------
// Processes
// ---------------------------------------------------------------------------

process FETCH_GENOMES {
    tag   { ref_name }
    label 'tools'
    label 'process_low'
    storeDir "${params.refdir}/${ref_name}"

    input:
    tuple val(ref_name), path(reference_tsv)

    output:
    tuple val(ref_name), path("genomes"), emit: genomes

    script:
    """
    # Accessions come from the reference set, so the download is driven by the
    # same file that defines the analysis -- they cannot drift apart. Matching
    # on the accession pattern skips both the comment preamble and the header
    # row without depending on how many comment lines precede them.
    accessions=\$(awk -F'\\t' '\$2 ~ /^GC[AF]_/ {print \$2}' ${reference_tsv} | sort -u)
    if [ -z "\$accessions" ]; then
        echo "error: no GCF_/GCA_ accessions found in ${reference_tsv}" >&2; exit 1
    fi
    echo "fetching: \$accessions"

    mkdir -p genomes
    datasets download genome accession \$accessions \\
        --include genome --filename genomes.zip
    unzip -q -o genomes.zip -d genomes
    rm -f genomes.zip

    n=\$(find genomes -name '*.fna' | wc -l | tr -d ' ')
    expected=\$(echo "\$accessions" | wc -w | tr -d ' ')
    if [ "\$n" -lt "\$expected" ]; then
        echo "error: expected \$expected genomes, found \$n" >&2
        exit 1
    fi
    """
}

process BUILD_REFERENCE {
    tag   { ref_name }
    label 'analysis'
    label 'process_low'
    publishDir "${params.outdir}/references/${ref_name}", mode: 'copy'

    input:
    tuple val(ref_name), path(reference_tsv), path(genomes)
    path script

    output:
    tuple val(ref_name), path("combined.fasta"),      emit: fasta
    tuple val(ref_name), path("contig_map.tsv"),      emit: contig_map
    tuple val(ref_name), path("genome_sizes.tsv"),    emit: genome_sizes
    path "reference_provenance.json",                 emit: provenance

    script:
    """
    python3 ${script} \\
        --reference-tsv    ${reference_tsv} \\
        --genome-dir       ${genomes} \\
        --out-fasta        combined.fasta \\
        --out-contig-map   contig_map.tsv \\
        --out-genome-sizes genome_sizes.tsv \\
        --out-provenance   reference_provenance.json
    """
}

process FILTER_READS {
    tag   { meta.sample_id }
    label 'tools'
    label 'process_low'
    publishDir "${params.outdir}/${meta.sample_id}/qc", mode: 'copy', pattern: '*.qfilter.tsv'

    input:
    tuple val(meta), path(fastq)
    path script

    output:
    tuple val(meta), path("${meta.sample_id}.filtered.fastq"), emit: reads
    path "${meta.sample_id}.qfilter.tsv",                      emit: stats

    script:
    // awk cannot read gzip, and it fails silently on compressed input rather
    // than erroring, so decompress explicitly when needed. Published datasets
    // (e.g. the Raghavendra reads) arrive gzipped; this study's do not.
    """
    FQ="${fastq}"
    if [ "\$FQ" != "\${FQ%.gz}" ]; then
        gzip -dc "\$FQ" \\
          | awk -v MINQ=${params.min_qscore} -v STATS=${meta.sample_id}.qfilter.tsv \\
                -f ${script} - > ${meta.sample_id}.filtered.fastq
    else
        awk -v MINQ=${params.min_qscore} -v STATS=${meta.sample_id}.qfilter.tsv \\
            -f ${script} "\$FQ" > ${meta.sample_id}.filtered.fastq
    fi

    if [ ! -s ${meta.sample_id}.filtered.fastq ]; then
        echo "error: quality filter produced no reads from \$FQ" >&2; exit 1
    fi
    """
}

process MAP_COMPETITIVE {
    tag   { meta.sample_id }
    label 'tools'
    label 'process_high'
    publishDir "${params.outdir}/${meta.sample_id}/alignments", mode: 'copy',
               enabled: params.keep_bams, pattern: '*.bam'

    input:
    tuple val(meta), path(fastq), path(fasta)

    output:
    tuple val(meta), path("${meta.sample_id}.qname.bam"), emit: bam
    path "${meta.sample_id}.minimap2.log",                emit: log

    script:
    // -N 10 keeps enough secondary alignments to see the runner-up organism,
    // which is what the assignment margin is computed from. Output stays
    // qname-grouped (unsorted) so the assigner can stream it.
    """
    minimap2 -ax ${params.minimap2_preset} \\
        -t ${task.cpus} \\
        -N 10 --secondary=yes \\
        ${fasta} ${fastq} \\
        2> ${meta.sample_id}.minimap2.log \\
      | samtools view -@ 2 -b -o ${meta.sample_id}.qname.bam -

    grep -q 'ERROR' ${meta.sample_id}.minimap2.log && exit 1 || true
    """
}

process ASSIGN_READS {
    tag   { "${meta.sample_id}:${mode}" }
    label 'analysis'
    label 'process_medium'
    publishDir "${params.outdir}/${meta.sample_id}/${mode}", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(contig_map), val(mode), path(consensus_hits)
    path script

    output:
    tuple val(meta), val(mode), path("${meta.sample_id}.counts.tsv"),        emit: counts
    tuple val(meta), val(mode), path("${meta.sample_id}.readlengths.tsv.gz"), emit: readlengths
    path "${meta.sample_id}.assignments.tsv.gz",                              emit: assignments

    script:
    // Nextflow has no optional path inputs, so the no-consensus case stages a
    // placeholder. Its name is the switch: assign_reads.py sees the flag only
    // when a real hits file was produced upstream.
    def consensus_arg = consensus_hits.name != 'NO_CONSENSUS_HITS'
                      ? "--consensus-hits ${consensus_hits}" : ''
    """
    python3 ${script} ${bam} \\
        --contig-map ${contig_map} \\
        --prefix     ${meta.sample_id} \\
        --min-mapq   ${params.min_mapq} \\
        --mode       ${mode} ${consensus_arg}
    """
}

process EXTRACT_CONTAMINANT_READS {
    tag   { meta.sample_id }
    label 'tools'
    label 'process_medium'

    input:
    tuple val(meta), path(bam), path(contig_map), path(fasta)

    output:
    tuple val(meta), path("contaminant.fastq"), path("contaminant.fasta"), emit: seed

    script:
    // breseq needs two things: the contaminant reference on its own, and a read
    // set to build the consensus from.
    //
    // Feeding it the whole run would be absurd -- ~10M reads through bowtie2 --
    // and pointless, since only reads that resemble the contaminant inform its
    // consensus. Seeding from reads whose PRIMARY alignment lands on the
    // contaminant is the cheap equivalent. The deviation from the original
    // analysis, which ran breseq on everything, is that a read whose primary
    // alignment goes elsewhere but which carries a supplementary alignment to
    // the contaminant does not contribute. Those reads are chimeric or
    // repeat-spanning; they would add little consensus evidence and a lot of
    // noise.
    """
    awk -F'\\t' 'NR>1 && \$3=="contaminant" {print \$1}' ${contig_map} > contaminant_contigs.txt

    if [ ! -s contaminant_contigs.txt ]; then
        echo "error: no contaminant-role contigs in ${contig_map}; --breseq_consensus has nothing to build a consensus of" >&2
        exit 1
    fi

    samtools faidx ${fasta}
    xargs samtools faidx ${fasta} < contaminant_contigs.txt > contaminant.fasta

    # -F 0x900 keeps primary alignments only, so each read appears at most once.
    # `samtools fastq` writes reads with neither READ1 nor READ2 set -- i.e. all
    # of these, the data being single-end -- to stdout, so no -0/-s redirection
    # is given; supplying -0 here would send every read to that file instead.
    # -n keeps the qname verbatim so it still matches the assignment BAM.
    samtools view -h -F 0x900 ${bam} \\
      | awk '
            NR == FNR      { keep[\$1] = 1; next }
            /^@SQ/         { split(\$2, a, ":"); if (a[2] in keep) print; next }
            /^@/           { print; next }
            (\$3 in keep)  { print }
        ' contaminant_contigs.txt - \\
      | samtools fastq -n - > contaminant.fastq

    n=\$(( \$(wc -l < contaminant.fastq) / 4 ))
    if [ "\$n" -eq 0 ]; then
        echo "error: no reads aligned to the contaminant, so no consensus can be built" >&2
        exit 1
    fi

    # Depth gate. breseq calls missing coverage before it calls a consensus, so
    # below roughly 10x it declares the whole reference deleted and gdtools
    # APPLY then writes an empty FASTA -- a failure that surfaces as "empty
    # consensus" several steps downstream and says nothing about the cause.
    # Checking here turns that into a number. Read bases over reference bases
    # is an upper bound on true depth (it ignores soft-clipping), which is the
    # safe direction for a gate: it only ever lets a marginal case through to
    # breseq's own, stricter judgement.
    ref_bases=\$(grep -v '^>' contaminant.fasta | tr -d '\\n' | wc -c | tr -d ' ')
    read_bases=\$(awk 'NR % 4 == 2 { n += length(\$0) } END { printf "%.0f", n }' contaminant.fastq)
    depth=\$(awk -v r="\$read_bases" -v g="\$ref_bases" 'BEGIN { printf "%.2f", r / g }')
    echo "seeding breseq with \$n contaminant reads, \$read_bases bases over \$ref_bases reference bases (~\${depth}x)"

    if awk -v d="\$depth" -v m=${params.breseq_min_depth} 'BEGIN { exit !(d < m) }'; then
        echo "error: contaminant depth ~\${depth}x is below --breseq_min_depth ${params.breseq_min_depth}x." >&2
        echo "       breseq would call the entire reference deleted and the consensus would come back empty." >&2
        echo "       The full replicates in this study reach 20-56x; the bundled test profile reaches ~0.2x," >&2
        echo "       which is why -profile test cannot exercise --breseq_consensus." >&2
        echo "       Lower the gate only if you have a reason to believe breseq can work at this depth." >&2
        exit 1
    fi
    """
}

process BRESEQ_CONSENSUS {
    tag   { meta.sample_id }
    label 'breseq'
    label 'process_medium'
    publishDir "${params.outdir}/${meta.sample_id}/breseq", mode: 'copy',
               pattern: '{consensus.fasta,output.gd,breseq.log,consensus_summary.txt}'

    input:
    tuple val(meta), path(reads), path(reference)

    output:
    tuple val(meta), path("consensus.fasta"), emit: consensus
    path "output.gd",                         emit: gd
    path "consensus_summary.txt",             emit: summary
    path "breseq.log",                        emit: log

    script:
    // -x is breseq's nanopore mode (0.38+): it splits long reads into shorter
    // subsequences so its short-read mapping and mutation calling apply. The
    // documented cost is that indels in homopolymers of 4+ bases are not
    // called, which for ONT data is the right trade -- those are exactly the
    // positions where the basecaller is least trustworthy, and calling them
    // would write basecalling error into the consensus.
    //
    // --no-junction-prediction: we want a corrected consensus sequence, not a
    // structural-variant catalogue. Junction prediction is the expensive part
    // and gdtools APPLY does not use its output.
    """
    breseq -x \\
        --no-junction-prediction \\
        -j ${task.cpus} \\
        -n ${meta.sample_id} \\
        -r ${reference} \\
        -o breseq_out \\
        ${reads} > breseq.log 2>&1 || {
            echo "error: breseq failed; tail of log follows" >&2
            tail -50 breseq.log >&2
            exit 1
        }

    if [ ! -s breseq_out/output/output.gd ]; then
        echo "error: breseq produced no output.gd" >&2
        tail -50 breseq.log >&2
        exit 1
    fi
    cp breseq_out/output/output.gd output.gd

    # breseq predicts missing coverage before it predicts a consensus. Given too
    # little data it emits a DEL spanning the whole reference, and gdtools APPLY
    # then faithfully deletes the genome. Catching that here names the cause;
    # left to APPLY it surfaces only as a zero-byte FASTA.
    ref_bases=\$(grep -v '^>' ${reference} | tr -d '\\n' | wc -c | tr -d ' ')
    if awk -v g="\$ref_bases" '\$1 == "DEL" && \$6 >= g { found = 1 } END { exit !found }' output.gd; then
        echo "error: breseq called the entire contaminant reference deleted (a DEL spanning all \$ref_bases bp)." >&2
        echo "       That is its missing-coverage prediction, not a consensus: there was too little" >&2
        echo "       evidence for it to call one. breseq saw:" >&2
        grep -E '^#=(INPUT|CONVERTED|MAPPED)-' output.gd >&2
        exit 1
    fi

    gdtools APPLY -r ${reference} -f FASTA -o consensus.fasta output.gd

    if [ ! -s consensus.fasta ]; then
        echo "error: gdtools APPLY produced an empty consensus despite no whole-reference deletion. Inspect output.gd." >&2
        exit 1
    fi

    # A consensus identical to the reference means the whole step changed
    # nothing, which is worth stating plainly in the results rather than
    # leaving the reader to diff two FASTAs.
    {
        echo "sample\t${meta.sample_id}"
        echo "reference_bases\t\$ref_bases"
        echo "consensus_bases\t\$(grep -v '^>' consensus.fasta | tr -d '\\n' | wc -c | tr -d ' ')"
        echo "mutations_applied\t\$(awk '\$1 ~ /^(SNP|SUB|DEL|INS|MOB|AMP|CON|INV)\$/' output.gd | wc -l | tr -d ' ')"
        grep -E '^#=(INPUT|CONVERTED|MAPPED)-' output.gd | sed 's/^#=//;s/ /\t/'
    } > consensus_summary.txt
    """
}

process MAP_CONSENSUS {
    tag   { meta.sample_id }
    label 'tools'
    label 'process_high'
    publishDir "${params.outdir}/${meta.sample_id}/breseq", mode: 'copy',
               pattern: 'consensus_hits.tsv'

    input:
    tuple val(meta), path(fastq), path(consensus)

    output:
    tuple val(meta), path("consensus_hits.txt"), emit: hits
    path "consensus_hits.tsv",                   emit: stats

    script:
    // Every read in the run is tested against the consensus, not just the ones
    // that already aligned to the stock contaminant reference. That is the
    // point of the exercise: the consensus is supposed to catch reads the stock
    // reference misses, and restricting the test to reads it already caught
    // would guarantee it never does.
    //
    // -F 0x904 drops unmapped, secondary and supplementary records, so each
    // surviving read contributes its primary alignment once.
    """
    minimap2 -ax ${params.minimap2_preset} -t ${task.cpus} \\
        ${consensus} ${fastq} 2> minimap2.log \\
      | samtools view -F 0x904 -q ${params.min_mapq} - \\
      | cut -f1 \\
      | sort -u > consensus_hits.txt

    grep -q 'ERROR' minimap2.log && exit 1 || true

    if [ ! -s consensus_hits.txt ]; then
        echo "error: no reads matched the breseq consensus. An empty subtraction is indistinguishable from a broken step, so this is an error rather than a silent no-op." >&2
        exit 1
    fi

    {
        echo "sample\t${meta.sample_id}"
        echo "consensus_hits\t\$(wc -l < consensus_hits.txt | tr -d ' ')"
        echo "min_mapq\t${params.min_mapq}"
    } > consensus_hits.tsv
    """
}

process COMPUTE_METRICS {
    tag   { "${meta.sample_id}:${mode}" }
    label 'analysis'
    label 'process_low'
    publishDir "${params.outdir}/${meta.sample_id}/${mode}", mode: 'copy'

    input:
    tuple val(meta), val(mode), path(counts), path(reference_tsv)
    path script

    output:
    tuple val(meta), val(mode), path("${meta.sample_id}.metrics.tsv"), emit: metrics
    path "${meta.sample_id}.summary.json",                             emit: summary

    script:
    def lib = meta.library_dna_ng ?: ''
    """
    python3 ${script} \\
        --counts         ${counts} \\
        --reference-tsv  ${reference_tsv} \\
        --sample-id      ${meta.sample_id} \\
        --experiment     ${meta.experiment} \\
        --replicate      ${meta.replicate} \\
        --library-dna-ng '${lib}' \\
        --carrier-dna-ng '${meta.carrier_dna_ng}' \\
        --mode           ${mode} \\
        --out-tsv        ${meta.sample_id}.metrics.tsv \\
        --out-json       ${meta.sample_id}.summary.json
    """
}

process COVERAGE_PROFILE {
    tag   { meta.sample_id }
    label 'tools'
    label 'process_medium'
    publishDir "${params.outdir}/${meta.sample_id}/coverage", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(contig_map)

    output:
    tuple val(meta), path("${meta.sample_id}.depth.tsv.gz"), emit: depth

    script:
    // Only community (role=sample) contigs are profiled: carrier depth is
    // uninformative and the carrier is >95% of the reads.
    //
    // The filter is applied to the SAM *stream* rather than by passing regions
    // to `samtools view`. Region arguments require a coordinate-sorted, indexed
    // BAM, but this BAM comes qname-grouped straight from minimap2, so the
    // region form fails and any fallback ends up sorting all ~10M reads --
    // which made this the single slowest step in the pipeline (43-58 min per
    // replicate, versus 7-55 min for the mapping that produced the file).
    // Filtering first means `samtools sort` only ever sees community reads.
    //
    // Non-sample @SQ lines are dropped from the header too, so `depth -a`
    // enumerates only community genomes instead of every reference present.
    """
    awk -F'\\t' 'NR>1 && \$3=="sample" {print \$1}' ${contig_map} > sample_contigs.txt

    if [ ! -s sample_contigs.txt ]; then
        echo "error: no sample-role contigs in ${contig_map}" >&2; exit 1
    fi

    samtools view -h -F 0x900 ${bam} \\
      | awk '
            NR == FNR                { keep[\$1] = 1; next }
            /^@SQ/                   { split(\$2, a, ":"); if (a[2] in keep) print; next }
            /^@/                     { print; next }
            (\$3 in keep)            { print }
        ' sample_contigs.txt - \\
      | samtools sort -@ ${task.cpus} -m 1G -o sorted.bam -

    samtools index sorted.bam
    samtools depth -a -@ ${task.cpus} sorted.bam | gzip -c > ${meta.sample_id}.depth.tsv.gz
    rm -f sorted.bam sorted.bam.bai
    """
}

process COVERAGE_SUMMARY {
    tag   { meta.sample_id }
    label 'analysis'
    label 'process_medium'
    publishDir "${params.outdir}/${meta.sample_id}/coverage", mode: 'copy'

    input:
    tuple val(meta), path(depth), path(contig_map)
    path script

    output:
    path "${meta.sample_id}.coverage_summary.tsv", emit: summary
    path "${meta.sample_id}.coverage_profile.tsv", emit: profile

    script:
    """
    python3 ${script} ${depth} \\
        --contig-map  ${contig_map} \\
        --sample-id   ${meta.sample_id} \\
        --window      ${params.coverage_window} \\
        --out-summary ${meta.sample_id}.coverage_summary.tsv \\
        --out-profile ${meta.sample_id}.coverage_profile.tsv
    """
}

process AGGREGATE {
    label 'analysis'
    label 'process_medium'
    publishDir "${params.outdir}/summary", mode: 'copy'

    input:
    // With --mode both, competitive and sequential emit identically named files
    // for the same sample (they are kept apart by publishDir, not by filename).
    // Collecting them into one process would therefore collide, so each file is
    // staged into its own numbered subdirectory. Basenames are preserved, which
    // matters because aggregate_results.py recovers the sample id from them.
    path metrics,     stageAs: 'metrics*/*'
    path summaries,   stageAs: 'summaries*/*'
    path readlengths, stageAs: 'readlengths*/*'
    path measurements
    path samplesheet
    path script

    output:
    path "per_organism.tsv",       emit: per_organism
    path "per_sample.tsv",         emit: per_sample
    path "experiment_summary.tsv", optional: true, emit: experiment_summary
    path "*.pdf",                  optional: true, emit: figures
    path "*.png",                  optional: true
    path "*.csv",                  optional: true
    path "*.json",                 optional: true

    script:
    """
    python3 ${script} \\
        --metrics     ${metrics} \\
        --summaries   ${summaries} \\
        --readlengths ${readlengths} \\
        --measurements ${measurements} \\
        --samplesheet ${samplesheet} \\
        --outdir      .
    """
}

// ---------------------------------------------------------------------------
// Samplesheet parsing
// ---------------------------------------------------------------------------

// Quote-aware split of one delimited line. Nextflow's splitCsv is a channel
// operator, and parsing inside an operator would push validation onto a
// dataflow thread where error messages are swallowed (see parseSamplesheet).
// Ten lines here buy eager parsing and diagnostics that actually reach the user.
def splitDelimited(String line, char delim) {
    def out = []
    def cur = new StringBuilder()
    boolean inQuotes = false
    for (int i = 0; i < line.length(); i++) {
        char c = line.charAt(i)
        if (c == '"' as char) {
            if (inQuotes && i + 1 < line.length() && line.charAt(i + 1) == ('"' as char)) {
                cur.append('"'); i++
            } else {
                inQuotes = !inQuotes
            }
        } else if (c == delim && !inQuotes) {
            out << cur.toString(); cur = new StringBuilder()
        } else {
            cur.append(c)
        }
    }
    out << cur.toString()
    return out
}

// Rows of a delimited file as maps, with `#` comments and blank lines dropped.
def readTable(path, char delim) {
    def lines = file(path).readLines()
                    .findAll { l -> l.trim() && !l.trim().startsWith('#') }
    if (!lines) return []
    def cols = splitDelimited(lines[0], delim)*.trim()
    return lines.drop(1).collect { l ->
        def vals = splitDelimited(l, delim)
        def row = [:]
        cols.eachWithIndex { c, i -> row[c] = (i < vals.size() ? vals[i].trim() : '') }
        return row
    }
}

def loadMeasurements(path) {
    // Experimental masses live in assets/measurements.tsv, not in the
    // samplesheets. A samplesheet mixes local file paths with experimental
    // facts, and keeping a second copy of a mass next to a path is how a stale
    // carrier value once survived beside a Methods section that contradicted
    // it. This file is the single source for every quantity a per-femtogram or
    // enrichment number divides by; `make measurements` checks it.
    if (!file(path).exists()) return [:]
    def out = [:]
    readTable(path, '\t' as char).each { row ->
        if (!row.sample_id) return
        out[row.sample_id] = [
            sample_dna_ng:       row.sample_dna_ng,
            carrier_dna_ng:      row.carrier_dna_ng,
            include_in_headline: row.include_in_headline == '1',
        ]
    }
    return out
}

def parseSamplesheet(path) {
    // The samplesheets carry a `#` comment preamble documenting every column;
    // readTable drops it, so the first non-comment line is the header.
    file(path, checkIfExists: true)
    def rows = readTable(path, ',' as char)
    if (!rows) error "${path} contains no data rows"

    if (params.fetch_from_sra) {
        error "--fetch_from_sra is not implemented yet; see docs/TODO.md"
    }

    def measured = loadMeasurements(params.measurements)

    // Validation happens HERE, in the function body, and not inside a channel
    // operator. `error "message"` raised from inside a `.map{}` closure runs on
    // a dataflow thread: the process does exit non-zero, but the message is
    // swallowed and the user gets a bare exit code with no diagnostic. Every
    // check below used to live in that closure and was therefore silent.
    def records = []
    rows.each { row ->
        if (!row.sample_id) return

        // measurements.tsv is authoritative for every sample it names. The
        // synthetic smoke-test sheet and the reanalysed prior-study sheet are not
        // experiments of ours and carry their masses inline; those fall through
        // to the samplesheet columns.
        def m = measured[row.sample_id]
        if (m) {
            [library_dna_ng: 'sample_dna_ng', carrier_dna_ng: 'carrier_dna_ng'].each { sheetCol, measCol ->
                def v = row[sheetCol]?.trim()
                if (v && v != m[measCol]) {
                    error "${row.sample_id}: ${path} says ${sheetCol}=${v} but " +
                            "${params.measurements} says ${measCol}=${m[measCol]}. Masses belong " +
                            "in measurements.tsv only; remove the column from the samplesheet."
                }
            }
        }

        def meta = [
            sample_id:           row.sample_id,
            experiment:          row.experiment,
            replicate:           row.replicate,
            reference_set:       row.reference_set,
            library_dna_ng:      m ? m.sample_dna_ng  : row.library_dna_ng?.trim(),
            carrier_dna_ng:      (m ? m.carrier_dna_ng : row.carrier_dna_ng?.trim()) ?: '0',
            include_in_headline: m ? m.include_in_headline
                                   : (row.include_in_headline?.trim() ?: '1') == '1',
            sra_accession:       row.sra_accession?.trim(),
        ]

        // A mass still marked PENDING is not a number; treat it as absent so
        // downstream code emits blanks rather than trying to parse "PENDING".
        if (meta.library_dna_ng == 'PENDING') meta.library_dna_ng = ''
        if (meta.carrier_dna_ng == 'PENDING') meta.carrier_dna_ng = '0'

        if (!row.fastq?.trim()) {
            error "${meta.sample_id} has no fastq path. Reads are not yet " +
                    "deposited, so a local path is required (see docs/TODO.md)."
        }
        def fq = file(row.fastq.trim())
        if (!fq.isAbsolute()) fq = file("${projectDir}/${row.fastq.trim()}")
        if (!fq.exists()) error "FASTQ not found for ${meta.sample_id}: ${fq}"

        if (!row.reference_set?.trim()) {
            error "${meta.sample_id} has no reference_set"
        }
        def ref = file(row.reference_set.trim())
        if (!ref.isAbsolute()) ref = file("${projectDir}/${row.reference_set.trim()}")
        if (!ref.exists()) error "reference set not found: ${ref}"

        records << tuple(meta, fq, ref)
    }

    if (!records) error "${path} yielded no usable samples"
    Channel.fromList(records)
}

// ---------------------------------------------------------------------------
// Workflow
// ---------------------------------------------------------------------------

workflow {

    log.info """
    ================================================================
     low-input-nanopore ${workflow.manifest.version}
     samplesheet : ${params.samplesheet}
     mode        : ${params.mode}
     outdir      : ${params.outdir}
    ================================================================
    """.stripIndent()

    // Analysis scripts are staged as explicit process inputs rather than relying
    // on Nextflow's bin/ PATH injection, which breaks when the project path
    // contains spaces, as it does on the reference platform.
    ch_script_build   = file("${projectDir}/bin/build_reference_set.py", checkIfExists: true)
    ch_script_assign  = file("${projectDir}/bin/assign_reads.py",        checkIfExists: true)
    ch_script_metrics = file("${projectDir}/bin/compute_metrics.py",     checkIfExists: true)
    ch_script_agg     = file("${projectDir}/bin/aggregate_results.py",   checkIfExists: true)
    ch_script_cov     = file("${projectDir}/bin/coverage_summary.py",    checkIfExists: true)
    ch_script_qfilt   = file("${projectDir}/bin/filter_by_qscore.awk",   checkIfExists: true)

    // Placeholder for the optional consensus-hits input. Nextflow process
    // inputs are positional and non-optional, so the absence of a real file has
    // to be represented by a real file whose name says so.
    ch_no_consensus   = file("${projectDir}/assets/NO_CONSENSUS_HITS", checkIfExists: true)

    ch_samples = parseSamplesheet(params.samplesheet)

    // One reference build per distinct reference set, shared by its samples.
    ch_refsets = ch_samples
        .map { meta, fq, ref -> tuple(ref.simpleName, ref) }
        .unique { it[0] }

    FETCH_GENOMES(ch_refsets)

    BUILD_REFERENCE(
        ch_refsets.join(FETCH_GENOMES.out.genomes).map { name, tsv, genomes ->
            tuple(name, tsv, genomes)
        },
        ch_script_build
    )

    // Re-key samples by reference-set name so each joins its built reference.
    ch_keyed = ch_samples.map { meta, fq, ref -> tuple(ref.simpleName, meta, fq, ref) }

    // Optional mean-read-quality filter. Off by default (min_qscore = 0) so the
    // primary analysis uses every basecalled read. Set --min_qscore 10 to match
    // the wf-metagenomics --min_read_qual 10 used when reanalysing the prior
    // studies, which is what makes a cross-study comparison quality-matched.
    if (params.min_qscore > 0) {
        FILTER_READS(ch_keyed.map { name, meta, fq, ref -> tuple(meta, fq) },
                     ch_script_qfilt)
        ch_reads = FILTER_READS.out.reads
    } else {
        ch_reads = ch_keyed.map { name, meta, fq, ref -> tuple(meta, fq) }
    }

    ch_to_map = ch_keyed
        .map { name, meta, fq, ref -> tuple(meta.sample_id, name, meta) }
        .join(ch_reads.map { meta, fq -> tuple(meta.sample_id, fq) })
        .map { sid, name, meta, fq -> tuple(name, meta, fq) }
        .combine(BUILD_REFERENCE.out.fasta, by: 0)
        .map { name, meta, fq, fasta -> tuple(meta, fq, fasta) }

    MAP_COMPETITIVE(ch_to_map)

    ch_base_assign = ch_keyed
        .combine(BUILD_REFERENCE.out.contig_map, by: 0)
        .map { name, meta, fq, ref, cmap -> tuple(meta.sample_id, meta, cmap) }
        .join(MAP_COMPETITIVE.out.bam.map { meta, bam -> tuple(meta.sample_id, bam) })
        .map { sid, meta, cmap, bam -> tuple(meta, bam, cmap) }
        // Both assignment rules read the same alignments, so `--mode both` costs
        // one extra cheap pass over the BAM rather than a second mapping run.
        .combine(Channel.fromList(
            params.mode == 'both' ? ['competitive', 'sequential'] : [params.mode]))

    // ---- optional breseq contaminant consensus -----------------------------
    // Subtract against the E. coli actually present in the carrier prep rather
    // than against the stock MG1655 reference, which is what the original
    // lowinput_s1 analysis did. Only sequential mode subtracts anything, so
    // competitive rows carry the placeholder even when this is on.
    if (params.breseq_consensus) {
        EXTRACT_CONTAMINANT_READS(
            ch_keyed
                .combine(BUILD_REFERENCE.out.contig_map, by: 0)
                .map { name, meta, fq, ref, cmap -> tuple(name, meta, cmap) }
                .combine(BUILD_REFERENCE.out.fasta, by: 0)
                .map { name, meta, cmap, fasta -> tuple(meta.sample_id, meta, cmap, fasta) }
                .join(MAP_COMPETITIVE.out.bam.map { meta, bam -> tuple(meta.sample_id, bam) })
                .map { sid, meta, cmap, fasta, bam -> tuple(meta, bam, cmap, fasta) }
        )

        BRESEQ_CONSENSUS(EXTRACT_CONTAMINANT_READS.out.seed)

        MAP_CONSENSUS(
            ch_reads
                .map { meta, fq -> tuple(meta.sample_id, meta, fq) }
                .join(BRESEQ_CONSENSUS.out.consensus
                          .map { meta, cons -> tuple(meta.sample_id, cons) })
                .map { sid, meta, fq, cons -> tuple(meta, fq, cons) }
        )

        ch_to_assign = ch_base_assign
            .map { meta, bam, cmap, mode -> tuple(meta.sample_id, meta, bam, cmap, mode) }
            .combine(MAP_CONSENSUS.out.hits
                         .map { meta, hits -> tuple(meta.sample_id, hits) }, by: 0)
            .map { sid, meta, bam, cmap, mode, hits ->
                   tuple(meta, bam, cmap, mode,
                         mode == 'sequential' ? hits : ch_no_consensus) }
    } else {
        ch_to_assign = ch_base_assign
            .map { meta, bam, cmap, mode -> tuple(meta, bam, cmap, mode, ch_no_consensus) }
    }

    ASSIGN_READS(ch_to_assign, ch_script_assign)

    // combine(by:0) rather than join(): with --mode both there are two rows per
    // sample_id, and join() does not handle duplicate keys on one side.
    ch_to_metrics = ASSIGN_READS.out.counts
        .map { meta, mode, counts -> tuple(meta.sample_id, meta, mode, counts) }
        .combine(ch_keyed.map { name, meta, fq, ref -> tuple(meta.sample_id, ref) }, by: 0)
        .map { sid, meta, mode, counts, ref -> tuple(meta, mode, counts, ref) }

    COMPUTE_METRICS(ch_to_metrics, ch_script_metrics)

    ch_to_coverage = ch_keyed
        .combine(BUILD_REFERENCE.out.contig_map, by: 0)
        .map { name, meta, fq, ref, cmap -> tuple(meta.sample_id, meta, cmap) }
        .join(MAP_COMPETITIVE.out.bam.map { meta, bam -> tuple(meta.sample_id, bam) })
        .map { sid, meta, cmap, bam -> tuple(meta, bam, cmap) }

    COVERAGE_PROFILE(ch_to_coverage)

    COVERAGE_SUMMARY(
        COVERAGE_PROFILE.out.depth
            .map { meta, depth -> tuple(meta.sample_id, meta, depth) }
            .combine(ch_keyed.map { name, meta, fq, ref -> tuple(meta.sample_id, name) }, by: 0)
            .map { sid, meta, depth, name -> tuple(name, meta, depth) }
            .combine(BUILD_REFERENCE.out.contig_map, by: 0)
            .map { name, meta, depth, cmap -> tuple(meta, depth, cmap) },
        ch_script_cov
    )

    // Study-level tables and display items, once every sample has finished.
    AGGREGATE(
        COMPUTE_METRICS.out.metrics.map { meta, mode, m -> m }.collect(),
        COMPUTE_METRICS.out.summary.collect(),
        // One mode only. With --mode both the same reads appear twice (relabelled
        // by the other rule), which would double-count every read in the
        // length distribution. Competitive is the primary rule.
        ASSIGN_READS.out.readlengths
            .filter { meta, mode, rl -> mode == (params.mode == 'sequential' ? 'sequential' : 'competitive') }
            .map { meta, mode, rl -> rl }
            .collect(),
        file(params.measurements, checkIfExists: true),
        file(params.samplesheet),
        ch_script_agg
    )
}

workflow.onComplete {
    log.info """
    ================================================================
     ${workflow.success ? 'Completed' : 'FAILED'} in ${workflow.duration}
     results: ${params.outdir}
    ================================================================
    """.stripIndent()
}
