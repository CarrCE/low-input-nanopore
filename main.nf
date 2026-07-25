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
    tuple val(meta), path(bam), path(contig_map), val(mode)
    path script

    output:
    tuple val(meta), val(mode), path("${meta.sample_id}.counts.tsv"),        emit: counts
    tuple val(meta), val(mode), path("${meta.sample_id}.readlengths.tsv.gz"), emit: readlengths
    path "${meta.sample_id}.assignments.tsv.gz",                              emit: assignments

    script:
    """
    python3 ${script} ${bam} \\
        --contig-map ${contig_map} \\
        --prefix     ${meta.sample_id} \\
        --min-mapq   ${params.min_mapq} \\
        --mode       ${mode}
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
    label 'process_low'
    publishDir "${params.outdir}/summary", mode: 'copy'

    input:
    path metrics
    path summaries
    path readlengths
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
        --samplesheet ${samplesheet} \\
        --outdir      .
    """
}

// ---------------------------------------------------------------------------
// Samplesheet parsing
// ---------------------------------------------------------------------------

def parseSamplesheet(path) {
    // The samplesheets carry a `#` comment preamble documenting every column.
    // splitCsv would take the first comment line as the header, so strip
    // comments and blank lines before parsing.
    def sheet = file(path, checkIfExists: true)
    def body = sheet.readLines()
                    .findAll { line -> line.trim() && !line.trim().startsWith('#') }
                    .join('\n')
    if (!body) exit 1, "error: ${path} contains no data rows"

    Channel
        .of(body)
        .splitCsv(header: true, strip: true, sep: ',')
        .filter { row -> row.sample_id }
        .map { row ->
            def meta = [
                sample_id:           row.sample_id,
                experiment:          row.experiment,
                replicate:           row.replicate,
                reference_set:       row.reference_set,
                library_dna_ng:      row.library_dna_ng?.trim(),
                carrier_dna_ng:      (row.carrier_dna_ng?.trim() ?: '0'),
                include_in_headline: (row.include_in_headline?.trim() ?: '1') == '1',
                sra_accession:       row.sra_accession?.trim(),
            ]

            if (params.fetch_from_sra) {
                exit 1, "error: --fetch_from_sra is not implemented yet; see docs/TODO.md"
            }
            if (!row.fastq?.trim()) {
                exit 1, "error: ${meta.sample_id} has no fastq path. Reads are not yet " +
                        "deposited, so a local path is required (see docs/TODO.md)."
            }

            def fq = file(row.fastq.trim())
            if (!fq.isAbsolute()) fq = file("${projectDir}/${row.fastq.trim()}")
            if (!fq.exists()) exit 1, "error: FASTQ not found for ${meta.sample_id}: ${fq}"

            def ref = file(row.reference_set.trim())
            if (!ref.isAbsolute()) ref = file("${projectDir}/${row.reference_set.trim()}")
            if (!ref.exists()) exit 1, "error: reference set not found: ${ref}"

            tuple(meta, fq, ref)
        }
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

    ch_to_map = ch_keyed
        .combine(BUILD_REFERENCE.out.fasta, by: 0)
        .map { name, meta, fq, ref, fasta -> tuple(meta, fq, fasta) }

    MAP_COMPETITIVE(ch_to_map)

    ch_to_assign = ch_keyed
        .combine(BUILD_REFERENCE.out.contig_map, by: 0)
        .map { name, meta, fq, ref, cmap -> tuple(meta.sample_id, meta, cmap) }
        .join(MAP_COMPETITIVE.out.bam.map { meta, bam -> tuple(meta.sample_id, bam) })
        .map { sid, meta, cmap, bam -> tuple(meta, bam, cmap) }
        // Both assignment rules read the same alignments, so `--mode both` costs
        // one extra cheap pass over the BAM rather than a second mapping run.
        .combine(Channel.fromList(
            params.mode == 'both' ? ['competitive', 'sequential'] : [params.mode]))

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
        ASSIGN_READS.out.readlengths.map { meta, mode, rl -> rl }.collect(),
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
