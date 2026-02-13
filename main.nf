#!/usr/bin/env nextflow

/*
 * BaitBench - Probe Capture Efficiency Testing Pipeline
 *
 * Simulates probe-based capture and evaluates detection performance
 * using TP/FP/FN/TN metrics.
 *
 * Usage:
 *   nextflow run main.nf \
 *     --targets <targets.fa> \
 *     --distractors <distractors.fa> \
 *     --probes <probes.fa> \
 *     --outdir <output_directory>
 */

nextflow.enable.dsl=2

// Validate required parameters
if (!params.targets) {
    error "Please specify --targets (FASTA file with target sequences)"
}
if (!params.distractors) {
    error "Please specify --distractors (FASTA file with distractor sequences)"
}
if (!params.probes) {
    error "Please specify --probes (FASTA file with probe sequences)"
}

// Generate run name if not provided
def run_name = params.run_name ?: "run_${new Date().format('yyyyMMdd_HHmmss')}"

log.info """
    =============================================
    BaitBench - Probe Capture Testing Pipeline
    =============================================
    Run name            : ${run_name}
    Targets FASTA       : ${params.targets}
    Distractors FASTA   : ${params.distractors}
    Probe FASTA         : ${params.probes}
    Capture method      : ${params.capture_method}
    ${params.capture_method == 'blast' ? "BLAST DB dir        : ${params.blast_db_dir}" : ""}
    Host FASTA          : ${params.host_fasta ?: 'none (skip host filtering)'}
    Num reads           : ${params.num_reads}
    Distractor fraction : ${params.distractor_fraction}
    Max mismatches      : ${params.max_mismatches}
    Min match bases     : ${params.min_match_bases}
    Seed                : ${params.seed ?: 'random'}
    Output dir          : ${params.outdir}
    =============================================
""".stripIndent()


/*
 * Step 1: Prepare combined reference and weights
 */
process PREPARE_REFERENCE {
    tag "${run_name}"
    publishDir "${params.outdir}/${run_name}", mode: 'copy', pattern: '*.txt'

    input:
    path targets
    path distractors
    val distractor_fraction

    output:
    path "combined_reference.fa", emit: reference
    path "weights.txt", emit: weights
    path "targets.txt", emit: targets_list
    path "distractors.txt", emit: distractors_list

    script:
    """
    python ${projectDir}/bin/prepare_reference.py \\
        --targets ${targets} \\
        --distractors ${distractors} \\
        --output-fasta combined_reference.fa \\
        --output-weights weights.txt \\
        --distractor-fraction ${distractor_fraction} \\
        --targets-list targets.txt \\
        --distractors-list distractors.txt
    """
}


/*
 * Step 2: Generate simulated reads from reference genomes
 */
process GENERATE_READS {
    tag "${run_name}"
    publishDir "${params.outdir}/${run_name}", mode: 'copy'

    input:
    path reference_fasta
    path weights_file
    val num_reads
    val seed

    output:
    path "reads.fa", emit: reads
    path "generate_reads.log", emit: log

    script:
    def seed_arg = seed ? "-s ${seed}" : ""
    """
    python ${projectDir}/bin/fasta_sampler.py \\
        -f ${reference_fasta} \\
        -w ${weights_file} \\
        -n ${num_reads} \\
        -o reads.fa \\
        ${seed_arg} \\
        2>&1 | tee generate_reads.log
    """
}


/*
 * Step 3a: In silico capture using minimap2
 */
process CAPTURE {
    tag "${run_name}"
    publishDir "${params.outdir}/${run_name}", mode: 'copy'

    input:
    path probe_fasta
    path reads
    val max_mismatches
    val min_match_bases

    output:
    path "captured.fa", emit: captured
    path "capture.log", emit: log

    script:
    """
    # Map reads to probes (PAF output with CIGAR)
    minimap2 -x sr -c --cs -A 4 -B 2 -O 12,32 --secondary=yes ${probe_fasta} ${reads} > mappings.paf 2> capture.log

    # Filter PAF: no indels, ≤max_mismatches mismatches, ≥min_match_bases matching
    awk -F'\\t' -v max_mm=${max_mismatches} -v min_match=${min_match_bases} '
    {
        # Extract NM (mismatches) from tags
        nm = 0
        for (i=13; i<=NF; i++) {
            if (\$i ~ /^NM:i:/) {
                split(\$i, a, ":")
                nm = a[3]
            }
        }

        # Check for indels in CIGAR (cg tag)
        has_indel = 0
        for (i=13; i<=NF; i++) {
            if (\$i ~ /^cg:Z:/) {
                if (\$i ~ /[ID]/) has_indel = 1
            }
        }

        # \$10 = number of matching bases
        matches = \$10

        # Apply filters
        if (has_indel == 0 && nm <= max_mm && matches >= min_match) {
            print \$1
        }
    }' mappings.paf | sort -u > passing_reads.txt

    # Extract passing reads from FASTA
    seqtk subseq ${reads} passing_reads.txt > captured.fa

    # Log statistics
    echo "Total mappings: \$(wc -l < mappings.paf)" >> capture.log
    echo "Reads passing filters: \$(wc -l < passing_reads.txt)" >> capture.log
    echo "Captured sequences: \$(grep -c '^>' captured.fa || echo 0)" >> capture.log
    """
}


/*
 * Step 3b: In silico capture using BLAST (alternative)
 */
process BLAST_CAPTURE {
    tag "${run_name}"
    publishDir "${params.outdir}/${run_name}", mode: 'copy'

    input:
    val blast_db
    path reads
    val min_match_bases

    output:
    path "captured.fa", emit: captured
    path "capture.log", emit: log

    script:
    """
    # Run BLAST alignment
    blastn -task blastn-short \\
        -db ${blast_db} \\
        -query ${reads} \\
        -out blast_results.tsv \\
        -outfmt "6 qseqid sseqid nident length mismatch gapopen" \\
        -max_hsps 1 \\
        -max_target_seqs 1 \\
        -evalue 1000 \\
        -word_size 7 \\
        -dust no \\
        -soft_masking false \\
        -num_threads ${task.cpus} \\
        2> capture.log

    # Filter BLAST results
    awk -F'\\t' -v min_match=${min_match_bases} '
        !seen[\$1] && \$6 == 0 && \$3 >= min_match {
            print \$1
            seen[\$1] = 1
        }
    ' blast_results.tsv > passing_reads.txt

    # Extract passing reads from FASTA
    seqtk subseq ${reads} passing_reads.txt > captured.fa

    # Log statistics
    echo "Total BLAST hits: \$(wc -l < blast_results.tsv)" >> capture.log
    echo "Reads passing filters: \$(wc -l < passing_reads.txt)" >> capture.log
    echo "Captured sequences: \$(grep -c '^>' captured.fa || echo 0)" >> capture.log
    """
}


/*
 * Step 4: Filter out host reads (optional)
 */
process FILTER_HOST {
    tag "${run_name}"
    publishDir "${params.outdir}/${run_name}", mode: 'copy'

    input:
    path host_fasta
    path captured_reads
    val minimap_preset

    output:
    path "filtered.fa", emit: filtered
    path "host_filter.log", emit: log

    script:
    """
    # Map captured reads to host genome
    minimap2 -ax ${minimap_preset} ${host_fasta} ${captured_reads} > host_mapped.sam 2> host_filter.log

    # Extract read IDs that map to host
    grep -v '^@' host_mapped.sam | awk '\$3 != "*"' | cut -f1 | sort -u > host_reads.txt

    # Get all captured read IDs
    grep '^>' ${captured_reads} | sed 's/^>//' | cut -d' ' -f1 > all_reads.txt

    # Get reads that don't map to host
    comm -23 <(sort all_reads.txt) <(sort host_reads.txt) > passing_reads.txt

    # Extract non-host reads
    seqtk subseq ${captured_reads} passing_reads.txt > filtered.fa

    # Log statistics
    echo "Input sequences: \$(grep -c '^>' ${captured_reads} || echo 0)" >> host_filter.log
    echo "Host matches: \$(wc -l < host_reads.txt)" >> host_filter.log
    echo "Filtered sequences: \$(grep -c '^>' filtered.fa || echo 0)" >> host_filter.log
    """
}


/*
 * Step 5: Map reads to reference genomes
 */
process MAP_READS {
    tag "${run_name}"
    publishDir "${params.outdir}/${run_name}", mode: 'copy'

    input:
    path reference_fasta
    path captured_reads
    val minimap_preset

    output:
    path "mapped.sam", emit: sam

    script:
    """
    minimap2 -ax ${minimap_preset} --secondary=no \\
        ${reference_fasta} \\
        ${captured_reads} \\
        > mapped.sam
    """
}


/*
 * Step 6: Generate detection list from SAM file
 */
process GENERATE_LIST {
    tag "${run_name}"
    publishDir "${params.outdir}/${run_name}", mode: 'copy'

    input:
    path sam_file

    output:
    path "detected.list", emit: detected_list

    script:
    """
    grep -v '^@' ${sam_file} | \\
        awk -F'\\t' '
            \$3 == "*" { next }
            { flag = \$2 + 0 }
            (int(flag / 256) % 2 == 1) { next }
            (int(flag / 2048) % 2 == 1) { next }
            !seen[\$1]++ { print \$3 }
        ' | \\
        sort | \\
        uniq -c | \\
        sort -n -k1 | \\
        awk '{print \$2 "\\t" \$1}' > detected.list
    """
}


/*
 * Step 7: Calculate metrics
 */
process CALCULATE_METRICS {
    tag "${run_name}"
    publishDir "${params.outdir}/${run_name}", mode: 'copy'

    input:
    path targets_list
    path distractors_list
    path detected_list
    path reads
    path captured
    val run_name
    val num_reads
    val seed

    output:
    path "results.tsv", emit: results
    path "detected_detail.tsv", emit: detail
    path "results.json", emit: json

    script:
    def seed_str = seed ?: "NA"
    """
    python ${projectDir}/bin/metrics.py \\
        --targets ${targets_list} \\
        --distractors ${distractors_list} \\
        --detected ${detected_list} \\
        --reads ${reads} \\
        --captured ${captured} \\
        --run-name "${run_name}" \\
        --num-reads ${num_reads} \\
        --seed "${seed_str}" \\
        --output-summary results.tsv \\
        --output-detail detected_detail.tsv \\
        --output-json results.json
    """
}


/*
 * Step 8: Generate HTML report
 */
process GENERATE_REPORT {
    tag "${run_name}"
    publishDir "${params.outdir}/${run_name}", mode: 'copy'

    input:
    path results_tsv
    path detail_tsv
    val run_name

    output:
    path "report.html", emit: report

    script:
    """
    python ${projectDir}/bin/generate_report.py \\
        --summary ${results_tsv} \\
        --detail ${detail_tsv} \\
        --output report.html \\
        --run-name "${run_name}"
    """
}


/*
 * Main workflow
 */
workflow {
    // Create channels from input files
    targets_ch = Channel.fromPath(params.targets, checkIfExists: true)
    distractors_ch = Channel.fromPath(params.distractors, checkIfExists: true)
    probe_ch = Channel.fromPath(params.probes, checkIfExists: true)

    // Step 1: Prepare combined reference
    PREPARE_REFERENCE(
        targets_ch,
        distractors_ch,
        params.distractor_fraction
    )

    // Step 2: Generate reads
    GENERATE_READS(
        PREPARE_REFERENCE.out.reference,
        PREPARE_REFERENCE.out.weights,
        params.num_reads,
        params.seed
    )

    // Step 3: Capture (conditional based on method)
    if (params.capture_method == 'blast') {
        def probe_basename = file(params.probes).getBaseName()
        def blast_db_dir_abs = file(params.blast_db_dir).toAbsolutePath().toString()
        def blast_db_path = "${blast_db_dir_abs}/${probe_basename}"

        BLAST_CAPTURE(
            blast_db_path,
            GENERATE_READS.out.reads,
            params.min_match_bases
        )
        captured_reads = BLAST_CAPTURE.out.captured
    } else {
        CAPTURE(
            probe_ch,
            GENERATE_READS.out.reads,
            params.max_mismatches,
            params.min_match_bases
        )
        captured_reads = CAPTURE.out.captured
    }

    // Step 4: Optional host filtering
    if (params.host_fasta) {
        host_ch = Channel.fromPath(params.host_fasta, checkIfExists: true)
        FILTER_HOST(
            host_ch,
            captured_reads,
            params.host_minimap_preset
        )
        reads_for_mapping = FILTER_HOST.out.filtered
    } else {
        reads_for_mapping = captured_reads
    }

    // Step 5: Map reads to reference
    MAP_READS(
        PREPARE_REFERENCE.out.reference,
        reads_for_mapping,
        params.minimap_preset
    )

    // Step 6: Generate detection list
    GENERATE_LIST(
        MAP_READS.out.sam
    )

    // Step 7: Calculate metrics
    CALCULATE_METRICS(
        PREPARE_REFERENCE.out.targets_list,
        PREPARE_REFERENCE.out.distractors_list,
        GENERATE_LIST.out.detected_list,
        GENERATE_READS.out.reads,
        captured_reads,
        run_name,
        params.num_reads,
        params.seed
    )

    // Step 8: Generate HTML report
    GENERATE_REPORT(
        CALCULATE_METRICS.out.results,
        CALCULATE_METRICS.out.detail,
        run_name
    )
}


workflow.onComplete {
    log.info """
    =============================================
    Pipeline completed!
    =============================================
    Run name    : ${run_name}
    Status      : ${workflow.success ? 'SUCCESS' : 'FAILED'}
    Duration    : ${workflow.duration}
    Output dir  : ${params.outdir}/${run_name}
    Report      : ${params.outdir}/${run_name}/report.html
    =============================================
    """.stripIndent()
}
