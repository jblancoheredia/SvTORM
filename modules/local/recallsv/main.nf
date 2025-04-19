process RECALL_SV {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://drgiovianco/gridss:2.13.2':
        'drgiovianco/gridss:2.13.2' }"

    input:
    tuple val(meta) , path(tbam), path(tbai)
    tuple val(meta1), path(fasta)
    tuple val(meta2), path(fasta_fai)
    tuple val(meta3), path(interval_list)
    tuple val(meta4), path(known_sites), path(known_sites_tbi)
    path(refflat)
    path(nbam)
    path(nbai)
    path(bed)
    path(blocklist)
    path(bwa_index)
    path(kraken2db)
    path(pon_dir)

    output:
    path "versions.yml",                         emit: versions
    tuple val(meta), path("*.unfiltered.vcf"),   emit: vcf

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '2.13.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    def bwa = bwa_index ? "cp -s ${bwa_index}/* ." : ""
    def CollectGridssMetrics = "java -cp /opt/gridss/gridss--gridss-jar-with-dependencies.jar gridss.analysis.CollectGridssMetrics"
    """
    samtools view -h -F 256 -o ${prefix}_N_filtered.bam ${nbam}
    samtools index ${prefix}_N_filtered.bam

    samtools view -h -F 256 -o ${prefix}_T_filtered.bam ${tbam}
    samtools index ${prefix}_T_filtered.bam

    rm ${fasta} ${fasta_fai}

    ${bwa}

    mkdir ${prefix}-N.bam.gridss.working

    ${CollectGridssMetrics} \\
        INPUT=${prefix}_N_filtered.bam \\
        PROGRAM=RnaSeqMetrics \\
        THRESHOLD_COVERAGE=5000 \\
        PROGRAM=MeanQualityByCycle \\
        PROGRAM=CollectGcBiasMetrics \\
        PROGRAM=CollectInsertSizeMetrics \\
        PROGRAM=QualityScoreDistribution  \\
        PROGRAM=CollectQualityYieldMetrics \\
        PROGRAM=CollectBaseDistributionByCycle \\
        PROGRAM=CollectAlignmentSummaryMetrics  \\
        PROGRAM=CollectSequencingArtifactMetrics \\
        OUTPUT=${prefix}-N.bam.gridss.working/${prefix}-N.bam \\
        INTERVALS=${interval_list} \\
        REF_FLAT=${refflat} \\
        REFERENCE_SEQUENCE=${fasta} \\
        DB_SNP=${known_sites}

    gridss_extract_overlapping_fragments \\
        --targetbed ${bed} \\
        -t ${task.cpus} \\
        -o ${prefix}-N.bam \\
        ${prefix}_N_filtered.bam

    mkdir ${prefix}-T.bam.gridss.working

    ${CollectGridssMetrics} \\
        INPUT=${prefix}_T_filtered.bam \\
        PROGRAM=RnaSeqMetrics \\
        THRESHOLD_COVERAGE=5000 \\
        PROGRAM=MeanQualityByCycle \\
        PROGRAM=CollectGcBiasMetrics \\
        PROGRAM=CollectInsertSizeMetrics \\
        PROGRAM=QualityScoreDistribution  \\
        PROGRAM=CollectQualityYieldMetrics \\
        PROGRAM=CollectBaseDistributionByCycle \\
        PROGRAM=CollectAlignmentSummaryMetrics  \\
        PROGRAM=CollectSequencingArtifactMetrics \\
        OUTPUT=${prefix}-T.bam.gridss.working/${prefix}-T.bam \\
        INTERVALS=${interval_list} \\
        VALIDATION_STRINGENCY=LENIENT \\
        REF_FLAT=${refflat} \\
        REFERENCE_SEQUENCE=${fasta} \\
        DB_SNP=${known_sites}

    gridss_extract_overlapping_fragments \\
        --targetbed ${bed} \\
        -t ${task.cpus} \\
        -o ${prefix}-T.bam \\
        ${prefix}_T_filtered.bam

    gridss \\
        --threads ${task.cpus} \\
        --labels "NORMAL",${prefix} \\
        --jvmheap ${task.memory.toGiga() - 1}g \\
        --otherjvmheap ${task.memory.toGiga() - 1}g \\
        --output ${prefix}_all_calls.vcf \\
        --reference ${fasta} \\
        --picardoptions VALIDATION_STRINGENCY=LENIENT \\
        -b ${blocklist} \\
        ${prefix}-N.bam \\
        ${prefix}-T.bam

    gridss_annotate_vcf_kraken2 \\
        -t ${task.cpus} \\
        -o ${prefix}.recall.unfiltered.vcf \\
        --kraken2db ${kraken2db} \\
        ${prefix}_all_calls.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gridss: ${VERSION}
    END_VERSIONS
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '2.13.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch ${prefix}.recall.unfiltered.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gridss: ${VERSION}
    END_VERSIONS
    """
}
