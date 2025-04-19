process GRIDSS {
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
    path(nbam)
    path(nbai)
    path(bed)
    path(blocklist)
    path(bwa_index)
    path(kraken2db)

    output:
    path "versions.yml",                                emit: versions
    tuple val(meta), path("*.gridss.unfiltered.vcf"),   emit: vcf

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '2.13.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    def bwa = bwa_index ? "cp -s ${bwa_index}/* ." : ""
    """
    samtools view -h -F 256 -o ${prefix}_N_filtered.bam ${nbam}
    samtools index ${prefix}_N_filtered.bam

    samtools view -h -F 256 -o ${prefix}_T_filtered.bam ${tbam}
    samtools index ${prefix}_T_filtered.bam

    rm ${fasta} ${fasta_fai}

    ${bwa}

    gridss_extract_overlapping_fragments \\
        --targetbed ${bed} \\
        -t ${task.cpus} \\
        -o ${prefix}_GRIDSS-N.bam \\
        ${prefix}_N_filtered.bam

    gridss_extract_overlapping_fragments \\
        --targetbed ${bed} \\
        -t ${task.cpus} \\
        -o ${prefix}_GRIDSS-T.bam \\
        ${prefix}_T_filtered.bam

    gridss \\
        --output ${prefix}.vcf.gz \\
        --reference ${fasta} \\
        --threads ${task.cpus} \\
        --jvmheap ${task.memory.toGiga() - 1}g \\
        --otherjvmheap ${task.memory.toGiga() - 1}g \\
        -b ${blocklist} \\
        ${prefix}_GRIDSS-N.bam \\
        ${prefix}_GRIDSS-T.bam

    gridss_annotate_vcf_kraken2 \\
        -t ${task.cpus} \\
        -o ${prefix}_AVK.vcf \\
        --kraken2db ${kraken2db} \\
        ${prefix}.vcf.gz

    awk 'BEGIN {FS=OFS=\"\\\\t\"}  /^#/ {print}' ${prefix}_AVK.vcf > ${prefix}.gridss.unfiltered.vcf
    awk 'BEGIN {FS=OFS=\"\\\\t\"}  \$1 ~ /^(1?[0-9]|2[0-2]|X|Y)\$/ {print}' ${prefix}_AVK.vcf >> ${prefix}.gridss.unfiltered.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gridss: ${VERSION}
    END_VERSIONS
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '2.13.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch ${prefix}.gridss.unfiltered.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gridss: ${VERSION}
    END_VERSIONS
    """
}
