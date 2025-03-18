process GRIDSS {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gridss:2.13.2--h270b39a_0':
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

    output:
    tuple val(meta), path("*.vcf.gz")       , emit: vcf
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '2.13.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    def bwa = bwa_index ? "cp -s ${bwa_index}/* ." : ""
    """
    rm ${fasta} ${fasta_fai}

    ${bwa}

    gridss_extract_overlapping_fragments \\
        --targetbed ${bed} \\
        -t ${task.cpus} \\
        -o ${prefix}_GRIDSS-N.bam \\
        ${nbam}

    gridss_extract_overlapping_fragments \\
        --targetbed ${bed} \\
        -t ${task.cpus} \\
        -o ${prefix}_GRIDSS-T.bam \\
        ${tbam}

    gridss \\
        --output ${prefix}.vcf.gz \\
        --reference ${fasta} \\
        --threads ${task.cpus} \\
        --jvmheap ${task.memory.toGiga() - 1}g \\
        --otherjvmheap ${task.memory.toGiga() - 1}g \\
        -b ${blocklist} \\
        ${prefix}_GRIDSS-N.bam \\
        ${prefix}_GRIDSS-T.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gridss: ${VERSION}
    END_VERSIONS
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '2.13.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gridss: ${VERSION}
    END_VERSIONS
    """
}
