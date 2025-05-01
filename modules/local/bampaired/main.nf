process BAM_PAIRED {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::samtools=1.16.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.16.1--h6899075_0' :
        'biocontainers/samtools:1.16.1--h6899075_0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path(bam), path("is_singleend.txt"), emit:reads
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    BAM=\$(ls *.bam)

    if [ \$(samtools view \$BAM -@ ${task.cpus} | head -n 1000 | awk '{ if(and(\$2, 1)) count++ } END { print count+0 }') -gt 0 ]; then
        echo false > is_singleend.txt
    else
        echo true > is_singleend.txt
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}