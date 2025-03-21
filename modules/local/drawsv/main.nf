process DRAWSV {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/arriba:2.4.0--hdbdd923_4':
        'biocontainers/arriba:2.4.0--hdbdd923_4' }"

    input:
    tuple val(meta) , path(bam), path(bai)
    tuple val(meta1), path(tsv)
    path(gtf)
    path(cytobands)
    path(protein_domains)

    output:
    tuple val(meta), path("*.pdf"), emit: pdf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python PRE_DRAW_SV.py \\
        ${prefix} \\
        --input ${tsv} \\
        --annotations ${gtf} \\
        ${args}

    Rscript DrawSVs.R \\
        --fusions=${tsv} \\
        --alignments=${bam} \\
        --annotation=${gtf}   \\
        --cytobands=${cytobands} \\
        --output=${prefix}_DrawSVs.pdf \\
        --transcriptSelection=canonical \\
        --minConfidenceForCircosPlot=High \\
        --proteinDomains=${protein_domains} \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        drawsv: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_DrawSVs.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        drawsv: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """
}
