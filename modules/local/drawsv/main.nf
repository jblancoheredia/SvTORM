process DRAWSV {
    tag "$meta_tumour.patient_id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/drawsv:1.0.1':
        'blancojmskcc/drawsv:1.0.1' }"

    input:
    tuple val(patient_id), 
          val(meta_normal), path(normal_bam), path(normal_bai),
          val(meta_tumour), path(tumour_bam), path(tumour_bai)
    tuple val(meta1), path(tsv)
    path(gtf)
    path(cytobands)
    path(protein_domains)

    output:
    tuple val(meta_tumour), path("*.pdf"), emit: pdf
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta_tumour.patient_id}"
    def bam = "${tumour_bam}"
    """
    cp "${workflow.projectDir}/bin/PRE_DRAW_SV.py" .
    cp "${workflow.projectDir}/bin/DrawSVs.R" .

    python3 PRE_DRAW_SV.py \\
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
    def prefix = task.ext.prefix ?: "${meta_tumour.patient_id}"
    """
    touch ${prefix}_DrawSVs.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        drawsv: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """
}
