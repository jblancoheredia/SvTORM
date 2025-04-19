process SVABA {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/svaba:1.2.0' :
        'blancojmskcc/svaba:1.2.0' }"

    input:
    tuple val(meta),  path(tbam), path(tbai)
    path(bwa)
    path(nbam)
    path(nbai)
    path(dbsnp)
    path(dbsnp_tbi)
    path(bed)

    output:
    tuple val(meta), path("*.svaba.unfiltered.vcf"), emit: vcf
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def args2  = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    svaba run \\
        -t ${tbam} \\
        -n ${nbam} \\
        $args \\
        -a ${prefix} \\
        -D ${dbsnp} \\
        -k ${bed} \\
        -p $task.cpus \\
        $args2

    awk 'BEGIN {FS=OFS=\"\\\\t\"}  /^#/ {print}' ${prefix}.svaba.unfiltered.somatic.sv.vcf > ${prefix}.svaba.unfiltered.vcf
    awk 'BEGIN {FS=OFS=\"\\\\t\"}  \$1 ~ /^(1?[0-9]|2[0-2]|X|Y)\$/ {print}' ${prefix}.svaba.unfiltered.somatic.sv.vcf  >> ${prefix}.svaba.unfiltered.vcf
    awk 'BEGIN {FS=OFS=\"\\\\t\"}  \$1 ~ /^(1?[0-9]|2[0-2]|X|Y)\$/ {print}' ${prefix}.svaba.unfiltered.germline.sv.vcf >> ${prefix}.svaba.unfiltered.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        svaba: \$(svaba --version |& sed '1!d ; s/svaba //')
    END_VERSIONS
    """
    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.svaba.unfiltered.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        svaba: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """
}
