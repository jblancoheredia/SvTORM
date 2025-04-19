process DELLY {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/delly:1.3.1' :
        'blancojmskcc/delly:1.3.1' }"

    input:
    tuple val(meta) , path(input)
    tuple val(meta0), path(input_index)
    tuple val(meta1), path(fasta)
    tuple val(meta2), path(fai)
    tuple val(meta3), path(normal)
    tuple val(meta4), path(normal_index)
    path(exclude_bed)

    output:
    tuple val(meta), path("*.{csi,tbi}")            , emit: csi
    tuple val(meta), path("*.{bcf,delly.vcf.gz}")   , emit: bcf
    tuple val(meta), path("*.delly.unfiltered.vcf") , emit: vcf
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def normal = "NORMAL"
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = task.ext.suffix ?: "bcf"
    def exclude = exclude_bed ? "--exclude ${exclude_bed}" : ""
    def bcf_output = suffix == "bcf" ? "--outfile ${prefix}.bcf" : ""
    def bcf_filter = suffix == "bcf" ? "--outfile ${prefix}.filtered.bcf" : ""
    """
    echo -e "${meta.id}\\\\ttumor\\\\n${normal}\\\\tcontrol" > sample_file.tsv

    delly \\
        call \\
        ${args} \\
        ${bcf_output} \\
        --genome ${fasta} \\
        ${exclude} \\
        ${input} \\
        ${normal}.bam

    delly filter \\
        -t \\
        -v 3 \\
        -m 50 \\
        -a 0.01 \\
        -f somatic \\
        -s sample_file.tsv \\
        ${bcf_filter} \\
        ${prefix}.bcf

    bcftools convert -O v -o ${prefix}.delly.vcf ${prefix}.filtered.bcf

    awk 'BEGIN {FS=OFS=\"\\\\t\"}  /^#/ {print}' ${prefix}.delly.vcf > ${prefix}.delly.unfiltered.vcf
    awk 'BEGIN {FS=OFS=\"\\\\t\"}  \$1 ~ /^(1?[0-9]|2[0-2]|X|Y)\$/ {print}' ${prefix}.delly.vcf >> ${prefix}.delly.unfiltered.vcf

    bgzip ${prefix}.delly.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        delly: \$( echo \$(delly --version 2>&1) | sed 's/^.*Delly version: v//; s/ using.*\$//')
    END_VERSIONS
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.delly.vcf.gz
    touch ${prefix}.delly.unfiltered.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        delly: \$( echo \$(delly --version 2>&1) | sed 's/^.*Delly version: v//; s/ using.*\$//')
    END_VERSIONS
    """
}
