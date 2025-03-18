process SURVIVOR_MERGE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/survivor:1.0.7':
        'blancojmskcc/survivor:1.0.7' }"

    input:
    tuple val(meta),  path(delly_vcf)
    tuple val(meta1), path(svaba_vcf)
    tuple val(meta2), path(manta_vcf)
    tuple val(meta4), path(gridss_vcf)
    path(chr_length)
    val(max_distance_breakpoints)
    val(min_supporting_callers)
    val(account_for_type)
    val(account_for_sv_strands)
    val(estimate_distanced_by_sv_size)
    val(min_sv_size)

    output:
    tuple val(meta), path("*.bed")   , emit: bed
    tuple val(meta), path("*.vcf")   , emit: vcf
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    SURVIVOR merge \\
        <(ls *.vcf) \\
        ${max_distance_breakpoints} \\
        ${min_supporting_callers} \\
        ${account_for_type} \\
        ${account_for_sv_strands} \\
        ${estimate_distanced_by_sv_size} \\
        ${min_sv_size} \\
        ${prefix}.survivor_sv_sup.vcf

    SURVIVOR vcftobed ${prefix}.survivor_sv_sup.vcf 0 -1 ${prefix}.survivor_sv_sup.bed.tmp

    awk '{ 
        new_start = (\$2 - 500 < 0) ? 0 : \$2 - 500;
        new_end = \$2 + 500;
        print \$1 \"\t\" new_start \"\t\" new_end 
    }' ${prefix}.survivor_sv_sup.bed.tmp >> ${prefix}.survivor_sv_sup.bed

    awk '{ 
        new_start = (\$5 - 500 < 0) ? 0 : \$5 - 500;
        new_end = \$5 + 500;
        print \$4 \"\t\" new_start \"\t\" new_end 
    }' ${prefix}.survivor_sv_sup.bed.tmp >> ${prefix}.survivor_sv_sup.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        survivor: \$(echo \$(SURVIVOR 2>&1 | grep "Version" | sed 's/^Version: //'))
    END_VERSIONS
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.survivor_sv_sup.vcf
    touch ${prefix}.survivor_sv_sup.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        survivor: \$(echo \$(SURVIVOR 2>&1 | grep "Version" | sed 's/^Version: //'))
    END_VERSIONS
    """
}
