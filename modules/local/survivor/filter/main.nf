process SURVIVOR_FILTER {
    tag "$meta.patient_id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/survivor_filter:1.0.7':
        'blancojmskcc/survivor_filter:1.0.7' }"

    input:
    tuple val(meta), 
          val(meta2), path('delly.vcf'),
          val(meta3), path('svaba.vcf'),
          val(meta4), path('manta.vcf'),
          val(meta5), path('gridss.vcf'),
          val(meta6), path('recall.vcf')
    val(max_distance_breakpoints)
    val(min_supporting_callers)
    val(account_for_type)
    val(account_for_sv_strands)
    val(estimate_distanced_by_sv_size)
    val(min_sv_size)

    output:
    tuple val(meta), path("*_SURVOR_SV_FIL.sort.vcf.gz"), path("*_SURVOR_SV_FIL.sort.vcf.gz.tbi"), emit: filtered_vcf
    tuple val(meta), path("*_SURVOR_SV_FIL.tsv")                                                 , emit: filtered_tsv
    tuple val(meta), path("*_ANNOTE_SV_INN.tsv")                                                 , emit: annote_input
    path "versions.yml"                                                                          , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.patient_id}"
    """
    awk 'BEGIN {FS=OFS="\\t"} {for (i=1; i<=NF; i++) if (i != 11) printf "%s%s", \$i, (i==NF || i==10 ? "\\n" : OFS)}' delly.vcf > ${prefix}_DELLY__SV_FIL.vcf
    awk 'BEGIN {FS=OFS="\\t"} {for (i=1; i<=NF; i++) if (i != 10) printf "%s%s", \$i, (i == NF ? "\\n" : OFS)}' svaba.vcf > ${prefix}_SVABA__SV_FIL.vcf
    awk 'BEGIN {FS=OFS="\\t"} {print}' manta.vcf  > ${prefix}_MANTA__SV_FIL.vcf
    awk 'BEGIN {FS=OFS="\\t"} {for (i=1; i<=NF; i++) if (i != 10) printf "%s%s", \$i, (i == NF ? "\\n" : OFS)}' gridss.vcf > ${prefix}_GRIDSS_SV_FIL.vcf
    awk 'BEGIN {FS=OFS="\\t"} {for (i=1; i<=NF; i++) if (i != 10) printf "%s%s", \$i, (i == NF ? "\\n" : OFS)}' recall.vcf > ${prefix}_RECALL_SV_FIL.vcf

    echo ${prefix}_DELLY__SV_FIL.vcf >> Original_VCFs_List.txt
    echo ${prefix}_SVABA__SV_FIL.vcf >> Original_VCFs_List.txt
    echo ${prefix}_MANTA__SV_FIL.vcf >> Original_VCFs_List.txt
    echo ${prefix}_GRIDSS_SV_FIL.vcf >> Original_VCFs_List.txt
    echo ${prefix}_RECALL_SV_FIL.vcf >> Original_VCFs_List.txt

    SURVIVOR merge \\
        <(ls *_SV_FIL.vcf) \\
        ${max_distance_breakpoints} \\
        ${min_supporting_callers} \\
        ${account_for_type} \\
        ${account_for_sv_strands} \\
        ${estimate_distanced_by_sv_size} \\
        ${min_sv_size} \\
        ${prefix}_SURVOR_SV_UNF.vcf || true

    grep \"#\" ${prefix}_SURVOR_SV_UNF.vcf > ${prefix}_SURVOR_SV_FIL.vcf

    awk 'BEGIN {FS=OFS=\"\\t\"} \$8 ~ /CHR2=(1?[0-9]|2[0-2])([^0-9]|\$)/ {print}' ${prefix}_SURVOR_SV_UNF.vcf >> ${prefix}_SURVOR_SV_FIL.vcf || true

    grep -v \"#\" ${prefix}_SURVOR_SV_UNF.vcf | grep \"CHR2=X\" >> ${prefix}_SURVOR_SV_FIL.vcf || true

    grep -v \"#\" ${prefix}_SURVOR_SV_UNF.vcf | grep \"CHR2=Y\" >> ${prefix}_SURVOR_SV_FIL.vcf || true
    
    SVRVOR2TSV \\
        --merged_vcf ${prefix}_SURVOR_SV_FIL.vcf \\
        --original_vcf_list Original_VCFs_List.txt \\
        --tsv ${prefix}_SURVOR_SV_FIL.tsv

    VCF2iANN \\
        ${prefix}
    
    grep \"#\" ${prefix}_SURVOR_SV_FIL.vcf > ${prefix}_SURVOR_SV_FIL.sort.vcf

    grep -v \"#\" ${prefix}_SURVOR_SV_FIL.vcf | sort -k1,1V -k2,2n >> ${prefix}_SURVOR_SV_FIL.sort.vcf

    bgzip ${prefix}_SURVOR_SV_FIL.sort.vcf

    tabix -p vcf ${prefix}_SURVOR_SV_FIL.sort.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        survivor: \$(echo \$(SURVIVOR 2>&1 | grep "Version" | sed 's/^Version: //'))
    END_VERSIONS
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.patient_id}"
    """
    touch ${prefix}_SURVOR_SV_FIL.vcf.gz
    touch ${prefix}_SURVOR_SV_FIL.vcf.gz.tbi
    touch ${prefix}_ANNOTE_SV_INN.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        survivor: \$(echo \$(SURVIVOR 2>&1 | grep "Version" | sed 's/^Version: //'))
    END_VERSIONS
    """
}
