/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                    IMPORT PLUGINS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryMap                                                              } from 'plugin/nf-schema'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                    IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { DELLY                                                                         } from '../modules/local/delly/main'
include { SVABA                                                                         } from '../modules/local/svaba/main'
include { DRAWSV                                                                        } from '../modules/local/drawsv/main'
include { GRIDSS                                                                        } from '../modules/local/gridss/main'
include { MULTIQC                                                                       } from '../modules/nf-core/multiqc/main'
include { RECALL_SV                                                                     } from '../modules/local/recallsv/main'
include { MANTA_SOMATIC                                                                 } from '../modules/nf-core/manta/somatic/main'
include { SURVIVOR_MERGE                                                                } from '../modules/local/survivor/merge/main'
include { SURVIVOR_FILTER                                                               } from '../modules/local/survivor/filter/main'
include { GATK4_BEDTOINTERVALLIST                                                       } from '../modules/nf-core/gatk4/bedtointervallist/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                 IMPORT SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryMultiqc                                                          } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                                                        } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                                                        } from '../subworkflows/local/utils_nfcore_svtorm_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                CREATE CUSTOM CHANNELS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                 RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SVTORM {

    take:
    ch_bam
    ch_bai
    ch_fai
    ch_fasta
    ch_normal_bam
    ch_normal_bai
    ch_known_sites
    ch_targets_bed
    ch_targets_bed_tbi

    main:
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // MODULE: Run Manta in Somatic Mode
    //
    ch_final_input = ch_normal_bam
        .combine(ch_normal_bai)
        .combine(ch_bam)
        .combine(ch_bai)
        .combine(ch_targets_bed)
        .combine(ch_targets_bed_tbi)
        .map { normal_bam, normal_bai, tumor_bam, tumor_bai, target_bed, target_bed_tbi -> 
            tuple(normal_bam[0], normal_bam[1], normal_bai[1], tumor_bam[1], tumor_bai[1], target_bed[1], target_bed_tbi[1])
        }
    MANTA_SOMATIC(ch_bam, ch_intervals_gunzip, ch_intervals_gunzip_index, ch_fasta, ch_fai, [])
    ch_versions = ch_versions.mix(MANTA_SOMATIC.out.versions)
    ch_manta_vcf = MANTA_SOMATIC.out.vcf

    //
    // MODULE: Run Delly Call
    //
    DELLY(ch_bam, ch_bai, ch_fasta, ch_fai, ch_normal_bam, ch_normal_bai, params.exclude_bed)
    ch_versions = ch_versions.mix(DELLY.out.versions)
    ch_delly_vcf = DELLY.out.vcf

    //
    // MODULE: Run SvABA Note: version 1.2.0
    //
    SVABA(ch_bam, params.bwa, params.normal_bam, params.normal_bai, params.known_sites, params.known_sites_tbi, params.intervals)
    ch_versions = ch_versions.mix(SVABA.out.versions)
    ch_svaba_vcf = SVABA.out.vcf

    //
    // MODULE: Run Gridds (Extract overlapping fragments & calling)
    //
    GRIDSS(ch_bam, ch_fasta, ch_fai, params.normal_bam, params.normal_bai, params.intervals, params.blocklist_bed, params.bwa, params.kraken2db)
    ch_versions = ch_versions.mix(GRIDSS.out.versions)
    ch_gridss_vcf = GRIDSS.out.vcf

    //
    // MODULE: Run Survivor to merge Unfiltered VCFs
    //
    SURVIVOR_MERGE(ch_delly_vcf, ch_svaba_vcf, ch_manta_vcf, ch_gridss_vcf, params.chromosomes, 1000, 2, 0, 0, 0, 1000)
    ch_versions = ch_versions.mix(SURVIVOR_MERGE.out.versions)
    ch_merged_bed = SURVIVOR_MERGE.out.bed
    ch_merged_vcf = SURVIVOR_MERGE.out.vcf

    //
    // MODULE: Run GATK4 to Convert BED into Interval List
    //
    GATK4_BEDTOINTERVALLIST(ch_merged_bed, ch_dict)
    ch_versions = ch_versions.mix(GATK4_BEDTOINTERVALLIST.out.versions)
    ch_merged_int_list = GATK4_BEDTOINTERVALLIST.out.interval_list

    //
    // MODULE: Run Gridds in ReCall mode
    //
    RECALL_SV(ch_bam, ch_fasta, ch_fai, ch_merged_int_list, ch_known_sites, params.refflat, params.normal_bam, params.normal_bai, params.intervals,
              params.blocklist_bed, params.bwa, params.kraken2db, params.pon_directory)
    ch_versions = ch_versions.mix(RECALL_SV.out.versions)
    ch_recall_vcf = RECALL_SV.out.vcf

    //
    // MODULE: Run Survivor to filter Unfiltered VCFs
    //
    SURVIVOR_FILTER(ch_delly_vcf, ch_svaba_vcf, ch_manta_vcf, ch_gridss_vcf, ch_recall_vcf, 1000, 3, 0, 0, 0, 1000)
    ch_versions = ch_versions.mix(SURVIVOR_FILTER.out.versions)
    ch_filtered_vcf = SURVIVOR_FILTER.out.filtered_vcf

    //
    // MODULE: Run iAnnotateSV 
    //
    IANNOTATESV(ch_delly_vcf, ch_svaba_vcf, ch_manta_vcf, ch_gridss_vcf, ch_recall_vcf, 1000, 3, 0, 0, 0, 1000)
    ch_versions = ch_versions.mix(IANNOTATESV.out.versions)
    ch_filtered_vcf = IANNOTATESV.out.filtered_vcf

    //
    // MODULE: Run DrawSV
    //
    DRAWSV(ch_delly_vcf, ch_svaba_vcf, ch_manta_vcf, ch_gridss_vcf, ch_recall_vcf, 1000, 3, 0, 0, 0, 1000)
    ch_versions = ch_versions.mix(DRAWSV.out.versions)
    ch_filtered_vcf = DRAWSV.out.filtered_vcf

    //
    // MODULE: Run Check-Up and Clean-Up
    //
    CHECKUPNCLEANUP(ch_delly_vcf, ch_svaba_vcf, ch_manta_vcf, ch_gridss_vcf, ch_recall_vcf, 1000, 3, 0, 0, 0, 1000)
    ch_versions = ch_versions.mix(CHECKUPNCLEANUP.out.versions)
    ch_filtered_vcf = CHECKUPNCLEANUP.out.filtered_vcf

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'svtorm_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:multiqc_report = MULTIQC.out.report.toList()
    versions       = ch_versions

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                       THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
