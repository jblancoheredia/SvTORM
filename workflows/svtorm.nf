#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                    IMPORT PLUGINS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryMap                                                              } from 'plugin/nf-schema'
include { samplesheetToList                                                             } from 'plugin/nf-schema'

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
include { BAM_MATCH                                                                     } from '../modules/local/bammatch/main'
include { CAT_FASTQ                                                                     } from '../modules/nf-core/cat/fastq/main'  
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

ch_fai                                          = Channel.fromPath(params.fai).map                              { it -> [[id:it.Name], it] }.collect()
ch_fasta                                        = Channel.fromPath(params.fasta).map                            { it -> [[id:it.Name], it] }.collect()
ch_known_sites                                  = Channel.fromPath(params.known_sites).map                      { it -> [[id:it.Name], it] }.collect()
ch_targets_bed                                  = Channel.fromPath(params.intervals_bed_gunzip).map             { it -> [[id:it.Name], it] }.collect()
ch_targets_bed_tbi                              = Channel.fromPath(params.intervals_bed_gunzip_index).map       { it -> [[id:it.Name], it] }.collect()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                 RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SVTORM {

    take:
    ch_samplesheet

    main:
    ch_versions         = Channel.empty()
    ch_multiqc_files    = Channel.empty()

    // Create input channel depending on it's format 
    ch_samplesheet
        .branch { meta, files ->
            bam : files[0].getExtension() == "bam"
            bai : files[1].getExtension() == "bai"
            fastq_multiple :
                (meta.single_end && files.size() > 1) ||
                (!meta.single_end && files.size() > 2)
            fastq_single : true
        }
        .set { ch_input_files }

    //
    // MODULE: Concatenate FastQ files from same sample if required
    //
    CAT_FASTQ(ch_input_files.fastq_multiple).reads
    .set { ch_cat_fastq }
    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions.first())

    // determine BAM pairedness for fastq conversion
    BAM_MATCH (ch_input_files.bam )
    BAM_MATCH.out.reads
        .map {meta, reads, single_end ->
            meta["single_end"] = single_end.text.toBoolean()
            [meta, reads]
        }
        .set { ch_bam_pe_corrected }
    ch_versions = ch_versions.mix(BAM_MATCH.out.versions)

//    //
//    // Create channel from input file provided through params.input
//    //
//    Channel
//        .fromList(samplesheetToList(params.input, "${workflow.projectDir}/assets/schema_input.json"))
//        .map {
//            meta, bam_n, bai_n, bam_t, bai_t ->
//                return [ meta.id, meta + [ bam_n ] ]
//        }
//        .groupTuple()
//        .map {
//            validateInputSamplesheet(it)
//        }
//        .set { ch_nbam }
//
//    Channel
//        .fromList(samplesheetToList(params.input, "${workflow.projectDir}/assets/schema_input.json"))
//        .map {
//            meta, bam_n, bai_n, bam_t, bai_t ->
//                return [ meta.id, meta + [ bai_n ] ]
//        }
//        .groupTuple()
//        .map {
//            validateInputSamplesheet(it)
//        }
//        .set { ch_nbai }
//
//    Channel
//        .fromList(samplesheetToList(params.input, "${workflow.projectDir}/assets/schema_input.json"))
//        .map {
//            meta, bam_n, bai_n, bam_t, bai_t ->
//                return [ meta.id, meta + [ bam_t ] ]
//        }
//        .groupTuple()
//        .map {
//            validateInputSamplesheet(it)
//        }
//        .set { ch_tbam }
//
//    Channel
//        .fromList(samplesheetToList(params.input, "${workflow.projectDir}/assets/schema_input.json"))
//        .map {
//            meta, bam_n, bai_n, bam_t, bai_t ->
//                return [ meta.id, meta + [ bai_t ] ]
//        }
//        .groupTuple()
//        .map {
//            validateInputSamplesheet(it)
//        }
//        .set { ch_tbai }

//    //
//    // MODULE: Run Manta in Somatic Mode
//    //
//    ch_manta_input = ch_input_files
//        .combine(ch_targets_bed)
//        .combine(ch_targets_bed_tbi)
//        .map { input_files, target_bed, target_bed_tbi -> 
//            tuple(input_files[0], input_files[1], input_files[2], input_files[3], input_files[4], target_bed[1], target_bed_tbi[1])
//        }
//    MANTA_SOMATIC(ch_manta_input, ch_fasta, ch_fai, [])
//    ch_versions = ch_versions.mix(MANTA_SOMATIC.out.versions)
//    ch_manta_vcf = MANTA_SOMATIC.out.somatic_sv_vcf
//    ch_manta_vcf_tbi = MANTA_SOMATIC.out.somatic_sv_vcf_tbi
//    ch_manta_candidate_small_indels_vcf = MANTA_SOMATIC.out.candidate_small_indels_vcf
//    ch_manta_candidate_small_indels_vcf_tbi = MANTA_SOMATIC.out.candidate_small_indels_vcf_tbi
//
//    //
//    // MODULE: Run Delly Call
//    //
//    DELLY(ch_bam, ch_bai, ch_fasta, ch_fai, ch_normal_bam, ch_normal_bai, params.exclude_bed)
//    ch_versions = ch_versions.mix(DELLY.out.versions)
//    ch_delly_vcf = DELLY.out.vcf
//
//    //
//    // MODULE: Run SvABA Note: version 1.2.0
//    //
//    SVABA(ch_bam, params.bwa, params.normal_bam, params.normal_bai, params.known_sites, params.known_sites_tbi, params.intervals)
//    ch_versions = ch_versions.mix(SVABA.out.versions)
//    ch_svaba_vcf = SVABA.out.vcf
//
//    //
//    // MODULE: Run Gridds (Extract overlapping fragments & calling)
//    //
//    GRIDSS(ch_bam, ch_fasta, ch_fai, params.normal_bam, params.normal_bai, params.intervals, params.blocklist_bed, params.bwa, params.kraken2db)
//    ch_versions = ch_versions.mix(GRIDSS.out.versions)
//    ch_gridss_vcf = GRIDSS.out.vcf
//
//    //
//    // MODULE: Run Survivor to merge Unfiltered VCFs
//    //
//    SURVIVOR_MERGE(ch_delly_vcf, ch_svaba_vcf, ch_manta_vcf, ch_gridss_vcf, params.chromosomes, 1000, 2, 0, 0, 0, 1000)
//    ch_versions = ch_versions.mix(SURVIVOR_MERGE.out.versions)
//    ch_merged_bed = SURVIVOR_MERGE.out.bed
//    ch_merged_vcf = SURVIVOR_MERGE.out.vcf
//
//    //
//    // MODULE: Run GATK4 to Convert BED into Interval List
//    //
//    GATK4_BEDTOINTERVALLIST(ch_merged_bed, ch_dict)
//    ch_versions = ch_versions.mix(GATK4_BEDTOINTERVALLIST.out.versions)
//    ch_merged_int_list = GATK4_BEDTOINTERVALLIST.out.interval_list
//
//    //
//    // MODULE: Run Gridds in ReCall mode
//    //
//    RECALL_SV(ch_bam, ch_fasta, ch_fai, ch_merged_int_list, ch_known_sites, params.refflat, params.normal_bam, params.normal_bai, params.intervals,
//              params.blocklist_bed, params.bwa, params.kraken2db, params.pon_directory)
//    ch_versions = ch_versions.mix(RECALL_SV.out.versions)
//    ch_recall_vcf = RECALL_SV.out.vcf
//
//    //
//    // MODULE: Run Survivor to filter Unfiltered VCFs
//    //
//    SURVIVOR_FILTER(ch_delly_vcf, ch_svaba_vcf, ch_manta_vcf, ch_gridss_vcf, ch_recall_vcf, 1000, 3, 0, 0, 0, 1000)
//    ch_versions = ch_versions.mix(SURVIVOR_FILTER.out.versions)
//    ch_filtered_vcf = SURVIVOR_FILTER.out.filtered_vcf
//
//    //
//    // MODULE: Run iAnnotateSV 
//    //
//    IANNOTATESV(ch_delly_vcf, ch_svaba_vcf, ch_manta_vcf, ch_gridss_vcf, ch_recall_vcf, 1000, 3, 0, 0, 0, 1000)
//    ch_versions = ch_versions.mix(IANNOTATESV.out.versions)
//    ch_annotated_tsv = IANNOTATESV.out.tsv

//    //
//    // MODULE: Run DrawSV
//    //
//    DRAWSV(ch_bam, ch_annotated_tsv, params.annotations, params.cytobands, params.protein_domains)
//    ch_versions = ch_versions.mix(DRAWSV.out.versions)
//    ch_drawsv_pdf = DRAWSV.out.pdf

//    //
//    // MODULE: Run Check-Up and Clean-Up
//    //
//    CHECKUPNCLEANUP(ch_annotated_tsv)
//    ch_versions = ch_versions.mix(CHECKUPNCLEANUP.out.versions)

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

    emit:
    versions       = ch_versions
    multiqc_report = MULTIQC.out.report.toList()

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                       THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
