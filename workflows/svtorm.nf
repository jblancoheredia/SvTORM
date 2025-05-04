#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                        IMPORT PLUGINS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryMap                                                              } from 'plugin/nf-schema'
include { samplesheetToList                                                             } from 'plugin/nf-schema'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                        IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { DELLY                                                                         } from '../modules/local/delly/main'
include { SVABA                                                                         } from '../modules/local/svaba/main'
include { DRAWSV                                                                        } from '../modules/local/drawsv/main'
include { FASTQC                                                                        } from '../modules/nf-core/fastqc/main'
include { GRIDSS                                                                        } from '../modules/local/gridss/main'
include { MULTIQC                                                                       } from '../modules/nf-core/multiqc/main'
include { CAT_FASTQ                            as CAT_FASTQ_NORMAL                      } from '../modules/nf-core/cat/fastq/main'  
include { CAT_FASTQ                            as CAT_FASTQ_TUMOUR                      } from '../modules/nf-core/cat/fastq/main'  
include { RECALL_SV                                                                     } from '../modules/local/recallsv/main'
include { BAM_PAIRED                                                                    } from '../modules/local/bampaired/main'
include { IANNOTATESV                                                                   } from '../modules/local/iannotatesv/main'
include { MANTA_SOMATIC                                                                 } from '../modules/local/manta/somatic/main'
include { SURVIVOR_MERGE                                                                } from '../modules/local/survivor/merge/main'
include { SURVIVOR_STATS                                                                } from '../modules/local/survivor/stats/main'
include { SURVIVOR_FILTER                                                               } from '../modules/local/survivor/filter/main'
include { GATK4_BEDTOINTERVALLIST                                                       } from '../modules/nf-core/gatk4/bedtointervallist/main'
include { PICARD_COLLECTMULTIPLEMETRICS                                                 } from '../modules/nf-core/picard/collectmultiplemetrics/main'
include { VCF_ANNOTATE_ENSEMBLVEP_SNPEFF                                                } from '../subworkflows/nf-core/vcf_annotate_ensemblvep_snpeff/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                     IMPORT SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryMultiqc                                                          } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                                                        } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                                                        } from '../subworkflows/local/utils_nfcore_svtorm_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

                                                                    CREATE CUSTOM CHANNELS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_fai                                          = Channel.fromPath(params.fai).map                              { it -> [[id:it.Name], it] }.collect()
ch_dict                                         = Channel.fromPath(params.dict).map                             { it -> [[id:it.Name], it] }.collect()
ch_fasta                                        = Channel.fromPath(params.fasta).map                            { it -> [[id:it.Name], it] }.collect()
ch_known_sites                                  = Channel.fromPath(params.known_sites).map                      { it -> [[id:it.Name], it] }.collect()
ch_targets_bed                                  = Channel.fromPath(params.intervals_bed_gunzip).map             { it -> [[id:it.Name], it] }.collect()
ch_known_sites_tbi                              = Channel.fromPath(params.known_sites_tbi).map                  { it -> [[id:it.Name], it] }.collect()
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
    ch_reports  = Channel.empty()
    ch_versions = Channel.empty()

    //
    // STEP 1 — Normalize and branch input by type
    //
    ch_samplesheet
        .branch { meta, files ->
            is_bam_input : files[0].getExtension() == "bam"
            is_fastq     : true
        }
        .set { ch_input_branch }
    ch_input_branch.is_fastq.set { ch_fastq_input }
    ch_input_branch.is_bam_input.set { ch_bam_input }

    //
    // STEP 2 — Get single_end from BAM using BAM_PAIRED
    //
    // MODULE: Determine BAM pairedness for fastq conversion
    //
    BAM_PAIRED(ch_bam_input)
    ch_bam_with_pairedness = BAM_PAIRED.out.reads.map { meta, bam_tuple, single_end ->
        meta["single_end"] = single_end.text.trim().toBoolean()
        [meta, *bam_tuple]
    }
    ch_versions = ch_versions.mix(BAM_PAIRED.out.versions)

    //
    // STEP 3 — Prepare FASTQ input with placeholder `single_end` and lane count
    //
    ch_fastq_prepped = ch_fastq_input.map { meta, fastq_1, fastq_2 ->
        def lane_count = fastq_1 instanceof List ? fastq_1.size() : 1
        meta["single_end"] = fastq_2 == null
        meta["lane_count"] = lane_count
        [meta, fastq_1, fastq_2, null, null]
    }

    //
    // STEP 4 — Combine both BAM and FASTQ inputs into one unified stream
    //
    ch_all_inputs = ch_bam_with_pairedness.map { meta, bam, bai ->
        [meta, null, null, bam, bai]
    }.mix(ch_fastq_prepped)

    //
    // STEP 5 — Restructure by patient for grouping
    //
    ch_grouped_by_patient = ch_all_inputs
        .map { meta, fastq_1, fastq_2, bam, bai ->
            def sample_id = "${meta.patient}_${meta.sample}"
            [meta.patient, [meta, fastq_1, fastq_2, bam, bai]]
        }
        .groupTuple()
        .map { patient_id, samples ->

            def normal = null
            def tumour = null

            for (sample in samples) {
                def (meta, fastq_1, fastq_2, bam, bai) = sample

                def sample_data = [
                    meta       : meta,
                    fastq_1    : fastq_1,
                    fastq_2    : fastq_2,
                    bam        : bam,
                    bai        : bai,
                    single_end : meta.single_end,
                    lane_count : meta.lane_count ?: (fastq_1 instanceof List ? fastq_1.size() : 1)
                ]

                if (meta.status == 0) {
                    normal = sample_data
                } else if (meta.status == 1) {
                    tumour = sample_data
                }
            }

            if (!normal) {
                if (!params.normal_bam || !params.normal_bai) {
                    error("Missing matched normal and no fallback BAM provided for patient: ${patient_id}")
                }
                normal = [
                    meta       : [ id: "external_normal", input_type: "bam", status: 0 ],
                    bam        : file(params.normal_bam),
                    bai        : file(params.normal_bai),
                    fastq_1    : null,
                    fastq_2    : null,
                    single_end : false,
                    lane_count : 1
                ]
            }

            if (!tumour) {
                error("Tumour sample is required for patient: ${patient_id}")
            }

            def patient_meta = [
                patient_id           : patient_id,

                normal_meta          : normal.meta,
                tumour_meta          : tumour.meta,

                normal_bam           : normal.bam,
                normal_bai           : normal.bai,
                tumour_bam           : tumour.bam,
                tumour_bai           : tumour.bai,

                normal_fastq1        : normal.fastq_1,
                normal_fastq2        : normal.fastq_2,
                tumour_fastq1        : tumour.fastq_1,
                tumour_fastq2        : tumour.fastq_2,

                normal_is_single_end : normal.single_end,
                tumour_is_single_end : tumour.single_end,

                normal_lane_count    : normal.lane_count,
                tumour_lane_count    : tumour.lane_count
            ]

            return tuple(patient_id, patient_meta)
        }
        .set { ch_patient_inputs }

    //
    // STEP 6 — Create Normal & Tumour BAM input channels
    //
    ch_bam_normal = ch_patient_inputs.map { patient_id, meta ->
        def normal_meta = [
            patient_id : meta.patient_id,
            sample_id  : meta.normal_meta.id,
            status     : 'normal'
        ]
        tuple(normal_meta, meta.normal_bam, meta.normal_bai)
    }
    
    ch_bam_tumour = ch_patient_inputs.map { patient_id, meta ->
        def tumour_meta = [
            patient_id : meta.patient_id,
            sample_id  : meta.tumour_meta.id,
            status     : 'tumour'
        ]
        tuple(tumour_meta, meta.tumour_bam, meta.tumour_bai)
    }

    //
    // STEP 7 — Branch for FASTQ that needs concatenation
    //
    ch_patient_inputs
        .branch {
            tumour_fastq_multiple :
                (!it[1].tumour_is_single_end && it[1].tumour_lane_count > 2) ||
                ( it[1].tumour_is_single_end && it[1].tumour_lane_count > 1 )
            normal_fastq_multiple :
                (!it[1].normal_is_single_end && it[1].normal_lane_count > 2) ||
                ( it[1].normal_is_single_end && it[1].normal_lane_count > 1 )
            tumour_fastq_single : true
            normal_fastq_single : true
        }
        .set { ch_fastq_branch }
    ch_fastq_branch.tumour_fastq_single.set { ch_tumour_fastq_single }
    ch_fastq_branch.normal_fastq_single.set { ch_normal_fastq_single }
    ch_fastq_branch.tumour_fastq_multiple.set { ch_tumour_fastq_multi }
    ch_fastq_branch.normal_fastq_multiple.set { ch_normal_fastq_multi }

//    if (ch_tumour_fastq_multi  || ch_normal_fastq_multi  || 
//        ch_tumour_fastq_single || ch_normal_fastq_single) {
//    
//        //
//        // MODULE: Run CAT_FASTQ if there's FASTQ files
//        //
//        CAT_FASTQ_TUMOUR(ch_tumour_fastq_multi.map { _, meta -> tuple(meta.tumour_meta, meta.tumour_fastq1, meta.tumour_fastq2) })
//        ch_cat_tumour_fastq = CAT_FASTQ_TUMOUR.out.reads
//        ch_versions = ch_versions.mix(CAT_FASTQ_TUMOUR.out.versions.first())
//        
//        CAT_FASTQ_NORMAL(ch_normal_fastq_multi.map { _, meta -> tuple(meta.normal_meta, meta.normal_fastq1, meta.normal_fastq2) })
//        ch_cat_normal_fastq = CAT_FASTQ_NORMAL.out.reads  // Adjust this name!
//        ch_versions = ch_versions.mix(CAT_FASTQ_NORMAL.out.versions.first())
//
//        ch_all_fastq_for_fastqc = Channel.empty()
//            .mix(ch_cat_tumour_fastq)
//            .mix(ch_cat_normal_fastq)
//            .mix(ch_tumour_fastq_single)
//            .mix(ch_normal_fastq_single)
//
//        //
//        // MODULE: Run FastQC if there's FASTQ files
//        //
//        FASTQC(ch_all_fastq_for_fastqc)
//        ch_reports = ch_reports.mix(FASTQC.out.zip.collect{it[1]})
//        ch_versions = ch_versions.mix(FASTQC.out.versions)
//    }

    //
    // MODULE: Run Picard Collect Multiple Metrics
    //
    PICARD_COLLECTMULTIPLEMETRICS(ch_bam_tumour, ch_fasta, ch_fai)
    ch_versions = ch_versions.mix(PICARD_COLLECTMULTIPLEMETRICS.out.versions)
    ch_reports  = ch_reports.mix(PICARD_COLLECTMULTIPLEMETRICS.out.metrics.collect{it[1]}.ifEmpty([]))

    //
    // STEP 8: Pairing the BAM files by patient_id
    //
    ch_bam_normal = ch_bam_normal.map { meta, bam, bai -> tuple(meta.patient_id, meta, bam, bai) }
    ch_bam_tumour = ch_bam_tumour.map { meta, bam, bai -> tuple(meta.patient_id, meta, bam, bai) }

    ch_bam_pairs = ch_bam_normal
        .join(ch_bam_tumour)
        .map { patient_id, meta_n, bam_n, bai_n, meta_t, bam_t, bai_t ->
            tuple(patient_id, meta_n, bam_n, bai_n, meta_t, bam_t, bai_t)
        }

    //
    // MODULE: Run Delly Call
    //
    DELLY(ch_bam_pairs, ch_fasta, ch_fai, params.exclude_bed)
    ch_versions = ch_versions.mix(DELLY.out.versions)
    ch_delly_vcf = DELLY.out.vcf
    ch_delly_vcf = ch_delly_vcf.map { meta, vcf -> tuple(meta.patient_id, meta, vcf) }

    //
    // MODULE: Run SvABA Note: version 1.2.0
    //
    SVABA(ch_bam_pairs, ch_fasta, ch_fai, params.known_sites_tbi, params.known_sites, params.intervals, params.bwa)
    ch_versions = ch_versions.mix(SVABA.out.versions)
    ch_svaba_vcf = SVABA.out.vcf
    ch_svaba_vcf = ch_svaba_vcf.map { meta, vcf -> tuple(meta.patient_id, meta, vcf) }

    //
    // MODULE: Run Manta in Somatic Mode
    //
    MANTA_SOMATIC(ch_bam_pairs, ch_targets_bed, ch_targets_bed_tbi, ch_fasta, ch_fai, [])
    ch_versions = ch_versions.mix(MANTA_SOMATIC.out.versions)
    ch_manta_vcf = MANTA_SOMATIC.out.vcf
    ch_manta_candidate_small_indels_vcf = MANTA_SOMATIC.out.candidate_small_indels_vcf
    ch_manta_candidate_small_indels_vcf_tbi = MANTA_SOMATIC.out.candidate_small_indels_vcf_tbi
    ch_manta_vcf = ch_manta_vcf.map { meta, vcf -> tuple(meta.patient_id, meta, vcf) }

    //
    // MODULE: Run Gridds (Extract overlapping fragments & calling)
    //
    GRIDSS(ch_bam_pairs, ch_fasta, ch_fai, params.intervals, params.blocklist_bed, params.bwa, params.kraken2db)
    ch_versions = ch_versions.mix(GRIDSS.out.versions)
    ch_gridss_vcf = GRIDSS.out.vcf
    ch_gridss_vcf = ch_gridss_vcf.map { meta, vcf -> tuple(meta.patient_id, meta, vcf) }

    //
    // Combine the vcf by patient_id
    //
    ch_vcf_merged = ch_delly_vcf
        .join(ch_svaba_vcf)
        .join(ch_manta_vcf)
        .join(ch_gridss_vcf)
        .map { patient_id, meta_delly, delly_vcf, meta_svaba, svaba_vcf, meta_manta, manta_vcf, meta_gridss, gridss_vcf ->
            tuple(
                meta_delly, 
                meta_delly, delly_vcf,
                meta_svaba, svaba_vcf,
                meta_manta, manta_vcf,
                meta_gridss, gridss_vcf
            )
        }
    
    //
    // MODULE: Run Survivor to merge Unfiltered VCFs
    //
    SURVIVOR_MERGE(ch_vcf_merged, params.chromosomes, 1000, 2, 0, 0, 0, 30)
    ch_versions = ch_versions.mix(SURVIVOR_MERGE.out.versions)
    ch_merged_bed = SURVIVOR_MERGE.out.bed
    ch_merged_vcf = SURVIVOR_MERGE.out.vcf

    //
    // MODULE: Run GATK4 to Convert BED into Interval List
    //
    GATK4_BEDTOINTERVALLIST(ch_merged_bed, ch_dict)
    ch_versions = ch_versions.mix(GATK4_BEDTOINTERVALLIST.out.versions)
    ch_merged_int_list = GATK4_BEDTOINTERVALLIST.out.interval_list
    ch_merged_int_list = ch_merged_int_list.map { meta, interval_list -> tuple(meta.patient_id, meta, interval_list) }

    //
    // Join interval lists with BAM pairs based on patient_id
    //
    ch_recall_input = ch_bam_pairs
        .join(ch_merged_int_list)
        .map { patient_id, meta_n, bam_n, bai_n, meta_t, bam_t, bai_t, meta_i, interval_list ->
            tuple(
                meta_t, 
                meta_n, bam_n, bai_n, 
                meta_t, bam_t, bai_t, 
                meta_i, interval_list
            )
        }

    //
    // MODULE: Run Gridds in ReCall mode
    //
    RECALL_SV(ch_recall_input, ch_fasta, ch_fai, ch_known_sites, ch_known_sites_tbi, params.refflat, params.intervals, params.blocklist_bed, params.bwa, params.kraken2db, params.pon_directory)
    ch_versions = ch_versions.mix(RECALL_SV.out.versions)
    ch_recall_vcf = RECALL_SV.out.vcf
    ch_recall_vcf = ch_recall_vcf.map { meta, vcf -> tuple(meta.patient_id, meta, vcf) }

    //
    // Combine the vcf by patient_id
    //
    ch_vcfs_merged = ch_delly_vcf
        .join(ch_svaba_vcf)
        .join(ch_manta_vcf)
        .join(ch_gridss_vcf)
        .join(ch_recall_vcf)
        .map { patient_id, meta_delly, delly_vcf, meta_svaba, svaba_vcf, meta_manta, manta_vcf, meta_gridss, gridss_vcf, meta_recall, recall_vcf ->
            tuple(
                meta_delly, 
                meta_delly, delly_vcf,
                meta_svaba, svaba_vcf,
                meta_manta, manta_vcf,
                meta_gridss, gridss_vcf,
                meta_recall, recall_vcf
            )
        }

    //
    // MODULE: Run Survivor to filter Unfiltered VCFs
    //
    SURVIVOR_FILTER(ch_vcfs_merged, 10000, 3, 1, 1, 0, 50)
    ch_versions = ch_versions.mix(SURVIVOR_FILTER.out.versions)
    ch_filtered_vcf = SURVIVOR_FILTER.out.filtered_vcf
    ch_filtered_tsv = SURVIVOR_FILTER.out.filtered_tsv
    ch_annote_input = SURVIVOR_FILTER.out.annote_input

    //
    // MODULE: Run Survivor Stats
    //
    SURVIVOR_STATS(ch_filtered_vcf, -1, -1, -1)
    ch_versions = ch_versions.mix(SURVIVOR_STATS.out.versions)
    ch_reports  = ch_reports.mix(SURVIVOR_STATS.out.stats.collect{it[1]}.ifEmpty([]))

    //
    // MODULE: Run iAnnotateSV 
    //
    IANNOTATESV(ch_filtered_vcf, ch_filtered_tsv, ch_annote_input)
    ch_versions = ch_versions.mix(IANNOTATESV.out.versions)
    ch_annotated_tsv = IANNOTATESV.out.tsv

//    //
//    // MODULE: Run 
//    //
//    VCF_ANNOTATE_ENSEMBLVEP_SNPEFF(ch_filtered_vcf, ch_fasta, )
//    ch_versions = ch_versions.mix(VCF_ANNOTATE_ENSEMBLVEP_SNPEFF.out.versions)
//
//    //
//    // MODULE: Run DrawSV
//    //
//    DRAWSV(ch_bam_pairs, ch_filtered_tsv, params.annotations, params.cytobands, params.protein_domains)
//    ch_versions = ch_versions.mix(DRAWSV.out.versions)
//    ch_drawsv_pdf = DRAWSV.out.pdf
//    
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
    versions = ch_collated_versions

    //
    // MODULE: MultiQC
    //
    val_multiqc_report = Channel.empty()
    if (!params.skip_multiqc){
        multiqc_files = Channel.empty()
        multiqc_config = Channel.fromPath("${projectDir}/assets/multiqc_config.yml", checkIfExists: true)
        multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
        multiqc_logo = params.multiqc_logo ? Channel.fromPath(params.multiqc_logo, checkIfExists: true) : Channel.empty()
        summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
        workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
        multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
        methods_description = Channel.value(methodsDescriptionText(multiqc_custom_methods_description))
        multiqc_files = multiqc_files.mix(ch_reports)
        multiqc_files = multiqc_files.mix(workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        multiqc_files = multiqc_files.mix(versions)
        multiqc_files = multiqc_files.mix(methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))

        MULTIQC (
            multiqc_files.collect(),
            multiqc_config.toList(),
            multiqc_custom_config.toList(),
            multiqc_logo.toList(),
            [],
            []
        )
        val_multiqc_report = MULTIQC.out.report.toList()
        ch_versions = ch_versions.mix(MULTIQC.out.versions)
    }

    emit:
    multiqc_report = val_multiqc_report
    versions = ch_versions
}
//    }
//
//
//
//
//    ch_multiqc_config                     = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
//    ch_multiqc_custom_config              = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
//    ch_multiqc_logo                       = params.multiqc_logo ? Channel.fromPath(params.multiqc_logo, checkIfExists: true) : Channel.empty()
//    summary_params                        = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
//    ch_workflow_summary                   = Channel.value(paramsSummaryMultiqc(summary_params))
//    ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
//    ch_methods_description                = Channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
//    ch_multiqc_files                      = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
//    ch_multiqc_files                      = ch_multiqc_files.mix(ch_collated_versions)
//    ch_multiqc_files                      = ch_multiqc_files.mix(ch_reports)
//    ch_multiqc_files                      = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))
//
//    MULTIQC (
//        ch_multiqc_files.collect(),
//        ch_multiqc_config.toList(),
//        ch_multiqc_custom_config.toList(),
//        ch_multiqc_logo.toList(),
//        [],
//        []
//    )
//    multiqc_report = MULTIQC.out.report.toList()
//
//
//    emit:
//    multiqc_report // channel: /path/to/multiqc_report.html
//    versions       // channel: [ path(versions.yml) ]
//}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                            THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

*/
