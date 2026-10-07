/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FQ_SUBSAMPLE           } from '../modules/nf-core/fq/subsample/main'
include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { FASTP                  } from '../modules/nf-core/fastp/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { HISAT2_EXTRACTSPLICESITES } from '../modules/nf-core/hisat2/extractsplicesites/main'
include { HISAT2_BUILD              } from '../modules/nf-core/hisat2/build/main'
include { HISAT2_ALIGN              } from '../modules/nf-core/hisat2/align/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_analysisrnaseq_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ANALYSISRNASEQ {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir

    main:

    def ch_versions = channel.empty()
    def ch_multiqc_files = channel.empty()

    //
    // MODULE: Subsample reads (optional, only if --subsample is set)
    // All following steps use ch_reads (subsampled reads or all reads)
    //

    def ch_reads = ch_samplesheet
    if (params.subsample) {
        FQ_SUBSAMPLE(ch_samplesheet)
        ch_reads = FQ_SUBSAMPLE.out.fastq
    }

    //
    // MODULE: Run FastQC (quality control of raw reads)
    //
    FASTQC(ch_reads)
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.map{ _meta, file -> file })

    //
    // MODULE: Run fastp (adapter and quality trimming)
    //
    FASTP(
        ch_reads.map { meta, reads -> [ meta, reads, [] ] }, // [] = no adapter file, fastp detects adapters automatically
        false,  // discard_trimmed_pass: keep the trimmed reads
        false,  // save_trimmed_fail: do not save reads that fail filtering
        false   // save_merged: do not merge paired-end reads
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json.map{ _meta, file -> file })

    // trimmed reads for the next step (alignment)
    def ch_trimmed_reads = FASTP.out.reads


        //
    // Reference files (genome and annotation)
    //
    def ch_fasta = channel.value([ [id: 'genome'], file(params.fasta, checkIfExists: true) ])
    def ch_gtf   = channel.value([ [id: 'genome'], file(params.gtf,   checkIfExists: true) ])

    //
    // MODULE: Extract splice sites from the GTF (needed for spliced RNA-seq reads)
    //
    HISAT2_EXTRACTSPLICESITES(ch_gtf)
    def ch_splicesites = HISAT2_EXTRACTSPLICESITES.out.txt.first()

    //
    // MODULE: Build HISAT2 index (only if no prebuilt index is given)
    //
    def ch_hisat2_index = channel.empty()
    if (params.hisat2_index) {
        ch_hisat2_index = channel.value([ [id: 'genome'], file(params.hisat2_index, checkIfExists: true) ])
    } else {
        HISAT2_BUILD(
            ch_fasta
                .combine(ch_gtf)
                .combine(ch_splicesites)
                .map { meta, fasta, _meta2, gtf, _meta3, splicesites -> [ meta, fasta, gtf, splicesites ] },
            '200.GB' 
        )
        ch_hisat2_index = HISAT2_BUILD.out.index.first()
    }

    //
    // MODULE: Align trimmed reads with HISAT2
    //
    HISAT2_ALIGN(
        ch_trimmed_reads,   // trimmed reads from fastp
        ch_hisat2_index,    // HISAT2 index
        ch_splicesites,     // known splice sites
        false               // save_unaligned: do not save unmapped reads
    )
    ch_multiqc_files = ch_multiqc_files.mix(HISAT2_ALIGN.out.summary.map{ _meta, file -> file })

    // aligned reads (BAM) for the next step (mark duplicates)
    def ch_bam = HISAT2_ALIGN.out.bam


    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name:  'analysisrnaseq_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    //
    // MODULE: MultiQC
    //
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    def ch_multiqc_custom_methods_description = multiqc_methods_description
        ? file(multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    def ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))
    MULTIQC(
        ch_multiqc_files.flatten().collect().map { files ->
            [
                [id: 'analysisrnaseq'],
                files,
                multiqc_config
                    ? file(multiqc_config, checkIfExists: true)
                    : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
                [],
                [],
            ]
        }
    )

    emit:
    multiqc_report = MULTIQC.out.report.map { _meta, report -> [report] }.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                                                    // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/