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
include { GUNZIP as GUNZIP_FASTA; GUNZIP as GUNZIP_GTF } from '../modules/nf-core/gunzip/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_analysisrnaseq_pipeline'
include { SAMTOOLS_SORT             } from '../modules/nf-core/samtools/sort/main'
include { PICARD_MARKDUPLICATES     } from '../modules/nf-core/picard/markduplicates/main'
include { SUBREAD_FEATURECOUNTS     } from '../modules/nf-core/subread/featurecounts/main'
include { TPM_TABLE                 } from '../modules/local/tpm_table/main'
include { PLOT_TPM                  } from '../modules/local/plot_tpm/main' 

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
    if (params.fasta.endsWith('.gz')) {
        GUNZIP_FASTA(ch_fasta)
        ch_fasta = GUNZIP_FASTA.out.gunzip
    }
    def ch_gtf = channel.value([ [id: 'genome'], file(params.gtf, checkIfExists: true) ])
    if (params.gtf.endsWith('.gz')) {
        GUNZIP_GTF(ch_gtf)
        ch_gtf = GUNZIP_GTF.out.gunzip
    }

    //
    // MODULE: Extract splice sites from the GTF (needed for spliced RNA-seq reads)
    //
    HISAT2_EXTRACTSPLICESITES(ch_gtf)
    def ch_splicesites = HISAT2_EXTRACTSPLICESITES.out.txt

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
            '200.GB'  // min memory for building the index with splice sites
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
    // MODULE: Sort BAM files by genomic coordinate (required by Picard)
    //
    SAMTOOLS_SORT(
        ch_bam,
        [ [:], [], [] ],   // no reference needed for BAM output
        ''                 // no index here, Picard creates it after marking
    )
    def ch_bam_sorted = SAMTOOLS_SORT.out.bam

    //
    // MODULE: Mark duplicate reads (duplicates are flagged, not removed)
    //
    PICARD_MARKDUPLICATES(
        ch_bam_sorted,
        [ [:], [], [] ]    // no reference needed for BAM input
    )
    ch_multiqc_files = ch_multiqc_files.mix(PICARD_MARKDUPLICATES.out.metrics.map{ _meta, file -> file })

    // marked BAM (+ index) for the next step (quantification)
    def ch_bam_markdup = PICARD_MARKDUPLICATES.out.bam

    //
    // MODULE: Count reads per gene (featureCounts)
    // combine() attaches the single GTF file to every BAM -> [ meta, bam, gtf ]
    // Strandedness and single/paired-end are read from meta by the module
    //
    SUBREAD_FEATURECOUNTS(
        ch_bam_markdup.combine( ch_gtf.map { _meta, gtf -> gtf } )
    )
    ch_multiqc_files = ch_multiqc_files.mix(SUBREAD_FEATURECOUNTS.out.summary.map{ _meta, file -> file })

        //
    // MODULE (local): Convert counts to TPM and merge all samples into one table
    // collect() waits until all samples are counted and passes all files at once
    //
    TPM_TABLE(
        SUBREAD_FEATURECOUNTS.out.counts.map { _meta, file -> file }.collect()
    )

    //
    // MODULE (local): Plot TPM overview (top genes, distribution, sample correlation)
    // Files ending in _mqc.png are picked up by MultiQC and shown in the report
    //
    PLOT_TPM(TPM_TABLE.out.tpm)
    ch_multiqc_files = ch_multiqc_files.mix(PLOT_TPM.out.plots)
    
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