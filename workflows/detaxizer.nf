/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap; fromSamplesheet } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

WorkflowDetaxizer.initialise(params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

nextflow.enable.dsl=2

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { KRAKEN2PREPARATION } from '../modules/local/kraken2preparation'
include { ISOLATE_IDS_FROM_KRAKEN2_TO_BLASTN } from '../modules/local/isolate_ids_from_kraken2_to_blastn'
include { PREPARE_FASTA4BLASTN } from '../modules/local/prepare_fasta4blastn'
include { FILTER_BLASTN_IDENTCOV } from '../modules/local/filter_blastn_identcov'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                      } from '../modules/nf-core/fastqc/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { FASTP } from '../modules/nf-core/fastp/main'
include { KRAKEN2_KRAKEN2 } from '../modules/nf-core/kraken2/kraken2/main'
include { BLAST_BLASTN } from '../modules/nf-core/blast/blastn/main'
include { BLAST_MAKEBLASTDB } from '../modules/nf-core/blast/makeblastdb'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow DETAXIZER {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    //INPUT_CHECK (
    //    file(params.input)
    //)
    // TODO: OPTIONAL, you can use nf-validation plugin to create an input channel from the samplesheet with Channel.fromSamplesheet("input")
    // See the documentation https://nextflow-io.github.io/nf-validation/samplesheets/fromSamplesheet/
    // ! There is currently no tooling to help you write a sample sheet schema
    // TODO: Make the map for single ended 
    ch_input = Channel.fromSamplesheet('input')
        .map{meta, fastq_1, fastq_2 -> [meta,[fastq_1, fastq_2]]}
    
    //
    // MODULE: Run FastQC
    //
    FASTQC (
        ch_input
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    //
    // MODULE: Run fastp
    //

    FASTP (
        ch_input,
        [],
        params.fastp_save_trimmed_fail,
        []
    )
    ch_versions = ch_versions.mix(FASTP.out.versions)

    //
    // MODULE: Prepare Kraken2 Database
    //

    KRAKEN2PREPARATION (
        params.kraken2db
    )
    ch_versions = ch_versions.mix(KRAKEN2PREPARATION.out.versions)

    // 
    // MODULE: Run Kraken2
    //

    KRAKEN2_KRAKEN2 (
        FASTP.out.reads,
        KRAKEN2PREPARATION.out.db,
        params.save_output_fastqs,
        params.save_reads_assignment
    )
    ch_versions = ch_versions.mix(KRAKEN2_KRAKEN2.out.versions)

    //    
    // MODULE: Isolate the hits for a certain taxa 
    //

    ISOLATE_IDS_FROM_KRAKEN2_TO_BLASTN (
        KRAKEN2_KRAKEN2.out.classified_reads_assignment
    )
    ch_versions = ch_versions.mix(ISOLATE_IDS_FROM_KRAKEN2_TO_BLASTN.out.versions)
    
    //
    // MODULE: Extract the hits to fasta format
    //

    ch_combined = FASTP.out.reads
            .join(
                ISOLATE_IDS_FROM_KRAKEN2_TO_BLASTN.out.classified, by: 0
                ) 
    // TODO: REPLACE awk AND seqtk WITH TOOLS THAT CAN DISPLAY THEIR VERSION NUMBER
    PREPARE_FASTA4BLASTN(
        ch_combined
    )
    ch_versions = ch_versions.mix(PREPARE_FASTA4BLASTN.out.versions)    

    // 
    // MODULE: Run BLASTN if --skip_blastn = true
    //
    ch_reference_fasta = Channel.empty()
    if (!params.skip_blastn) {  // If skip_blastn is false, then execute the process
        ch_reference_fasta =  file( params.fasta )
        }
    
    BLAST_MAKEBLASTDB (
            ch_reference_fasta
        )
    ch_versions = ch_versions.mix(BLAST_MAKEBLASTDB.out.versions)

    ch_fasta4blastn = PREPARE_FASTA4BLASTN.out.fasta
        .flatMap { meta, fastaList ->
            [ 
                [['id': "${meta.id}_R1"], fastaList[0]],
                [['id': "${meta.id}_R2"], fastaList[1]]
            ]
        }


    BLAST_BLASTN (
            ch_fasta4blastn,
            BLAST_MAKEBLASTDB.out.db
        )
    ch_versions = ch_versions.mix(BLAST_BLASTN.out.versions)
    
    FILTER_BLASTN_IDENTCOV (
        BLAST_BLASTN.out.txt
    )
    ch_versions = ch_versions.mix(FILTER_BLASTN_IDENTCOV.out.versions)
    

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowDetaxizer.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowDetaxizer.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
