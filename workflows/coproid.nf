/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowCoproid.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.fasta ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input       ) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }
if (params.genomes     ) { ch_genomes = Channel.fromPath(params.genomes) } else { exit 1, 'Genomes sheet not specified!' }
if (params.sam2lca_db  ) { ch_sam2lca_db = Channel.fromPath(params.sam2lca_db).first() } else { exit 1, 'SAM2LCA database path not specified!' }
if (params.kraken2_db  ) { ch_kraken2_db = Channel.fromPath(params.kraken2_db).first() } else { exit 1, 'Kraken2 database path not specified!' }
if (params.sp_sources  ) { ch_sp_sources = Channel.fromPath(params.sp_sources).first() } else { exit 1, 'SourcePredict sources file not specified!' }
if (params.ch_sp_labels) { ch_sp_labels = Channel.fromPath(params.ch_sp_labels).first() } else { exit 1, 'SourcePredict labels file not specified!' }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK               } from '../subworkflows/local/input_check'
include { PREPARE_GENOMES           } from '../subworkflows/local/prepare_genomes_indices'
include { ALIGN_INDEX               } from '../subworkflows/local/align_index'
include { MERGE_SORT_INDEX_SAMTOOLS } from '../subworkflows/local/merge_sort_index_samtools'
include { KRAKEN2_WF                } from '../subworkflows/local/kraken2_wf'

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                      } from '../modules/nf-core/fastqc/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { FASTP                       } from '../modules/nf-core/fastp/main'
include { BOWTIE2_BUILD               } from '../modules/nf-core/bowtie2/build/main'
include { SAMTOOLS_MERGE              } from '../modules/nf-core/samtools/merge/main'
include { SAM2LCA                     } from '../modules/nf-core/sam2lca/main'
include { SOURCEPREDICT               } from '../modules/local/sourcepredict'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow COPROID {

    ch_versions = Channel.empty()

    /*
    SUBWORKFLOW: Read in samplesheet, validate and stage input files
    */
    INPUT_CHECK (
        ch_input
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions.first())

    /*
    SUBWORKFLOW: Genomes in genome sheet, validate and stage genome fasta, and indices
    */
    PREPARE_GENOMES (
        genomes
    )
    ch_versions = ch_versions.mix(PREPARE_GENOMES.out.versions.first())

    FASTQC (
        INPUT_CHECK.out.reads
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    FASTP (
        INPUT_CHECK.out.reads
    )
    ch_versions = ch_versions.mix(FASTP.out.versions.first())

    /*
    Align all samples to all genomes
    */

    FASTP.out.reads_merged // [meta[ID, single_end], merged_reads]
        .combine(PREPARE_GENOMES.out.genomes) //[meta[genome_name], fasta, index]
        .map {
            meta_reads, reads, meta_genome, genome_fasta, genome_index ->
            [
                [
                    'id': meta_reads.id + '_' + meta_genome.genome_name,
                    'genome_name': meta_genome.genome_name,
                    'genome_taxid': meta_genome.taxid,
                    'genome_size': meta.genome_size,
                    'sample_name': meta_reads.id,
                    'single_end': meta_reads.single_end
                ],
                reads,
                genome_index
            ]
        }
        .set { ch_reads_genomes_index }

    /*
    SUBWORKFLOW: Align reads to genomes, and index alignments
    */

    ALIGN_INDEX (
        ch_reads_genomes_index
    )
    ch_versions = ch_versions.mix(ALIGN_INDEX.out.versions.first())


    ALIGN_BOWTIE2.out.bam.join(
        ALIGN_BOWTIE2.out.bai
    ).map {
        meta, bam, bai -> [['id':meta.sample_name], bam] // meta.id, bam
    }.groupTuple()
    .set { bams_synced }

    MERGE_SORT_INDEX_SAMTOOLS (
        bams_synced
    )

    SAM2LCA (
        MERGE_SORT_INDEX_SAMTOOLS.out.bam.join(
            MERGE_SORT_INDEX_SAMTOOLS.out.bai
        ),
        ch_sam2lca_db
    )

    KRAKEN2_WF(
        FASTP.out.reads_merged,
        ch_kraken2_db
    )

    SOURCEPREDICT (
        KRAKEN2_WF.out.kraken2_report,
        ch_sp_sources,
        ch_sp_labels
    )


    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowCoproid.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowCoproid.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
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
