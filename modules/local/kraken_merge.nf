process KRAKEN_PARSE {
    tag "${meta.id}"
    label 'process_single'

    conda "conda-forge::pandas=1.4.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.4.3' :
        'quay.io/biocontainers/pandas:1.4.3' }"

    input:
    path kraken2_reports

    output:
    path("*.csv") into kraken_merged_report

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    kraken_merge.py \\
        -or ${prefix}.kraken2_merged_report.csv \\
    """

}
