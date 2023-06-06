process KRAKEN_PARSE {
    tag "${meta.id}"
    label 'process_single'

    conda "conda-forge::pandas=1.4.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.4.3' :
        'quay.io/biocontainers/pandas:1.4.3' }"

    input:
    path kraken2_report

    output:
    path("*.read_kraken_parsed.csv") into kraken_read_count
    path("*.kmer_kraken_parsed.csv") into kraken_kmer_count

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    kraken_parse.py $kraken2_report
    """

}
