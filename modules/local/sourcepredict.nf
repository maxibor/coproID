process SOURCEPREDICT {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::sourcepredict=0.5" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'quay.io/biocontainers/YOUR-TOOL-HERE' }"

    input:
    tuple val(meta), path(sink)
    path(sources)
    path(labels)

    output:
    tuple val(meta), path("*.sourcepredict.csv"), emit: prediction,
    tuple val(meta), path("*_embedding.csv"),     emit: embedding
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    sourcepredict \\
        $args \\
        -s $sources
        -l $labels
        -t $task.cpus \\
        -o $prefix \\
        -e ${prefix}_embedding.csv \\
        $sink

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sourcepredict: \$(echo \$(sourcepredict --help | grep "^SourcePredict v") | sed 's/SourcePredict v//g' ))
    END_VERSIONS
    """
}
