process DEEPBGC_PIPELINE {
    tag "${meta.id}"

    container 'quay.io/biocontainers/deepbgc:0.1.31--pyhca03a8a_0'

    input:
    tuple val(meta), path(fasta)
    path db

    output:
    tuple val(meta), path("${meta.id}/LOG.txt"), emit: log
    tuple val(meta), path("${meta.id}/${meta.id}.full.gbk"), emit: full_gbk
    tuple val(meta), path("${meta.id}/${meta.id}.bgc.gbk"), optional: true, emit: bgc_gbk
    tuple val(meta), path("${meta.id}/${meta.id}.bgc.tsv"), optional: true, emit: bgc_tsv
    tuple val(meta), path("${meta.id}/${meta.id}.pfam.tsv"), optional: true, emit: pfam_tsv
    tuple val(meta), path("${meta.id}/${meta.id}.antismash.json"), optional: true, emit: antismash_json
    path 'versions.yml', emit: versions

    script:
    def args = task.ext.args ?: ''

    """
    export DEEPBGC_DOWNLOADS_DIR="\$PWD/${db}"

    deepbgc pipeline \\
        ${args} \\
        --output "${meta.id}" \\
        "${fasta}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deepbgc: \$(python -c "from pkg_resources import get_distribution; print(get_distribution('deepbgc').version)")
    END_VERSIONS
    """
}