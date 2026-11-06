process DEEPBGC_DOWNLOAD {
    label 'process_low'
    container 'quay.io/biocontainers/deepbgc:0.1.31--pyhca03a8a_0'

    output:
    path 'deepbgc_db', emit: db
    path 'versions.yml', emit: versions

    script:
    """
    export DEEPBGC_DOWNLOADS_DIR="\$PWD/deepbgc_db"

    deepbgc download

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deepbgc: \$(python -c "from pkg_resources import get_distribution; print(get_distribution('deepbgc').version)")
    END_VERSIONS
    """
}