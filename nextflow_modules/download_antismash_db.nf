#!/usr/bin/env nextflow

process DOWNLOAD_ANTISMASH_DB {
    label 'download_antismash_db'
    storeDir "${params.antismash_db}"

    output:
    path "db_antismash", emit: db_antismash

    script: 
    """
    if [ ! -f ${params.antismash_db}/.done ]; then
        mkdir -p ${params.antismash_db}
        download-antismash-databases --database-dir ${params.antismash_db}
        touch ${params.antismash_db}/.done
    fi
    touch db_antismash
    """
}