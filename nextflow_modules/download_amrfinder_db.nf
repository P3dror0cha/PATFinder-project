#!/usr/bin/env nextflow

process DOWNLOAD_AMRFINDER_DB {
    label 'download_amrfinder_db'
    publishDir "results/amrfinder_db", mode: 'copy'

    output:
    path "amrfinder_db_downloaded", emit: db_ready

    script:
    """
    amrfinder -U
    touch amrfinder_db_downloaded
    """
}