#!/usr/bin/env nextflow

process DOWNLOAD_AMRFINDER_DB {

    label 'download_amrfinder_db'

    publishDir "results/amrfinder_db", mode: 'copy'

    input:
    val db_preinstalled

    output:
    path "amrfinder_db_downloaded", emit: db_ready

    script:

    if (db_preinstalled) {
        """
        echo "AMRFinderPlus database already installed. Skipping download."
        touch amrfinder_db_downloaded
        """
    } else {
        """
        echo "Downloading AMRFinderPlus database..."
        amrfinder -U
        touch amrfinder_db_downloaded
        """
    }
}