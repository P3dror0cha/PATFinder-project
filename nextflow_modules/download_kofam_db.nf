#!/usr/bin/env nextflow

process DOWNLOAD_KOFAM_DB {
    label 'download_kofam_db'
    storeDir "${params.kofam_db}"

    output:
    path "ko_list", emit: ko_list
    path "profiles", emit: profiles

    script: 
    """
    wget ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz
    gunzip ko_list.gz

    wget ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz
    tar -xzf profiles.tar.gz
    rm profiles.tar.gz
    """
}
