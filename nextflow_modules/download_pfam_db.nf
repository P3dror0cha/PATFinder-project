#!/usr/bin/env nextflow

process DOWNLOAD_PFAM_DB {
    label 'download_pfam_db'
    storeDir "${params.pfam_db}"

    output:
    path "Pfam-A.hmm"   , emit: pfam_hmm

    script: 
    """
    wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
    gunzip Pfam-A.hmm.gz
    """
}
