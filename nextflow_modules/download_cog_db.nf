#!/usr/bin/env nextflow

process DOWNLOAD_COG_DB {
    label 'download_cog_db'
    storeDir "${params.cog_db}"

    output:
    path "cog", emit: cog_db

    script: 
    """
    mkdir -p cog/cog2014

    echo "Download COG database"
    wget -q -c ftp://ftp.ncbi.nih.gov/pub/COG/COG2014/data/*2014* -P cog/

    echo "Unzip cog fasta and format the file by diamond"
    gunzip cog/prot2003-2014.fa.gz
    
    diamond makedb \\
        --in cog/prot2003-2014.fa \\
        -d cog/prot2003-2014.fa
    """
}