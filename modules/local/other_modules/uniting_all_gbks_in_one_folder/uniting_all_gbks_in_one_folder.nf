#!/usr/bin/env nextflow

process UNITING_ALL_GBKS {

    publishDir "results", mode: 'copy'

    input:
    path gbk_files

    output:
    path "all_BGCs/*.gbk", emit: bgc_dir

    script:
    """
    mkdir -p all_BGCs

    for f in ${gbk_files}; do
        if [[ "\$f" == *_*.gbk ]]; then
            cp "\$f" all_BGCs/
        fi
    done
    """
}


