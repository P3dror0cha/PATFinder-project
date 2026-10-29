process ANTISMASH_REGION_VERIFICATION {
    container 'antismash/standalone:8.0.4'
    label 'antismash_regions'
    publishDir "results/antismash_regions", mode: 'copy'

    input:
    path gbk_files 

    output:
    path "output_gbks/*.gbk", emit: renamed_gbks

    script: 
    '''
    #!/bin/bash
    # This script prevents problems from the same sample with different regions
    # (multiple files from antismash).

    mkdir -p output_gbks

    for file in *.gbk; do

        [ -e "$file" ] || continue

        # Usando [.] para representar o ponto literal, evitando problemas de escape
        if [[ "$file" =~ ^sample_[0-9]+_[0-9]+[.]region[0-9]+[.]gbk$ ]]; then

            new_name="${file/.region/_region}"
            
            cp -L "$file" "output_gbks/$new_name"
        else
            cp -L "$file" "output_gbks/$file" 
        fi
    done
    '''
}