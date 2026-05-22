#!/usr/bin/env nextflow

process DEEPSEA {
    label 'deepsea'
    publishDir "results/deepsea", mode: 'copy'

    input:
    path cds_files  

    output:
    path "deepsea.table", emit: deepsea_table
    path "merged_deepsea_results.tsv", emit: deepsea_merged_table

script: 
    """
    ln -s ${projectDir}/DeepSEA-project/class-encoder .
    ln -s ${projectDir}/DeepSEA-project/models .
    ln -s ${projectDir}/DeepSEA-project/deepsea .

    cat *.faa > merged_metagenome_cds_raw.faa
    
    python -c "
import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
import tensorflow as tf
from Bio import SeqIO

model = tf.keras.models.load_model('${projectDir}/DeepSEA-project/models/cnn-model')

alphabet = 'ACDEFGHIKLMNPQRSTVWY'
valid_chars = set()

for char in alphabet:
    try:
        model.predict(tf.convert_to_tensor([char]), verbose=0)
        valid_chars.add(char)
    except Exception:
        pass

banned = set(alphabet) - valid_chars
valid_str = ''.join(sorted(valid_chars))
banned_str = ''.join(sorted(banned))

with open('merged_metagenome_cds.faa', 'w') as out:
    for r in SeqIO.parse('merged_metagenome_cds_raw.faa', 'fasta'):
        seq = str(r.seq).upper().replace('*', '')
        clean_seq = ''.join([c if c in valid_chars else 'A' for c in seq])
        if len(clean_seq) > 0:
            out.write('>' + str(r.id) + chr(10) + clean_seq + chr(10))
"

    python ${projectDir}/DeepSEA-project/DeepSEA.py run \\
        --input merged_metagenome_cds.faa \\
        --outname deepsea.table
    
    cp deepsea.table.tsv deepsea.table

    tsv_files=\$(ls *.tsv | grep -v merged_deepsea_results.tsv || true)
    
    if [ -n "\$tsv_files" ]; then
        cat \$tsv_files > merged_deepsea_results.tsv
    else
        touch merged_deepsea_results.tsv
    fi
    """
}
