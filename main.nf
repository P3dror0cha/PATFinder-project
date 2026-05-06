include { MGNIFY_DOWNLOAD } from './nextflow_modules/metagenomes_download.nf'
include { ANTISMASH } from './nextflow_modules/antismash.nf'
include { BIGSCAPE } from './nextflow_modules/bigscape.nf'
include { UNITING_ALL_GBKS } from './nextflow_modules/uniting_all_gbks_in_one_folder.nf'
include { FILTERING_BIGSCAPE_RESULTS } from './nextflow_modules/filtering_bigscape_results.nf'
include { DEEPSEA } from './nextflow_modules/deepsea.nf'

workflow {

    results = MGNIFY_DOWNLOAD(params.max_pages, params.output_prefix)

    deepsea_output = DEEPSEA(results.faa.collect())

    antismash_output = ANTISMASH(results.fna.flatten())

    all_BGCs = UNITING_ALL_GBKS(antismash_output.collect())
    all_BGCs.bgc_dir.view { "Collected GBK files: ${it.size()}" }

    bigscape_output = BIGSCAPE(all_BGCs.bgc_dir, file(params.pfam_db))

    FILTERING_BIGSCAPE_RESULTS(bigscape_output.bigscape_fullnetwork)
}
