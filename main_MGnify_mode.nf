include { CREATE_CONDA_ENVS } from './modules/local/other_modules/create_conda_envs/create_conda_envs.nf'
include { MGNIFY_DOWNLOAD } from './modules/local/other_modules/metagenome_download/metagenomes_download.nf'
include { MGNIFY_PLOTS } from './modules/local/plots_and_graphs/MGnify_plot/MGnify_plot.nf'
include { DOWNLOAD_ANTISMASH_DB } from './modules/local/databases_download/antismash_db/main.nf'
include { DOWNLOAD_PFAM_DB } from './modules/local/databases_download/pfam_db/main.nf'
include { DOWNLOAD_KOFAM_DB } from './modules/local/databases_download/kofam_db/main.nf'
include { ANTISMASH } from './modules/local/antismash/main.nf'
include { BGC_PLOTS } from './modules/local/plots_and_graphs/bgc_plot/bgc_plots.nf'
include { BIGSCAPE } from './modules/local/bigscape/main.nf'
include { UNITING_ALL_GBKS } from './modules/local/other_modules/uniting_all_gbks_in_one_folder/uniting_all_gbks_in_one_folder.nf'
include { FILTERING_BIGSCAPE_RESULTS } from './modules/local/filtering_results/bigscape_results/filtering_bigscape_results.nf'
include { DEEPSEA } from './modules/local/deepsea/main.nf'
include { FILTERING_DEEPSEA_RESULTS } from './modules/local/filtering_results/deepsea_results/filtering_deepsea_results.nf'
include { RENAME_FASTA } from './modules/local/other_modules/rename_fasta/rename_fasta.nf'
include { KOFAM } from './modules/local/kofam/main.nf'
include { FILTERING_KOFAM_RESULTS } from './modules/local/filtering_results/kofam_results/filtering_kofam_results.nf'
include { RENAME_FAA } from './modules/local/other_modules/rename_faa/rename_faa.nf'
include { EXTRACTING_GBKS_SEQUENCES } from './modules/local/extracting_gbks_information/extracting_gbks_sequences/main.nf'
include { EXTRACTING_GBKS_CDS } from './modules/local/extracting_gbks_information/extracting_gbks_cds/main.nf'
include { DIAMOND } from './modules/local/diamond/main.nf'
include { FILTERING_DIAMOND_RESULTS } from './modules/local/filtering_results/diamond_results/filtering_diamond_results.nf'
include { CONCAT_ALL_RESULTS } from './modules/local/other_modules/concat_all_results/concat_all_results.nf'
include { DOWNLOAD_AMRFINDER_DB } from './modules/local/databases_download/amrfinder_db/main.nf'
include { AMRFINDER } from './modules/local/amrfinder/main.nf'
include { FILTERING_AMRFINDER_RESULTS } from './modules/local/filtering_results/amrfinder_results/filtering_amrfinder_results.nf'

workflow {

    if (params.create_conda_envs) {
        CREATE_CONDA_ENVS()
        ch_envs_done = CREATE_CONDA_ENVS.out.done
    } else {
        ch_envs_done = Channel.of('skip')
    }   
    // ========================================================================
    // 1. DOWNLOAD METAGENOMES
    // ========================================================================
    results = MGNIFY_DOWNLOAD(params.max_pages, params.output_prefix, ch_envs_done.first())
    MGNIFY_PLOTS(results.metadata, params.output_prefix)
    

    // ========================================================================
    // 2. PREPARING CHANNELS AND METADATA
    // ========================================================================
    ch_fna_by_name = results.fna.flatten().map { tuple(it.baseName, it) }
    ch_faa_by_name = results.faa.flatten().map { tuple(it.baseName, it) }
    
    ch_meta_files = ch_fna_by_name
        .join(ch_faa_by_name) 
        .toList() 
        .flatMap { it.withIndex() } 
        .map { item, index ->

            def (original_name, fna_file, faa_file) = item
            
            def clean_name = String.format("sample_%08d", index + 1)
            def meta = [ id: clean_name, original_id: original_name ]
            
            return tuple(meta, fna_file, faa_file)
        }

    ch_fna_concat = ch_meta_files.map { meta, fna, faa -> tuple(meta, fna) }
    ch_faa_concat     = ch_meta_files.map { meta, fna, faa -> tuple(meta, faa) }

    ids_correlation = ch_fna_concat
        .map { meta, file -> "${meta.id}\t${meta.original_id}\n" }
        .collectFile(
            name: 'ids_correlation.tsv', 
            storeDir: 'results', 
            seed: "Padronized_ID\tOriginal_ID\n"
        )

    // ========================================================================
    // 3. INSTALLING DATABASES
    // ========================================================================
    DOWNLOAD_ANTISMASH_DB()
    DOWNLOAD_PFAM_DB()
    DOWNLOAD_KOFAM_DB()
    DOWNLOAD_AMRFINDER_DB(params.amrfinder_db_preinstalled)

    // ========================================================================
    // 4. ANTISMASH PIPELINE
    // ========================================================================
    ch_db_antismash = DOWNLOAD_ANTISMASH_DB.out.db_antismash
    ch_renamed_fasta = RENAME_FASTA(ch_fna_concat, ch_envs_done.first())

    ch_antismash_input = ch_renamed_fasta.combine(ch_db_antismash)

    antismash_output = ANTISMASH(ch_antismash_input)
    
    ch_gbks = antismash_output.gbk.map { meta, gbk -> gbk }.collect()

    bgc_plots = BGC_PLOTS(ch_gbks)

    // ========================================================================
    // 5. DEEPSEA PIPELINE
    // ========================================================================
    ch_renamed_faa = RENAME_FAA(ch_faa_concat)
    ch_all_renamed_faa = ch_renamed_faa.map { meta, faa -> faa }.collect()

    deepsea = DEEPSEA(ch_all_renamed_faa)
    deepsea_results = FILTERING_DEEPSEA_RESULTS(
        deepsea.deepsea_merged_table,
        ch_all_renamed_faa,
        ch_gbks
    )

    // ========================================================================
    // 6. BiG-SCAPE ANALYSIS
    // ========================================================================
    all_BGCs = UNITING_ALL_GBKS(ch_gbks)

    bigscape_output = BIGSCAPE(all_BGCs.bgc_dir, file(params.pfam_db))
    filtered_bigscape_results = FILTERING_BIGSCAPE_RESULTS(bigscape_output.bigscape_fullnetwork, ids_correlation)

    // ========================================================================
    // 7. KOFAM PIPELINE
    // ========================================================================
    cds_from_all_gbks = EXTRACTING_GBKS_CDS(ch_gbks)
    kofam = KOFAM(cds_from_all_gbks)
    kofam_results = FILTERING_KOFAM_RESULTS(kofam.kofam_output, filtered_bigscape_results.filtered_bigscape_results)

    // ========================================================================
    // 8. DIAMOND ALIGNMENT (RESISTANCE PROTEINS AND ANTIBIOTIC BIOSYNTHESIS PROTEINS)
    // ========================================================================
    diamond_results = DIAMOND(cds_from_all_gbks)
    filtered_diamond_results = FILTERING_DIAMOND_RESULTS(diamond_results.diamond_result)

    // ========================================================================
    // 9. AMRFINDER PIPELINE
    // ========================================================================

    amrfinder_results = AMRFINDER(cds_from_all_gbks, DOWNLOAD_AMRFINDER_DB.out.db_ready)
    filtered_amrfinder_results = FILTERING_AMRFINDER_RESULTS(amrfinder_results.amrfinder_result)

    // ========================================================================
    // 10. UNITING ALL RESULTS 
    // ========================================================================
    concat = CONCAT_ALL_RESULTS(
        all_BGCs,
        deepsea_results.deepsea_csv,
        kofam_results.kofam_bigscape_filtered_results,
        filtered_diamond_results.diamond_filtered_result,
        ids_correlation,
        filtered_amrfinder_results.amrfinder_filtered_result
        )

}

