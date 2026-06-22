include { MGNIFY_DOWNLOAD } from './nextflow_modules/metagenomes_download.nf'
include { DOWNLOAD_ANTISMASH_DB } from './nextflow_modules/download_antismash_db.nf'
include { DOWNLOAD_PFAM_DB } from './nextflow_modules/download_pfam_db.nf'
include { DOWNLOAD_KOFAM_DB } from './nextflow_modules/download_kofam_db.nf'
include { ANTISMASH } from './nextflow_modules/antismash.nf'
include { BIGSCAPE } from './nextflow_modules/bigscape.nf'
include { UNITING_ALL_GBKS } from './nextflow_modules/uniting_all_gbks_in_one_folder.nf'
include { FILTERING_BIGSCAPE_RESULTS } from './nextflow_modules/filtering_bigscape_results.nf'
include { DEEPSEA } from './nextflow_modules/deepsea.nf'
include { FILTERING_DEEPSEA_RESULTS } from './nextflow_modules/filtering_deepsea_results.nf'
include { RENAME_FASTA } from './nextflow_modules/rename_fasta.nf'
include { KOFAM } from './nextflow_modules/kofam.nf'
include { FILTERING_KOFAM_RESULTS } from './nextflow_modules/filtering_kofam_results.nf'
include { RENAME_FAA } from './nextflow_modules/rename_faa.nf'
include { EXTRACTING_GBKS_SEQUENCES } from './nextflow_modules/extracting_gbks_sequences.nf'
include { EXTRACTING_GBKS_CDS } from './nextflow_modules/extracting_gbks_cds.nf'
include { POEM } from './nextflow_modules/POEM_pipeline.nf'
include { FILTERING_POEM_RESULTS } from './nextflow_modules/filtering_POEM_results.nf'
include { DIAMOND } from './nextflow_modules/diamond.nf'
include { FILTERING_DIAMOND_RESULTS } from './nextflow_modules/filtering_diamond_results.nf'

workflow {

    // ========================================================================
    // 1. DOWNLOAD METAGENOMES
    // ========================================================================
    results = MGNIFY_DOWNLOAD(params.max_pages, params.output_prefix)
    

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

    ch_fna_concat
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

    // ========================================================================
    // 4. ANTISMASH PIPELINE
    // ========================================================================
    ch_db_antismash = DOWNLOAD_ANTISMASH_DB.out.db_antismash
    ch_renamed_fasta = RENAME_FASTA(ch_fna_concat)

    ch_antismash_input = ch_renamed_fasta.combine(ch_db_antismash)

    antismash_output = ANTISMASH(ch_antismash_input)
    
    ch_gbks = antismash_output.gbk.map { meta, gbk -> gbk }.collect()

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
    filtered_bigscape_results = FILTERING_BIGSCAPE_RESULTS(bigscape_output.bigscape_fullnetwork)

    // ========================================================================
    // 7. KOFAM PIPELINE
    // ========================================================================
    cds_from_all_gbks = EXTRACTING_GBKS_CDS(ch_gbks)
    kofam = KOFAM(cds_from_all_gbks)
    kofam_results = FILTERING_KOFAM_RESULTS(kofam.kofam_output, filtered_bigscape_results.filtered_bigscape_results)

    // ========================================================================
    // 8. POEM PIPELINE (OPERON)
    // ========================================================================
    sequences_from_all_gbks = EXTRACTING_GBKS_SEQUENCES(ch_gbks)
    poem_output = POEM(sequences_from_all_gbks.bgc_sequences)
    poem_results = FILTERING_POEM_RESULTS(poem_output.operon_file)

    // ========================================================================
    // 9. DIAMOND ALIGNMENT (RESISTANCE PROTEINS AND ANTIBIOTIC BIOSYNTHESIS PROTEINS)
    // ========================================================================
    diamond_results = DIAMOND(cds_from_all_gbks)
    filtered_diamond_results = FILTERING_DIAMOND_RESULTS(diamond_results.diamond_result)

    // ========================================================================
    // 10. UNITING ALL RESULTS 
    // ========================================================================

}

