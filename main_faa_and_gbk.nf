include { MGNIFY_DOWNLOAD } from './nextflow_modules/metagenomes_download.nf'
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

params.gbk_files = null
params.faa_files = null

workflow {

    // ========================================================================
    // 1 & 2. READING INPUTS, PAIRING AND PREPARING METADATA
    // ========================================================================

    if (!params.gbk_files || !params.faa_files) {
        error "ERRO: Please, provide the path to your gbk and faa files. Use --gbk_files and --faa_files."
    }

    def ch_gbk_by_name = Channel.fromPath(params.gbk_files)
        .map { file ->
            def id = file.name.replaceAll(/_[0-9]+\.region[0-9]+\.gbk$/, "")
            tuple(id, file)
        }
        .groupTuple() 

    def ch_faa_by_name = Channel.fromPath(params.faa_files)
        .map { file -> tuple(file.baseName, file) }

    def ch_paired_inputs = ch_gbk_by_name
        .join(ch_faa_by_name)
        .toList()                  
        .flatMap { it.withIndex() } 
        
    ch_faa_concat = ch_paired_inputs.map { item, index ->
        def (original_name, gbk_files_list, faa_file) = item
        def clean_name = String.format("sample_%08d", index + 1)
        def meta = [ id: clean_name, original_id: original_name ]
        
        return tuple(meta, faa_file) 
    }

    ch_gbks = ch_paired_inputs.flatMap { item, index ->
        def (original_name, gbk_files_list, faa_file) = item
        return gbk_files_list
    }.collect()

    ch_faa_concat
        .map { meta, file -> "${meta.id}\t${meta.original_id}\n" }
        .collectFile(
            name: 'ids_correlation.tsv', 
            storeDir: 'results', 
            seed: "Padronized_ID\tOriginal_ID\n"
        )

    // ========================================================================
    // 4. DEEPSEA PIPELINE (Pula a etapa 3 do AntiSMASH)
    // ========================================================================
    ch_renamed_faa = RENAME_FAA(ch_faa_concat)
    ch_all_renamed_faa = ch_renamed_faa.map { meta, faa -> faa }.collect()

    deepsea = DEEPSEA(ch_all_renamed_faa)
    deepsea_results = FILTERING_DEEPSEA_RESULTS(deepsea.deepsea_merged_table, ch_all_renamed_faa, ch_gbks)

    // ========================================================================
    // 5. BiG-SCAPE ANALYSIS
    // ========================================================================
    all_BGCs = UNITING_ALL_GBKS(ch_gbks)

    bigscape_output = BIGSCAPE(all_BGCs.bgc_dir, file(params.pfam_db))
    filtered_bigscape_results = FILTERING_BIGSCAPE_RESULTS(bigscape_output.bigscape_fullnetwork)

    // ========================================================================
    // 6. KOFAM PIPELINE
    // ========================================================================
    cds_from_all_gbks = EXTRACTING_GBKS_CDS(ch_gbks)
    kofam = KOFAM(cds_from_all_gbks)
    kofam_results = FILTERING_KOFAM_RESULTS(kofam.kofam_output, filtered_bigscape_results.filtered_bigscape_results)

    // ========================================================================
    // 6. POEM PIPELINE (OPERON)
    // ========================================================================
    sequences_from_all_gbks = EXTRACTING_GBKS_SEQUENCES(ch_gbks)
    poem_output = POEM(sequences_from_all_gbks.bgc_sequences)
    poem_results = FILTERING_POEM_RESULTS(poem_output.operon_file)

    // ========================================================================
    // 7. DIAMOND ALIGNMENT (RESISTANCE PROTEINS AND ANTIBIOTIC BIOSYNTHESIS PROTEINS)
    // ========================================================================
    diamond_results = DIAMOND(cds_from_all_gbks)
}