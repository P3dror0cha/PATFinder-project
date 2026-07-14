process CONCAT_ALL_RESULTS {
    label 'final_results'
    publishDir "results/final_results", mode: 'copy'

    input:
    path folder_with_gbks
    path poem_final_results
    path deepsea_final_results
    path kofam_final_results
    path diamond_final_results
    path ids_correlation

    output:
    path "concat_all_results.csv", emit: final_result
    path "bgc_class_correlation.csv", emit: bgc_class_correlation

script:
"""
python3 ${projectDir}/bin/Concat_all_results_process.py \\
    --gbk ${folder_with_gbks} \\
    --poem ${poem_final_results} \\
    --deepsea ${deepsea_final_results} \\
    --kofam ${kofam_final_results} \\
    --diamond ${diamond_final_results} \\
    --output concat_all_results.csv \\
    --ids-correlation ${ids_correlation}
"""
}
