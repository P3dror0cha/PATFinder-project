# Notebooks Descriptions:

## 1. BGC_architecture:

### Name:
bgc_arquitecture_analysis.ipynb

### Objective:
- Notebook dedicated to the analysis of the arquitecture from BGC's metagenomes.
- Makes the following graphs: Frequency, length, CDS Count, Gene Kind Count, Number of Domains of BGC by classes.
- Makes a PCA analysis of the domains composition of the BGC population
   

### Input:
- A folder with all .gbk from antiSMASH. 

### Processament:
- Separated in seven parts. Most parts produce descritive graphs of BGC's population characteristics 

### Output:
- Generate six images about the BGC's population.

### Observations:
- The PCA part (last one) do not function. But the graphs still can be useful.


## 2. BIG-SCAPE:

### Name:
BIG-SCAPE_exploratory_notebook.ipynb

### Objective:
- This notebook associates each BGC in the sample with the nearest MIBIG BGC.
- Makes the following graphs: boxplot of distance, jaccard, adjacency and dss.
- Finally, it returns a dataframe containing the MIBiG database BGC that most closely matches each BGC identified in the sample.
   
### Input:
- A file named XXX_full.network produced by BIG-SCAPE, providing similarity metrics between BGC pairs.

### Processament:
- Process the input data by applying the specified filtering criteria.

### Output:
- Generate a Boxplot graph and a table with BGCs similarity (from sample and MIBIG).

### Observations:
- The command used for BIG-SCAPE was:
- python3 bigscape.py cluster  -i /file/with/all_BGCs/together -o /your/directory/output --include-singletons -p /path/to/pfam_database/Pfam-A.hmm --mibig-version 4.0

## 3. DeepSEA:

### Name:
deepsea_results.ipynb

### Objective:
- This notebook is designed to construct a dataframe comprising BGCs associated with resistance proteins.
- Makes the following graphs: "Distribution by resistance class", "Number of resistance proteins in each metagenome", "Confusion matrix with BGC class and Resistance class".
   

### Input:
- All DeepSEA .tsv outputs. The file will have the following columns: "Name", "Class", "Prob".
Obs: DeepSEA will always produce one tsv for each metagenomic sample.
- Files with annotated CDS from each metagenome sample (.faa format).
Obs: All .faa files in one folder. 

### Processament:
- Take the raw result from DeepSEA and filters the dataframe.

### Output:
- Generate the graphs above and the filtered dataframe.

### Observations:
- The command used with DeepSEA was: python DeepSEA.py run --input test/test.fasta --outname test.table
