# Notebooks Descriptions:

## 1. BGC_architecture:

### Name:
**bgc_arquitecture_analysis.ipynb**

#### Objective:
- Notebook dedicated to the analysis of the arquitecture from BGC's metagenomes.
- Makes the following graphs: Frequency, length, CDS Count, Gene Kind Count, Number of Domains of BGC by classes.
- Makes a PCA analysis of the domains composition of the BGC population
   

#### Input:
- A folder with all .gbk from antiSMASH. 

#### Output:
- Generate six images about the BGC's population.

#### Observations:
- The PCA part (last one) do not function. But the graphs still can be useful.


## 2. BIG-SCAPE:

### Name:
**BIG-SCAPE_exploratory_notebook.ipynb**

#### Objective:
- This notebook associates each BGC in the sample with the nearest MIBIG BGC.
- Makes the following graphs: boxplot of distance, jaccard, adjacency and dss.
- Finally, it returns a dataframe containing the MIBiG database BGC that most closely matches each BGC identified in the sample.
   
#### Input:
- A file named XXX_full.network produced by BIG-SCAPE, providing similarity metrics between BGC pairs.

#### Output:
- Generate a Boxplot graph and a table with BGCs similarity (from sample and MIBIG).

#### Observations:
- The command used for BIG-SCAPE was:
- python3 bigscape.py cluster  -i /file/with/all_BGCs/together -o /your/directory/output --include-singletons -p /path/to/pfam_database/Pfam-A.hmm --mibig-version 4.0

## 3. DeepSEA:

### Name:
**deepsea_results.ipynb**

#### Objective:
- This notebook is designed to construct a dataframe comprising BGCs associated with resistance proteins.
- Makes the following graphs: "Distribution by resistance class", "Number of resistance proteins in each metagenome", "Confusion matrix with BGC class and Resistance class".
   
#### Input:
- All DeepSEA .tsv outputs. The file will have the following columns: "Name", "Class", "Prob".
Obs: DeepSEA will always produce one tsv for each metagenomic sample.
- Files with annotated CDS from each metagenome sample (.faa format).
Obs: All .faa files in one folder. 

#### Output:
- Generate the graphs above and the filtered dataframe.

#### Observations:
- The command used with DeepSEA was: python DeepSEA.py run --input test/test.fasta --outname test.table

## 4. Operon Finder:

### 1 - Name:
**processing_POEM_results.ipynb**

#### Objective:
- This notebook takes the output from POEM, specificaly the output named "input.fsa.operon" that has the information about the operon genes in each BGC.
   
#### Input:
-  input.fsa.operon. This file has all genes in sequence for each operon detected.

#### Output:
- POEM_operon_information.csv. Dataframe containing the following informations: "metagenome_ids", "strand" and "genes_for_each_operon"

#### Observations:
- The operon finder used was POEM pipeline (https://github.com/Rinoahu/POEM_py3k).
- The command used was: bash ./bin/run_poem.sh  -f /home/pedro/resultados_bgc/fna_dos_bgcs/BGCS.fasta -a n -p pka -l y

### 2 - Name:
**Uniting_POEM_with_BGC.ipynb**

#### Objective:
- Here we take the processed result from POEM pipeline and merge with the BGCs results (BIG-SCAPE, antiSMASH and DeepSEA results).
   
#### Input:
-  POEM_operon_information.csv (Table made in the "processing_POEM_results.ipynb" notebook)
-  A folder path that contains all .gbk files from BGCs.
-  Dataframe with BIG-SCAPE informations (BIG-SCAPE_exploratory_notebook.ipynb)
-  Dataframe with DeepSEA informations (deepsea_results.ipynb)

#### Output:
- updated_BGCs_information.csv (csv with all informations - BIG-SCAPE, antiSMASH, DeepSEA and POEM)

#### Observations:

## 5. Uniprot Alignment - Antibiotic biosynthesis:

### 1 - Name:
**gbk_cds_sequence_extraction.ipynb**

#### Objective:
- Here we take the output from antismash and create a dataframe with sequence information for future alignment.
   
#### Input:
- A folder with all .gbk files

#### Output:
- antismash_cds.faa. A file with all CDS from different samples.

### Observations:
- It is necessary for the BlastP alignment.

### 2 - Name:
**Uniprot_protein_download.ipynb**

#### Objective:
- This notebook uses the output from the UniProt API to perform alignments between CDSs from sample BGCs and proteins associated with “antibiotic biosynthesis.”
   
#### Input:
- None

#### Output:
- uniprot_table.csv. Makes a dataframe with parsed informations, such as: "accession", "submitted_name", "gene_id", "organism", "lineage", "sequence" and others.

### Observations:
- 

### 3 - Name:
**Uniprot_protein_alignment.ipynb**

#### Objective:
- This notebook takes the output from DIAMOND  and analise the alignment between the CDS from samples BGCs and sequences of antibiotic biosynthesis associated proteins (from Uniprot).

   
#### Input:
- Dataframe from **Uniprot_protein_download.ipynb**
- Aligment result from DIAMOND (output.tsv)

#### Output:
- Is incomplete.

### Observations:
- This notebook is incomplete.
- The command used for DIAMOND alignment was: "diamond blastp   -q antismash_cds.faa   -d uniprot_db   -o cds_vs_uniprot.tsv   -f 6   qseqid sseqid pident length qstart qend sstart send evalue bitscore   --evalue 1e-5   --max-target-seqs 5"


## 6. antiSMASH:

### Name:
**download_of_metagenomes.ipynb**

#### Objective:
- Notebook dedicated to characterize the metagenomes pool, using the python functions
   
#### Input:
- None

#### Output:
- Download the .fna (with DNA sequences) of each metagenome searched

### Observations:
- 
