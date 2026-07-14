# PATFinder-project
Repository dedicated to PATFinder, a bioinformatics tool designed to identify and explore potential antibiotic compounds within biosynthetic gene clusters (BGCs).

# BGC pipeline project
This project focuses on the development of a bioinformatics pipeline for the analysis of Biosynthetic Gene Clusters (BGCs), with an emphasis on identifying and associating them with antimicrobial potential.

The pipeline integrates the following tools:
- antiSMASH (https://github.com/antismash);
- BIG-SCAPE (https://github.com/medema-group/BiG-SCAPE);
- POEM-pipeline (https://github.com/Rinoahu/POEM_py3k);
- DeepSEA (https://github.com/tiagocabralborelli/DeepSEA-code);
- KOfam (https://github.com/takaram/kofam_scan);
- DIAMOND (https://github.com/bbuchfink/diamond);

Below is the current pipeline workflow.

<img width="4617" height="3054" alt="Imagem meu github" src="https://github.com/user-attachments/assets/0c7d1c8f-aeae-4574-a79e-f6e9a03c2cb6" />


# Installation

Follow the steps below to set up the environment and run the pipeline.

## 1. Clone the repository and activate the submodules

```bash
git clone https://github.com/P3dror0cha/Segunda_IC.git](https://github.com/P3dror0cha/PATFinder-project.git
cd PATFinder-project
git submodule update --init --recursive
```

## 2. Installing Nextflow in your ambient 

```bash
# If you are using conda
conda install -c conda-forge -c bioconda nextflow

# If you are using micromamba
micromamba install -c conda-forge -c bioconda nextflow

```
Note that nextflow requires Java version 17 or higher.

## 3. Running PATFinder

The workflow has two different pipelines. 

### Option 1 – Public metagenomes (`main.nf`):
This workflow automatically retrieves publicly available metagenomes from the MGnify API. The search criteria can be customized by modifying the URL defined in:
```text
bin/MGnify_download_functions.py
```

If you want to increase the number of results, increase the number of pages queried from the API. Use the max_pages param in nextflow.config

For usage of this mode, run the following commands:
```bash
# If you are using conda
nextflow run main.nf -resume -profile conda

# If you are using micromamba
nextflow run main.nf -resume -profile micromamba
```

### Option 2 – User-provided GBK and FAA files (`main_faa_and_gbk.nf`):
In this mode PATFinder takes your .gbk and .faa files as input. Before usage of the pipeline, the .gbk files must have a name similar to antiSMASH output and be present in the same folder. The name convention are: {sample_name}_{contig_number}.region{region_number}.gbk (ex. MGYG000296008_22.region001.gbk). The .faa files must have corresponding names with the .gbk files.

An example of accepted inputs are described below:
```text
PATFinder-project/
├── gbk_files/
│   ├── MGYG000296008_2.region001.gbk
│   ├── MGYG000296008_6.region001.gbk
│   ├── MGYG000296008_22.region001.gbk
│   ├── MGYG000296006_405.region001.gbk
│   ├── MGYG000296006_71.region001.gbk
│   ├── MGYG000296009_1.region001.gbk
│   └── MGYG000296014_17.region001.gbk
├── faa_files/
│   ├── MGYG000296008.faa
│   ├── MGYG000296006.faa
│   ├── MGYG000296009.faa
│   └── MGYG000296014.faa
```
For running this mode, use:
```bash
# If you are using conda
nextflow run main_faa_and_gbk.nf -resume -profile conda --gbk_files "path/to/your/gbk_files/*.gbk" --faa_files "path/to/your/faa_files/*.faa"

# If you are using micromamba
nextflow run main_faa_and_gbk.nf -resume -profile micromamba --gbk_files "path/to/your/gbk_files/*.gbk" --faa_files "path/to/your/faa_files/*.faa"
```
