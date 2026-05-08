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


## Installation

Follow the steps below to set up the environment and run the pipeline.

### 1. Clone the repository

```bash
git clone https://github.com/P3dror0cha/Segunda_IC.git
cd Segunda_IC
```
### 2. Installing dependencies and creating ambients

This section shows the commands used to make the necessary ambients for this pipeline to work. For instance, the ambients create are: 1 - antismash_bigscape, 2 - deepsea-project, 3 - DIAMOND_ambient and 4 - POEM_pipeline_ambient.

## 2.1 antiSMASH and BIG-SCAPE
Create a conda environment for antismash_bigscape.yml file in envs folder:

```bash
conda env create -f envs/antismash_bigscape.yml
```
For BIG-SCAPE to work, you will need the Pfam database. In the repository folder, run the following commands:

```bash
mkdir databases
cd databases
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
gunzip Pfam-A.hmm.gz
hmmpress Pfam-A.hmm
cd ..
```
## 2.2 DeepSEA
In the repository folder (PATFinder_project), use the following commands:

```bash
git clone https://github.com/computational-chemical-biology/DeepSEA-project.git
cd DeepSEA-project
conda env create -f environment-gpu.yml 
```
## 2.3 DIAMOND and KOfam
In the repository folder, use the following commands to make the env and install KOfam database:

```bash
conda env create -f envs/diamond_kofam.yml

cd databases
wget ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz
wget ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz

gunzip ko_list.gz
tar -xzf profiles.tar.gz
```
## 2.4 POEM-pipeline
In the repository folder, use the following commands to make the env

```bash
conda activate diamond_kofam
git clone https://github.com/Rinoahu/POEM_py3k
cd ./POEM_py3k
bash ./install.sh
```
Obs: Note that if you do not use conda, changes in install.sh are necessary!

### 3. Installing important repositories

To ensure this workflow functions properly, you must first install BIG-SCAPE, KOfam, POEM-pipeline, and DeepSEA by following their respective installation instructions.
