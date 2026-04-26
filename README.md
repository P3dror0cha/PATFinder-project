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
### 2. Installing dependencies

Create a conda environment using the following comands:

```bash
conda env create -f environment.yml
conda activate BGC_workflow
```
### 3. Installing important repositories

To ensure this workflow functions properly, you must first install BIG-SCAPE, KOfam, POEM-pipeline, and DeepSEA by following their respective installation instructions.
