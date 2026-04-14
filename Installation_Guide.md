This file has the purpose to create an ambient where you can run the nextflow pipeline.
# Creating a conda/mamba ambient
blabla
# Antismash installation
Following the installation guide (https://docs.antismash.secondarymetabolites.org/install/) make the following comands inside your ambient:
```bash
conda install -c conda-forge -c bioconda -c defaults antismash=8.0.4
download-antismash-databases (verificar se é necessário dps)
```
# BIG-SCAPE installation
For BIG-SCAPE to work we must first install Pfam database
```bash
wget https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
gunzip Pfam-A.hmm.gz
hmmpress Pfam-A.hmm
```
After installing the database, run the command below and install bigscape. Remember to change the path to Pfam database in BIGSCAPE.nf
```bash
conda install -c conda-forge -c bioconda bigscape=2.0.0
```
# DeepSEA installation
Following the installation guide described in (https://github.com/computational-chemical-biology/DeepSEA-project), make the following steps inside the PATFinder-project repository.
```bash
git clone https://github.com/computational-chemical-biology/DeepSEA-project.git
cd DeepSEA-project
conda env update --file environment-gpu.yml 
```
