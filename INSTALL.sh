#!/bin/bash

############################## INSTALL REQUIRED PIPELINES ##############################

if [[ ! -d $(git root)/pipelines ]]; 
then 
    mkdir -p $(git root)/pipelines

    ## FastQC 0.11.9
    wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip -P $(git root)/pipelines/fastqc
    unzip  $(git root)/pipelines/fastqc/fastqc_v0.11.9.zip -d $(git root)/pipelines/fastqc
    rm $(git root)/pipelines/fastqc/fastqc_v0.11.9.zip
    mv $(git root)/pipelines/fastqc/FastQC $(git root)/pipelines/fastqc/FastQC-0.11.9
    chmod 755 $(git root)/pipelines/fastqc/FastQC-0.11.9/fastqc

    ## Cutadapt 3.2
    mkdir -p $(git root)/pipelines/cutadapt
    python3 -m venv $(git root)/pipelines/cutadapt/cutadapt-3.2
    $(git root)/pipelines/cutadapt/cutadapt-3.2/bin/pip install --upgrade pip
    $(git root)/pipelines/cutadapt/cutadapt-3.2/bin/pip install cutadapt==3.2

    ## Salmon 1.4.0
    wget https://github.com/COMBINE-lab/salmon/releases/download/v1.4.0/salmon-1.4.0_linux_x86_64.tar.gz -P $(git root)/pipelines/salmon 
    tar -xvf $(git root)/pipelines/salmon/salmon-1.4.0_linux_x86_64.tar.gz -C $(git root)/pipelines/salmon
    rm $(git root)/pipelines/salmon/salmon-1.4.0_linux_x86_64.tar.gz
    mv $(git root)/pipelines/salmon/salmon-latest_linux_x86_64/ $(git root)/pipelines/salmon/salmon-1.4.0
else
    echo "Pipelines already installed."
fi



############################## DOWNLOAD REFERENCE GENOME AND TRANSCRIPTOME ##############################

ref=human/assembly__GRCh38-hg38/annotation__gencode/gencode_34

if [[ ! -d $(git root)/references ]]; 
then 
    mkdir -p $(git root)/references

    wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/gencode.v34.pc_transcripts.fa.gz -P $(git root)/references/$ref/transcriptome
    wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/GRCh38.primary_assembly.genome.fa.gz -P $(git root)/references/$ref/genome
else
    echo "References already available."
fi


############################## MAKE ALL SCRIPTS EXECUTABLE ##############################

chmod u+x $(find $(git root) -name '*.sh')
