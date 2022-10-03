#!/bin/bash


############################## GIT CONFIG ##############################

git config --global alias.root 'rev-parse --show-toplevel'


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

    ## SRA Toolkit
    wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.0.0/sratoolkit.3.0.0-centos_linux64.tar.gz -P $(git root)/pipelines/sratoolkit
    tar -xvf $(git root)/pipelines/sratoolkit/sratoolkit.3.0.0-centos_linux64.tar.gz -C $(git root)/pipelines/sratoolkit
    rm $(git root)/pipelines/sratoolkit/sratoolkit.3.0.0-centos_linux64.tar.gz

    ## Samtools
    wget https://github.com/samtools/samtools/releases/download/1.16/samtools-1.16.tar.bz2 -P $(git root)/pipelines/samtools
    tar -vxjf $(git root)/pipelines/samtools/samtools-1.16.tar.bz2 -C $(git root)/pipelines/samtools
    rm $(git root)/pipelines/samtools/samtools-1.16.tar.bz2

    ## STAR
    wget https://github.com/alexdobin/STAR/archive/2.7.10a.tar.gz -P $(git root)/pipelines/star
    tar -xzf $(git root)/pipelines/star/2.7.10a.tar.gz -C $(git root)/pipelines/star
    rm $(git root)/pipelines/star/2.7.10a.tar.gz

    ## RSEM
    wget https://github.com/deweylab/RSEM/archive/refs/tags/v1.3.3.tar.gz -P $(git root)/pipelines/rsem
    tar -xzf $(git root)/pipelines/rsem/v1.3.3.tar.gz -C $(git root)/pipelines/rsem
    rm $(git root)/pipelines/rsem/v1.3.3.tar.gz
    make -C $(git root)/pipelines/rsem/RSEM-1.3.3

else
    echo "Pipelines already installed."
fi



############################## CREATE PYTHON ENVIRONMENT AND INSTALL PACKAGES ##############################

python3 -m venv $(git root)/pipelines/python368_venv
source $(git root)/pipelines/python368_venv/bin/activate
pip install --upgrade pip
# python3 -m pip install -r $(git root)/pipelines/requirements.txt

