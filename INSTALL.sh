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
else
    echo "Pipelines already installed."
fi



############################## DOWNLOAD REFERENCE GENOME AND TRANSCRIPTOME ##############################

assembly=human/assembly__GRCh38-hg38
annot=annotation__gencode/gencode_34

if [[ ! -d $(git root)/references ]]; 
then 
    mkdir -p $(git root)/references

    wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/gencode.v34.transcripts.fa.gz -P $(git root)/references/$assembly/$annot
    wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/GRCh38.primary_assembly.genome.fa.gz -P $(git root)/references/$assembly
else
    echo "References already available."
fi

echo ">EGFP
ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAA
GTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGC
TGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAG
CAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTA
CAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGG
ACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAAC
GGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACAC
CCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACG
AGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAA" > $references_path/exogenous/EGFP/genome.EGFP.fa
gzip $references_path/exogenous/EGFP/genome.EGFP.fa
cat $references_path/$assembly/$annot/gencode.v34.transcripts.fa.gz $references_path/exogenous/EGFP/genome.EGFP.fa.gz > $references_path/$assembly/$annot/gencode.v34.EGFP_transcripts.fa.gz



######################## GET TRANSCRIPTS TO GENES MAPPING ########################


zcat $references_path/$assembly/$annot/gencode.v34.transcripts.fa.gz | grep '>' | cut -d"|" -f1,2 > $references_path/$assembly/$annot/gencode.v34.transcripts_txp2gene.tsv
sed -i 's/|/\t/g' $references_path/$assembly/$annot/gencode.v34.transcripts_txp2gene.tsv
sed -i 's/>//g' $references_path/$assembly/$annot/gencode.v34.transcripts_txp2gene.tsv

echo -en 'EGFP\tEGFP' > $references_path/$assembly/$annot/tmp.tsv
cat $references_path/$assembly/$annot/gencode.v34.transcripts_txp2gene.tsv $references_path/$assembly/$annot/tmp.tsv > $references_path/$assembly/$annot/gencode.v34.EGFP_transcripts_txp2gene.tsv
rm $references_path/$assembly/$annot/tmp.tsv