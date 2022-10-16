#!/bin/bash

torque=${1:-"1"}
path=${2:-"$(git root)"}
key=${3:-"2212-e9fcdc80cfd658683c786d89297a28fd"}






#################################################################################
######################### PREPARE WORKING ENVIRONNEMENT #########################
#################################################################################


export PATH="$path/scripts/bash/utils/:$PATH" 
export PATH="$path/scripts/bash/quantification/:$PATH"


if [[ $torque = 0 ]]
then
    submit=run-all-fastqs.sh
else
    module load torque
    submit=qsub-all-fastqs.sh
fi


assembly=human/assembly__GRCh38-hg38
annot=annotation__gencode/gencode_34

arrayGet() { 
    local array=$1 index=$2
    local i="${array}_$index"
    printf '%s' "${!i}"
}








##################################################################################
############################## DOWNLOAD IRIC FASTQS ##############################  
##################################################################################

project_IDs=(762 779 992 1090)


for p in ${project_IDs[@]}; do
    echo $p
    if [[ ! -d $path/data/iric/sc/dsp$p/downloaded-fastqs ]]; then
        wget --no-check-certificate -O - \
        "https://genomique.iric.ca/FastQList?key=${key}&projectID=$p&wget=1" \
        | wget --no-check-certificate -P $path/data/iric/sc/dsp$p/downloaded-fastqs -cri -
    fi
    merge-fastq-iric.sh $path/data/iric/sc/dsp$p/downloaded-fastqs $path/data/iric/sc/dsp$p/merged-fastqs
done








#########################################################################################################
############################## DOWNLOAD REFERENCE GENOME AND TRANSCRIPTOME ##############################
#########################################################################################################


if [[ ! -f $path/data/references/$assembly/GRCh38.primary_assembly.genome.fa.gz ]]; then 
    mkdir -p $path/data/references
    wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/GRCh38.primary_assembly.genome.fa.gz -P $path/data/references/$assembly
else
    echo "Reference $assembly already available."
fi


if [[ ! -f $path/data/references/$assembly/$annot/gencode.v34.transcripts.fa.gz ]]; then 
    mkdir -p $path/data/references
    wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/gencode.v34.transcripts.fa.gz -P $path/data/references/$assembly/$annot
else
    echo "Reference $annot already available."
fi


if [[ ! -f $path/data/references/exogenous/EGFP/genome.EGFP.fa.gz ]]; then 
    mkdir -p $path/data/references/exogenous/EGFP/
    echo ">EGFP
    ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAA
    GTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGC
    TGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAG
    CAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTA
    CAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGG
    ACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAAC
    GGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACAC
    CCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACG
    AGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAA" > $path/data/references/exogenous/EGFP/genome.EGFP.fa
    gzip $path/data/references/exogenous/EGFP/genome.EGFP.fa
else
    echo "Reference EGFP genome already available."
fi


if [[ ! -f $path/data/references/$assembly/$annot/gencode.v34.EGFP_transcripts.fa.gz ]]; then 
    cat $path/data/references/$assembly/$annot/gencode.v34.transcripts.fa.gz $path/data/references/exogenous/EGFP/genome.EGFP.fa.gz > $path/data/references/$assembly/$annot/gencode.v34.EGFP_transcripts.fa.gz
else
    echo "Reference $annot + EGFP already available."
fi









##################################################################################
######################## GET TRANSCRIPTS TO GENES MAPPING ########################
##################################################################################

if [[ ! -f $path/data/references/$assembly/$annot/gencode.v34.transcripts_txp2gene.tsv ]]; then
    zcat $path/data/references/$assembly/$annot/gencode.v34.transcripts.fa.gz | grep '>' | cut -d"|" -f1,2 > $path/data/references/$assembly/$annot/gencode.v34.transcripts_txp2gene.tsv
    sed -i 's/|/\t/g' $path/data/references/$assembly/$annot/gencode.v34.transcripts_txp2gene.tsv
    sed -i 's/>//g' $path/data/references/$assembly/$annot/gencode.v34.transcripts_txp2gene.tsv
fi

if [[ ! -f $path/data/references/$assembly/$annot/gencode.v34.EGFP_transcripts_txp2gene.tsv ]]; then
    echo -en 'EGFP\tEGFP' > $path/data/references/$assembly/$annot/tmp.tsv
    cat $path/data/references/$assembly/$annot/gencode.v34.transcripts_txp2gene.tsv $path/data/references/$assembly/$annot/tmp.tsv > $path/data/references/$assembly/$annot/gencode.v34.EGFP_transcripts_txp2gene.tsv
    rm $path/data/references/$assembly/$annot/tmp.tsv
fi







#############################################################################
############################## SALMON INDEXING ##############################    
#############################################################################


output_folder_index=$path/pipelines/salmon/salmon-1.4.0/files/$assembly/$annot
mkdir -p $output_folder_index/decoys
mkdir -p $output_folder_index/no-decoys



# Get decoys with gencode 34 #
###############################


if [[ ! -f $output_folder_index/decoys/decoys.txt && ! -f $output_folder_index/decoys/gentrome.fa ]] 
then
    salmon-get-decoys.sh -g $path/data/references/$assembly/GRCh38.primary_assembly.genome.fa.gz \
    -t $path/data/references/$assembly/$annot/gencode.v34.transcripts.fa.gz -o $output_folder_index/decoys
fi



# Indexing with decoys (several kmers sizes) #
##############################################


kmers=(17 19 23 27 29 31) 
for k in ${kmers[@]} 
do
    if [[ -d $output_folder_index/decoys/index_k$k ]] 
    then
        echo "decoys/index_k$k already generated.
        
        "
    else
        if [[ $torque = 0 ]]
        then
            export PATH="$path/pipelines/salmon/salmon-1.4.0/bin:$PATH"
            salmon index -k $k -t $output_folder_index/decoys/gentrome.fa.gz -d $output_folder_index/decoys/decoys.txt -i $output_folder_index/decoys/index_k$k --gencode

        else
            qsub-salmon-indexing.sh -w 5:00:00 -m 150gb -t $output_folder_index/decoys/gentrome.fa.gz -d $output_folder_index/decoys/decoys.txt -o $output_folder_index/decoys -k $k -p $path/pipelines/salmon/salmon-1.4.0 -args '--gencode'
            declare "pids_index_decoys_k$k=$(qstat -u $USER | tail -n 1 | awk '{print $1}' | cut -d"." -f1)"

        fi
    fi
done




# Indexing without decoys (k31) #
#################################


if [[ -d $output_folder_index/no-decoys/index_k31 ]] 
then
    echo "no-decoys/index_k31 already generated.
    
    "
else
    if [[ $torque = 0 ]]
    then
        export PATH="$path/pipelines/salmon/salmon-1.4.0/bin:$PATH"
        salmon index -k 31 -t $path/data/references/$assembly/$annot/gencode.v34.transcripts.fa.gz -i $output_folder_index/no-decoys/index_k31 --gencode
    else
        qsub-salmon-indexing.sh -w 1:00:00 -m 10gb -t $path/data/references/$assembly/$annot/gencode.v34.transcripts.fa.gz -o $output_folder_index/no-decoys -k 31 -p $path/pipelines/salmon/salmon-1.4.0 -args '--gencode'
        pids_index_no_decoys_k31=$(qstat -u $USER | tail -n 1 | awk '{print $1}' | cut -d"." -f1)
    fi
fi




# Get decoys with GFP seq & gencode 34 #
#########################################


if [[ ! -f $output_folder_index/EGFP-decoys/decoys.txt && ! -f $output_folder_index/EGFP-decoys/gentrome.fa ]]; then
    salmon-get-decoys.sh -g $path/data/references/$assembly/GRCh38.primary_assembly.genome.fa.gz \
    -t $path/data/references/$assembly/$annot/gencode.v34.EGFP_transcripts.fa.gz -o $output_folder_index/EGFP-decoys
fi


# Indexing with decoys + GFP seq & gencode 34 #
###############################################


if [[ -d $output_folder_index/EGFP-decoys/index_k19 ]]; then
    echo "EGFP-decoys/index_k19 already generated.
    
    "
else
    if [[ $torque = 0 ]]
    then
        export PATH="$path/pipelines/salmon/salmon-1.4.0/bin:$PATH"
        salmon index -k 19 -t $output_folder_index/EGFP-decoys/gentrome.fa.gz -d $output_folder_index/EGFP-decoys/decoys.txt -i $output_folder_index/EGFP-decoys/index_k19 --gencode
    else
        qsub-salmon-indexing.sh -w 15:00:00 -m 150gb -t $output_folder_index/EGFP-decoys/gentrome.fa.gz -d $output_folder_index/EGFP-decoys/decoys.txt -o $output_folder_index/EGFP-decoys -k 19 -p $path/pipelines/salmon/salmon-1.4.0 -args '--gencode'
        pids_index_decoys_gfp_k19=$(qstat -u $USER | tail -n 1 | awk '{print $1}' | cut -d"." -f1)
    fi
fi








######################################################################
############################## TRIMMING ############################## 
######################################################################



for p in ${project_IDs[@]}
do
    $submit -w 03:00:00 -m 20gb -r $path/scripts/bash/quality/run-fastqc.sh -f $path/data/iric/sc/dsp$p/merged-fastqs -o $path/data/iric/sc/dsp$p/merged-fastqs -l "PE" -s "all" -p $path/pipelines/fastqc/FastQC-0.11.9

    $submit -w 5:00:00 -m 10gb -r $path/scripts/bash/trimming/run-cutadapt.sh -f $path/data/iric/sc/dsp$p/merged-fastqs -o $path/data/iric/sc/dsp$p/trimmed-fastqs/cutadapt \
    -l "PE" -s "all" -p $path/pipelines/cutadapt/cutadapt-3.2 -savepids $path/tmp/qsub_pids \
    -args '-m 20:20 -a CTGTCTCTTATACACATCTC;min_overlap=6 -A A{100};min_overlap=6 -A N{20}GTACTCTGCGTTGATACCACTGCTTCCGCGGACAGGCGTGTAGATCTCGGTGGTCGCCGTATCATT;min_overlap=26'
    
    if [[ $torque = 1  && -s "$path/tmp/qsub_pids/pids.txt" ]]
    then
        sample=$(cut -f1 -d ' ' $path/tmp/qsub_pids/pids.txt)
        pid=$(cut -f2 -d ' ' $path/tmp/qsub_pids/pids.txt)
        declare "pids_trim_$p=$(echo $sample , $pid)"
    fi
    
    wait_pid="$(arrayGet pids_trim $p)"
    $submit -w 03:00:00 -m 20gb -h "$wait_pid" -r $path/scripts/bash/quality/run-fastqc.sh -f $path/data/iric/sc/dsp$p/trimmed-fastqs/cutadapt -o $path/data/iric/sc/dsp$p/trimmed-fastqs/cutadapt -l "PE" -s "all" -p $path/pipelines/fastqc/FastQC-0.11.9
done












#################################################################
##################### PARAMETERS COMPARISON #####################
#################################################################



# Decoys, trimming & default parameters + mt_genes + rrna genes #
#################################################################

for p in ${project_IDs[@]}
do

    if [[ $p = 1090 ]]
    then
        wait_pid="$(arrayGet pids_trim $p) ; $(arrayGet pids_index  decoys_gfp_k19)"
        mapping_params="EGFP-decoys-k19-1.4.0"
        index="EGFP-decoys"
        txp2gene="gencode.v34.EGFP_transcripts_txp2gene.tsv"

    else
        wait_pid="$(arrayGet pids_trim $p) ; $(arrayGet pids_index  decoys_k19)" 
        mapping_params="decoys-k19-1.4.0"
        index="decoys"
        txp2gene="gencode.v34.transcripts_txp2gene.tsv"
    fi


    args="-l ISR -i $output_folder_index/$index/index_k19 --tgMap $path/data/references/$assembly/$annot/$txp2gene --dropseq --dumpMtx --dumpFeatures --mrna $path/data/references/$assembly/$annot/biomart_ens100/mt_genes_biomart_ens100.tsv --rrna $path/data/references/$assembly/$annot/biomart_ens100/rRNA_genes_biomart_ens100.tsv"

    $submit -w 6:00:00 -m 60gb -h "$wait_pid" -r $path/scripts/bash/quantification/run-alevin.sh -f $path/data/iric/sc/dsp$p/trimmed-fastqs/cutadapt \
    -o $path/data/iric/sc/dsp$p/quant/alevin/$assembly/$annot/trimmed-reads-cutadapt/${mapping_params}/default -l "PE" -s "all" -p $path/pipelines/salmon/salmon-1.4.0 -args "$args"

done





# Decoys, trimming, forceCells 3000 & different kmer sizes #
############################################################


for i in ${!kmers[@]}
do
    k=${kmers[$i]}
    pid="$(arrayGet pids_trim 779) ; $(arrayGet pids_index decoys_k$k)"
    $submit -w 6:00:00  -h "$pid" -m 50gb -r $path/scripts/bash/quantification/run-alevin.sh -f $path/data/iric/sc/dsp779/trimmed-fastqs/cutadapt \
    -o $path/data/iric/sc/dsp779/quant/alevin/$assembly/$annot/trimmed-reads-cutadapt/decoys-k$k-1.4.0/forceCells-3000 -l "PE" -s "all" -p $path/pipelines/salmon/salmon-1.4.0 \
    -args "-l ISR -i $output_folder_index/decoys/index_k$k --tgMap $path/data/references/$assembly/$annot/gencode.v34.transcripts_txp2gene.tsv \
    --dropseq --dumpMtx --dumpFeatures --forceCells 3000"
done






# # Decoys, no trimming, forceCells 3000 & different kmer sizes #
# ###############################################################

# for i in ${!kmers[@]}
# do
#     k=${kmers[$i]}
#     pid="$(arrayGet pids_trim 779) ; $(arrayGet pids_index decoys_k$k)"
#     $submit -w 6:00:00  -h "$pid" -m 50gb -r $path/scripts/bash/quantification/run-alevin.sh -f $path/data/iric/sc/dsp779/merged-fastqs \
#     -o $path/data/iric/sc/dsp779/quant/alevin/$assembly/$annot/raw-reads/decoys-k$k-1.4.0/forceCells-3000 -l "PE" -s "all" -p $path/pipelines/salmon/salmon-1.4.0 \
#     -args "-l ISR -i $output_folder_index/decoys/index_k$k --tgMap $path/data/references/$assembly/$annot/gencode.v34.transcripts_txp2gene.tsv \
#     --dropseq --dumpMtx --dumpFeatures --forceCells 3000 --noWhitelist"
# done





# No decoys, k31 & forceCells 3000 #
####################################

pid="$(arrayGet pids_trim 779) ; $(arrayGet pids_index  no_decoys_k31)"
$submit -w 6:00:00 -m 50gb -h "$pid" -r $path/scripts/bash/quantification/run-alevin.sh -f $path/data/iric/sc/dsp779/trimmed-fastqs/cutadapt \
-o $path/data/iric/sc/dsp779/quant/alevin/$assembly/$annot/trimmed-reads-cutadapt/no-decoys-k31-1.4.0/forceCells-3000 -l "PE" -s "all" -p $path/pipelines/salmon/salmon-1.4.0 \
-args "-l ISR -i $output_folder_index/no-decoys/index_k31 --tgMap $path/data/references/$assembly/$annot/gencode.v34.transcripts_txp2gene.tsv \
--dropseq --dumpMtx --dumpFeatures --forceCells 3000"





# Decoys, k19 & forceCells 400 ER1 sample #
###########################################

$submit -w 6:00:00  -h "$pid" -m 50gb -r $path/scripts/bash/quantification/run-alevin.sh -f $path/data/iric/sc/dsp779/trimmed-fastqs/cutadapt \
    -o $path/data/iric/sc/dsp779/quant/alevin/$assembly/$annot/trimmed-reads-cutadapt/decoys-k19-1.4.0/forceCells-400 -l "PE" -s "Sample_N705_-_ER1" -p $path/pipelines/salmon/salmon-1.4.0 \
    -args "-l ISR -i $output_folder_index/decoys/index_k19 --tgMap $path/data/references/$assembly/$annot/gencode.v34.transcripts_txp2gene.tsv \
    --dropseq --dumpMtx --dumpFeatures --forceCells 400"








###############################################################
##################### SELECTED PARAMETERS #####################
###############################################################



for p in ${project_IDs[@]}
do
    if [[ $p = 1090 ]]
    then
        pid="$(arrayGet pids_trim $p) ; $(arrayGet pids_index  decoys_gfp_k19)"
        mapping_params="EGFP-decoys-k19-1.4.0"
        index="EGFP-decoys"
        txp2gene="gencode.v34.EGFP_transcripts_txp2gene.tsv"
    else
        pid="$(arrayGet pids_trim $p) ; $(arrayGet pids_index  decoys_k19)"
        mapping_params="decoys-k19-1.4.0"
        index="decoys"
        txp2gene="gencode.v34.transcripts_txp2gene.tsv"
    fi

    $submit -w 6:00:00 -m 50gb -h "$pid" -r $path/scripts/bash/quantification/run-alevin.sh -f $path/data/iric/sc/dsp$p/trimmed-fastqs/cutadapt \
    -o $path/data/iric/sc/dsp$p/quant/alevin/$assembly/$annot/trimmed-reads-cutadapt/${mapping_params}/forceCells-3000 -l "PE" -s "all" -p $path/pipelines/salmon/salmon-1.4.0 \
    -args "-l ISR -i $output_folder_index/$index/index_k19 --tgMap $path/data/references/$assembly/$annot/$txp2gene \
    --dropseq --dumpMtx --dumpFeatures --forceCells 3000"
done



