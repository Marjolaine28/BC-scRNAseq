#!/bin/bash

torque=${1:-"0"}
data_path=${2:-"$(git root)/BC-Cell-Lines_project/data"}
scripts_path=${3:-"$(git root)/scripts"}
references_path=${4:-"$(git root)/references"}
pipelines_path=${5:-"$(git root)/pipelines"}
key=${6:-"2212-e9fcdc80cfd658683c786d89297a28fd"}







#################################################################################
######################### PREPARE WORKING ENVIRONNEMENT #########################
#################################################################################


export PATH="$scripts_path/bash/utils/:$PATH" 
export PATH="$scripts_path/bash/quantification/:$PATH" 

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


# for p in ${project_IDs[@]}
# do
#     # wget --no-check-certificate -O - \
#     # "https://genomique.iric.ca/FastQList?key=${key}&projectID=$p&wget=1" \
#     # | wget --no-check-certificate -P $data_path/dsp$p/downloaded -cri -

#     merge-fastq-iric.sh $data_path/dsp$p/downloaded $data_path/dsp$p/raw-fastqs

# done









#########################################################################################################
############################## DOWNLOAD REFERENCE GENOME AND TRANSCRIPTOME ##############################
#########################################################################################################

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









##################################################################################
######################## GET TRANSCRIPTS TO GENES MAPPING ########################
##################################################################################

zcat $references_path/$assembly/$annot/gencode.v34.transcripts.fa.gz | grep '>' | cut -d"|" -f1,2 > $references_path/$assembly/$annot/gencode.v34.transcripts_txp2gene.tsv
sed -i 's/|/\t/g' $references_path/$assembly/$annot/gencode.v34.transcripts_txp2gene.tsv
sed -i 's/>//g' $references_path/$assembly/$annot/gencode.v34.transcripts_txp2gene.tsv

echo -en 'EGFP\tEGFP' > $references_path/$assembly/$annot/tmp.tsv
cat $references_path/$assembly/$annot/gencode.v34.transcripts_txp2gene.tsv $references_path/$assembly/$annot/tmp.tsv > $references_path/$assembly/$annot/gencode.v34.EGFP_transcripts_txp2gene.tsv
rm $references_path/$assembly/$annot/tmp.tsv









#############################################################################
############################## SALMON INDEXING ##############################    
#############################################################################


output_folder_index=$pipelines_path/salmon/salmon-1.4.0/files/$assembly/$annot 




################# Indexing with decoys (several kmers sizes)

mkdir -p $output_folder_index/decoys


if [[ ! -f $output_folder_index/decoys/decoys.txt && ! -f $output_folder_index/decoys/gentrome.fa ]] 
then
    salmon-get-decoys.sh -g $references_path/$assembly/GRCh38.primary_assembly.genome.fa.gz \
    -t $references_path/$assembly/$annot/gencode.v34.transcripts.fa.gz -o $output_folder_index/decoys
fi 


kmers=(17 19 23 25 27 29 31) 
for k in ${kmers[@]} 
do
    if [[ -d $output_folder_index/decoys/index_k$k ]] 
    then
        echo "decoys/index_k$k already generated.
        
        "
    else
        if [[ $torque = 0 ]]
        then
            export PATH="$pipelines_path/salmon/salmon-1.4.0/bin:$PATH"
            salmon index -k $k -t $output_folder_index/decoys/gentrome.fa.gz -d $output_folder_index/decoys/decoys.txt -i $output_folder_index/decoys/index_k$k --gencode

        else
            qsub-salmon-indexing.sh -w 5:00:00 -m 150gb -t $output_folder_index/decoys/gentrome.fa.gz -d $output_folder_index/decoys/decoys.txt -o $output_folder_index/decoys -k $k -p $pipelines_path/salmon/salmon-1.4.0 -args '--gencode'
            declare "pids_index_decoys_k$k=$(qstat -u $USER | tail -n 1 | awk '{print $1}' | cut -d"." -f1)"

        fi
    fi
done




################# Indexing without decoys

mkdir -p $output_folder_index/no-decoys

if [[ -d $output_folder_index/no-decoys/index_k19 ]] 
then
    echo "no-decoys/index_k19 already generated.
    
    "
else
    if [[ $torque = 0 ]]
    then
        export PATH="$pipelines_path/salmon/salmon-1.4.0/bin:$PATH"
        salmon index -k 19 -t $references_path/$assembly/$annot/gencode.v34.transcripts.fa.gz -i $output_folder_index/no-decoys/index_k19 --gencode
    else
        qsub-salmon-indexing.sh -w 1:00:00 -m 10gb -t $references_path/$assembly/$annot/gencode.v34.transcripts.fa.gz -o $output_folder_index/decoys -k 19 -p $pipelines_path/salmon/salmon-1.4.0 -args '--gencode'
        pids_index_no_decoys_k19=$(qstat -u $USER | tail -n 1 | awk '{print $1}' | cut -d"." -f1)
    fi
fi




################# Indexing using decoys + GFP seq

mkdir -p $output_folder_index/EGFP-decoys

#### Generate decoys

if [[ ! -f $output_folder_index/EGFP-decoys/decoys.txt && ! -f $output_folder_index/EGFP-decoys/gentrome.fa ]] 
then
    salmon-get-decoys.sh -g $references_path/$assembly/GRCh38.primary_assembly.genome.fa.gz \
    -t $references_path/$assembly/$annot/gencode.v34.EGFP_transcripts.fa.gz -o $output_folder_index/EGFP-decoys
 fi 


if [[ -d $output_folder_index/EGFP-decoys/index_k19 ]] 
then
    echo "EGFP-decoys/index_k19 already generated.
    
    "
else
    if [[ $torque = 0 ]]
    then
        export PATH="$pipelines_path/salmon/salmon-1.4.0/bin:$PATH"
        salmon index -k 19 -t $output_folder_index/EGFP-decoys/gentrome.fa.gz -d $output_folder_index/EGFP-decoys/decoys.txt -i $output_folder_index/EGFP-decoys/index_k19 --gencode
    else
        qsub-salmon-indexing.sh -w 15:00:00 -m 150gb -t $output_folder_index/EGFP-decoys/gentrome.fa.gz -d $output_folder_index/EGFP-decoys/decoys.txt -o $output_folder_index/EGFP-decoys -k 19 -p $pipelines_path/salmon/salmon-1.4.0 -args '--gencode'
        pids_index_decoys_gfp_k19=$(qstat -u $USER | tail -n 1 | awk '{print $1}' | cut -d"." -f1)
    fi
fi









######################################################################
############################## TRIMMING ############################## 
######################################################################



for p in ${project_IDs[@]}
do
    $submit -w 03:00:00 -m 20gb -r $scripts_path/bash/quality/run-fastqc.sh -f $data_path/dsp$p/raw-fastqs -o $data_path/dsp$p/raw-fastqs -l "PE" -s "all" -p $pipelines_path/fastqc/FastQC-0.11.9

    $submit -w 5:00:00 -m 10gb -r $scripts_path/bash/trimming/run-cutadapt.sh -f $data_path/dsp$p/raw-fastqs -o $data_path/dsp$p/trimmed-fastqs/cutadapt-all \
    -l "PE" -s "all" -p $pipelines_path/cutadapt/cutadapt-3.2 -savepids 1 \
    -args '-m 20:20 -a CTGTCTCTTATACACATCTC;min_overlap=6 -A A{100};min_overlap=6 -A N{20}GTACTCTGCGTTGATACCACTGCTTCCGCGGACAGGCGTGTAGATCTCGGTGGTCGCCGTATCATT;min_overlap=26'
    
    if [[ $torque = 1  && -s "~/tmp/qsub_pids/pids.txt" ]]
    then
        sample=$(cut -f1 -d ' ' ~/tmp/qsub_pids/pids.txt)
        pid=$(cut -f2 -d ' ' ~/tmp/qsub_pids/pids.txt)
        declare "pids_trim_$p=$(echo $sample , $pid)"
    fi
    
    wait_pid="$(arrayGet pids_trim $p)"
    $submit -w 03:00:00 -m 20gb -h "$wait_pid" -r $scripts_path/bash/quality/run-fastqc.sh -f $data_path/dsp$p/trimmed-fastqs/cutadapt-all -o $data_path/dsp$p/trimmed-fastqs/cutadapt-all -l "PE" -s "all" -p $pipelines_path/fastqc/FastQC-0.11.9
done








##########################################################################
############################## WHITELISTING ##############################
##########################################################################



for p in ${project_IDs[@]}
do
    wait_pid="$(arrayGet pids_trim $p)"
    $submit -w 00:30:00 -m 5gb -h "$wait_pid" -r $scripts_path/bash/whitelisting/get_cb_whitelist.sh -f $data_path/dsp$p/trimmed-fastqs/cutadapt-all \
    -o $data_path/dsp$p/trimmed-fastqs/cutadapt-all  -s "all" -p $scripts_path/python/rnaseq/whitelisting.py -savepids 1 -args "--top=3000"
    
    if [[ $torque = 1 && -s "~/tmp/qsub_pids/pids.txt" ]]
    then
        sample=$(cut -f1 -d ' ' ~/tmp/qsub_pids/pids.txt)
        pid=$(cut -f2 -d ' ' ~/tmp/qsub_pids/pids.txt)
        echo $pid
        declare "pids_wh_$p=$(echo $sample , $pid)"
    fi
done







############################################################################
############################## QUANTIFICATION ##############################
############################################################################



for p in ${project_IDs[@]}
do

    if [[ $p = 1090 ]]
    then
        wait_pid="$(arrayGet pids_wh $p) ; $(arrayGet pids_index  decoys_gfp_k19)"
        mapping_params="EGFP-decoys-k19-1.4.0"
        index="EGFP-decoys"
        txp2gene="biomart_ens100/txp2gene_biomart_ens100_EGFP.tsv"   #"gencode.v34.EGFP_transcripts_txp2gene.tsv"

    else
        wait_pid="$(arrayGet pids_wh $p) ; $(arrayGet pids_index  decoys_k19)"
        echo $wait_pid
        mapping_params="decoys-k19-1.4.0"
        index="decoys"
        txp2gene="biomart_ens100/txp2gene_biomart_ens100.tsv"
    fi

    args="-l ISR -i $output_folder_index/$index/index_k19 --tgMap $references_path/$assembly/$annot/$txp2gene --dropseq --dumpMtx --dumpFeatures --customWhitelist" # --writeMappings=mapping.sam

    $submit -w 6:00:00 -m 60gb -h "$wait_pid" -r $scripts_path/bash/quantification/run-alevin.sh -f $data_path/dsp$p/trimmed-fastqs/cutadapt-all \
    -o $data_path/dsp$p/quant/alevin/$assembly/$annot/trimmed-reads-cutadapt-all/${mapping_params}/customWh-top-3000 -l "PE" -s "all" -p $pipelines_path/salmon/salmon-1.4.0 -args "$args"

done












######################################################
################# Compare parameters #################
######################################################


# Decoys (recommanded), trimming & default parameters + mt_genes + rrna genes  #### EDIT : run on all samples

for p in ${project_IDs[@]}
do

    if [[ $p = 1090 ]]
    then
        wait_pid="$(arrayGet pids_wh $p) ; $(arrayGet pids_index  decoys_gfp_k19)"
        mapping_params="EGFP-decoys-k19-1.4.0"
        index="EGFP-decoys"
        txp2gene="biomart_ens100/txp2gene_biomart_ens100_EGFP.tsv"   #"gencode.v34.EGFP_transcripts_txp2gene.tsv"

    else
        wait_pid="$(arrayGet pids_wh $p) ; $(arrayGet pids_index  decoys_k19)" 
        mapping_params="decoys-k19-1.4.0"
        index="decoys"
        txp2gene="biomart_ens100/txp2gene_biomart_ens100.tsv"
    fi


    args="-l ISR -i $output_folder_index/$index/index_k19 --tgMap $references_path/$assembly/$annot/$txp2gene --dropseq --dumpMtx --dumpFeatures --mrna $references_path/$assembly/$annot/biomart_ens100/mt_genes_biomart_ens100.tsv --rrna $references_path/$assembly/$annot/biomart_ens100/rRNA_genes_biomart_ens100.tsv"

    $submit -w 6:00:00 -m 60gb -h "$wait_pid" -r $scripts_path/bash/quantification/run-alevin.sh -f $data_path/dsp$p/trimmed-fastqs/cutadapt-all \
    -o $data_path/dsp$p/quant/alevin/$assembly/$annot/trimmed-reads-cutadapt-all/${mapping_params}/default -l "PE" -s "all" -p $pipelines_path/salmon/salmon-1.4.0 -args "$args"

done





# Decoys, trimming, forceCells 3000 & different kmer sizes 


for i in ${!kmers[@]}
do
    k=${kmers[$i]}
    pid="$(arrayGet pids_trim 779) ; $(arrayGet pids_index decoys_k$k)"
    $submit -w 6:00:00  -h "$pid" -m 50gb -r $scripts_path/bash/quantification/run-alevin.sh -f $data_path/dsp779/trimmed-fastqs/cutadapt-all \
    -o $data_path/dsp779/quant/alevin/$assembly/$annot/trimmed-reads-cutadapt-all/decoys-k$k-1.4.0/forceCells-3000 -l "PE" -s "all" -p $pipelines_path/salmon/salmon-1.4.0 \
    -args "-l ISR -i $output_folder_index/decoys/index_k$k --tgMap $references_path/$assembly/$annot/biomart_ens100/txp2gene_biomart_ens100.tsv \
    --dropseq --dumpMtx --dumpFeatures --forceCells 3000 --mrna $references_path/$assembly/$annot/biomart_ens100/mt_genes_biomart_ens100.tsv --rrna $references_path/$assembly/$annot/biomart_ens100/rRNA_genes_biomart_ens100.tsv"
done



# Decoys, no trimming, forceCells 3000 & different kmer sizes 

for i in ${!kmers[@]}
do
    k=${kmers[$i]}
    pid="$(arrayGet pids_trim 779) ; $(arrayGet pids_index decoys_k$k)"
    $submit -w 6:00:00  -h "$pid" -m 50gb -r $scripts_path/bash/quantification/run-alevin.sh -f $data_path/dsp779/raw-fastqs \
    -o $data_path/dsp779/quant/alevin/$assembly/$annot/raw-reads/decoys-k$k-1.4.0/forceCells-3000 -l "PE" -s "all" -p $pipelines_path/salmon/salmon-1.4.0 \
    -args "-l ISR -i $output_folder_index/decoys/index_k$k --tgMap $references_path/$assembly/$annot/biomart_ens100/txp2gene_biomart_ens100.tsv \
    --dropseq --dumpMtx --dumpFeatures --forceCells 3000 --noWhitelist"
done



# No decoys, k19 & forceCells 3000

pid="$(arrayGet pids_trim 779) ; $(arrayGet pids_index  no_decoys_k19)"
$submit -w 6:00:00 -m 50gb -h "$pid" -r $scripts_path/bash/quantification/run-alevin.sh -f $data_path/dsp779/trimmed-fastqs/cutadapt-all \
-o $data_path/dsp779/quant/alevin/$assembly/$annot/trimmed-reads-cutadapt-all/no-decoys-k19-1.4.0/forceCells-3000 -l "PE" -s "all" -p $pipelines_path/salmon/salmon-1.4.0 \
-args "-l ISR -i $output_folder_index/no-decoys/index_k19 --tgMap $references_path/$assembly/$annot/biomart_ens100/txp2gene_biomart_ens100.tsv \
--dropseq --dumpMtx --dumpFeatures --forceCells 3000 --mrna $references_path/$assembly/$annot/biomart_ens100/mt_genes_biomart_ens100.tsv --rrna $references_path/$assembly/$annot/biomart_ens100/rRNA_genes_biomart_ens100.tsv"


# Decoys, k19 & forceCells 400 ER1 sample

$submit -w 6:00:00  -h "$pid" -m 50gb -r $scripts_path/bash/quantification/run-alevin.sh -f $data_path/dsp779/trimmed-fastqs/cutadapt-all \
    -o $data_path/dsp779/quant/alevin/$assembly/$annot/trimmed-reads-cutadapt-all/decoys-k19-1.4.0/forceCells-400 -l "PE" -s "Sample_N705_-_ER1" -p $pipelines_path/salmon/salmon-1.4.0 \
    -args "-l ISR -i $output_folder_index/decoys/index_k19 --tgMap $references_path/$assembly/$annot/biomart_ens100/txp2gene_biomart_ens100.tsv \
    --dropseq --dumpMtx --dumpFeatures --forceCells 400 --mrna $references_path/$assembly/$annot/biomart_ens100/mt_genes_biomart_ens100.tsv --rrna $references_path/$assembly/$annot/biomart_ens100/rRNA_genes_biomart_ens100.tsv"


# ################# Run Alevin on all other projects using selected parameters (kmer 19, forceCells 3000, decoys)             BEFORE WHITELISTING !
# ###############################################################################################################


# for p in ${project_IDs[@]}
# do
#     if [[ $p = 1090 ]]
#     then
#         # pid="n${pid_trim[$i]} ${pid_decoys_gfp}"    # Use also EGFP sequence to quantify dsp1090
#         pid="$(arrayGet pids_trim $p) ; $(arrayGet pids_index  k19-decoys-gfp)"
#         mapping_params="EGFP-decoys-k19-1.4.0"
#         index="EGFP-decoys"
#         txp2gene="biomart_ens100/txp2gene_biomart_ens100_EGFP.tsv"   #"gencode.v34.EGFP_transcripts_txp2gene.tsv"
#     else
#         # pid="n${pid_trim[$i]} ${pid_decoys[3]}"
#         pid="$(arrayGet pids_trim $p) ; $(arrayGet pids_index  k19-decoys)"
#         mapping_params="decoys-k19-1.4.0"
#         index="decoys"
#         txp2gene="biomart_ens100/txp2gene_biomart_ens100.tsv"   #"gencode.v34.transcripts_txp2gene.tsv"
#     fi

#     $submit -w 6:00:00 -m 50gb -h "$pid" -r $scripts_path/bash/quantification/run-alevin.sh -f $data_path/dsp$p/trimmed-fastqs/cutadapt-all \
#     -o $data_path/dsp$p/quant/alevin/$assembly/$annot/trimmed-reads-cutadapt-all/${mapping_params}/forceCells-3000 -l "PE" -s "all" -p $pipelines_path/salmon/salmon-1.4.0 \
#     -args "-l ISR -i $output_folder_index/$index/index_k19 --tgMap $txp2gene \
#     --dropseq --dumpMtx --dumpFeatures --forceCells 3000"
# done





# # Align T47D in bulk mode with STAR

# $submit -w 30:00:00 -m 100gb -r $scripts_path/bash/quantification/run-star.sh -f $data_path/dsp762/raw-fastqs \
# -o $data_path/dsp762/quant/star/$assembly/$annot/raw-reads -l "PEr2" -s "Sample_T47D" \
# -args "--genomeDir $pipelines_path/star/STAR-2.7.9a/files/human/assembly__GRCh38-hg38/annotation__gencode/gencode_34/index/ \
# --runThreadN 8 \
# --outSAMtype BAM SortedByCoordinate \
# --quantMode TranscriptomeSAM \
# --outSAMunmapped Within \
# --outSAMattributes Standard"