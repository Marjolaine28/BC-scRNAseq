#!/bin/bash

torque=${1:-"0"}
data_path=${2:-"$(git root)/data"}
scripts_path=${3:-"$(git root)/scripts"}
pipelines_path=${4:-"$(git root)/pipelines"}
key=${5:-"2212-e9fcdc80cfd658683c786d89297a28fd"}







#################################################################################
######################### PREPARE WORKING ENVIRONNEMENT #########################
#################################################################################


export PATH="$scripts_path/bash/utils/:$PATH" 
export PATH="$scripts_path/bash/quantification/:$PATH"
export PATH="$pipelines_path/sratoolkit/sratoolkit.3.0.0-centos_linux64/bin/:$PATH"
# ajouter path pour samtools et star

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
    wget --no-check-certificate -O - \
    "https://genomique.iric.ca/FastQList?key=${key}&projectID=$p&wget=1" \
    | wget --no-check-certificate -P $data_path/iric/dsp$p/downloaded -cri -

    merge-fastq-iric.sh $data_path/iric/dsp$p/downloaded $data_path/iric/dsp$p/raw-fastqs
done









#########################################################################################################
############################## DOWNLOAD REFERENCE GENOME AND TRANSCRIPTOME ##############################
#########################################################################################################

assembly=human/assembly__GRCh38-hg38
annot=annotation__gencode/gencode_34

if [[ ! -d $data_path/references ]]; 
then 
    mkdir -p $data_path/references

    wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/gencode.v34.transcripts.fa.gz -P $data_path/references/$assembly/$annot
    wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/GRCh38.primary_assembly.genome.fa.gz -P $data_path/references/$assembly
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
AGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAA" > $data_path/references/exogenous/EGFP/genome.EGFP.fa
gzip $data_path/references/exogenous/EGFP/genome.EGFP.fa
cat $data_path/references/$assembly/$annot/gencode.v34.transcripts.fa.gz $data_path/references/exogenous/EGFP/genome.EGFP.fa.gz > $data_path/references/$assembly/$annot/gencode.v34.EGFP_transcripts.fa.gz









##################################################################################
######################## GET TRANSCRIPTS TO GENES MAPPING ########################
##################################################################################

zcat $data_path/references/$assembly/$annot/gencode.v34.transcripts.fa.gz | grep '>' | cut -d"|" -f1,2 > $data_path/references/$assembly/$annot/gencode.v34.transcripts_txp2gene.tsv
sed -i 's/|/\t/g' $data_path/references/$assembly/$annot/gencode.v34.transcripts_txp2gene.tsv
sed -i 's/>//g' $data_path/references/$assembly/$annot/gencode.v34.transcripts_txp2gene.tsv

echo -en 'EGFP\tEGFP' > $data_path/references/$assembly/$annot/tmp.tsv
cat $data_path/references/$assembly/$annot/gencode.v34.transcripts_txp2gene.tsv $data_path/references/$assembly/$annot/tmp.tsv > $data_path/references/$assembly/$annot/gencode.v34.EGFP_transcripts_txp2gene.tsv
rm $data_path/references/$assembly/$annot/tmp.tsv









#############################################################################
############################## SALMON INDEXING ##############################    
#############################################################################


output_folder_index=$pipelines_path/salmon/salmon-1.4.0/files/$assembly/$annot 




################# Indexing with decoys (several kmers sizes)

mkdir -p $output_folder_index/decoys


if [[ ! -f $output_folder_index/decoys/decoys.txt && ! -f $output_folder_index/decoys/gentrome.fa ]] 
then
    salmon-get-decoys.sh -g $data_path/references/$assembly/GRCh38.primary_assembly.genome.fa.gz \
    -t $data_path/references/$assembly/$annot/gencode.v34.transcripts.fa.gz -o $output_folder_index/decoys
fi 


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
        salmon index -k 19 -t $data_path/references/$assembly/$annot/gencode.v34.transcripts.fa.gz -i $output_folder_index/no-decoys/index_k19 --gencode
    else
        qsub-salmon-indexing.sh -w 1:00:00 -m 10gb -t $data_path/references/$assembly/$annot/gencode.v34.transcripts.fa.gz -o $output_folder_index/decoys -k 19 -p $pipelines_path/salmon/salmon-1.4.0 -args '--gencode'
        pids_index_no_decoys_k19=$(qstat -u $USER | tail -n 1 | awk '{print $1}' | cut -d"." -f1)
    fi
fi




################# Indexing using decoys + GFP seq

mkdir -p $output_folder_index/EGFP-decoys

#### Generate decoys

if [[ ! -f $output_folder_index/EGFP-decoys/decoys.txt && ! -f $output_folder_index/EGFP-decoys/gentrome.fa ]] 
then
    salmon-get-decoys.sh -g $data_path/references/$assembly/GRCh38.primary_assembly.genome.fa.gz \
    -t $data_path/references/$assembly/$annot/gencode.v34.EGFP_transcripts.fa.gz -o $output_folder_index/EGFP-decoys
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
    $submit -w 03:00:00 -m 20gb -r $scripts_path/bash/quality/run-fastqc.sh -f $data_path/iric/dsp$p/raw-fastqs -o $data_path/iric/dsp$p/raw-fastqs -l "PE" -s "all" -p $pipelines_path/fastqc/FastQC-0.11.9

    $submit -w 5:00:00 -m 10gb -r $scripts_path/bash/trimming/run-cutadapt.sh -f $data_path/iric/dsp$p/raw-fastqs -o $data_path/iric/dsp$p/trimmed-fastqs/cutadapt-all \
    -l "PE" -s "all" -p $pipelines_path/cutadapt/cutadapt-3.2 -savepids 1 \
    -args '-m 20:20 -a CTGTCTCTTATACACATCTC;min_overlap=6 -A A{100};min_overlap=6 -A N{20}GTACTCTGCGTTGATACCACTGCTTCCGCGGACAGGCGTGTAGATCTCGGTGGTCGCCGTATCATT;min_overlap=26'
    
    if [[ $torque = 1  && -s "~/tmp/qsub_pids/pids.txt" ]]
    then
        sample=$(cut -f1 -d ' ' ~/tmp/qsub_pids/pids.txt)
        pid=$(cut -f2 -d ' ' ~/tmp/qsub_pids/pids.txt)
        declare "pids_trim_$p=$(echo $sample , $pid)"
    fi
    
    wait_pid="$(arrayGet pids_trim $p)"
    $submit -w 03:00:00 -m 20gb -h "$wait_pid" -r $scripts_path/bash/quality/run-fastqc.sh -f $data_path/iric/dsp$p/trimmed-fastqs/cutadapt-all -o $data_path/iric/dsp$p/trimmed-fastqs/cutadapt-all -l "PE" -s "all" -p $pipelines_path/fastqc/FastQC-0.11.9
done












#################################################################
##################### PARAMETERS COMPARISON #####################
#################################################################



# Decoys, trimming & default parameters + mt_genes + rrna genes  #### EDIT : run on all samples

for p in ${project_IDs[@]}
do

    if [[ $p = 1090 ]]
    then
        wait_pid="$(arrayGet pids_wh $p) ; $(arrayGet pids_index  decoys_gfp_k19)"
        mapping_params="EGFP-decoys-k19-1.4.0"
        index="EGFP-decoys"
        txp2gene="biomart_ens100/gencode.v34.EGFP_transcripts_txp2gene.tsv"

    else
        wait_pid="$(arrayGet pids_wh $p) ; $(arrayGet pids_index  decoys_k19)" 
        mapping_params="decoys-k19-1.4.0"
        index="decoys"
        txp2gene="biomart_ens100/gencode.v34.transcripts_txp2gene.tsv"
    fi


    args="-l ISR -i $output_folder_index/$index/index_k19 --tgMap $data_path/references/$assembly/$annot/$txp2gene --dropseq --dumpMtx --dumpFeatures --mrna $data_path/references/$assembly/$annot/biomart_ens100/mt_genes_biomart_ens100.tsv --rrna $data_path/references/$assembly/$annot/biomart_ens100/rRNA_genes_biomart_ens100.tsv"

    $submit -w 6:00:00 -m 60gb -h "$wait_pid" -r $scripts_path/bash/quantification/run-alevin.sh -f $data_path/iric/dsp$p/trimmed-fastqs/cutadapt-all \
    -o $data_path/iric/dsp$p/quant/alevin/$assembly/$annot/trimmed-reads-cutadapt-all/${mapping_params}/default -l "PE" -s "all" -p $pipelines_path/salmon/salmon-1.4.0 -args "$args"

done




# Decoys, trimming, forceCells 3000 & different kmer sizes 


for i in ${!kmers[@]}
do
    k=${kmers[$i]}
    pid="$(arrayGet pids_trim 779) ; $(arrayGet pids_index decoys_k$k)"
    $submit -w 6:00:00  -h "$pid" -m 50gb -r $scripts_path/bash/quantification/run-alevin.sh -f $data_path/iric/dsp779/trimmed-fastqs/cutadapt-all \
    -o $data_path/iric/dsp779/quant/alevin/$assembly/$annot/trimmed-reads-cutadapt-all/decoys-k$k-1.4.0/forceCells-3000 -l "PE" -s "all" -p $pipelines_path/salmon/salmon-1.4.0 \
    -args "-l ISR -i $output_folder_index/decoys/index_k$k --tgMap $data_path/references/$assembly/$annot/biomart_ens100/gencode.v34.transcripts_txp2gene.tsv \
    --dropseq --dumpMtx --dumpFeatures --forceCells 3000 --mrna $data_path/references/$assembly/$annot/biomart_ens100/mt_genes_biomart_ens100.tsv --rrna $data_path/references/$assembly/$annot/biomart_ens100/rRNA_genes_biomart_ens100.tsv"
done




# # Decoys, no trimming, forceCells 3000 & different kmer sizes 

# for i in ${!kmers[@]}
# do
#     k=${kmers[$i]}
#     pid="$(arrayGet pids_trim 779) ; $(arrayGet pids_index decoys_k$k)"
#     $submit -w 6:00:00  -h "$pid" -m 50gb -r $scripts_path/bash/quantification/run-alevin.sh -f $data_path/iric/dsp779/raw-fastqs \
#     -o $data_path/iric/dsp779/quant/alevin/$assembly/$annot/raw-reads/decoys-k$k-1.4.0/forceCells-3000 -l "PE" -s "all" -p $pipelines_path/salmon/salmon-1.4.0 \
#     -args "-l ISR -i $output_folder_index/decoys/index_k$k --tgMap $data_path/references/$assembly/$annot/biomart_ens100/gencode.v34.transcripts_txp2gene.tsv \
#     --dropseq --dumpMtx --dumpFeatures --forceCells 3000 --noWhitelist"
# done




# No decoys, k19 & forceCells 3000

pid="$(arrayGet pids_trim 779) ; $(arrayGet pids_index  no_decoys_k19)"
$submit -w 6:00:00 -m 50gb -h "$pid" -r $scripts_path/bash/quantification/run-alevin.sh -f $data_path/iric/dsp779/trimmed-fastqs/cutadapt-all \
-o $data_path/iric/dsp779/quant/alevin/$assembly/$annot/trimmed-reads-cutadapt-all/no-decoys-k19-1.4.0/forceCells-3000 -l "PE" -s "all" -p $pipelines_path/salmon/salmon-1.4.0 \
-args "-l ISR -i $output_folder_index/no-decoys/index_k19 --tgMap $data_path/references/$assembly/$annot/biomart_ens100/gencode.v34.transcripts_txp2gene.tsv \
--dropseq --dumpMtx --dumpFeatures --forceCells 3000 --mrna $data_path/references/$assembly/$annot/biomart_ens100/mt_genes_biomart_ens100.tsv --rrna $data_path/references/$assembly/$annot/biomart_ens100/rRNA_genes_biomart_ens100.tsv"




# Decoys, k19 & forceCells 400 ER1 sample

$submit -w 6:00:00  -h "$pid" -m 50gb -r $scripts_path/bash/quantification/run-alevin.sh -f $data_path/iric/dsp779/trimmed-fastqs/cutadapt-all \
    -o $data_path/iric/dsp779/quant/alevin/$assembly/$annot/trimmed-reads-cutadapt-all/decoys-k19-1.4.0/forceCells-400 -l "PE" -s "Sample_N705_-_ER1" -p $pipelines_path/salmon/salmon-1.4.0 \
    -args "-l ISR -i $output_folder_index/decoys/index_k19 --tgMap $data_path/references/$assembly/$annot/biomart_ens100/gencode.v34.transcripts_txp2gene.tsv \
    --dropseq --dumpMtx --dumpFeatures --forceCells 400 --mrna $data_path/references/$assembly/$annot/biomart_ens100/mt_genes_biomart_ens100.tsv --rrna $data_path/references/$assembly/$annot/biomart_ens100/rRNA_genes_biomart_ens100.tsv"





################# Run Alevin on all other projects using selected parameters (kmer 19, forceCells 3000, decoys)
###############################################################################################################


for p in ${project_IDs[@]}
do
    if [[ $p = 1090 ]]
    then
        # pid="n${pid_trim[$i]} ${pid_decoys_gfp}"    # Use also EGFP sequence to quantify dsp1090
        pid="$(arrayGet pids_trim $p) ; $(arrayGet pids_index  k19-decoys-gfp)"
        mapping_params="EGFP-decoys-k19-1.4.0"
        index="EGFP-decoys"
        txp2gene="biomart_ens100/gencode.v34.EGFP_transcripts_txp2gene.tsv"
    else
        # pid="n${pid_trim[$i]} ${pid_decoys[3]}"
        pid="$(arrayGet pids_trim $p) ; $(arrayGet pids_index  k19-decoys)"
        mapping_params="decoys-k19-1.4.0"
        index="decoys"
        txp2gene="biomart_ens100/gencode.v34.transcripts_txp2gene.tsv"
    fi

    $submit -w 6:00:00 -m 50gb -h "$pid" -r $scripts_path/bash/quantification/run-alevin.sh -f $data_path/iric/dsp$p/trimmed-fastqs/cutadapt-all \
    -o $data_path/iric/dsp$p/quant/alevin/$assembly/$annot/trimmed-reads-cutadapt-all/${mapping_params}/forceCells-3000 -l "PE" -s "all" -p $pipelines_path/salmon/salmon-1.4.0 \
    -args "-l ISR -i $output_folder_index/$index/index_k19 --tgMap $txp2gene \
    --dropseq --dumpMtx --dumpFeatures --forceCells 3000"
done










##############################################################################
############################## NEW WHITELISTING ##############################
##############################################################################



for p in ${project_IDs[@]}
do
    wait_pid="$(arrayGet pids_trim $p)"
    $submit -w 00:30:00 -m 5gb -h "$wait_pid" -r $scripts_path/bash/whitelisting/get_cb_whitelist.sh -f $data_path/iric/dsp$p/trimmed-fastqs/cutadapt-all \
    -o $data_path/iric/dsp$p/trimmed-fastqs/cutadapt-all  -s "all" -p $scripts_path/python/rnaseq/whitelisting.py -savepids 1 -args "--top=3000"
    
    if [[ $torque = 1 && -s "~/tmp/qsub_pids/pids.txt" ]]; then
        sample=$(cut -f1 -d ' ' ~/tmp/qsub_pids/pids.txt)
        pid=$(cut -f2 -d ' ' ~/tmp/qsub_pids/pids.txt)
        echo $pid
        declare "pids_wh_$p=$(echo $sample , $pid)"
    fi
done








#################################################################################################################################
############################## FINAL MAPPING & QUANT (with selected parameters and new whitelist) ###############################
#################################################################################################################################



for p in ${project_IDs[@]}
do

    if [[ $p = 1090 ]]; then
        wait_pid="$(arrayGet pids_wh $p) ; $(arrayGet pids_index  decoys_gfp_k19)"
        mapping_params="EGFP-decoys-k19-1.4.0"
        index="EGFP-decoys"
        txp2gene="biomart_ens100/gencode.v34.EGFP_transcripts_txp2gene.tsv"

    else
        wait_pid="$(arrayGet pids_wh $p) ; $(arrayGet pids_index  decoys_k19)"
        echo $wait_pid
        mapping_params="decoys-k19-1.4.0"
        index="decoys"
        txp2gene="biomart_ens100/gencode.v34.transcripts_txp2gene.tsv"
    fi

    args="-l ISR -i $output_folder_index/$index/index_k19 --tgMap $data_path/references/$assembly/$annot/$txp2gene --dropseq --dumpMtx --dumpFeatures --customWhitelist" # --writeMappings=mapping.sam

    $submit -w 6:00:00 -m 60gb -h "$wait_pid" -r $scripts_path/bash/quantification/run-alevin.sh -f $data_path/iric/dsp$p/trimmed-fastqs/cutadapt-all \
    -o $data_path/iric/dsp$p/quant/alevin/$assembly/$annot/trimmed-reads-cutadapt-all/${mapping_params}/customWh-top-3000 -l "PE" -s "all" -p $pipelines_path/salmon/salmon-1.4.0 -args "$args"

done














################################################################################
##############################  MAPPING BULK DATA ##############################
################################################################################


# Map 56 cells lines
####################

mkdir $data_path/public/bulk/56_cell_lines/downloaded

 
x=$(cat SRR_Acc_List.txt) ; 
for d in $x; do 
    prefetch $d ; 
    fasterq-dump --split-files $d; 
    rm -r $d; 
    gzip ${d}_*.fastq; 
done

sed -i -e 's/T47D Kbluc/T47D_Kbluc/g' SraRunTable.txt
sed -i -e 's/MDAMB175/MDAMB175VII/g' SraRunTable.txt
sed -i -e 's/SUM1315/SUM1315MO2/g' SraRunTable.txt

x=($(echo $(cut -f1 -d',' SraRunTable.txt)))
y=($(echo $(cut -f8 -d',' SraRunTable.txt)))

echo ${#x[@]}
echo ${#y[@]}

for i in "${!x[@]}"; do if [[ $i != 0 ]] ; then mkdir -p ${y[$i]} ; mv ${x[$i]}*1.fastq.gz ${y[$i]}/${x[$i]}*R1.fastq.gz ; mv ${x[$i]}*2.fastq.gz ${y[$i]}/${x[$i]}*R2.fastq.gz ; fi ; done
for i in "${!x[@]}"; do if [[ $i != 0 ]] ; then echo ${y[$i]} ; echo ${x[$i]} ; fi ; done


merge-fastq-iric.sh /home/arion/davidm/Data/datasets/public/RNA-seq/bulk/BC-Cell-Lines-Panel_GSE48213/raw-fastqs



# Map bulk MCF7 (DSP356) with Salmon
####################################

# No decoys k31 dsp356

# $submit -w 10:00:00 -m 60gb -r $scripts_path/bash/quantification/run-salmon.sh -f $data_path/dsp356/raw-fastqs \
# -o $data_path/dsp356/quant/salmon/$assembly/$annot/raw-reads/no-decoys-k31-1.4.0 -l "PE" -s "Sample_E-1  Sample_E-2  Sample_E-3 Sample_V-1  Sample_V-2  Sample_V-3" \
# -p $pipelines_path/salmon/salmon-1.4.0 -args "-l A -i $output_folder_index/no-decoys/index_k31 \
# -g $references_path/$assembly/$annot/biomart_ens100/txp2gene_biomart_ens100.tsv --validateMappings"


# Decoys k31 dsp356

# $submit -w 10:00:00 -m 60gb -r $scripts_path/bash/quantification/run-salmon.sh -f $data_path/dsp356/raw-fastqs \
# -o $data_path/dsp356/quant/salmon/$assembly/$annot/raw-reads/decoys-k31-1.4.0 -l "PE" -s "Sample_E-1  Sample_E-2  Sample_E-3 Sample_V-1  Sample_V-2  Sample_V-3" \
# -p $pipelines_path/salmon/salmon-1.4.0 -args "-l A -i $output_folder_index/decoys/index_k31 \
# -g $references_path/$assembly/$annot/biomart_ens100/txp2gene_biomart_ens100.tsv --validateMappings"







####################################################################################################
############################## STAR ALIGNMENTS FOR CORRELATIONS & IGV ##############################
####################################################################################################



# Create STAR index with GRCh38-hg38 assembly and gencode 37 annotation
#######################################################################

echo "module load star; \
STAR --runThreadN 8 \
--runMode genomeGenerate \
--genomeDir $sw/RNA-seq/star/star-2.7.1a/files/human/assembly__GRCh38-hg38/annotation__gencode/gencode_37/index \
--genomeFastaFiles $references/human/assembly__GRCh38-hg38/GRCh38.primary_assembly.genome.fa \
--sjdbGTFfile $references/human/assembly__GRCh38-hg38/annotation__gencode/gencode_37/gencode.v37.annotation.gtf" \
| qsub -V -l nodes=1,mem=200gb,vmem=200gb,walltime=48:00:00 -j oe -d $sw/RNA-seq/star/star-2.7.1a/files/human/assembly__GRCh38-hg38/annotation__gencode/gencode_37/index/logs/ -N "STAR-index-hg38-v37"



# Align sc T47D (DSP762) in bulk mode with STAR (just for IGV, no need to quantify)
###############################################

# change to gencode 37 !!
$submit -w 30:00:00 -m 100gb -r $scripts_path/bash/quantification/run-star.sh -f $data_path/iric/dsp762/raw-fastqs \
-o $data_path/iric/dsp762/quant/star/$assembly/$annot/raw-reads -l "PEr2" -s "Sample_T47D" \
-args "--genomeDir $pipelines_path/star/STAR-2.7.9a/files/human/assembly__GRCh38-hg38/annotation__gencode/gencode_34/index/ \
--runThreadN 8 \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM \
--outSAMunmapped Within \
--outSAMattributes Standard"



# Create RSEM reference with Gencode 37
#######################################

mkdir -p $sw/RNA-seq/rsem/rsem-1.2.28/files/human/assembly__GRCh38-hg38/annotation__gencode/gencode_37/index/

rsem-prepare-reference --gtf /home/arion/davidm/Data/references/human/assembly__GRCh38-hg38/annotation__gencode/gencode_37/gencode.v37.annotation.gtf \
/home/arion/davidm/Data/references/human/assembly__GRCh38-hg38/GRCh38.primary_assembly.genome.fa \
/home/arion/davidm/Softwares/RNA-seq/rsem/rsem-1.2.28/files/human/assembly__GRCh38-hg38/annotation__gencode/gencode_37/index/




# Align bulk T47D (DSP280) with Star
####################################

data_path=$dsp280
scripts_path=$scripts/RNA-seq
references_path=$references
pipelines_path=$sw/RNA-seq
assembly=human/assembly__GRCh38-hg38
annot=annotation__gencode/gencode_37

export PATH="$PATH:$scripts_path/bash/utils/" 
export PATH="$PATH:$scripts_path/bash/quantification/" 

qsub-all-fastqs.sh -w 30:00:00 -m 100gb -h 0 -r $scripts_path/bash/quantification/run-star.sh -f $data_path/raw-fastqs \
-o $data_path/quant/star/$assembly/$annot/raw -l "PE" -s "all" \
-args "--genomeDir $sw/RNA-seq/star/star-2.7.1a/files/human/assembly__GRCh38-hg38/annotation__gencode/gencode_37/index/ \
--runThreadN 8 \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM \
--outSAMunmapped Within \
--outSAMattributes Standard"


x=($(echo $(find $dsp280/quant/star/human/assembly__GRCh38-hg38/annotation__gencode/gencode_37/raw/ -mindepth 1 -maxdepth 1 -type d)))

for s in ${x[@]};  
do  
if [[ $(basename $s) != 'logs' ]]; 
then  
echo "rsem-calculate-expression -p 8 --paired-end --bam --no-bam-output \
$dsp280/quant/star/human/assembly__GRCh38-hg38/annotation__gencode/gencode_37/raw/ $(basename $s)/Aligned.toTranscriptome.out.bam \
$sw/RNA-seq/rsem/rsem-1.2.28/files/human/assembly__GRCh38-hg38/annotation__gencode/gencode_37/index/ \
$dsp280/quant/star/human/assembly__GRCh38-hg38/annotation__gencode/gencode_37/raw/ $(basename $s)/rsem" run-star.sh_| qsub -V -l nodes=1,mem=35gb,vmem=35gb,walltime=2:00:00 -j oe -d $dsp280/quant/star/human/assembly__GRCh38-hg38/annotation__gencode/gencode_37/raw/$(basename $s)/logs/ -N RSEM-bulk-T47D-$(basename $s);
  fi






###### VOIR AVEC PATRICK COMMENT TELECHARGER

# Get alignments already launched by the plateforme

# project_IDs=(550 589 1111)

# for p in ${project_IDs[@]}; do
#     wget --no-check-certificate -O - \
#     "https://genomique.iric.ca/BamList?key=${key}&projectID=$p&wget=1" \
#     | wget --no-check-certificate -P $data_path/iric/dsp$p/downloaded -cri -
# done


# mkdir -p $sw/RNA-seq/rsem/rsem-1.2.28/files/human/assembly__GRCh38-hg38/annotation__gencode/gencode_34/index/

# rsem-prepare-reference --gtf $references/human/assembly__GRCh38-hg38/annotation__gencode/gencode_34/gencode.v34.basic.annotation.gtf \
# /home/arion/davidm/Data/references/human/assembly__GRCh38-hg38/GRCh38.primary_assembly.genome.fa \
# $sw/RNA-seq/rsem/rsem-1.2.28/files/human/assembly__GRCh38-hg38/annotation__gencode/gencode_34/index/







########################################################################################
############################## INDEXING BAM FILES FOR IGV ##############################
########################################################################################


# sc T47D DSP762
################

samtools index $dsp762/quant/star/human/assembly__GRCh38-hg38/annotation__gencode/gencode_34/trimmed-reads-cutadapt/T47D/T47D_Aligned.sortedByCoord.out.bam

# bulk T47D DSP280 (NI-E2)
##########################

samtools index $dsp762/quant/star/human/assembly__GRCh38-hg38/annotation__gencode/gencode_34/trimmed-reads-cutadapt/T47D/T47D_Aligned.sortedByCoord.out.bam


# sc T47D DSP550 (pBabe n1)
###########################

samtools index $dsp762/quant/star/human/assembly__GRCh38-hg38/annotation__gencode/gencode_34/trimmed-reads-cutadapt/T47D/T47D_Aligned.sortedByCoord.out.bam

# sc ZR75 DSP550 (pBabe n1)
###########################

samtools index $dsp762/quant/star/human/assembly__GRCh38-hg38/annotation__gencode/gencode_34/trimmed-reads-cutadapt/T47D/T47D_Aligned.sortedByCoord.out.bam


# sc T47D DSP1111 (pMIG n1)
###########################

samtools index $dsp762/quant/star/human/assembly__GRCh38-hg38/annotation__gencode/gencode_34/trimmed-reads-cutadapt/T47D/T47D_Aligned.sortedByCoord.out.bam
