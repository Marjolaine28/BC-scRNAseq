#!/u/davidm/VirtualEns/python_3.6.8/bin/python3

from Bio.SeqIO.QualityIO import FastqGeneralIterator
import gzip
import pandas as pd
import csv
import sys

# sample = "Sample_N705_-_ER1"
# project_path = "/home/arion/davidm/Data/datasets/raw/private/RNA-seq/sc/dsp779/trimmed-fastqs/cutadapt"

project_path = sys.argv[1]
sample = sys.argv[2]
# cbs = set(line.strip() for line in open(f"{project_path}/trimmed-fastqs/cutadapt/{sample}/raw_cb_freq.txt"))
cbs = pd.read_table(f"{project_path}/trimmed-fastqs/cutadapt/{sample}/raw_cb_freq.tsv", header=None)[1].values
file = f"{project_path}/trimmed-fastqs/cutadapt/{sample}/trimmed-{sample}_R1.fastq.gz"
phreds = {cb : [] for cb in cbs}
umis = {cb : [] for cb in cbs}
i=0

with open(f"{project_path}/trimmed-fastqs/cutadapt/{sample}/new_phreds.csv", "w") as outfile:
    writer = csv.writer(outfile)
    writer.writerow(phreds.keys())
with open(f"{project_path}/trimmed-fastqs/cutadapt/{sample}/new_umis.csv", "w") as outfile:
    writer = csv.writer(outfile)
    writer.writerow(umis.keys())
    
    
for _, seq, qual in FastqGeneralIterator(gzip.open(file, "rt")) :

    cb = seq[:12]
    if cb in cbs :
        phreds[cb].append(qual)
        umis[cb].append(seq[12:])
    i +=1
    if i > 1000000 :
        with open(f"{project_path}/trimmed-fastqs/cutadapt/{sample}/new_phreds.csv", "a") as outfile:
            writer = csv.writer(outfile)
            writer.writerow(phreds.values())
        with open(f"{project_path}/trimmed-fastqs/cutadapt/{sample}/new_umis.csv", "a") as outfile:
            writer = csv.writer(outfile)
            writer.writerow(umis.values())

        phreds = {cb : [] for cb in cbs}
        umis = {cb : [] for cb in cbs}
        i = 0


with open(f"{project_path}/trimmed-fastqs/cutadapt/{sample}/new_phreds.csv", "a") as outfile:
    writer = csv.writer(outfile)
    writer.writerow(phreds.values())
with open(f"{project_path}/trimmed-fastqs/cutadapt/{sample}/new_umis.csv", "a") as outfile:
    writer = csv.writer(outfile)
    writer.writerow(umis.values())