#!/u/davidm/VirtualEns/python_3.6.8/bin/python3

from Bio.SeqIO.QualityIO import FastqGeneralIterator
import csv
import sys

# sample = "Sample_N705_-_ER1"
# project_path = "/home/arion/davidm/Data/datasets/raw/private/RNA-seq/sc/dsp779/raw-fastqs"

project_path = sys.argv[1]
sample = sys.argv[2]
cbs = set(line.strip() for line in open(f"{project_path}/raw-fastqs/{sample}/CBs.txt"))
file = f"{project_path}/raw-fastqs/{sample}/{sample}_R1.fastq"
phreds = {cb : [] for cb in cbs}
umis = {cb : [] for cb in cbs}
i=0

with open(f"{project_path}/raw-fastqs/{sample}/new_phreds.csv", "w") as outfile:
    writer = csv.writer(outfile)
    writer.writerow(phreds.keys())
with open(f"{project_path}/raw-fastqs/{sample}/new_umis.csv", "w") as outfile:
    writer = csv.writer(outfile)
    writer.writerow(umis.keys())
    
    
for _, seq, qual in FastqGeneralIterator(open(file)) :

    cb = seq[:12]
    if cb in cbs :
        phreds[cb].append(qual)
        umis[cb].append(seq[12:])
    i +=1
    if i > 1000000 :
        with open(f"{project_path}/raw-fastqs/{sample}/new_phreds.csv", "a") as outfile:
            writer = csv.writer(outfile)
            writer.writerow(phreds.values())
        with open(f"{project_path}/raw-fastqs/{sample}/new_umis.csv", "a") as outfile:
            writer = csv.writer(outfile)
            writer.writerow(umis.values())

        phreds = {cb : [] for cb in cbs}
        umis = {cb : [] for cb in cbs}
        i = 0


with open(f"{project_path}/raw-fastqs/{sample}/new_phreds.csv", "a") as outfile:
    writer = csv.writer(outfile)
    writer.writerow(phreds.values())
with open(f"{project_path}/raw-fastqs/{sample}/new_umis.csv", "a") as outfile:
    writer = csv.writer(outfile)
    writer.writerow(umis.values())