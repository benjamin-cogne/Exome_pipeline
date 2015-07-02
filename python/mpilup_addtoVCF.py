#!/bin/python

### Usage: python mpilup_addtoVCF.py Samtools_mpilup_file VCF_from_haplotypecaller 

import sys
samtools = sys.argv[1]
HC=sys.argv[2]

file=open(samtools,"r")

bib_sam={}
for line in file:
    if line[0] != '#':
        field=line.split('\t')	
        bib_sam[field[0]+field[1]+field[3]+field[4]]=field[5]

OUT=open(HC+".QUALSAM","w")

file2=open(HC,"r")

for line in file2:
    if line[0] == '#':
        OUT.write(line)
    else:
        field=line.split('\t')
        if field[0]+field[1]+field[3]+field[4] in bib_sam.keys():
            OUT.write(field[0]+'\t'+field[1]+'\t'+field[2]+'\t'+field[3]+'\t'+field[4]+'\t'+field[5]+'\t'+field[6]+'\t'+field[7]+';QUALSAM='+bib_sam[field[0]+field[1]+field[3]+field[4]]+'\t'+field[8]+'\t'+field[9])
        else:
            OUT.write(line)

OUT.close()

OUT=open(HC+".QUALSAM","r")
OUT2=open(HC+".QUALSAMHEAD","w")
cnt=0
for line in OUT:
    if cnt==0:
        if line.startswith('##contig'):
            OUT2.write("##INFO=<ID=QUALSAM,Number=1,Type=Float,Description=\"Quality variant calling by samtools\">\n")
            cnt=1
            OUT2.write(line)
        else:
            OUT2.write(line)
    else:
        OUT2.write(line)
OUT2.close()
