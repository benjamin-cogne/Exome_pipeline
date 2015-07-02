#!/bin/python

import sys

####### Programme pour mise en forme des gènes ciblés pour l'exome depuis ucsc table browser #######
##### Nécessite 3 arguments: 1) les transcripts canoniques obtenus à partir du table browser de ucsc 
######                       2) L'intersection entre le bed obtenu sur ucsc table browser pour les gènes ciblés (avec tous les transcrits) et le design de capture choisi (ex Sureselect V4)
#####                        3) Le bed refgene contenant les coordonnées des gènes ciblés (avec tous les transcrits)

cano = sys.argv[1]
Intersect_SurV4 = sys.argv[2]
Refgene = sys.argv[3]

file=open(cano,"r")

### création dictionnaire avec le nom des transcrits canoniques: [NM_X:nom_gène] pour les gènes ciblés (ex: MLH1,PMS2...)
bib_cano={}
for line in file:
    if line[0] != '#':
        field=line.split('\t')
        bib_cano[field[4][:-1]]=field[3]
print bib_cano.keys()

### Ecrit un nouveau bed à partir du bed des refgenes (coordonnées exons) des gènes ciblés en extrayant uniquement les exons des gènes canoniques
file_refgene=open(Refgene,'r')
refgene_cano=open(Refgene+'.cano','w')
for line in file_refgene:
    if line[0] != '#':
        field=line.split('\t')
        field_NM=field[3].split('_')
        if field_NM[0]+'_'+field_NM[1] in bib_cano.keys():
            refgene_cano.write(line)
refgene_cano.close()

### création distionnaire à partir du refgene canonique [NM_X/numexon:line]
bib_refgene=open(Refgene+'.cano','r')
bib_refgene_cano={}
for line in bib_refgene:
    field=line.split('\t')
    field_NM=field[3].split('_')
    bib_refgene_cano[field_NM[0]+'_'+field_NM[1]+'/'+field_NM[3]]=line

### A partir du fichier d'intersection entre les refgenes et le design de capture: itération ligne par ligne et écriture du fichier de sortie sous un 
### format tab délimité contenant les coordonnées start et stop des transcrits canoniques et les start et stop capturés dans le design + le nom du gène et transcrit 
Intersect=open(Intersect_SurV4,'r')
OUT=open(Intersect_SurV4+".cano.bed","w")
OUT.write("#CHR\tStart_ref\tStop_ref\tStart_capt\tStop_capt\tNM\tExon\tGene_name\tFor/Rev\n")
done=[]
for line in Intersect:
    if line[0] != '#':
        field=line.split('\t')
        field_NM=field[3].split('_')
        print(field_NM[0]+'_'+field_NM[1]+'/'+field_NM[3])
        done.append(field_NM[0]+'_'+field_NM[1]+'/'+field_NM[3])
        if field_NM[0]+'_'+field_NM[1] in bib_cano.keys():
            start_ref=bib_refgene_cano[field_NM[0]+'_'+field_NM[1]+'/'+field_NM[3]].split('\t')[1]
            stop_ref=bib_refgene_cano[field_NM[0]+'_'+field_NM[1]+'/'+field_NM[3]].split('\t')[2]
            OUT.write(field[0]+'\t'+start_ref+'\t'+stop_ref+'\t'+field[1]+'\t'+field[2]+'\t'+field_NM[0]+'_'+field_NM[1]+'\t'+str(int(field_NM[3])+1)+'\t'+bib_cano[field_NM[0]+'_'+field_NM[1]]+'\t'+field[5][:-1]+'\n')

### Création d'une liste contenant les exons non couverts par le design de capture dans le refgene canonique
rest=list(set(bib_refgene_cano.keys())-set(done))

### Ecriture dans le fichier de sortie des exons non couverts en ajoutant des . au niveau des strat et stop du design de capture
for i in rest:
    field=bib_refgene_cano[i].split('\t')
    field_NM=field[3].split('_')
    OUT.write(field[0]+'\t'+field[1]+'\t'+field[2]+'\t.\t.\t'+field_NM[0]+'_'+field_NM[1]+'\t'+str(int(field_NM[3])+1)+'\t'+bib_cano[field_NM[0]+'_'+field_NM[1]]+'\t'+field[5][:-1]+'\n')

OUT.close()
