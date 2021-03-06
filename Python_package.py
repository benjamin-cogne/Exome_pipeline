#!/bin/python

import subprocess
from subprocess import call
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import sys
import os
import MySQLdb
from datetime import datetime

if len(sys.argv) == 1:
    print("\n####################################################################################\n\n\tI am a package of python scripts to be used in an exome analysis pipeline\n\n\tUsage: python Python_package.py [Name Scripts] [Arguments]\n\n\tName Scripts : ExomePlot, Addmpilup, VariantsPrioritize, ReduceMultiallelicVCF, CheckCoverage, IGVSnapshots\n")    
    exit()

if sys.argv[1] == 'ReduceMultiallelicVCF':
    if len(sys.argv) != 5 :
        print('\nUsage: Python_package.py ReduceMultiallelicVCF [folder] [Name] [VCF]\n')
        exit() 
    Path = sys.argv[2]
    Name = sys.argv[3]
    VCF = open(sys.argv[4],'r')
    VCF_biall = open(Path+'variants_VQSR.'+Name+'.biallelic.vcf','w')

    it = 0
    print('Starting of python script: ReduceMultiallelicVCF')
    for line in VCF:
        it+=1
        if it % 10000 == 0:
            print(str(it)+' sites processed')
        if line[0] == '#':
            VCF_biall.write(line)
        else:
            field = line.split('\t')
            ### Check if multiallelic site
            if len(field[4].split(',')) == 1:
                VCF_biall.write(line)
            else:
                refpos = ''
                altpos = ''
                info = field[8].split(':')
                info_data = field[9].split(':')

            ### Recuperate genotype
                alleles = info_data[info.index('GT')].split('/')
                all1 = int(alleles[0])
                all2 = int(alleles[1])
                
            ### New AD and alt values corresponding to genotype
                AD_before = info_data[info.index('AD')].split(',')
                alt_before = field[4].split(',')
                if '<*:DEL>' in alt_before:
                    continue
                if all1 == 0 or (all1 > 0 and all1 == all2):
                    AD_after = AD_before[0]+','+AD_before[all2]
                    alt_after = alt_before[all2-1]

                else:
                    AD_after = AD_before[0]+','+AD_before[all1]+','+AD_before[all2]
                    alt_after = alt_before[all1-1]+','+alt_before[all2-1]
                
                if all1 == 0 and all2 > 1:
                    GT_after = '0/1'
                elif all1 > 1 and all2 > 1 and all1 == all2:
                    GT_after = '1/1'
                elif all1 > 1 and all2 > 1 and all1 != all2:
                    GT_after = '1/2'
                else:
                    GT_after = info_data[info.index('GT')]
            ### Retrieving constant sample info data fields
                if not 'DP' in info:
                    continue            
                DP = info_data[info.index('DP')]
                GQ = info_data[info.index('GQ')]
                PL =info_data[info.index('PL')]
                
            ### write new info data field
                if 'PID' in info:
                    PGT = info_data[info.index('PGT')]
                    PID = info_data[info.index('PID')]
                    info_data_new = GT_after+':'+AD_after+':'+DP+':'+GQ+':'+PGT+':'+PID+':'+PL
                else:
                    info_data_new = GT_after+':'+AD_after+':'+DP+':'+GQ+':'+PL
            
            ### write new bi-allelic line
                VCF_biall.write('\t'.join(field[:4])+'\t'+alt_after+'\t'+'\t'.join(field[5:9])+'\t'+info_data_new)

    VCF_biall.close()

    print('End of python script: ReduceMultiallelicVCF')    
    
if sys.argv[1] == 'ExomePlot':
    
    if len(sys.argv) != 3:
        print("\n\tExomePlot script need one argument: the path of the folder with hist.txt files generated by bedtools coverage ... | grep \'all\'\n\n\tUsage: python /.../Python_package.py ExomePlot pathToQualityControlFolder\n")
        exit()
    
    fig = plt.figure(1)
    ax = plt.subplot(111)

    for file in os.listdir(sys.argv[2]):
        if file.endswith("hist.txt"):
            f=open(sys.argv[2]+file,"r")
            x=[]
            y=[]
            cumul = 0
            for line in f:
                field=line.split("\t")
                if int(field[1]) < 201:
                    x.append(field[1])
                    cumul += float(field[4][:-2])
                    y.append(1-cumul)
                else:
                    break
            ax.plot(x,y,label=file[:-9])
            f.close()

    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.1,
                     box.width, box.height * 0.9])

    plt.title('Exome coverage')
    plt.xlabel('Sequencing Depth')
    plt.ylabel('Cumulative exome coverage')
    plt.xlim([0,200])
    plt.xticks([10,20,30,50,100,150,200])
    plt.yticks([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1])
    for i in [10,20,30]:
        plt.axvline(x=i,linewidth=0.5, color='0.2')
    plt.axhline(y=0.5, linewidth=0.5, color='0.2')
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1), ncol=5, fancybox=True)

    plt.savefig(sys.argv[2]+"Exome_cumulative_coverage.pdf")
    exit()
    
if sys.argv[1] == 'Addmpilup':
    
    ### Usage: python /.../Python_package.py Addmpilup Samtools_mpilup_file VCF_from_haplotypecaller 
    samtools = sys.argv[2]
    HC=sys.argv[3]

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

    exit()

if sys.argv[1] == 'VariantsPrioritize':

    if len(sys.argv) != 9 :
        print('\nUsage: Python_package.py VariantsPrioritize [Variants file] [OMIM file] [HPO file] [HPO file genes to phenotype][GDI file] [Name_analysis] [HPO list]\n')
        exit()

    Variants = sys.argv[2]
    OMIM = sys.argv[3]
    HPO_file_IN = sys.argv[4]
    HPO_file_genesTopheno = sys.argv[5]
    GDI_file_IN = sys.argv[6]
    Name = sys.argv[7]
    list_HP = sys.argv[8].split(',')
    
    # Create dictionary with genes related to phenotype
    file=open(HPO_file_IN,'r')
    dic_genes={}
    dic_desc={}
    for line in file:
        if line[0] != '#':
            field=line.split('\t')
            HP=field[0]
            desc=field[1]
            gene=field[3][:-1]
            if HP in list_HP:

                if not HP in dic_desc.keys():
                    dic_desc[HP]=desc
                if not HP in dic_genes.keys():
                    dic_genes[HP]=[gene]
                else:
                    dic_genes[HP].append(gene)

    file.close()
        
    # Create dictionary with OMIM gene to disease database [gene:disease]
    file=open(OMIM,'r')
    OMIM_lib={}
    for line in file:
        field=line.split('\t')
        OMIM_lib[field[0]]=field[1][:-1]
    file.close()

    # Create dictionary with the GDI phred of genes
    GDI_file=open(GDI_file_IN,'r')
    GDI_lib={}
    l=0
    for line in GDI_file:
        l+=1
        if l == 1:
            continue
        else:
            field=line.split('\t')
            GDI_lib[field[0]]=field[2]
    GDI_file.close()
    
    # Writing new file with new MAF filters and OMIM annotation
    file=open(Variants,'r')
    OUT=open(Variants+'.python_filtered.tab','w')
    Gene_before=''

    for line in file:
        if line[0] == '#':
            field=line.split('\t')
            OUT.write('\t'.join(field[:10])+'\t'+'GDI phred'+'\t'+'Biological process'+'\t'+'OMIM'+'\t'+'HPO'+'\t'+'Trans'+'\t'+'Clinic'+'\t'+'Internal DB'+'\t'+'\t'.join(field[10:]))
        else:
            field=line.split('\t')
            # Filter MAF<0.005 in EXAC and ESP
            if (field[10]=='.' or float(field[10])< 0.005) and (field[15]=='.' or float(field[15])< 0.005):

            # Add OMIM annotation and genes HP terms correspondind to phenotype query (HP terms)
                OMIM_disease = '.'
                gene=field[6]
                if len(gene.split(';')) > 1:
                    Gene=gene.split(';')[0]
                elif len(gene.split(',')) > 1:
                    Gene=gene.split(',')[0]
                else:
                    Gene=gene  
                    
                if Gene in OMIM_lib.keys():
                        OMIM_disease = OMIM_lib[Gene]
                
                #call(["wget", "http://www.genecards.org/cgi-bin/carddisp.pl?gene="+Gene,"-O" ])
                list_HP=[]
                for key in dic_genes.keys():
                    if Gene in dic_genes[key]:
                        list_HP.append(dic_desc[key])
                
            # Iteration through HPO gene to phenotype and creation of transmission and clinical symptoms lists
                list_clinic=[]
                transmission=[]
                dis_to_pheno=open(HPO_file_genesTopheno,'r')
                for i in dis_to_pheno:
                    if i[0] != '#':
                        field_HPO=i.split('\t')
                        if field_HPO[1] == Gene:
                            if field_HPO[3] in ['HP:0000006','HP:0000007','HP:0001419']:
                                if field_HPO[3] == 'HP:0000006':
                                    transmission.append('AD')
                                if field_HPO[3] == 'HP:0000007':
                                    transmission.append('AR')
                                if field_HPO[3] == 'HP:0001419':
                                    transmission.append('X')
                            else:
                                list_clinic.append(field_HPO[4][:-1])
                list_clinic = sorted(set(list_clinic))
                transmission = sorted(set(transmission))
            
            #Get Internal database frequency
                db = MySQLdb.connect(host="127.0.0.1",user="bcogne",passwd="*******",db='Exome_db')
                cur = db.cursor()
                command = "SELECT name FROM exome_id WHERE NOT name LIKE '%_E'"
                cur.execute(command)
                results = cur.fetchall()
                nb_variants = 0
                nb_exomes = 0
                for row in results:
                    nb_exomes+=1
                    name = row[0]
                    cur_exome = db.cursor()
                    command = "SELECT * FROM "+name+" WHERE chr='"+field[53]+"' AND pos='"+field[54]+"' AND ref='"+field[56]+"' AND alt='"+field[57]+"' LIMIT 1;"
                    cur_exome.execute(command)
                    if cur_exome.rowcount == 1:
                        nb_variants+=1
            
            # Get GO from Ensembl BioMart file giving gene to GO biological process relationship
                if Gene != Gene_before:
                    #call(["wget", "-O","/d2/Exome_analysis_BC/tmp/results.txt",'http://www.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "1" count = "" datasetConfigVersion = "0.6" ><Dataset name = "hsapiens_gene_ensembl" interface = "default" ><Filter name = "hgnc_symbol" value = "'+Gene+'"/><Filter name = "go_parent_name" value = "regulation of biological process"/><Attribute name = "name_1006" /><Attribute name = "hgnc_symbol" /></Dataset></Query>'])
                    GO_file=open('/d2/Exome_analysis_BC/Exome_pipeline_datas/GeneToGO_ensembl_biomart_11-2015.txt','r')
                    GO=[]
                    for line in GO_file:
                        field_GO=line.split('\t')
                        if field_GO[1] == Gene:
                            GO.append(field_GO[0])
                Gene_before=Gene
            
            # Get GDI phred:
                if Gene in GDI_lib.keys():
                    GDI=GDI_lib[Gene]
                else:
                    GDI='.'
                    
            #Add hypertext links to alignment igv snapshot and genecard and write variants lines

                start = int(field[1]) - 110
                stop = int(field[1]) + 100

                OUT.write('=LIEN_HYPERTEXTE("./snapshots_'+Name+'/'+Gene+'_'+str(start)+'-'+str(stop)+'.png";"'+field[0]+'")\t'+'\t'.join(field[1:6])+'\t'+'=LIEN_HYPERTEXTE("http://www.genecards.org/cgi-bin/carddisp.pl?gene='+Gene+'";"'+Gene+'")\t'+'\t'.join(field[7:10])+'\t'+GDI+'\t'+','.join(GO)+'\t'+OMIM_disease+'\t'+','.join(list_HP)+'\t'+','.join(transmission)+'\t'+','.join(list_clinic)+'\t'+str(nb_variants)+'/'+str(nb_exomes)+'\t'+'\t'.join(field[10:]))
                #OUT.write('=LIEN_HYPERTEXTE("./snapshots_'+Name+'/'+Gene+'_'+str(start)+'-'+str(stop)+'.png";"'+field[0]+'")\t'+'\t'.join(field[1:6])+'\t'+'=LIEN_HYPERTEXTE("http://www.genecards.org/cgi-bin/carddisp.pl?gene='+Gene+'";"'+Gene+'")\t'+'\t'.join(field[7:10])+'\t'+Gene_function+'\t'+OMIM_disease+'\t'+','.join(list_HP)+'\t'+','.join(transmission)+'\t'+','.join(list_clinic)+'\t'+str(nb_variants)+'/'+str(nb_exomes)+'\t'+'\t'.join(field[10:]))

                Gene_before=Gene

    db.commit()
    file.close()
    OUT.close()
    
    IN=open(Variants+'.python_filtered.tab','r')
    ### get the number of lines in file
    list_AR_genes={}
    nb_line=0
    for line in IN:
        if nb_line > 0:
            field=line.split('\t')
            gene=field[6]
            if len(gene.split(';')) > 1:
                Gene=gene.split(';')[0]
            elif len(gene.split(',')) > 1:
                Gene=gene.split(',')[0]
            else:
                Gene=gene
            if Gene in list_AR_genes.keys():
                list_AR_genes[Gene]=list_AR_genes[Gene]+1
            else:
                list_AR_genes[Gene]=1
        nb_line+=1
    IN.close()
    
    IN=open(Variants+'.python_filtered.tab','r')
    OUT_AR_X=open(Variants+'.python_filtered_AR-X.tab','w')
    OUT_AD=open(Variants+'.python_filtered_AD.tab','w')
    a=0
    var=0
    gene_before=''
    writed = False
    for line in IN:
        a+=1
        if a == 1:
            OUT_AR_X.write(line)
            OUT_AD.write(line)
        else:
            var+=1
            field=line.split('\t')
            gene=field[6]
            
            if len(gene.split(';')) > 1:
                Gene=gene.split(';')[0]
            elif len(gene.split(',')) > 1:
                Gene=gene.split(',')[0]
            else:
                Gene=gene
            if list_AR_genes[Gene] > 1 or field[57] == 'hom':
                OUT_AR_X.write(line)
            else:
                OUT_AD.write(line)
    
    IN.close()
    OUT_AD.close()
    OUT_AR_X.close()

if sys.argv[1] == 'IGVSnapshots':
    ### Create text file with genome coordinates of mutations to allow snapshots with igv
    ### Arguments : Python_package.py IGVSnapshots variants_IN.txt OUT.txt PATH_folder_images BAM_PATH
    IN = open(sys.argv[2], 'r')

    OUT = open(sys.argv[3],'w')
    path_img = sys.argv[4]
    OUT.write("new\ngenome hg19\nload "+sys.argv[5]+"\n")

    for line in IN:
        if line[0] != '#':
            field = line.split('\t')
            if (field[10]=='.' or float(field[10])< 0.005) and (field[15]=='.' or float(field[15])< 0.005) and field[59] == 'PASS':
                gene=field[6]
                if len(gene.split(';')) > 1:
                    Gene=gene.split(';')[0]
                elif len(gene.split(',')) > 1:
                    Gene=gene.split(',')[0]
                else:
                    Gene=gene
                pos=field[54]
                ref=field[56]
                alt=field[57]
                if len(ref) > 1 and len (alt) > 1:
                    if len(ref) > len(alt):
                        pos=str(int(pos)+(len(alt)-1))
                    elif len(ref) < len(alt):
                        pos=str(int(pos)+(len(ref)-1))
                start = int(pos) - 110
                stop = int(pos) + 100
                OUT.write("goto "+field[0]+":"+str(start)+"-"+str(stop)+"\n")
                OUT.write("snapshot "+path_img+Gene+"_"+str(start)+"-"+str(stop)+".png\n")

    OUT.write("exit")
    IN.close()
    OUT.close()
   
    call(["/opt/IGV_Snapshot/igv.sh", "--batch", sys.argv[3]])

if sys.argv[1] == 'CheckCoverage':

    ### Program to add each base coverage information in the exome analyzed on targeted genes:
    ### Argument 1: tab delimited file generated by program Gene_exon_coord.py
    ### Argument 2: bed file generated by GATK DepthOfCoverage on targeted genes
    
    if len(sys.argv) != 4 :
        print('\n\tUsage: Python_package.py CheckCoverage [tab delimited file generated by script Gene_exon_coord.py] [bed file generated by GATK DepthOfCoverage on targeted genes]\n')
        exit()

    # Write file with the coverage of each base in the canonical transcripts of genes covered
    cano_tab=open(sys.argv[2],'r')

    OUT=open(sys.argv[3]+'.targetedGenesCovered','w')

    for line in cano_tab:
        coverage_bed=open(sys.argv[3],'r')
        if line[0] != '#':
            field=line.split('\t')
            if field[3] == '.':
                OUT.write(line[:-1]+'\t.\n')
            else:
                list_coord=[]
                list_cov=[]
                for i in range(int(field[1])+1,int(field[2])):
                    list_coord.append(field[0]+':'+str(int(i)))
                for line2 in coverage_bed:
                    field2=line2.split('\t')
                    if field2[0] in list_coord:
                        list_cov.append(int(field2[1]))
                OUT.write(line[:-1]+'\t'+','.join(str(x) for x in list_cov)+'\t'+str(min(list_cov))+'\t'+str(sum(list_cov)/len(list_cov))+'\n')
        else:
            OUT.write(line[:-1]+'\tCoverage_per_base\tmin_cov\tmoy_cov\n')

    OUT.close()

    OUT=open(sys.argv[3]+'.targetedGenesCovered','r')

    # Write file with only low covered exons
    FILTERED=open(sys.argv[3]+'.targetedGenesCovered.lowCov.txt','w')

    for line in OUT:
        if line[0] != '#':
            field=line.split('\t')
            if field[3]=='.':
                FILTERED.write(line)
            elif int(field[10]) < 30:
                FILTERED.write(line)

    FILTERED.close()
    
    # Plotting of low covered exons
    
    path_fields = sys.argv[3].split('/')
    path = '/'+'/'.join(path_fields[:-1])+'/'
    pp = PdfPages(path+'Low_covered_genes_plots.pdf')
    FILTERED=open(sys.argv[3]+'.targetedGenesCovered.lowCov.txt','r')
    gene_before = ''
    gene_now = ''
    gene_exon_before = ''
    for line in FILTERED:
        field = line.split('\t')
        if not '.' in field[9]:
            gene_exon_now = field[7]+field[6]
            if gene_exon_now == gene_exon_before:
                continue
                
            gene_now = field[7]
             
            # if new gene: retrieving of codons start and sopt coordinates
            if not gene_before or gene_before != gene_now:
                refseq = open('/d2/Exome_analysis_BC/MMR_bed/MMR.refseq','r')
                for line2 in refseq:
                    field2 = line2.split('\t')
                    if field2[1] == field[5]:
                        cds_start = int(field2[6])
                        cds_stop = int(field2[7])
      
            fig = plt.figure(1)
            ax = plt.subplot(111)

            exon_start = int(field[1])
            exon_stop = int(field[2])
            exon_num = field[6]
            position = field[8]
            x = range(exon_start+1,exon_stop)
            y = [int(i) for i in field[9].split(',')]

            ax.plot(x,y, color='purple')
            
            print gene_now, exon_num
            #box = ax.get_position()
            #ax.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9]) 
                    
            plt.title(gene_now+' ('+position+') '+'exon '+exon_num+' coverage')
            #plt.xlabel('Exon genomic coordinates')
            plt.ylabel('Coverage')
            
            plt.xlim([exon_start-len(x)/20,exon_stop+len(x)/20])
            plt.ylim(0,max(y)+10)
            if position == '+':
                plt.xticks([exon_start+20,exon_stop-20],['start','end'])
            if position == '-':
                plt.xticks([exon_start+20,exon_stop-20],['end','start'])
            plt.yticks([10,20,30,max(y)])
            
            for i in [10,20,30]:
                plt.axhline(y=i,linewidth=0.5, color='0.1', linestyle='--')
                
            plt.axvline(x=exon_start, linewidth=0.5, color='g')
            plt.axvline(x=exon_stop, linewidth=0.5, color='g')
            plt.axvline(x=exon_start+20, linewidth=0.5, color='b', label='start')
            plt.axvline(x=exon_stop-20, linewidth=0.5, color='b', label='stop')
            
            if cds_start in x:
                plt.axvline(x=cds_start, linewidth=0.5, color='r')
            if cds_stop in x:
                plt.axvline(x=cds_stop, linewidth=0.5, color='r')
            #ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1), ncol=5, fancybox=True)

            plt.savefig(pp, format='pdf')
            
            gene_before = field[7]
            gene_exon_before = gene_before+exon_num
    pp.close()
    FILTERED.close()
    exit()
