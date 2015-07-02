#!/bin/bash

fastq1=$1
fastq2=$2
Name=$3
hg19="/d2/Exome_analysis_BC/ucsc.hg19.fasta"
Chemin="/d2/Exome_analysis_BC/Resultats_"$3

mkdir $Chemin || { echo 'mkdir failed' ; exit 1; }

######### Quality control of fastq files #############

#/opt/FastQC/fastqc $1 $2

############### Trimming des adaptateurs et sur la qualité de séquençage en 3' ###

python /opt/cutadapt-1.8.1/bin/cutadapt -q 30 -m 50 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -o $Chemin/trimmed.$3.1.fastq -p $Chemin/trimmed.$3.2.fastq $1 $2 || { echo 'cutadapt failed' ; exit 1; }

############## Alignement des pared-reads par BWA mem ###

mkdir $Chemin/temp || { echo 'mkdir temp failed' ; exit 1; }

#### Pour compatibilité avec GATK: nécessité d'utiliser l'option -R '@RG\tID:group1\tSM:sample1\tPL:illumina\tLB:lib1\tPU:unit1'

/opt/bwa mem -M -R "@RG\tID:exome\tSM:${Name}\tPL:illumina\tLB:lib1" -t 10 $hg19 $Chemin/trimmed.$3.1.fastq $Chemin/trimmed.$3.2.fastq > $Chemin/temp/$3.sam || { echo 'bwa mem failed' ; exit 1; }

############## Sam -> Bam avec Samtools - Marquage des duplicats de PCR avec Picard ##########

/opt/samtools view -hbS -@ 10 $Chemin/temp/$3.sam > $Chemin/temp/$3.bam
/opt/samtools sort -@ 10 $Chemin/temp/$3.bam $Chemin/temp/sorted_$3

java -jar /opt/picard-tools-1.119/MarkDuplicates.jar \
    INPUT=$Chemin/temp/sorted_$3.bam \
    OUTPUT=$Chemin/temp/sorted_dedup_$3.bam \
    METRICS_FILE=$Chemin/metrics_markduplicates.txt || { echo 'MarkDuplicates failed' ; exit 1; }

/opt/samtools index $Chemin/temp/sorted_dedup_$3.bam

############## Realignement des reads autour des indels avec GATK ###

java -Xmx4g -jar /opt/GenomeAnalysisTK.jar \
	-T RealignerTargetCreator \
	-L /d2/Exome_analysis_BC/Sureselect_V4_Covered.bed \
	-R $hg19 \
	-nt 10 \
	-known /d2/Exome_analysis_BC/1000G_phase1.indels.hg19.sites.vcf \
	-o $Chemin/temp/$3.bam.list \
	-I $Chemin/temp/sorted_dedup_$3.bam || { echo 'RealignerTargetCreator failed' ; exit 1; }

java -jar /opt/GenomeAnalysisTK.jar \
    -T IndelRealigner \
    -R $hg19 \
    -I $Chemin/temp/sorted_dedup_$3.bam \
	-known /d2/Exome_analysis_BC/1000G_phase1.indels.hg19.sites.vcf \
    -targetIntervals $Chemin/temp/$3.bam.list \
    -o $Chemin/temp/realigned_$3.bam || { echo 'IndelRealigner failed' ; exit 1; }
	
############### Base quality score recalibration with GATK #######

## Analyze patterns of covariation in the sequence dataset ##

java -jar /opt/GenomeAnalysisTK.jar \
    -T BaseRecalibrator \
	-L /d2/Exome_analysis_BC/Sureselect_V4_Covered.bed \
    -R $hg19 \
    -nct 8 \
    -I $Chemin/temp/realigned_$3.bam \
    -knownSites /d2/Exome_analysis_BC/dbsnp_138.hg19.vcf \
    -o $Chemin/temp/recal_data.table || { echo 'BaseRecalibrator failed' ; exit 1; }

## Do a second pass to analyze covariation remaining after recalibration ## FACULTATIF -> Verification sur plot recalibration

#java -jar /opt/GenomeAnalysisTK.jar \
 #   -T BaseRecalibrator \
  #  -R $hg19 \
   # -I $Chemin/realigned_$3.bam \
    #-knownSites /d2/hg19/00-All.vcf \
    #-BQSR $Chemin/Base_recalibration/recal_data.table \
    #-o $Chemin/Base_recalibration/post_recal_data.table

## Generate before/after plots ## FACULTATIF -> Verification sur plot recalibration

#java -jar /opt/GenomeAnalysisTK.jar \
 #   -T AnalyzeCovariates \
  #  -R $hg19 \
   # -before $Chemin/Base_recalibration/recal_data.table \
   # -after $Chemin/Base_recalibration/post_recal_data.table \
   # -plots $Chemin/Base_recalibration/recalibration_plots.pdf
	
## Apply recalibration to data ##

java -jar /opt/GenomeAnalysisTK.jar \
    -T PrintReads \
    -nct 8 \
    -R $hg19 \
    -I $Chemin/temp/realigned_$3.bam \
    -BQSR $Chemin/temp/recal_data.table \
    -o $Chemin/realigned_recal_$3.bam || { echo 'GATK PrintReads failed' ; exit 1; }
	
################## Check Bam quality ####################################

mkdir $Chemin/Quality_control

/opt/samtools flagstat $Chemin/realigned_recal_$3.bam > $Chemin/Stat_mapping_$3.txt

/opt/bedtools/bedtools coverage -abam $Chemin/realigned_recal_$3.bam -b /d2/Exome_analysis_BC/Sureselect_V4_Covered.bed -hist | grep ^all > $Chemin/Quality_control/realigned_recal_coverage.hist.txt

python /home/bcogne/Multiple_exome_plot.py $Chemin/Quality_control/ $Name

java -jar /opt/GenomeAnalysisTK.jar \
   -T DepthOfCoverage \
   -R $hg19 \
   -o $Chemin/Quality_control/coverage_$3 \
   -I $Chemin/realigned_recal_$3.bam \
   -geneList /d2/Exome_analysis_BC/refseq_gatk.txt \
   -ct 10 -ct 20 -ct 30 \
   -L /d2/Exome_analysis_BC/MMR_20bp_intersect_sureselectv4.bed

################ Variant Calling ########################################

mkdir $Chemin/Variants

#java -jar /opt/GenomeAnalysisTK.jar \
 #   -T HaplotypeCaller \
  #  -R $hg19 \
   # -I $Chemin/realigned_recal_$3.bam \
    #-L /d2/Exome_analysis_BC/Sureselect_V4_Covered.bed \
    #--genotyping_mode DISCOVERY \
    #-stand_emit_conf 10 \
    #-stand_call_conf 30 \
    #-o $Chemin/Variants/temp/variants_HC_$3.vcf || { echo 'GATK HaplotypeCaller failed' ; exit 1; }

#-A StrandAlleleCountsBySample
	
java -jar /opt/GenomeAnalysisTK.jar \
        -T HaplotypeCaller \
		-ERC GVCF \
		--variant_index_type LINEAR \
		--variant_index_parameter 128000 \
        -R $hg19 \
        -nct 8 \
        -A StrandAlleleCountsBySample \
        -I $Chemin/realigned_recal_$3.bam \
	    -L /d2/Exome_analysis_BC/Sureselect_V4_Covered.bed \
        --genotyping_mode DISCOVERY \
        -stand_emit_conf 10 \
        -stand_call_conf 30 \
        -o $Chemin/Variants/variants_HC_$3.g.vcf

java -jar /opt/GenomeAnalysisTK.jar \
        -T GenotypeGVCFs \
        -nt 10 \
        -L /d2/Exome_analysis_BC/Sureselect_V4_Covered.bed \
        -R $hg19 \
        -V $Chemin/Variants/variants_HC_$3.g.vcf \
        -o $Chemin/Variants/variants_HC_$3.vcf

/opt/samtools mpileup -Bug -f $hg19 -l /d2/Exome_analysis_BC/Sureselect_V4_Covered.bed $Chemin/realigned_recal_$3.bam | /opt/bcftools view -cvg - > $Chemin/Variants/$3.samtools 

################  Variant Filtration ######################################

java -jar /opt/GenomeAnalysisTK.jar \
    -T SelectVariants \
    -R $hg19 \
    -V $Chemin/Variants/variants_HC_$3.vcf \
    -L /d2/Exome_analysis_BC/Sureselect_V4_Covered.bed \
    -selectType SNP \
    -o $Chemin/Variants/raw_snps_$3.vcf

java -jar /opt/GenomeAnalysisTK.jar \
    -T VariantFiltration \
    -R $hg19 \
    -V $Chemin/Variants/raw_snps_$3.vcf \
    --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
    --filterName "snp_filter" \
    -o $Chemin/Variants/filtered_snps_$3.vcf

java -jar /opt/GenomeAnalysisTK.jar \
    -T SelectVariants \
    -R $hg19 \
    -V $Chemin/Variants/variants_HC_$3.vcf \
    -selectType INDEL \
    -o $Chemin/Variants/raw_indels_$3.vcf

java -jar /opt/GenomeAnalysisTK.jar \
    -T VariantFiltration \
    -R $hg19 \
    -V $Chemin/Variants/raw_indels_$3.vcf \
    --filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" \
    --filterName "indel_filter" \
    -o $Chemin/Variants/filtered_indels_$3.vcf


/opt/vcftools_0.1.11/bin/vcf-concat $Chemin/Variants/filtered_snps_$3.vcf $Chemin/Variants/filtered_indels_$3.vcf > $Chemin/Variants/vcf_concat_variants_$3.vcf  

cat $Chemin/Variants/vcf_concat_variants_$3.vcf | /opt/vcftools_0.1.11/bin/vcf-sort > $Chemin/Variants/HC_vcf_sorted_$3.vcf

python /home/bcogne/mpilup_add.py $Chemin/Variants/$3.samtools $Chemin/Variants/HC_vcf_sorted_$3.vcf

mv $Chemin/Variants/HC_vcf_sorted_$3.vcf.QUALSAMHEAD $Chemin/Variants/Final_variants_$3.vcf

rm -r $Chemin/temp
rm $Chemin/Variants/raw_snps_$3.vcf* $Chemin/Variants/filtered_snps_$3.vcf* $Chemin/Variants/raw_indels_$3.vcf* $Chemin/Variants/filtered_indels_$3.vcf* $Chemin/Variants/vcf_concat_variants_$3.vcf $Chemin/Variants/HC_vcf_sorted_$3.vcf.QUALSAM $Chemin/Variants/HC_vcf_sorted_$3.vcf

############################# Filtration sur MAF et annotation ###########################

/opt/annovar/convert2annovar.pl -format vcf4 -includeinfo -withfreq $Chemin/Variants/Final_variants_$3.vcf -outfile $Chemin/Variants/Final_variants_$3.annovar

/opt/annovar/annotate_variation.pl -filter -dbtype 1000g2014oct_all -buildver hg19 -out $Chemin/Variants/Final_variants_${3}_1000g_0.01 $Chemin/Variants/Final_variants.annovar /opt/annovar/humandb/ -score_threshold 0.01

/opt/annovar/table_annovar.pl $Chemin/Variants/Final_variants_$3*filtered /opt/annovar/humandb/ -buildver hg19 -out $Chemin/Variants/Variants_annotated_$3 -remove -protocol refGene,esp6500siv2_all,1000g2014oct_all,1000g2014oct_afr,1000g2014oct_eas,1000g2014oct_eur,exac03,snp138,ljb26_all -operation g,f,f,f,f,f,f,f,f -nastring . -otherinfo -arg '-splicing 10',,'-score_threshold 0.01',,,,,,

head -n 1 $Chemin/Variants/Variants_annotated_$3.hg19_multianno.txt > $Chemin/Variants/Final_variants_${3}_exon_splice_filtered.txt

awk '{if (($6~/^exonic/ || $6~/^splicing/) && !($9~/^synonymous/)) {print $0}}' $Chemin/Variants/Variants_annotated_$3.hg19_multianno.txt | grep 'PASS' >> $Chemin/Variants/Final_variants_${3}_exon_splice_filtered.txt
