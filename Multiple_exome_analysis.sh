#!/bin/bash


#### Benjamin COGNE 06/2015
##### Programme pour analyser en parallèle plusieurs exomes #######
### Nécessite un dossier avec tous les fastq.gz R1 et R2 des exomes ##
#### Nom des analyses sera le premier élément séparé par un '_' ###



function FastqToBAM {
    
	mkdir $4/Resultats_$3
	Path="$4/Resultats_$3"
	hg19="/d2/Exome_analysis_BC/ucsc.hg19.fasta"	
	
	######### Quality control of fastq files #############

   #/opt/FastQC/fastqc $1 $2

    ############### Trimming des adaptateurs et sur la qualité de séquençage en 3' ###

    python /opt/cutadapt-1.8.1/bin/cutadapt -q 30 -m 50 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -o $Path/trimmed.$3.R1.fastq -p $Path/trimmed.$3.R2.fastq $1 $2 || { echo 'cutadapt failed' ; exit 1; }

    ############## Alignement des pared-reads par BWA mem ###

    mkdir $Path/temp || { echo 'mkdir temp failed' ; exit 1; }

    # Pour compatibilité avec GATK: nécessité d'utiliser l'option -R '@RG\tID:group1\tSM:sample1\tPL:illumina\tLB:lib1\tPU:unit1'

    /opt/bwa mem -M -R "@RG\tID:exome\tSM:${3}\tPL:illumina\tLB:lib1" $hg19 $Path/trimmed.$3.R1.fastq $Path/trimmed.$3.R2.fastq > $Path/temp/$3.sam || { echo 'bwa mem failed' ; exit 1; }

    ############## Sam -> Bam avec Samtools - Marquage des duplicats de PCR avec Picard ##########

    /opt/samtools view -hbS $Path/temp/$3.sam > $Path/temp/$3.bam
    /opt/samtools sort $Path/temp/$3.bam $Path/temp/sorted_$3

    java -jar /opt/picard-tools-1.119/MarkDuplicates.jar \
        INPUT=$Path/temp/sorted_$3.bam \
        OUTPUT=$Path/temp/sorted_dedup_$3.bam \
        METRICS_FILE=$Path/metrics_markduplicates.txt || { echo 'MarkDuplicates failed' ; exit 1; }

    /opt/samtools index $Path/temp/sorted_dedup_$3.bam

    ############## Realignement des reads autour des indels avec GATK ###

    java -Xmx4g -jar /opt/GenomeAnalysisTK.jar \
	    -T RealignerTargetCreator \
	    -L /d2/Exome_analysis_BC/Sureselect_V4_Covered.bed \
	    -R $hg19 \
	    -known /d2/Exome_analysis_BC/1000G_phase1.indels.hg19.sites.vcf \
	    -o $Path/temp/$3.bam.list \
	    -I $Path/temp/sorted_dedup_$3.bam || { echo 'RealignerTargetCreator failed' ; exit 1; }

    java -jar /opt/GenomeAnalysisTK.jar \
        -T IndelRealigner \
        -R $hg19 \
        -I $Path/temp/sorted_dedup_$3.bam \
	    -known /d2/Exome_analysis_BC/1000G_phase1.indels.hg19.sites.vcf \
        -targetIntervals $Path/temp/$3.bam.list \
        -o $Path/temp/realigned_$3.bam || { echo 'IndelRealigner failed' ; exit 1; }
	
    ############### Base quality score recalibration with GATK #######

    ## Analyze patterns of covariation in the sequence dataset ##

    java -jar /opt/GenomeAnalysisTK.jar \
        -T BaseRecalibrator \
	    -L /d2/Exome_analysis_BC/Sureselect_V4_Covered.bed \
        -R $hg19 \
        -I $Path/temp/realigned_$3.bam \
        -knownSites /d2/Exome_analysis_BC/dbsnp_138.hg19.vcf \
        -o $Path/temp/recal_data.table || { echo 'BaseRecalibrator failed' ; exit 1; }

    ## Apply recalibration to data ##

    java -jar /opt/GenomeAnalysisTK.jar \
        -T PrintReads \
        -R $hg19 \
        -I $Path/temp/realigned_$3.bam \
        -BQSR $Path/temp/recal_data.table \
        -o $Path/realigned_recal_$3.bam || { echo 'GATK PrintReads failed' ; exit 1; }
	
    ################ Check Bam quality ####################################

    mkdir $Path/Quality_control

    /opt/samtools flagstat $Path/realigned_recal_$3.bam > $Path/Quality_control/Stat_mapping_$3.txt

    /opt/bedtools/bedtools coverage -abam $Path/realigned_recal_$3.bam -b /d2/Exome_analysis_BC/Sureselect_V4_Covered.bed -hist | grep ^all > $Path/Quality_control/$3.hist.txt

}

function gVCFToVCFfiltered {

    echo "Début de gVCFToVCFfiltered"

    Chemin=$1
    
	hg19="/d2/Exome_analysis_BC/ucsc.hg19.fasta"
	list_gvcf=($Chemin/*/*.g.vcf)
	
	touch $Chemin/gvcf.list
	
	for i in ${list_gvcf[@]}; 
	do
	    echo $i >> $Chemin/gvcf.list
	done
    
	##### Merge variants in one vcf file ##########
	
    java -jar /opt/GenomeAnalysisTK.jar \
        -T GenotypeGVCFs \
        -nt 10 \
        -L /d2/Exome_analysis_BC/MMR_20bp_intersect_sureselectv4.bed \
        -R $hg19 \
        -V $Chemin/gvcf.list \
        -o $Chemin/variants_all.vcf
		
	###### Variant Quality Score Recalibration: add -an InbreedingCoeff if samples > 10 ########
    
	java -jar /opt/GenomeAnalysisTK.jar \
        -T VariantRecalibrator \
        -R $hg19 \
        -input $Chemin/variants_all.vcf \
        -resource:hapmap,known=false,training=true,truth=true,prior=15.0 /d2/Exome_analysis_BC/hapmap_3.3.hg19.sites.vcf \
        -resource:omni,known=false,training=true,truth=true,prior=12.0 /d2/Exome_analysis_BC/1000G_omni2.5.hg19.sites.vcf \
        -resource:1000G,known=false,training=true,truth=false,prior=10.0 /d2/Exome_analysis_BC/1000G_phase1.snps.high_confidence.hg19.sites.vcf \
        -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /d2/Exome_analysis_BC/dbsnp_138.hg19.vcf \
        -an DP -an QD -an FS -an SOR -an MQ -an MQRankSum -an ReadPosRankSum \
        -mode SNP \
        -nt 10 \
        -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
        -recalFile $Chemin/recalibrate_SNP.recal \
        -tranchesFile $Chemin/recalibrate_SNP.tranches \
        -rscriptFile $Chemin/recalibrate_SNP_plots.R

    java -jar /opt/GenomeAnalysisTK.jar \
        -T ApplyRecalibration \
        -R $hg19 \
        -input $Chemin/variants_all.vcf \
        -mode SNP \
        --ts_filter_level 99.0 \
        -recalFile $Chemin/recalibrate_SNP.recal \
        -tranchesFile $Chemin/recalibrate_SNP.tranches \
        -o $Chemin/recalibrated_snps_raw_indels.vcf 
	
	java -jar /opt/GenomeAnalysisTK.jar \
        -T VariantRecalibrator \
        -R $hg19 \
        -input $Chemin/recalibrated_snps_raw_indels.vcf \
        -resource:mills,known=true,training=true,truth=true,prior=12.0  /d2/Exome_analysis_BC/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf \
        -an QD -an DP -an FS -an SOR -an MQRankSum -an ReadPosRankSum \
        -mode INDEL \
        -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
        --maxGaussians 4 \
        -recalFile $Chemin/recalibrate_INDEL.recal \
        -tranchesFile $Chemin/recalibrate_INDEL.tranches \
        -rscriptFile $Chemin/recalibrate_INDEL_plots.R
		
	java -jar /opt/GenomeAnalysisTK.jar \
        -T ApplyRecalibration \
        -R $hg19 \
        -input $Chemin/recalibrated_snps_raw_indels.vcf \
        -mode INDEL \
        --ts_filter_level 99.0 \
        -recalFile $Chemin/recalibrate_INDEL.recal \
        -tranchesFile $Chemin/recalibrate_INDEL.tranches \
        -o $Chemin/recalibrated_variants_all.vcf;
}

hg19="/d2/Exome_analysis_BC/ucsc.hg19.fasta"
Chemin="/d2/Exome_analysis_BC/Analyse_$(date +%Y_%m_%d_%H-%M)"

mkdir -p $Chemin


####### Création de la liste des patients #########

list_fastq=($1*.fastq.gz)
files=()

for i in ${list_fastq[@]};
do
    a=(${i//'/'/ })
    files+=(${a[@]:(-1)})
done

Name=()

for i in ${files[@]};
do
    a=(${i//_/ })
	Name+=($a)
done

list_Name="$(echo ${Name[@]} | tr " " "\n" | sort | uniq )"

####### Itération sur tous les patients et génération des BAM ############

#for i in ${list_Name[@]};
#do 
#    FastqToBAM "$1*$i*R1*fastq.gz" "$1*$i*R2*fastq.gz" "$i" "$Chemin"
#done

#echo -e ${list_Name[@]} | parallel -j 4 FastqToBAM "$1*{}*R1*fastq.gz" "$1*{}*R2*fastq.gz" "$i" "$Chemin"

echo -e "Itération sur tous les fastq.gz du dossier spécifié \n\n\t\t--> FastqToBAM en parallèle sur 15 CPUs"
export -f FastqToBAM
parallel -j 15 --eta --bibtex FastqToBAM "$1*{}*R1*fastq.gz" "$1*{}*R2*fastq.gz" "{}" "$Chemin" ::: ${list_Name[@]}

rm -R $Chemin/*/temp

################ Variant Calling ########################################

mkdir $Chemin/VCF

echo -e "Début variant calling avec HaplotypeCaller en parallèle sur 15 CPUs"

find $Chemin/*/*.bam | parallel -j 15 --bibtex --eta "java -Xmx4g -jar /opt/GenomeAnalysisTK.jar \
        -T HaplotypeCaller \
		-ERC GVCF \
		--variant_index_type LINEAR \
		--variant_index_parameter 128000 \
        -R $hg19 \
        -I {} \
	    -L /d2/Exome_analysis_BC/Sureselect_V4_Covered.bed \
        --genotyping_mode DISCOVERY \
        -stand_emit_conf 10 \
        -stand_call_conf 30 \
        -o {}.g.vcf"
	
#opt/samtools mpileup -Bug -f $hg19 -l /d2/Exome_analysis_BC/Sureselect_V4_Covered.bed $Path/realigned_recal_$3.bam | /opt/bcftools view -cvg - > $Path/Variants/$3.samtools 

################# Coverage analysis ######################

mkdir $Chemin/Coverage

mv $Chemin/*/Quality_control/*.hist.txt $Chemin/Coverage/

python /home/sdumont/Multiple_exome_plot.py $Chemin/Coverage/

################# GVCF merge into one VCF file ##################

gVCFToVCFfiltered "$Chemin"

exit 0
