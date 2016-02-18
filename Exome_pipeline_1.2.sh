#!/bin/bash

fastq1=$1
fastq2=$2
Name=$3
HPO=$4
hg19="/d2/Exome_analysis_BC/Exome_pipeline_datas/ucsc.hg19.fasta"
Chemin="/d2/Exome_analysis_BC/Results_"$Name

threads=18
Design='/d2/Exome_analysis_BC/Exome_pipeline_datas/SureSelect_CRE/SureSelect_CRE_Covered.bed'
Design_capture='/d2/Exome_analysis_BC/Exome_pipeline_datas/SureSelect_CRE/SureSelect_CRE_Covered_50bp_chrM.bed'
Design_capture_by_chr='/d2/Exome_analysis_BC/Exome_pipeline_datas/SureSelect_CRE/chr'
Genes_coverage_analysis='/d2/Exome_analysis_BC/Exome_pipeline_datas/MMR_bed/MMR_20_cano.bed'
OMIM_gene2patho='/d2/Exome_analysis_BC/Exome_pipeline_datas/OMIM/OMIM_gene_patho.txt'
GDI_file='/d2/Exome_analysis_BC/Exome_pipeline_datas/GDI_full_10282015.txt'
HPO_file='/d2/Exome_analysis_BC/Exome_pipeline_datas/HPO/ALL_SOURCES_ALL_FREQUENCIES_phenotype_to_genes.txt'
HPO_genesToPhenotype='/d2/Exome_analysis_BC/Exome_pipeline_datas/HPO/ALL_SOURCES_ALL_FREQUENCIES_diseases_to_genes_to_phenotypes.txt'
GO_file='/d2/Exome_analysis_BC/Exome_pipeline_datas/GeneToGO_ensembl_biomart_11-2015.txt'

# Chemins des programmes utilisés:
FastQC='/opt/FastQC/fastqc'
cutadapt='/opt/cutadapt-1.8.1/bin/cutadapt'
bwa='/opt/bwa'
samtools='/opt/samtools-1.3/bin/samtools-1.3'
bcftools='/opt/bcftools-1.3/bin/bcftools-1.3'
MarkDuplicates='/opt/picard-tools-1.119/MarkDuplicates.jar'
GATK='/opt/GenomeAnalysisTK.jar'
bedtools='/opt/bedtools/bedtools'
vcf_concat='/opt/vcftools_0.1.11/bin/vcf-concat'
annovar='/opt/annovar'
annovar_db='/d2/Exome_analysis_BC/Exome_pipeline_datas/humandb/'
Python_package='/home/bcogne/Python_package.py'
VariantsMysql='/home/bcogne/VariantsTomysql.py'

logfile="$Chemin/Quality_control/log.txt"

version='Exome_pipeline_nts_1.2'

### Check if correct number of arguments

if [ "$#" -ne 4 ]; then
    echo -e "\n\n\t\tYou didn't provide the correct number of arguments.\n\n\t\tUsage: ./final_exome_pipeline.sh /R1.fastq /R2.fastq Name_analysis\n"
	exit 0
fi


### Create files and folders
[ -d $Chemin ] || mkdir $Chemin || { echo "Creation of the main folder failed during ${Name} analysis - Check user authorization" ; exit 1; }
[ -d $Chemin/Quality_control ] || mkdir $Chemin/Quality_control $Chemin/temp $Chemin/Variants $Chemin/Datas || { echo "Creation of subfolders failed during ${Name} analysis" ; exit 1; }
[ -d $Chemin/Variants/samtools ] || mkdir $Chemin/Variants/samtools || { echo "Creation of subfolders failed during ${Name} analysis" ; exit 1; }
touch $Chemin/Quality_control/log.txt

### Create &3 to allow output on the screen 
exec 3>&1 4>&2
#exec 1>$Chemin/Quality_control/log.txt 2>&1

{
date_start=$(date +"%s")

echo -e "\n\t\t####################################################################\n" >&3
echo -e "\n\t\tWelcome in the ${version} from the CHU of Nantes genetics lab\n" >&3
echo -e "\n\t\tThe analysis of exome ${Name} begins at $(date)" >&3
echo -e "\n\t\t####################################################################\n\n" >&3




############### Quality control of fastq files ##################################

file=($Chemin/Quality_control/*fastqc.html)
[[ -f "$file" ]] || $FastQC $fastq1 $fastq2 -o $Chemin/Quality_control || { echo "FastQC failed during ${Name} analysis" ; exit 1; } & 

############### Trimming of paired-end reads with cutadapt: removal of Illumina adapters and reads with sequencing quality under Q30 - keeping only reads >40bases ##### 

date1=$(date +"%s")
echo -e -n "\n\t\tStart of cutadapt..." >&3

file=($Chemin/temp/trimmed.$Name*)
[[ -f "$file" ]] || python $cutadapt -q 30 -m 40 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -o $Chemin/temp/trimmed.$Name.1.fastq -p $Chemin/temp/trimmed.$Name.2.fastq $fastq1 $fastq2 || { echo "cutadapt failed during ${Name} analysis" ; exit 1; }
date2=$(date +"%s")
time_analysis=$(($date2-$date1))

echo -e " end of cutadapt in $(($time_analysis / 60)) min" >&3


############### Paired-end reads alignment by BWA mem ######

# For GATK compatibility, the bwa mem option -R is required: -R '@RG\tID:group1\tSM:sample1\tPL:illumina\tLB:lib1\tPU:unit1'
date1=$(date +"%s")
echo -e -n "\n\t\tStart of the alignment with bwa mem..." >&3

[ -f $Chemin/temp/$Name.sam ] || $bwa mem -M -R "@RG\tID:exome\tSM:${Name}\tPL:illumina\tLB:lib1" -t $threads $hg19 $Chemin/temp/trimmed.$Name.1.fastq $Chemin/temp/trimmed.$Name.2.fastq > $Chemin/temp/$Name.sam || { echo "bwa mem failed during ${Name} analysis" ; exit 1; } 

date2=$(date +"%s")
time_analysis=$(($date2-$date1))
echo -e " end of the alignment with bwa mem in $(($time_analysis / 60)) min" >&3


############## Sam -> sorted Bam conversion with samtools  - Duplicates PCR reads marking woth picard ##########
echo -e -n "\n\t\tStart of the Sam to Bam conversion..." >&3
date1=$(date +"%s")
[ -f $Chemin/temp/$Name.bam ] || $samtools view -hbS -@ $threads $Chemin/temp/$Name.sam > $Chemin/temp/$Name.bam
date2=$(date +"%s")
time_analysis=$(($date2-$date1))
echo -e " end of the Sam to Bam conversion in $(($time_analysis / 60)) min" >&3


echo -e -n "\n\t\tStart of the Bam sorting..." >&3
date1=$(date +"%s")
[ -f $Chemin/temp/sorted_$Name* ] || $samtools sort -@ $threads $Chemin/temp/$Name.bam -o $Chemin/temp/sorted_$Name.bam
date2=$(date +"%s")
time_analysis=$(($date2-$date1))
echo -e " end of samtools sort in $(($time_analysis / 60)) min" >&3

echo -e -n "\n\t\tStart of duplicates marking..." >&3
date1=$(date +"%s")
[ -f $Chemin/temp/sorted_dedup_$Name.bam ] || java -Xmx15g -jar $MarkDuplicates \
    INPUT=$Chemin/temp/sorted_$Name.bam \
    OUTPUT=$Chemin/temp/sorted_dedup_$Name.bam \
	REMOVE_DUPLICATES=true \
    METRICS_FILE=$Chemin/Quality_control/metrics_markduplicates.txt || { echo "Picard MarkDuplicates failed during ${Name} analysis" ; exit 1; }
date2=$(date +"%s")
time_analysis=$(($date2-$date1))
echo -e " end of duplicates marking in $(($time_analysis / 60)) min" >&3

echo -e -n "\n\t\tStart of samtools indexing..." >&3
date1=$(date +"%s")
[ -f $Chemin/temp/sorted_dedup_$Name*bai ] || $samtools index $Chemin/temp/sorted_dedup_$Name.bam
date2=$(date +"%s")
time_analysis=$(($date2-$date1))
echo -e " end of samtools indexing in $(($time_analysis / 60)) min" >&3


############## Realignment of reads on known indels with GATK #################


echo -e -n "\n\t\tStart of indels realignment with GATK..." >&3
date1=$(date +"%s")

[ -f $Chemin/temp/$Name.bam.list ] || java -Xmx15g -jar $GATK \
	-T RealignerTargetCreator \
	-L $Design_capture \
	-R $hg19 \
	-nt $threads \
	-known /d2/Exome_analysis_BC/Exome_pipeline_datas/1000G_phase1.indels.hg19.sites.vcf \
	-o $Chemin/temp/$Name.bam.list \
	-I $Chemin/temp/sorted_dedup_$Name.bam || { echo "GATK RealignerTargetCreator failed during ${Name} analysis" ; exit 1; }

[ -f $Chemin/temp/realigned_$Name.bam ] || java -Xmx15g -jar $GATK \
    -T IndelRealigner \
    -R $hg19 \
    -I $Chemin/temp/sorted_dedup_$Name.bam \
	-known /d2/Exome_analysis_BC/Exome_pipeline_datas/1000G_phase1.indels.hg19.sites.vcf \
    -targetIntervals $Chemin/temp/$Name.bam.list \
    -o $Chemin/temp/realigned_$Name.bam || { echo "GATK IndelRealigner failed during ${Name} analysis" ; exit 1; }

date2=$(date +"%s")
time_analysis=$(($date2-$date1))
	
echo -e " end of indels realignment with GATK in $(($time_analysis / 60)) min" >&3


############### Base quality score recalibration with GATK ##################

echo -e -n "\n\t\tStart of Base recalibration with GATK..." >&3
date1=$(date +"%s")

## Analyze patterns of covariation in the sequence dataset ##

[ -f $Chemin/temp/recal_data.table ] || java -Xmx15g -jar $GATK \
    -T BaseRecalibrator \
	-L $Design_capture \
    -R $hg19 \
    -nct 2 \
    -I $Chemin/temp/realigned_$Name.bam \
    -knownSites /d2/Exome_analysis_BC/Exome_pipeline_datas/dbsnp_138.hg19.vcf \
    -o $Chemin/temp/recal_data.table || { echo "GATK BaseRecalibrator failed during ${Name} analysis" ; exit 1; }

date2=$(date +"%s")
time_analysis=$(($date2-$date1))
echo -e " end of Base recalibration with GATK in $(($time_analysis / 60)) min" >&3

# Do a second pass to analyze covariation remaining after recalibration ## FACULTATIF -> Verification sur plot recalibration
	
## Apply recalibration to dataset independently for each chromosome - multithreading mode ##

echo -e -n "\n\t\tPrinting of base recalibrated reads with GATK per chromosome..." >&3
date1=$(date +"%s")

function chrPR {
	chr_path=$1
	Name=$2
	hg19="/d2/Exome_analysis_BC/Exome_pipeline_datas/ucsc.hg19.fasta"
	Chemin="/d2/Exome_analysis_BC/Results_"$Name
	GATK='/opt/GenomeAnalysisTK.jar'
	
	Chem_split=(${chr_path//\// })
	Name_bed="${Chem_split[${#Chem_split[@]} - 1]}"
	Bed_split=(${Name_bed//\./ })
	Name_chr="${Bed_split[0]}"

	java -Xmx4g -jar $GATK \
		-T PrintReads \
		-R $hg19 \
		-I $Chemin/temp/realigned_$Name.bam \
		-BQSR $Chemin/temp/recal_data.table \
		-L $chr_path \
		-o $Chemin/Datas/realigned_recal_${Name}.${Name_chr}.bam
}

export -f chrPR
list=($Design_capture_by_chr/*.bed)
file=($Chemin/Datas/*.bam)
[[ -f "$file" ]] || parallel -j 5 --halt now,fail=1 --eta "chrPR {} ${Name}" ::: ${list[@]}

date2=$(date +"%s")
time_analysis=$(($date2-$date1))
echo -e " end of reads base recalibration with GATK in $(($time_analysis / 60)) min" >&3



################ Variant Calling : Generation of gVCF files with Haplotype Caller - multithreaded by chromosome  ########################################

echo -e -n "\n\t\tStart of variant calling with haplotype caller in gVCF mode..." >&3
date1=$(date +"%s")

function chrHC {
	
	chr_path=$1
	Name=$2
	hg19="/d2/Exome_analysis_BC/Exome_pipeline_datas/ucsc.hg19.fasta"
	Chemin="/d2/Exome_analysis_BC/Results_"$Name
	GATK='/opt/GenomeAnalysisTK.jar'
	
	Chem_split=(${chr_path//\// })
	Name_bed="${Chem_split[${#Chem_split[@]} - 1]}"
	Bed_split=(${Name_bed//\./ })
	Name_chr="${Bed_split[0]}"

	java -Xmx4g -jar $GATK \
        -T HaplotypeCaller \
		-ERC GVCF \
		--variant_index_type LINEAR \
		--variant_index_parameter 128000 \
        -R $hg19 \
        -I $Chemin/Datas/realigned_recal_$Name.${Name_chr}.bam \
	    -L $chr_path \
        --genotyping_mode DISCOVERY \
        -stand_emit_conf 10 \
        -stand_call_conf 30 \
        -o $Chemin/Variants/variants_HC_$Name.${Name_chr}.g.vcf || { echo "GATK HaplotypeCaller failed to generate gVCF during ${Name} analysis of ${Name_chr}"  ; exit 1; }
}

export -f chrHC

list=($Design_capture_by_chr/*.bed)
[ -f $Chemin/Variants/variants_HC_$Name.allchr.g.vcf ] || parallel -j 5 --halt now,fail=1 --eta "chrHC {} ${Name}" ::: ${list[@]}

#launch of samtools variant calling
$samtools mpileup -I -ugf $hg19 $Chemin/Datas/realigned_recal_$Name.bam | $bcftools call -V indels -vmO z -o $Chemin/Variants/samtools/variants_samtools.${Name}.vcf.gz &

date2=$(date +"%s")
time_analysis=$(($date2-$date1))
echo -e " end of variant calling in $(($time_analysis / 60)) min" >&3



echo -e -n "\n\t\tStart of per chromosome gVCFs merging with CatVariants..." >&3
date1=$(date +"%s")

[ -f $Chemin/Variants/variants_HC_$Name.allchr.g.vcf ] || java -cp $GATK org.broadinstitute.gatk.tools.CatVariants \
    -R $hg19 \
	-V $Chemin/Variants/variants_HC_$Name.chrM.g.vcf \
    -V $Chemin/Variants/variants_HC_$Name.chr1.g.vcf \
	-V $Chemin/Variants/variants_HC_$Name.chr2.g.vcf \
	-V $Chemin/Variants/variants_HC_$Name.chr3.g.vcf \
	-V $Chemin/Variants/variants_HC_$Name.chr4.g.vcf \
	-V $Chemin/Variants/variants_HC_$Name.chr5.g.vcf \
	-V $Chemin/Variants/variants_HC_$Name.chr6.g.vcf \
	-V $Chemin/Variants/variants_HC_$Name.chr7.g.vcf \
	-V $Chemin/Variants/variants_HC_$Name.chr8.g.vcf \
	-V $Chemin/Variants/variants_HC_$Name.chr9.g.vcf \
	-V $Chemin/Variants/variants_HC_$Name.chr10.g.vcf \
	-V $Chemin/Variants/variants_HC_$Name.chr11.g.vcf \
	-V $Chemin/Variants/variants_HC_$Name.chr12.g.vcf \
	-V $Chemin/Variants/variants_HC_$Name.chr13.g.vcf \
	-V $Chemin/Variants/variants_HC_$Name.chr14.g.vcf \
	-V $Chemin/Variants/variants_HC_$Name.chr15.g.vcf \
	-V $Chemin/Variants/variants_HC_$Name.chr16.g.vcf \
	-V $Chemin/Variants/variants_HC_$Name.chr17.g.vcf \
	-V $Chemin/Variants/variants_HC_$Name.chr18.g.vcf \
	-V $Chemin/Variants/variants_HC_$Name.chr19.g.vcf \
	-V $Chemin/Variants/variants_HC_$Name.chr20.g.vcf \
	-V $Chemin/Variants/variants_HC_$Name.chr21.g.vcf \
	-V $Chemin/Variants/variants_HC_$Name.chr22.g.vcf \
	-V $Chemin/Variants/variants_HC_$Name.chrX.g.vcf \
	-V $Chemin/Variants/variants_HC_$Name.chrY.g.vcf \
    -out $Chemin/Variants/variants_HC_$Name.allchr.g.vcf \
    -assumeSorted

file=($Chemin/Variants/variants_HC_$Name.chr*.g.vcf*)
[[ -f "$file" ]] && rm $Chemin/Variants/variants_HC_$Name.chr*.g.vcf*;

date2=$(date +"%s")
time_analysis=$(($date2-$date1))
echo -e " end of gVCF merging in $(($time_analysis / 60)) min" >&3

############### Genotyping of gVCF file together with other gVCFs in order to perform VQSR #####################################################################

echo -e -n "\n\t\tStart of the genotyping of gVCF with GATK GenotypeGVCFs..." >&3
date1=$(date +"%s")


touch $Chemin/temp/gvcf.list
echo -e "${Chemin}/Variants/variants_HC_${Name}.allchr.g.vcf" > $Chemin/temp/gvcf.list;
list_gvcf=(/d2/Exome_analysis_BC/Research-exomes/Hugodims/RUN1-5/Results_*/Variants/*.g.vcf)
	
for i in ${list_gvcf[@]}; 
	do
	    echo $i >> $Chemin/temp/gvcf.list;
	done;
    
[ -f $Chemin/temp/variants_VQSR.combined.vcf ] || java -Xmx15g -jar $GATK \
    -T GenotypeGVCFs \
    -nt $threads \
    -L $Design_capture \
    -R $hg19 \
    -V $Chemin/temp/gvcf.list \
    -o $Chemin/temp/variants_combined.vcf || { echo "GATK GenotypeGVCFs failed to generate VCF during ${Name}" ; exit 1; }
	
date2=$(date +"%s")
time_analysis=$(($date2-$date1))
echo -e " end of gVCF genotyping in $(($time_analysis / 60)) min" >&3

		
############### Variant Quality Score Recalibration: add -an InbreedingCoeff if samples > 10 ########
   
echo -e -n "\n\t\tStart of Variant Quality Score Recalibration..." >&3
date1=$(date +"%s")

[ -f $Chemin/temp/variants_VQSR.combined.vcf ] || java -Xmx10g -jar $GATK \
    -T VariantRecalibrator \
    -R $hg19 \
    -input $Chemin/temp/variants_combined.vcf \
    -resource:hapmap,known=false,training=true,truth=true,prior=15.0 /d2/Exome_analysis_BC/Exome_pipeline_datas/hapmap_3.3.hg19.sites.vcf \
    -resource:omni,known=false,training=true,truth=true,prior=12.0 /d2/Exome_analysis_BC/Exome_pipeline_datas/1000G_omni2.5.hg19.sites.vcf \
    -resource:1000G,known=false,training=true,truth=false,prior=10.0 /d2/Exome_analysis_BC/Exome_pipeline_datas/1000G_phase1.snps.high_confidence.hg19.sites.vcf \
    -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /d2/Exome_analysis_BC/Exome_pipeline_datas/dbsnp_138.hg19.vcf \
    -an QD -an FS -an SOR -an MQ -an MQRankSum -an ReadPosRankSum \
    -mode SNP \
    -nt $threads \
    -tranche 100.0 -tranche 99.9  -tranche 99.5 -tranche 99.0 -tranche 90.0 \
    -recalFile $Chemin/temp/recalibrate_SNP.recal \
    -tranchesFile $Chemin/temp/recalibrate_SNP.tranches \
    -rscriptFile $Chemin/temp/recalibrate_SNP_plots.R

[ -f $Chemin/temp/variants_VQSR.combined.vcf ] || java -Xmx10g -jar $GATK \
    -T ApplyRecalibration \
    -R $hg19 \
    -input $Chemin/temp/variants_combined.vcf \
    -mode SNP \
    --ts_filter_level 99.5 \
    -recalFile $Chemin/temp/recalibrate_SNP.recal \
    -tranchesFile $Chemin/temp/recalibrate_SNP.tranches \
    -o $Chemin/temp/recalibrated_snps_raw_indels.vcf 
	
[ -f $Chemin/temp/variants_VQSR.combined.vcf ] || java -Xmx10g -jar $GATK \
    -T VariantRecalibrator \
    -R $hg19 \
    -input $Chemin/temp/recalibrated_snps_raw_indels.vcf \
    -resource:mills,known=true,training=true,truth=true,prior=12.0  /d2/Exome_analysis_BC/Exome_pipeline_datas/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf \
    -an QD -an FS -an SOR -an MQRankSum -an ReadPosRankSum \
    -mode INDEL \
    -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
    --maxGaussians 4 \
    -recalFile $Chemin/temp/recalibrate_INDEL.recal \
    -tranchesFile $Chemin/temp/recalibrate_INDEL.tranches \
    -rscriptFile $Chemin/temp/recalibrate_INDEL_plots.R
		
[ -f $Chemin/temp/variants_VQSR.combined.vcf ] || java -Xmx10g -jar $GATK \
    -T ApplyRecalibration \
    -R $hg19 \
    -input $Chemin/temp/recalibrated_snps_raw_indels.vcf \
    -mode INDEL \
    --ts_filter_level 99.0 \
    -recalFile $Chemin/temp/recalibrate_INDEL.recal \
    -tranchesFile $Chemin/temp/recalibrate_INDEL.tranches \
    -o $Chemin/temp/variants_VQSR.combined.vcf
	
[[ -f $Chemin/temp/recalibrate_SNP.tranches.pdf ]] && mv $Chemin/temp/recalibrate_SNP.tranches.pdf $Chemin/Quality_control/recalibrate_SNP.tranches.${Name}.pdf

date2=$(date +"%s")
time_analysis=$(($date2-$date1))
echo -e " end of VQSR in $(($time_analysis / 60)) min" >&3


#$samtools mpileup -Bug -f $hg19 -l $Design_capture $Chemin/realigned_recal_$Name.bam | /opt/bcftools view -cvg - > $Chemin/Variants/$Name.samtools 

################ Extraction of the sample of interest with GATK Select Variant ##############################

echo -e -n "\n\t\tStart of the post-processing of candidate variants..." >&3
date1=$(date +"%s")

## -env option excludes non-variant sites

[ -f $Chemin/Variants/variants_VQSR.${Name}.vcf ] || java -Xmx4g -jar $GATK \
   -T SelectVariants \
   -R $hg19 \
   -V $Chemin/temp/variants_VQSR.combined.vcf \
   -o $Chemin/Variants/variants_VQSR.${Name}.vcf \
   -env \
   -sn ${Name}

################ Resolve multiallelic sites with one variant not in the sample ##########################

[ -f $Chemin/Variants/variants_VQSR.${Name}.biallelic.vcf ] || python $Python_package ReduceMultiallelicVCF $Chemin/Variants/ ${Name} $Chemin/Variants/variants_VQSR.${Name}.vcf

#python $Python_package Addmpilup $Chemin/Variants/$Name.samtools $Chemin/Variants/HC_vcf_sorted_$Name.vcf
#mv $Chemin/Variants/HC_vcf_sorted_$Name.vcf.QUALSAMHEAD $Chemin/Variants/Final_variants_$Name.vcf
#rm -r $Chemin/temp
#rm $Chemin/Variants/raw_snps_$Name.vcf* $Chemin/Variants/filtered_snps_$Name.vcf* $Chemin/Variants/raw_indels_$Name.vcf* $Chemin/Variants/filtered_indels_$Name.vcf* $Chemin/Variants/vcf_concat_variants_$Name.vcf $Chemin/Variants/HC_vcf_sorted_$Name.vcf.QUALSAM $Chemin/Variants/HC_vcf_sorted_$Name.vcf $Chemin/Variants/*.samtools


############################# MAF and mutation effect filtration and annotation using ANNOVAR and Python home-made scripts ###########################

# All databases used must be downloaded from annovar server and put in humandb folder of annovar

[ -f $Chemin/Variants/variants_VQSR.${Name}.biallelic.annovar ] || $annovar/convert2annovar.pl -format vcf4 -includeinfo -withzyg $Chemin/Variants/variants_VQSR.${Name}.biallelic.vcf -outfile $Chemin/Variants/variants_VQSR.${Name}.biallelic.annovar

# Annotation of all variants with annovar without filtration to put them in the exome database
file=($Chemin/Variants/variants_VQSR.${Name}.allforDB*)
[[ -f "$file" ]] || $annovar/table_annovar.pl $Chemin/Variants/variants_VQSR.${Name}.biallelic.annovar $annovar_db -buildver hg19 -out $Chemin/Variants/variants_VQSR.${Name}.allforDB -remove -protocol refGene,esp6500siv2_all,1000g2015aug_all,1000g2015aug_afr,1000g2015aug_eas,1000g2015aug_eur,exac03,avsnp144,clinvar_20150330,ljb26_all -operation g,f,f,f,f,f,f,f,f,f -nastring . -otherinfo -arg '-exonicsplicing -hgvs',,,,,,,,,

# Filtration and annotation 
file=($Chemin/Variants/variants_VQSR.${Name}.biallelic.1000g-0.01*)
[[ -f "$file" ]] || $annovar/annotate_variation.pl -filter -dbtype 1000g2015aug_all -buildver hg19 -out $Chemin/Variants/variants_VQSR.${Name}.biallelic.1000g-0.01 $Chemin/Variants/variants_VQSR.${Name}.biallelic.annovar $annovar_db -score_threshold 0.01

# Annotation avec annovar: splicing annotation à +- 2bp (changement avec argument --splicing_threshold X), si intronic annovar annote: "splicing", si exonic annovar annote: "exonic;splicing"

file=($Chemin/Variants/variants.${Name}.1000g-0.01.annotated*)
[[ -f "$file" ]] || $annovar/table_annovar.pl $Chemin/Variants/variants_VQSR.${Name}.biallelic.1000g-0.01*filtered $annovar_db -buildver hg19 -out $Chemin/Variants/variants.${Name}.1000g-0.01.annotated -remove -protocol refGene,esp6500siv2_all,1000g2015aug_all,1000g2015aug_afr,1000g2015aug_eas,1000g2015aug_eur,exac03,avsnp144,clinvar_20150330,ljb26_all -operation g,f,f,f,f,f,f,f,f,f -nastring . -otherinfo -arg '-exonicsplicing -hgvs',,,,,,,,,

# Keep only exonic nonsynonymous/indels and splicing intronic/exonic mutations

printf '#' > $Chemin/Variants/Final_variants_${Name}_exon_splice_filtered.txt
head -n 1 $Chemin/Variants/*.annotated.hg19_multianno.txt >> $Chemin/Variants/Final_variants_${Name}_exon_splice_filtered.txt
awk '{if (($6~/^exonic/ && $9!~/^synonymous/) || ($6~/splicing/)) {print $0}}' $Chemin/Variants/*.annotated.hg19_multianno.txt | grep 'PASS' | sort -k 7 >> $Chemin/Variants/Final_variants_${Name}_exon_splice_filtered.txt

### Variants prioritization and IGV Snapshots: generation of the final report

[ -f $Chemin/Variants/Final_variants_${Name}_exon_splice_filtered.txt.python_filtered.tab ] || python $Python_package VariantsPrioritize $Chemin/Variants/Final_variants_${Name}_exon_splice_filtered.txt $OMIM_gene2patho $HPO_file $HPO_genesToPhenotype $GDI_file $Name $HPO

date2=$(date +"%s")
time_analysis=$(($date2-$date1))
echo -e " end of variants post-processing in $(($time_analysis / 60)) min" >&3


################## Merging of per chromosome BAMs ####################################

echo -e -n "\n\t\tStart of per chromosome BAM merging with samtools..." >&3
date1=$(date +"%s")

[ -f $Chemin/Datas/realigned_recal_$Name.bam ] || $samtools merge -@ $threads $Chemin/Datas/realigned_recal_$Name.bam $Chemin/Datas/*.chrM.bam $Chemin/Datas/*.chr1.bam $Chemin/Datas/*.chr2.bam \
	$Chemin/Datas/*.chr3.bam $Chemin/Datas/*.chr4.bam $Chemin/Datas/*.chr5.bam $Chemin/Datas/*.chr6.bam $Chemin/Datas/*.chr7.bam \
	$Chemin/Datas/*.chr8.bam $Chemin/Datas/*.chr9.bam $Chemin/Datas/*.chr10.bam $Chemin/Datas/*.chr11.bam $Chemin/Datas/*.chr12.bam \
	$Chemin/Datas/*.chr13.bam $Chemin/Datas/*.chr14.bam $Chemin/Datas/*.chr15.bam $Chemin/Datas/*.chr16.bam $Chemin/Datas/*.chr17.bam \
	$Chemin/Datas/*.chr18.bam $Chemin/Datas/*.chr19.bam $Chemin/Datas/*.chr20.bam $Chemin/Datas/*.chr21.bam $Chemin/Datas/*.chr22.bam \
	$Chemin/Datas/*.chrX.bam $Chemin/Datas/*.chrY.bam 

[ -f $Chemin/Datas/realigned_recal_$Name.bam.bai ] || $samtools index $Chemin/Datas/realigned_recal_$Name.bam 

file=($Chemin/Datas/*.chr*.bam)
[[ -f "$file" ]] && rm $Chemin/Datas/*.chr*.ba*;

date2=$(date +"%s")
time_analysis=$(($date2-$date1))
echo -e " end of BAM merging in $(($time_analysis / 60)) min" >&3

####################### Pindel calling ###########################
'''
echo -e -n "\n\t\tStart pindel calling..." >&3
date1=$(date +"%s")

echo -e "${Chemin}/Datas/realigned_recal_${Name}.bam\t250\t${Name}" > ${Chemin}/temp/pindel.txt

function chrPindel {
	chr_path=$1
	Name=$2
    pindel_file=$3
	hg19="/d2/Exome_analysis_BC/Exome_pipeline_datas/ucsc.hg19.fasta"
	Chemin="/d2/Exome_analysis_BC/Results_"$Name
	pindel="/opt/pindel-master/pindel"
	
	Chem_split=(${chr_path//\// })
	Name_bed="${Chem_split[${#Chem_split[@]} - 1]}"
	Bed_split=(${Name_bed//\./ })
	Name_chr="${Bed_split[0]}"

	$pindel -f $hg19 -i $pindel_file -x 1 -j $chr_path -o ${Chemin}/temp/${Name}.${Name_chr}
	
}

export -f chrPindel
list=($Design_capture_by_chr/*.bed)
parallel -j 8 --halt now,fail=1 --eta "chrPindel {} ${Name} ${Chemin}/temp/pindel.txt" ::: ${list[@]} >/dev/null

cat ${Chemin}/temp/${Name}.${Name_chr}*_D > ${Chemin}/temp/${Name}..allchr_D
cat ${Chemin}/temp/${Name}.${Name_chr}*_SI > ${Chemin}/temp/${Name}..allchr_SI
cat ${Chemin}/temp/${Name}.${Name_chr}*_INV > ${Chemin}/temp/${Name}..allchr_INV
cat ${Chemin}/temp/${Name}.allchr_D ${Chemin}/temp/${Name}.allchr_SI ${Chemin}/temp/${Name}.allchr_INV > ${Chemin}/temp/${Name}.allchr_D_SI_INV

/opt/pindel-master/pindel2vcf -p ${Chemin}/temp/${Name}.allchr_D_SI_INV -e 5 -ir 2 -pr 2 -G -r $hg19 -R hg19 -d 20101123 -v ${Chemin}/temp/${Name}.allchr_D_SI_INV.vcf

date2=$(date +"%s")
time_analysis=$(($date2-$date1))
echo -e " end of pindel calling in $(($time_analysis / 60)) min" >&3
'''

################## Check Bam quality ####################################

echo -e -n "\n\t\tStart quality controls..." >&3
date1=$(date +"%s")

[ -f $Chemin/Quality_control/Stat_mapping_$Name.txt ] || $samtools flagstat $Chemin/Datas/realigned_recal_$Name.bam > $Chemin/Quality_control/Stat_mapping_$Name.txt &

file=($Chemin/Quality_control/GATK_DoC_MMR_$Name*)
[[ -f "$file" ]] || java -Xmx15g -jar $GATK \
   -T DepthOfCoverage \
   -R $hg19 \
   -o $Chemin/Quality_control/GATK_DoC_MMR_$Name \
   -I $Chemin/Datas/realigned_recal_$Name.bam \
   -geneList /d2/Exome_analysis_BC/Exome_pipeline_datas/refseq_gatk.txt \
   -ct 4 -ct 10 -ct 20 -ct 30 \
   -L $Genes_coverage_analysis &
   
[ -d $Chemin/Variants/snapshots_$Name ] || mkdir $Chemin/Variants/snapshots_$Name

file=($Chemin/Variants/snapshots_$Name/*.png)
[[ -f "$file" ]] || python $Python_package IGVSnapshots $Chemin/Variants/Final_variants_${Name}_exon_splice_filtered.txt $Chemin/temp/batch.txt $Chemin/Variants/snapshots_$Name/ $Chemin/Datas/realigned_recal_$Name.bam

[ -f $Chemin/Quality_control/$Name.hist.txt ] || $bedtools coverage -abam $Chemin/Datas/realigned_recal_$Name.bam -b $Design -hist | grep ^all > $Chemin/Quality_control/$Name.hist.txt 

[ -f $Chemin/Quality_control/Exome_cumulative_coverage.pdf ] || python $Python_package ExomePlot $Chemin/Quality_control/ 


date2=$(date +"%s")
time_analysis=$(($date2-$date1))
echo -e " end of quality controls in $(($time_analysis / 60)) min" >&3


############### Insert values into mysql table ############################

DATE=`date +%Y-%m-%d`

echo -e "USE Exome_db;" > $Chemin/temp/mysql_id.txt
echo -e "INSERT INTO exome_id VALUES (NULL,'${Name}','Clinical','Illumina','Integragen','Sureselect_CRE','${DATE}','','${HPO}','${version}','n');" >> $Chemin/temp/mysql_id.txt

mysql -u bcogne -p bencog15 < $Chemin/temp/mysql_id.txt

python $VariantsMysql $Chemin/Variants/variants_VQSR.${Name}.allforDB.hg19_multianno.txt ${Name}






############ Variant annotation of samtools calling#############################

echo -e -n "\n\t\tStart of Samtools Variant Calling annotation..." >&3
date1=$(date +"%s")

wait

gunzip $Chemin/Variants/samtools/variants_samtools.${Name}.vcf.gz

$annovar/convert2annovar.pl -format vcf4 -includeinfo -withzyg $Chemin/Variants/samtools/variants_samtools.${Name}.vcf -outfile $Chemin/temp/variants_samtools.${Name}.annovar

$annovar/annotate_variation.pl -filter -dbtype 1000g2015aug_all -buildver hg19 -out $Chemin/temp/variants_samtools.${Name}.1000g-0.01 $Chemin/temp/variants_samtools.${Name}.annovar $annovar_db -score_threshold 0.01

# Annotation avec annovar: splicing annotation à +- 2bp (changement avec argument --splicing_threshold X), si intronic annovar annote: "splicing", si exonic annovar annote: "exonic;splicing"

$annovar/table_annovar.pl $Chemin/temp/variants_samtools.${Name}.1000g-0.01*filtered $annovar_db -buildver hg19 -out $Chemin/temp/variants_samtools.${Name}.1000g-0.01.annotated -remove -protocol refGene,esp6500siv2_all,1000g2015aug_all,1000g2015aug_afr,1000g2015aug_eas,1000g2015aug_eur,exac03,avsnp144,clinvar_20150330,ljb26_all -operation g,f,f,f,f,f,f,f,f,f -nastring . -otherinfo -arg '-exonicsplicing -hgvs',,,,,,,,,

# Keep only exonic nonsynonymous/indels and splicing intronic/exonic mutations

printf '#' > $Chemin/temp/variants_samtools.${Name}_exon_splice_filtered.txt
head -n 1 $Chemin/temp/variants_samtools*.annotated.hg19_multianno.txt >> $Chemin/temp/variants_samtools.${Name}_exon_splice_filtered.txt
awk '{if (($6~/^exonic/) || ($6~/splicing/)) {print $0}}' $Chemin/temp/variants_samtools*.annotated.hg19_multianno.txt | sort -k 7 >> $$Chemin/temp/variants_samtools.${Name}_exon_splice_filtered.txt

python $Python_package VariantsPrioritize $Chemin/temp/variants_samtools.${Name}_exon_splice_filtered.txt $OMIM_gene2patho $HPO_file $HPO_genesToPhenotype $GDI_file $Name $HPO

mv $Chemin/temp/variants_samtools.${Name}_exon_splice_filtered.txt.python_filtered_A* $Chemin/Variants/samtools

date2=$(date +"%s")
time_analysis=$(($date2-$date1))
echo -e " end of Samtools Variant Calling annotation in $(($time_analysis / 60)) min" >&3


################ END OF PIPELINE #####################################

date_end=$(date +"%s")
diff_time=$(($date_end-$date_start))

echo -e "\n\t\tEnd of the analysis of exome ${Name} in $(($diff_time / 3600)) hour(s) and $((($diff_time / 60) %60)) minutes \n\n" >&3

} 1>$logfile 2>&1      
