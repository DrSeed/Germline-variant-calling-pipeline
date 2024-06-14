


# directories
ref="/home/drseed/VCF_analysis/supporting_files/hg38.fa"
known_sites="/home/drseed/VCF_analysis/supporting_files/Homo_sapiens_assembly38.dbsnp138.vcf"
aligned_reads="/home/drseed/VCF_analysis/aligned_reads"
reads="/home/drseed/VCF_analysis/reads"
results="/home/drseed/VCF_analysis/results"
data="/home/drseed/VCF_analysis/data"

#echo "STEP 1: QC - Run fastqc"
fastqc ${reads}/SRR062634_1.filt.fastq.gz -o ${reads}/
fastqc ${reads}/SRR062634_2.filt.fastq.gz -o ${reads}/ 

# --------------------------------------
# STEP 2: Map to reference using BWA-MEM
# --------------------------------------

echo "STEP 2: Map to reference using BWA-MEM"

# BWA index reference 
bwa index ${ref}


# BWA alignment
bwa mem -t 128 -R "@RG\tID:SRR062634\tPL:ILLUMINA\tSM:SRR062634" ${ref} ${reads}/SRR062634_1.filt.fastq.gz ${reads}/SRR062634_2.filt.fastq.gz > ${aligned_reads}/SRR062634.paired.sam


# -----------------------------------------
# STEP 3: Mark Duplicates and Sort - GATK4
# -----------------------------------------

echo "STEP 3: Mark Duplicates and Sort - GATK4"

gatk MarkDuplicatesSpark -I ${aligned_reads}/SRR062634.paired.sam -O ${aligned_reads}/SRR062634_sorted_dedup_reads.bam

# ----------------------------------
# STEP 4: Base quality recalibration
# ----------------------------------


echo "STEP 4: Base quality recalibration"

# 1. build the model
gatk BaseRecalibrator -I ${aligned_reads}/SRR062634_sorted_dedup_reads.bam -R ${ref} --known-sites ${known_sites} -O ${data}/recal_data.table 


# 2. Apply the model to adjust the base quality scores
gatk ApplyBQSR -I ${aligned_reads}/SRR062634_sorted_dedup_reads.bam -R ${ref} --bqsr-recal-file ${data}/recal_data.table -O ${aligned_reads}/SRR062634_sorted_dedup_bqsr_reads.bam 



# -----------------------------------------------
# STEP 5: Collect Alignment & Insert Size Metrics
# -----------------------------------------------


echo "STEP 5: Collect Alignment & Insert Size Metrics"

gatk CollectAlignmentSummaryMetrics R=${ref} I=${aligned_reads}/SRR062634_sorted_dedup_bqsr_reads.bam O=${aligned_reads}/alignment_metrics.txt
gatk CollectInsertSizeMetrics INPUT=${aligned_reads}/SRR062634_sorted_dedup_bqsr_reads.bam OUTPUT=${aligned_reads}/insert_size_metrics.txt HISTOGRAM_FILE=${aligned_reads}/insert_size_histogram.pdf


# ----------------------------------------------
# STEP 6: Call Variants - gatk haplotype caller
# ----------------------------------------------

echo "STEP 6: Call Variants - gatk haplotype caller"

gatk HaplotypeCaller -R ${ref} -I ${aligned_reads}/SRR062634_sorted_dedup_bqsr_reads.bam -O ${results}/raw_variants.vcf --native-pair-hmm-threads 128



# extract SNPs & INDELS

gatk SelectVariants -R ${ref} -V ${results}/raw_variants.vcf --select-type SNP -O ${results}/raw_snps.vcf
gatk SelectVariants -R ${ref} -V ${results}/raw_variants.vcf --select-type INDEL -O ${results}/raw_indels.vcf


# To filter varaints
# Filter SNPs
gatk VariantFiltration \
	-R ${ref} \
	-V ${results}/raw_snps.vcf \
	-O ${results}/filtered_snps.vcf \
	-filter-name "QD_filter" -filter "QD < 2.0" \
	-filter-name "FS_filter" -filter "FS > 60.0" \
	-filter-name "MQ_filter" -filter "MQ < 40.0" \
	-filter-name "SOR_filter" -filter "SOR > 4.0" \
	-filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" \
	-filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0" \
	-genotype-filter-expression "DP < 10" \
	-genotype-filter-name "DP_filter" \
	-genotype-filter-expression "GQ < 10" \
	-genotype-filter-name "GQ_filter"



# Filter INDELS
gatk VariantFiltration \
	-R ${ref} \
	-V ${results}/raw_indels.vcf \
	-O ${results}/filtered_indels.vcf \
	-filter-name "QD_filter" -filter "QD < 2.0" \
	-filter-name "FS_filter" -filter "FS > 200.0" \
	-filter-name "SOR_filter" -filter "SOR > 10.0" \
	-genotype-filter-expression "DP < 10" \
	-genotype-filter-name "DP_filter" \
	-genotype-filter-expression "GQ < 10" \
	-genotype-filter-name "GQ_filter"


# Select Variants that PASS filters
gatk SelectVariants \
	--exclude-filtered \
	-V ${results}/filtered_snps.vcf \
	-O ${results}/analysis-ready-snps.vcf


gatk SelectVariants \
	--exclude-filtered \
	-V ${results}/filtered_indels.vcf \
	-O ${results}/analysis-ready-indels.vcf



# Annotate using Funcotator
gatk Funcotator \
	--variant ${results}/analysis-ready-snps-filteredGT.vcf \
	--reference ${ref} \
	--ref-version hg38 \
	--data-sources-path /home/drseed/VCF_analysis/supporting_files/funcotator/funcotator_dataSources.v1.7.20200521g \
	--output ${results}/analysis-ready-snps-filteredGT-functotated.vcf \
	--output-file-format VCF

gatk Funcotator \
	--variant ${results}/analysis-ready-indels-filteredGT.vcf \
	--reference ${ref} \
	--ref-version hg38 \
	--data-sources-path /home/drseed/VCF_analysis/supporting_files/funcotator/funcotator_dataSources.v1.7.20200521g \
	--output ${results}/analysis-ready-indels-filteredGT-functotated.vcf \
	--output-file-format VCF

# Extract fields from a VCF file to a tab-delimited table

gatk VariantsToTable \
	-V ${results}/analysis-ready-snps-filteredGT-functotated.vcf -F AC -F AN -F DP -F AF -F FUNCOTATION \
	-O ${results}/output_snps.table



#GREPPING FOR GENE OF INTREST "CROCCP2"
# I am only intrested in looking at annotations in the gene "CROCCP2"
cat output_snps.table | less
cat analysis-ready-snps-filteredGT-functotated.vcf | less
cat analysis-ready-snps-filteredGT-functotated.vcf | grep " Funcotation fields are: " | less
cat analysis-ready-snps-filteredGT-functotated.vcf | grep " Funcotation fields are: " | sed 's/|/\t/g' | less
cat analysis-ready-snps-filteredGT-functotated.vcf | grep " Funcotation fields are: " | sed 's/|/\t/g' > output_curated_variants.txt
cat output_snps.table | less
cat output_snps.table | cut -f 5 | grep "CROCCP2" | less
cat output_snps.table | cut -f 5 | grep "CROCCP2" | sed 's/|/\t/g' | less
cat output_snps.table | cut -f 5 | grep "CROCCP2" | sed 's/|/\t/g' >> output_curated_variants.txt


# THE output_curated_variants.txt can then be looked at in excel for more information
