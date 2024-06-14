# README: Germline VCF Analysis Pipeline

## Introduction

This README provides a step-by-step guide for variant calling and analysis using the Genome Analysis Toolkit (GATK) and other bioinformatics tools. This pipeline is designed to handle raw sequencing data, align it to a reference genome, perform quality control, and call variants. The resulting variants are then filtered, annotated, and prepared for downstream analysis.

## Directories

I have organized my project into several directories:

- `ref`: Contains the reference genome file.
- `known_sites`: Contains known variant sites for base quality score recalibration.
- `aligned_reads`: Stores aligned read files.
- `reads`: Contains the raw and processed sequencing reads.
- `results`: Stores the final results of the variant calling process.
- `data`: Holds intermediate data files.

## Steps

### Step 1: Quality Control

I started by running FastQC on the filtered FASTQ files to assess the quality of my sequencing reads.

```bash
fastqc ${reads}/SRR062634_1.filt.fastq.gz -o ${reads}/
fastqc ${reads}/SRR062634_2.filt.fastq.gz -o ${reads}/
```

### Step 2: Mapping to Reference

Next, I mapped the reads to the reference genome using BWA-MEM.

```bash
# BWA index reference
bwa index ${ref}

# BWA alignment
bwa mem -t 128 -R "@RG\tID:SRR062634\tPL:ILLUMINA\tSM:SRR062634" ${ref} ${reads}/SRR062634_1.filt.fastq.gz ${reads}/SRR062634_2.filt.fastq.gz > ${aligned_reads}/SRR062634.paired.sam
```

### Step 3: Mark Duplicates and Sort

I used GATK MarkDuplicatesSpark to mark duplicates and sort the aligned reads.

```bash
gatk MarkDuplicatesSpark -I ${aligned_reads}/SRR062634.paired.sam -O ${aligned_reads}/SRR062634_sorted_dedup_reads.bam
```

### Step 4: Base Quality Recalibration

I performed base quality recalibration using GATK BaseRecalibrator and ApplyBQSR.

```bash
# Build the model
gatk BaseRecalibrator -I ${aligned_reads}/SRR062634_sorted_dedup_reads.bam -R ${ref} --known-sites ${known_sites} -O ${data}/recal_data.table

# Apply the model
gatk ApplyBQSR -I ${aligned_reads}/SRR062634_sorted_dedup_reads.bam -R ${ref} --bqsr-recal-file ${data}/recal_data.table -O ${aligned_reads}/SRR062634_sorted_dedup_bqsr_reads.bam
```

### Step 5: Collect Alignment & Insert Size Metrics

I collected alignment summary metrics and insert size metrics using GATK.

```bash
gatk CollectAlignmentSummaryMetrics R=${ref} I=${aligned_reads}/SRR062634_sorted_dedup_bqsr_reads.bam O=${aligned_reads}/alignment_metrics.txt
gatk CollectInsertSizeMetrics INPUT=${aligned_reads}/SRR062634_sorted_dedup_bqsr_reads.bam OUTPUT=${aligned_reads}/insert_size_metrics.txt HISTOGRAM_FILE=${aligned_reads}/insert_size_hist
```

### Step 6: Variant Calling

I called variants using GATK HaplotypeCaller, and then extracted SNPs and INDELs.

```bash
gatk HaplotypeCaller -R ${ref} -I ${aligned_reads}/SRR062634_sorted_dedup_bqsr_reads.bam -O ${results}/raw_variants.vcf --native-pair-hmm-threads 128

# Extract SNPs & INDELs
gatk SelectVariants -R ${ref} -V ${results}/raw_variants.vcf --select-type SNP -O ${results}/raw_snps.vcf
gatk SelectVariants -R ${ref} -V ${results}/raw_variants.vcf --select-type INDEL -O ${results}/raw_indels.vcf
```

### Variant Filtering

I applied filters to the variants to ensure high-quality results.

```bash
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
```

### Select Variants that PASS Filters

I selected variants that passed the filters.

```bash
gatk SelectVariants \
    --exclude-filtered \
    -V ${results}/filtered_snps.vcf \
    -O ${results}/analysis-ready-snps.vcf

gatk SelectVariants \
    --exclude-filtered \
    -V ${results}/filtered_indels.vcf \
    -O ${results}/analysis-ready-indels.vcf
```

### Variant Annotation

I used GATK Funcotator to annotate the variants.

```bash
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
```

### Extract Fields from VCF

I extracted fields from the VCF file to a tab-delimited table.

```bash
gatk VariantsToTable \
    -V ${results}/analysis-ready-snps-filteredGT-functotated.vcf -F AC -F AN -F DP -F AF -F FUNCOTATION \
    -O ${results}/output_snps.table
```

### Grepping for Gene of Interest "CROCCP2"

I filtered the annotations to focus on the gene of interest "CROCCP2".

```bash
cat output_snps.table | cut -f 5 | grep "CROCCP2" | sed 's/|/\t/g' >> output_curated_variants.txt
```

The resulting `output_curated_variants.txt` file can be viewed in Excel for further analysis.

## Conclusion

This workflow provides a comprehensive approach to variant calling and analysis, starting from raw sequencing data to annotated variants. The steps outlined ensure high-quality results that are ready for downstream analysis. For any questions or further assistance, feel free to reach out.
