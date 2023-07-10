#! /bin/bash 

# Assign file names 

deduplicated_bam_fles="deduplicated_files/B9569-1_S1_L001.sorted.sorted.RG.dedup.bam"
REF=MT.fasta
indexed_fasta=MT.fasta.fai
snp_file="merged_vcf/B9569-1_S1_L001_merged.vcf" 
basename="B9569-1_S1_L001"


#echo ${base_dir} 
#echo ${deduplicated_bam_fles}
#echo ${Fasta_file}
#echo ${indexed_fasta}
#echo ${snp_file}


# make results directory
#mkdir -p results
#results="${base_dir}/results"




#echo "Base quality recalibration"

# 1. build the model
#mkdir -p recalibration_tables

gatk BaseRecalibrator -I ${deduplicated_bam_fles} -R ${REF} --known-sites ${snp_file} -O ${recalibration_tables}/${basename}.recal_data.table



# 2. Apply the model to adjust the base quality scores

#mkdir -p deduplicated_bqsr

gatk ApplyBQSR -I ${deduplicated_bam_fles} -R ${REF} --bqsr-recal-file ${recalibration_tables}/${basename}.recal_data.table -O deduplicated_bqsr/${basename}_bqsr_reads.bam 


# assessing alignment and insert size metrics

# echo "Collecting alignment and insert size metrics"

mkdir -p   alignment_metrics  

gatk CollectAlignmentSummaryMetrics R=${REF} I=deduplicated_bqsr/${basename}_bqsr_reads.bam  O=alignment_metrics/${basename}_alignment_metrics.txt

gatk CollectInsertSizeMetrics INPUT=deduplicated_bqsr/${basename}_bqsr_reads.bam OUTPUT=alignment_metrics/${basename}_insert_size_metrics.txt HISTOGRAM_FILE=alignment_metrics/${basename}_insert_size_histogram.pdf


# 3. call variants using haplotype caller 
#mkdir -p haplotype_caller

gatk HaplotypeCaller -R ${REF} -I deduplicated_bqsr/${basename}_bqsr_reads.bam  -O haplotype_caller/${basename}_raw_variants.vcf 

# 4. Extract SNPs and INDELs
#mkdir -p snp_files indel_files

gatk SelectVariants -R ${REF} -V haplotype_caller/${basename}_raw_variants.vcf --select-type SNP -O snp_files/${basename}_raw_snps.vcf
gatk SelectVariants -R ${REF} -V haplotype_caller/${basename}_raw_variants.vcf --select-type INDEL -O indel_files/${basename}_raw_indels.vcf

# 5. Filter variants 
#mkdir -p filtered_snps filtered_indels

gatk VariantFiltration --variant snp_files/${basename}_raw_snps.vcf --filter-expression "QD < 2.0" --filter-name "QD2" --filter-expression "QUAL < 30.0" --filter-name "QUAL30" --filter-expression "SOR > 3.0" --filter-name "SOR3" --filter-expression "FS > 60.0" --filter-name "FS60" \
--filter-expression "MQ < 40.0" --filter-name "MQ40" \
--filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
--filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
--output filtered_snps/${basename}_SNP.filtered.vcf

gatk VariantFiltration \
--variant indel_files/${basename}_raw_indels.vcf \
--filter-expression "QD < 2.0" --filter-name "QD2" \
--filter-expression "QUAL < 30.0" --filter-name "QUAL30" \
--filter-expression "FS > 200.0" --filter-name "FS200" \
--filter-expression "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
--output filtered_indels/${basename}_INDEL.filtered.vcf


# 6.  select variants that passed the filters 
# mkdir snp_analysis_ready indels_analysis_ready

 gatk SelectVariants --exclude-filtered  -V filtered_snps/${basename}_SNP.filtered.vcf -O snp_analysis_ready/${basename}_analysis-ready-snps.vcf


 gatk SelectVariants --exclude-filtered -V filtered_indels/${basename}_INDEL.filtered.vcf -O indels_analysis_ready/${basename}_analysis-ready-indels.vcf

#7.  exclude variants that failed genotype filters
#mkdir -p passed_genotype_snps passed_genotype_indels

cat snp_analysis_ready/${basename}_analysis-ready-snps.vcf |grep -v -E "DP_filter|GQ_filter" > passed_genotype_snps/${basename}_analysis-ready-snps-filteredGT.vcf
cat indels_analysis_ready/${basename}_analysis-ready-indels.vcf|  grep -v -E "DP_filter|GQ_filter" > passed_genotype_indels/${basename}_analysis-ready-indels-filteredGT.vcf

# 8. merge snps and indel files
#mkdir -p merged_vcf

gatk MergeVcfs --INPUT passed_genotype_snps/${basename}_analysis-ready-snps-filteredGT.vcf --INPUT  passed_genotype_indels/${basename}_analysis-ready-indels-filteredGT.vcf --OUTPUT merged_vcf/${basename}_merged.vcf
