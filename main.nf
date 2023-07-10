#! /usr/bin/env nextflow 

params.reads="/home/bioinfo/bioinfo/Seby_2/raw_data/Fastq/B1195_S8_L001_{R1,R2}_001.fastq"
params.reference="/home/bioinfo/bioinfo/Seby_2/raw_data/Fastq/MT.fasta"
params.outdir= "/home/bioinfo/bioinfo/Seby_2/results/"
params.snp_file="/home/bioinfo/bioinfo/Seby_2/raw_data/Fastq/Homo_sapiens_assembly38.dbsnp138.vcf"
params.human_reference="/home/bioinfo/bioinfo/Seby_2/raw_data/Fastq/GCF_000001405.40_GRCh38.p14_genomic.fna.gz"


println "reads: $params.reads"
println "reference: $params.reference"
println "outdir: $params.outdir"
println "snp_file: $params.snp_file "
println "human_reference: $params.human_reference"


log.info """\
                        whole genome sequencing TB
                        ------------------------------------------
                        reads:       "$params.reads"
                        reference:   "$params.reference"
                        outdir:      "$params.outdir"
                        human_reference: "$params.human_reference"
                        ------------------------------------------
                        """
                .stripIndent()



// Quality control using fastqc

process fastqc {
                container 'biocontainers/fastqc:v0.11.9_cv8'
                publishDir "${params.outdir}", mode: 'copy'

                input:
                tuple val(sample_id), path(my_reads_1)

                output:
                tuple val(sample_id), 
                path("quality_control/${my_reads_1[0].baseName}_fastqc.html"), 
                path("quality_control/${my_reads_1[1].baseName}_fastqc.html"), 
                path("quality_control/${my_reads_1[0].baseName}_fastqc.zip"), 
                path("quality_control/${my_reads_1[1].baseName}_fastqc.zip")


                script:
                """
                mkdir  quality_control
                fastqc ${my_reads_1[0]} ${my_reads_1[1]}   -o quality_control

                """

}

// Removal of adapters using fastp

process fastp {
                    container 'nanozoo/fastp:latest'
                    publishDir "${params.outdir}", mode: 'copy'

                    input:
                    tuple val(sample_id), path(my_reads_1)

                    output:
                    tuple val(sample_id),  path("fastp_output/${my_reads_1[0].baseName}.trimmed.fastq"), path("fastp_output/${my_reads_1[1].baseName}.trimmed.fastq")
                    
                  

                    script:
                    """
                    mkdir -p fastp_output 
                    
                    fastp -i ${my_reads_1[0]} -I ${my_reads_1[1]} --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -Q 20 \
                    -o fastp_output/${my_reads_1[0].baseName}.trimmed.fastq -O fastp_output/${my_reads_1[1].baseName}.trimmed.fastq 
                
                    """
}

//  Trim reads with a quality score less than 20 using trimmomatic

process trimmomatic {
                    container 'staphb/trimmomatic:latest'
                    publishDir "${params.outdir}", mode: 'copy'

                    input:
                    tuple val(sample_id), path(my_reads_1), path(my_reads_2)

                    output:
                   
                    tuple val(sample_id), path("trimmomatic_output/${my_reads_1.baseName}.paired.fastq"), 
                    path("trimmomatic_output/${my_reads_1.baseName}.unpaired.fastq"),
                    path("trimmomatic_output/${my_reads_2.baseName}.paired.fastq"),
                    path("trimmomatic_output/${my_reads_2.baseName}.unpaired.fastq")
                  

                    script:
                    """
                    mkdir -p trimmomatic_output 
                    
                    trimmomatic-0.39.jar PE ${my_reads_1} ${my_reads_2}  trimmomatic_output/${my_reads_1.baseName}.paired.fastq trimmomatic_output/${my_reads_1.baseName}.unpaired.fastq \
                    trimmomatic_output/${my_reads_2.baseName}.paired.fastq trimmomatic_output/${my_reads_2.baseName}.unpaired.fastq SLIDINGWINDOW:4:20 
                    """
}

// Checking the quality of trimmed files 

process fastqc_trimmed{
                        container 'biocontainers/fastqc:v0.11.9_cv8'
                        publishDir "${params.outdir}", mode: 'copy'

                        input:
                        tuple val(sample_id), path(trimmed_reads_1), path(trimmed_reads_2), path(trimmed_reads_3), path(trimmed_reads_4)

                        output:
                        tuple val(sample_id), 
                        path("quality_control_trimmed/${trimmed_reads_1.baseName}_fastqc.html"), 
                        path("quality_control_trimmed/${trimmed_reads_3.baseName}_fastqc.html"), 
                        path("quality_control_trimmed/${trimmed_reads_1.baseName}_fastqc.zip"), 
                        path("quality_control_trimmed/${trimmed_reads_3.baseName}_fastqc.zip")


                        script:
                        """
                        mkdir  quality_control_trimmed
                        fastqc ${trimmed_reads_1}  ${trimmed_reads_3}    -o quality_control_trimmed

                        """
}



// Remove host reads 

// Creating an index for the reference genome 

process bowtie_human {
                        container 'biocontainers/bowtie2:v2.4.1_cv1'
                        publishDir "${params.outdir}", mode: 'copy'

                        input:
                        path reference

                        output:
                        tuple val("bowtie_human"), path("bowtie_human*")
                
                        script:
                        """
                        bowtie2-build ${reference}  bowtie_human
                        """
}

// map the indexed genome against the reads 

process bowtie_align_human {
                        container 'biocontainers/bowtie2:v2.4.1_cv1'
                        publishDir "${params.outdir}", mode: 'copy'

                        input:
                        tuple val(x), path(my_reads), val(index), path(align_input)
                    
                    
                        output:
                        tuple val(x), path("bowtie_output_human/${x}.sam")
                        
                        script:
                        """
                        mkdir bowtie_output_human
                        
                 
                        bowtie2 -x ${index} -1 ${my_reads[0]} -2 ${my_reads[1]} -S bowtie_output_human/${x}.sam                   
                        """
}

// Converting sam file to become Bam file followed by sorting and indexing the associated BAM file

process samtools_human {
                     container 'staphb/samtools:latest'
                     publishDir "${params.outdir}" , mode: 'copy'
                    

                     input:
                     tuple val(x),path(sam_file)

                     output:
                     tuple val(x), path("bam_files_human/${sam_file.baseName}.sorted.bam"), path("bam_files_human/${sam_file.baseName}.sorted.bam.bai")

                     script:
                     """
                     mkdir bam_files_human 
                    samtools view -O BAM ${sam_file} -o  bam_files_human/${sam_file.baseName}.bam
                    samtools sort  bam_files_human/${sam_file.baseName}.bam -o  bam_files_human/${sam_file.baseName}.sorted.bam -O BAM
                    samtools index  bam_files_human/${sam_file.baseName}.sorted.bam
                     """
}

// filter out unmapped reads 
process filter_unmapped_reads {
                     container 'staphb/samtools:latest'
                     publishDir "${params.outdir}" , mode: 'copy'
                    

                     input:
                     tuple val(x),path(bam_file), path(bam_file_index)

                     output:
                     tuple val(x), path("bam_files_unmapped/${bam_file.baseName}.unmapped.bam")

                     script:
                     """
                     mkdir bam_files_unmapped
                     samtools view -b -f 12 -F 256 ${bam_file} -o bam_files_unmapped/${bam_file.baseName}.unmapped.bam
                    
                     """
}

// split paired end reds into separated fastq files 

process split_file_sort  {
                     container 'staphb/samtools:latest'
                     publishDir "${params.outdir}" , mode: 'copy'
                    

                     input:
                     tuple val(x),path(bam_file)

                     output:
                     tuple val(x), path("bam_files_sort/${bam_file.baseName}.unmapped.sorted.bam")

                     script:
                     """
                     mkdir bam_files_sort
                     samtools sort ${bam_file} -o bam_files_sort/${bam_file.baseName}.unmapped.sorted.bam -O BAM
                           
                     """
}

// split files to fastq 

process split_file  {
                     container 'staphb/samtools:latest'
                     publishDir "${params.outdir}" , mode: 'copy'
                    

                     input:
                     tuple val(x),path(bam_file)

                     output:
                     tuple val(x), path("bam_files_split/${x}_R1_001.fastq"), path("bam_files_split/${x}_R2_001.fastq")

                     script:
                     """
                     mkdir bam_files_split
                     samtools fastq -1 bam_files_split/${x}_R1_001.fastq -2 bam_files_split/${x}_R2_001.fastq ${bam_file}
                           
                     """
}

// Mapping against the MTB reference genome 

// Creating an index for the reference genome 

process bwa_index {
                        container 'staphb/bwa:latest'
                        publishDir "${params.outdir}", mode: 'copy'

                        input:
                        path reference

                        output:
                        tuple val("MT.fasta"), 
                        path("${reference}.bwt"), 
                        path("${reference}.pac"), 
                        path("${reference}.ann"), 
                        path("${reference}.amb"), 
                        path("${reference}.sa")

                
                        script:
                        """
                        bwa index ${reference} 
                        """
} 

// map the indexed genome against the reads 

process bwa_align {
                        container 'staphb/bwa:latest'
                        publishDir "${params.outdir}", mode: 'copy'

                        input:
			 tuple val(x), path(reads_1), path(reads_2) 
                         tuple val(ref_id), path(index_bwt), path(index_pac), path(index_ann), path(index_amb), path(index_sa)
                         
                        output: 
                        
                        tuple val(x), path("sam_file_output/${x}.sam")
                    
                        
                        script:
                        """
                        mkdir -p indexed_files sam_file_output
                        mv ${index_bwt} ${index_pac} ${index_ann} ${index_amb} ${index_sa} indexed_files/
                        bwa mem indexed_files/${ref_id} ${reads_1} ${reads_2} > sam_file_output/${x}.sam 
                        """
}


// Converting sam file to become Bam file followed by sorting and indexing the associated BAM file

process samtools {
                     container 'staphb/samtools:latest'
                     publishDir "${params.outdir}" , mode: 'copy'
                    

                     input:
                     tuple val(x),path(sam_file)

                     output:
                     tuple val(x), path("bam_files/${sam_file.baseName}.sorted.bam"), path("bam_files/${sam_file.baseName}.sorted.bam.bai")

                     script:
                     """
                     mkdir bam_files 
                    samtools view -O BAM ${sam_file} -o  bam_files/${sam_file.baseName}.bam
                    samtools sort  bam_files/${sam_file.baseName}.bam -o  bam_files/${sam_file.baseName}.sorted.bam -O BAM
                    samtools index  bam_files/${sam_file.baseName}.sorted.bam
                     """
}

// Checking the number of reads which mapped to the SARS-CoV-2 reference genome

process samtools_flagstat {
                                        container 'staphb/samtools:latest'
                                        publishDir "${params.outdir}" , mode: 'copy'
                                        

                                        input:
                                        tuple val(x), path(bam_file), path(bam_file_1)

                                        output:
                                        path "mapped_stats/${bam_file.baseName}.stats"

                                        script:
                                        """
                                        mkdir mapped_stats
                                        samtools flagstat ${bam_file} > mapped_stats/${bam_file.baseName}.stats
                                        """
}




process samtools_index {
                                        container 'staphb/samtools:latest'
                                        publishDir "${params.outdir}" , mode: 'copy'
                                        

                                        input:
                                        path(reference)

                                        output:
                                        path ("${reference}.fai")

                                        script:
                                        """
                                        samtools faidx ${reference}
                                        """
}

process picard_sort {
                                        conda '/home/bioinfo/anaconda3/envs/bioinfo'
                                        publishDir "${params.outdir}" , mode: 'copy'
                                        

                                        input: 
                                        tuple val(x), path(bam_file), path(bam_file_1)

                                        output:
                                        tuple val(x), path("picard_sort/${bam_file.baseName}.sorted.bam")

                                        script:
                                        """
                                        mkdir picard_sort
                                        picard SortSam I=${bam_file} O=picard_sort/${bam_file.baseName}.sorted.bam SORT_ORDER=coordinate
                                        """

}

process picard_read_groups {
                                        conda '/home/bioinfo/anaconda3/envs/bioinfo'
                                        publishDir "${params.outdir}" , mode: 'copy'
                                        

                                        input: 
                                        tuple val(x), path(bam_file)

                                        output:
                                        tuple val(x), path("picard_RG/${bam_file.baseName}.RG.bam")

                                        script:
                                        """
                                        mkdir picard_RG
                                        picard AddOrReplaceReadGroups  I=${bam_file} O=picard_RG/${bam_file.baseName}.RG.bam RGLB=${x} RGPL=illumina RGPU=None RGSM=${x}
                                        """

}

process picard_mark_duplicates {
                                        conda '/home/bioinfo/anaconda3/envs/bioinfo'
                                        publishDir "${params.outdir}" , mode: 'copy'
                                        

                                        input: 
                                        tuple val(x), path(bam_file)

                                        output:
                                        tuple val(x), path("picard_duplicates/${bam_file.baseName}.dedup.bam"), path("picard_duplicates/${bam_file.baseName}.metrics.txt")

                                        script:
                                        """
                                        mkdir picard_duplicates
                                        picard MarkDuplicates  I=${bam_file} O=picard_duplicates/${bam_file.baseName}.dedup.bam METRICS_FILE=picard_duplicates/${bam_file.baseName}.metrics.txt REMOVE_DUPLICATES=TRUE \
                                        ASSUME_SORTED=TRUE MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000
                                        """

}

process lofreq_indels {
                                        container 'nanozoo/lofreq:2.1.5--229539a'
                                        publishDir "${params.outdir}" , mode: 'copy'
                                        

                                        input: 
                                        tuple val(x), path(bam_file), path(metric_file)
                                        path(reference)

                                        output:
                                        tuple val(x), path("lofreq_indels/${bam_file.baseName}.indels.bam")

                                        script:
                                        """
                                        mkdir lofreq_indels
                                        lofreq indelqual -f ${reference} --dindel  ${bam_file} -o lofreq_indels/${bam_file.baseName}.indels.bam
                                        """
}

process lofreq_call_variant {
                                        container 'nanozoo/lofreq:2.1.5--229539a'
                                        publishDir "${params.outdir}" , mode: 'copy'
                                        

                                        input: 
                                        tuple val(x), path(bam_file)
                                        path(reference)

                                        output:
                                        tuple val(x), path("vcf_file/${bam_file.baseName}.vcf")

                                        script:
                                        """
                                        mkdir vcf_file
                                        lofreq call -f ${reference} --call-indels -o vcf_file/${bam_file.baseName}.vcf ${bam_file}
                                        """
}

process lofreq_filter_variant {
                                        container 'nanozoo/lofreq:2.1.5--229539a'
                                        publishDir "${params.outdir}" , mode: 'copy'
                                        

                                        input: 
                                        tuple val(x), path(vcf_file)
                                        

                                        output:
                                        tuple val(x), path("lofreq_filter/${vcf_file.baseName}.vcf")

                                        script:
                                        """
                                        mkdir lofreq_filter
                                        lofreq filter -i ${vcf_file} -v 5 -a 0.3 -Q 20 -K 20  -o lofreq_filter/${vcf_file.baseName}.vcf
                                        """
}

process variant_calling_bcf{
                                        conda'/home/bioinfo/anaconda3/envs/bioinfo'
                                        publishDir "${params.outdir}" , mode: 'copy'
                                        

                                        input: 
                                        tuple val(x), path(bam_file), path(metric_file)
                                        path(reference)

                                        
                                        output:
                                        tuple val(x), path("bcftools_variant/${x}.vcf")

                                        script:
                                        """
                                        mkdir bcftools_variant 
                                        samtools  mpileup -u -f ${reference} ${bam_file}|bcftools call -v -c  --ploidy 1 > bcftools_variant/${x}.vcf
                                        """

}

process filter_bcf{
                                        conda'/home/bioinfo/anaconda3/envs/bioinfo'
                                        publishDir "${params.outdir}" , mode: 'copy'
                                        

                                        input: 
                                        tuple val(x), path(vcf_file)
                                       
                                        output:
                                        tuple val(x), path("bcftools_filter/${x}.vcf")

                                        script:
                                        """
                                        mkdir bcftools_filter
                                        bcftools filter -O v  ${vcf_file} -i 'TYPE="snp" && MIN(DP)>5 && QUAL>20' -o bcftools_filter/${x}.vcf

                                        """

}

process filter_bcf_indels{
                                        conda'/home/bioinfo/anaconda3/envs/bioinfo'
                                        publishDir "${params.outdir}" , mode: 'copy'


                                        input:
                                        tuple val(x), path(vcf_file)

                                        output:
                                        tuple val(x), path("bcftools_filter_indels/${x}.vcf")

                                        script:
                                        """
                                        mkdir bcftools_filter_indels
                                        bcftools filter -O v  ${vcf_file} -i 'TYPE="indels" && MIN(DP)>5 && QUAL>20' -o bcftools_filter_indels/${x}.vcf

                                        """

}


workflow {
my_reads=Channel.fromFilePairs("$params.reads")
my_reference=Channel.fromPath("$params.reference")
my_known_sites=Channel.fromPath("$params.snp_file")
human_reads=Channel.fromPath("$params.human_reference")
//my_reads.view()
fastqc(my_reads)
//fastp_ch=fastp(my_reads)
//fastp_ch.view()
//trimmomatic_ch=trimmomatic(fastp_ch)
//trimmomatic_ch.view()
//fastqc_trimmed(trimmomatic_ch)
bowtie_human_ch=bowtie_human(human_reads)
bowtie_align_human_ch=my_reads.combine(bowtie_human_ch)
//bowtie_align_human_ch.view()
bowtie_human_sam=bowtie_align_human(bowtie_align_human_ch)
samtools_human_ch=samtools_human(bowtie_human_sam)
unmapped_bam_ch =filter_unmapped_reads(samtools_human_ch)
split_file_ch=split_file_sort(unmapped_bam_ch)
split_bam_to_fastq=split_file(split_file_ch)
//split_bam_to_fastq.view()
bwa_ch=bwa_index(my_reference)
//bwa_ch.view()
//bwa_combined_ch=split_bam_to_fastq.combine(bwa_ch)
//bwa_combined_ch.view()
bwa_aligned_ch=bwa_align(split_bam_to_fastq, bwa_ch)
//bwa_aligned_ch.view()
samtools_ch=samtools(bwa_aligned_ch)
samtools_flagstat(samtools_ch)
samtools_index(my_reference)
picard_sort_ch=picard_sort(samtools_ch)
picard_RG_ch=picard_read_groups(picard_sort_ch)
picard_duplicates_ch=picard_mark_duplicates(picard_RG_ch)
lofreq_indels_ch=lofreq_indels(picard_duplicates_ch,my_reference)
lofreq_indels_ch.view()
lofreq_variant_ch=lofreq_call_variant(lofreq_indels_ch, my_reference)
lofreq_filter_variant(lofreq_variant_ch)
bcf_variant_ch=variant_calling_bcf(picard_duplicates_ch,my_reference)
//bcf_variant_ch.view()
filter_bcf(bcf_variant_ch)
filter_bcf_indels(bcf_variant_ch)
}
