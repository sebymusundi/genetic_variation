#! /usr/bin/env nextflow 

params.reads="/home/bioinfo/bioinfo/Seby_2/raw_data/Fastq/*_{R1,R2}_001.fastq"
params.reference="/home/bioinfo/bioinfo/Seby_2/raw_data/Fastq/MT.fasta"
params.outdir= "/home/bioinfo/bioinfo/Seby_2/results/"


println "reads: $params.reads"
println "reference: $params.reference"
println "outdir: $params.outdir"


log.info """\
                        whole genome sequencing TB
                        ------------------------------------------
                        reads:       "$params.reads"
                        reference:   "$params.reference"
                        outdir:      "$params.outdir"
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
                fastqc ${my_reads_1[0]} ${my_reads_1[1]}    -o quality_control

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

// Mapping against the MTB reference genome 

// Creating an index for the reference genome 

process bowtie_index {
                        container 'biocontainers/bowtie2:v2.4.1_cv1'
                        publishDir "${params.outdir}", mode: 'copy'

                        input:
                        path reference

                        output:
                        tuple val("bowtie_index"), path("bowtie_index*")
                
                        script:
                        """
                        bowtie2-build ${reference}  bowtie_index
                        """
}

// map the indexed genome against the reads 

process bowtie_align {
                        container 'biocontainers/bowtie2:v2.4.1_cv1'
                        publishDir "${params.outdir}", mode: 'copy'

                        input:
                        tuple val(x), path(my_reads), val(index), path(align_input)
                    
                    
                        output:
                        tuple val(x), path("bowtie_output/${x}.sam")
                        
                        script:
                        """
                        mkdir bowtie_output
                        
                 
                        bowtie2 -x ${index} -1 ${my_reads[0]} -2 ${my_reads[1]} -S bowtie_output/${x}.sam                   
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
                                        samtools faidx ${reference} -o ${reference}.fai
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


process samtools_index_2 {
                                        container 'staphb/samtools:latest'
                                        publishDir "${params.outdir}" , mode: 'copy'
                                        

                                        input:
                                        tuple val(x), path(bam_file), path(metrics_file)

                                        output:
                                        tuple val(x), path("index_duplicate/${bam_file.baseName}.index.bam.bai")

                                        script:
                                        """
                                        mkdir index_duplicate
                                        samtools index ${bam_file} -o index_duplicate/${bam_file.baseName}.index.bam.bai
                                        """
}


process GATK_indel_realign {
                                        container 'broadinstitute/gatk:latest'
                                        publishDir "${params.outdir}" , mode: 'copy'
                                        

                                        input:
                                        tuple val(x), path(bam_file), path(metrics_file)

                                        output:
                                        tuple val(x), path("index_duplicate/${bam_file.baseName}.index.bam.bai")

                                        script:
                                        """
                                        mkdir indel_realign
                                        gatk 
                                        """
}



workflow {
my_reads=Channel.fromFilePairs("$params.reads")
my_reference=Channel.fromPath("$params.reference")
//my_reads.view()
fastqc(my_reads)
//fastp_ch=fastp(my_reads)
//fastp_ch.view()
//trimmomatic_ch=trimmomatic(fastp_ch)
//trimmomatic_ch.view()
//fastqc_trimmed(trimmomatic_ch)
bowtie_ch=bowtie_index(my_reference)
//bowtie_ch.view()
bowtie_combined_ch=my_reads.combine(bowtie_ch)
//bowtie_combined_ch.view()
bowtie_aligned_ch=bowtie_align(bowtie_combined_ch)
bam_ch=samtools_ch=samtools(bowtie_aligned_ch)
samtools_flagstat(samtools_ch)
samtools_index(my_reference)
picard_sort_ch=picard_sort(bam_ch)
picard_RG_ch=picard_read_groups(picard_sort_ch)
picard_duplicates_ch=picard_mark_duplicates(picard_RG_ch)
indexed_2_ch=samtools_index_2(picard_duplicates_ch)
}