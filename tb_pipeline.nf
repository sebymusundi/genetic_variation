#! /usr/bin/env nextflow 

params.reads="/home/bioinfo/bioinfo/Seby_2/raw_data_pool_1/*L001_{R1,R2}_001.fastq.gz"
params.reference="/home/bioinfo/bioinfo/Seby_2/raw_data/Fastq/MT.fasta"
params.outdir= "/home/bioinfo/bioinfo/Seby_2/results/"
params.human_reference="/home/bioinfo/bioinfo/Seby_2/raw_data/Fastq/GRCh38_latest_genomic.fna.gz"


println "reads: $params.reads"
println "reference: $params.reference"
println "outdir: $params.outdir"
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
                container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'
                publishDir "${params.outdir}", mode: 'copy'

                input:
                tuple val(sample_id), path(my_reads)

                output:
                tuple val(sample_id), path("*")

                script:
                """
                mkdir -p quality_control
               fastqc ${my_reads[0]} ${my_reads[1]} -o quality_control

                """

}

// Removal of adapters using fastp

process fastp {
                    container 'quay.io/biocontainers/fastp:0.23.3--h5f740d0_0'
                    publishDir "${params.outdir}", mode: 'link'

                    input:
                    tuple val(sample_id), path(my_reads_1)

                    output:
                    tuple val(sample_id),  path("fastp_output/${my_reads_1[0].baseName}.trimmed.fastq"), path("fastp_output/${my_reads_1[1].baseName}.trimmed.fastq")
                    
                  

                    script:
                    """
                    mkdir -p fastp_output 
                    
                    fastp -i ${my_reads_1[0]} -I ${my_reads_1[1]} --adapter_sequence=CTGTCTCTTATACACATCT -Q 20 \
                    -o fastp_output/${my_reads_1[0].baseName}.trimmed.fastq -O fastp_output/${my_reads_1[1].baseName}.trimmed.fastq 
                
                    """
}

//  Trim reads with a quality score less than 20 using trimmomatic

process trimmomatic {
                    conda '/home/bioinfo/anaconda3'
                    publishDir "${params.outdir}", mode: 'link'

                    input:
                    tuple val(x), path(my_reads_1)

                    output:
                   
                    tuple val(x), path("trimmomatic_output/${my_reads_1[0].baseName}.paired.fastq"), 
                    path("trimmomatic_output/${my_reads_1[0].baseName}.unpaired.fastq"),
                    path("trimmomatic_output/${my_reads_1[1].baseName}.paired.fastq"),
                    path("trimmomatic_output/${my_reads_1[1].baseName}.unpaired.fastq")
                  

                    script:
                    """
                    mkdir -p trimmomatic_output 
                    
                    trimmomatic PE ${my_reads_1[0]} ${my_reads_1[1]}  trimmomatic_output/${my_reads_1[0].baseName}.paired.fastq trimmomatic_output/${my_reads_1[0].baseName}.unpaired.fastq \
                    trimmomatic_output/${my_reads_1[1].baseName}.paired.fastq trimmomatic_output/${my_reads_1[1].baseName}.unpaired.fastq SLIDINGWINDOW:4:20 MINLEN:20
                    """
}

// Checking the quality of trimmed files 

process fastqc_trimmed{
                        container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'
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
                        mkdir -p quality_control_trimmed
                        fastqc ${trimmed_reads_1}  ${trimmed_reads_3}    -o quality_control_trimmed

                        """
}



// Remove host reads

// Creating an index for the reference genome 

process bwa_human {
                        conda '/home/bioinfo/anaconda3'
                        publishDir "${params.outdir}", mode: 'link'
                        memory '28 GB'

                        input:
                        path reference

                        output:
                        tuple val("GRCh38_latest_genomic.fna.gz"), path("*")
                
                        script:
                        """
                        bwa index ${reference} 
                        """
}

// map the indexed genome against the reads 

process bwa_align_human {
                        conda '/home/bioinfo/anaconda3'
                        publishDir "${params.outdir}", mode: 'copy'
                         memory '28 GB'

                        input:
		        tuple val(ref_id), path(indexes), val(x), path(reads_1), path(reads_2), path(reads_3), path(read_4)
                         
                         
                        output: 
                        
                        tuple val(x), path("bam_files/${x}.bam")
                    
                        
                        script:
                        """
                        mkdir -p bam_files indexed_files
                        mv ${indexes}  indexed_files/
                       
                        bwa mem indexed_files/${ref_id} ${reads_1} ${reads_3} | samtools view -bS - > bam_files/${x}.bam
                       
                        """
}

// sort and index mapped reads to human genome 
process samtools_human {
                    conda '/home/bioinfo/anaconda3'
                     publishDir "${params.outdir}" , mode: 'link'
                    

                     input:
                     tuple val(x),path(sam_file)

                     output:
                     tuple val(x), path("bam_files_sorted/${sam_file.baseName}.sorted.bam"), path("bam_files_sorted/${sam_file.baseName}.sorted.bam.bai")

                     script:
                     """
                    mkdir -p bam_files_sorted 
                    samtools sort ${sam_file} -o  bam_files_sorted/${sam_file.baseName}.sorted.bam -O BAM
                    samtools index  bam_files_sorted/${sam_file.baseName}.sorted.bam
                    
                     """
}


// filter out unmapped reads 

// filter out unmapped reads 
process filter_unmapped_reads {
                    conda '/home/bioinfo/anaconda3'
                     publishDir "${params.outdir}" , mode: 'link'
                    

                     input:
                     tuple val(x),path(bam_file), path(bam_file_index)

                     output:
                     tuple val(x), path("bam_files_unmapped/${bam_file.baseName}.unmapped.bam"), path("bam_files_unmapped/${bam_file.baseName}.unmapped.sorted.bam")

                     script:
                     """
                     mkdir -p bam_files_unmapped
                     samtools view -b -f 4 ${bam_file} -o bam_files_unmapped/${bam_file.baseName}.unmapped.bam
                     samtools sort bam_files_unmapped/${bam_file.baseName}.unmapped.bam -o bam_files_unmapped/${bam_file.baseName}.unmapped.sorted.bam -O BAM
                    
                     """
}


// split files to fastq 

process split_file  {
                     conda '/home/bioinfo/anaconda3'
                     publishDir "${params.outdir}" , mode: 'copy'
                    

                     input:
                     tuple val(x),path(bam_file_1), path(bam_file_2)

                     output:
                     tuple val(x), path("bam_files_splitted/${x}_1.fastq"), path("bam_files_splitted/${x}_2.fastq")

                     script:
                     """
                     mkdir -p bam_files_splitted
                     samtools fastq ${bam_file_2} -1  bam_files_splitted/${x}_1.fastq -2 bam_files_splitted/${x}_2.fastq
                           
                     """
}


// Mapping against the MTB reference genome 

// Creating an index for the reference genome 

process bwa_index {
                        conda '/home/bioinfo/anaconda3'
                        publishDir "${params.outdir}", mode: 'link'

                        input:
                        path reference

                        output:
                        tuple val("MT.fasta"), path("*")

                
                        script:
                        """
                        bwa index ${reference} 
                        """
} 

// map the indexed genome against the reads 

process bwa_align {
                        conda '/home/bioinfo/anaconda3'
                        publishDir "${params.outdir}", mode: 'copy'

                        input:
			tuple val(ref_id), path(indexes), val(x), path(reads_1), path(reads_2)
                         
                        output: 
                        
                        tuple val(x), path("bam_file_output/${x}.bam")
                    
                        
                        script:
                        """
                        mkdir -p indexed_files bam_file_output
                        mv ${indexes} indexed_files/
                        bwa mem indexed_files/${ref_id} ${reads_1} ${reads_2} | samtools view -bS - > bam_file_output/${x}.bam
                        """
}



// sort and index mapped reads to TB genome
process samtools_TB {
                    conda '/home/bioinfo/anaconda3'
                     publishDir "${params.outdir}" , mode: 'link'
                    

                     input:
                     tuple val(x),path(sam_file)

                     output:
                     tuple val(x), path("bam_files_sorted/${sam_file.baseName}.sorted.bam"), path("bam_files_sorted/${sam_file.baseName}.sorted.bam.bai")

                     script:
                     """
                    mkdir -p bam_files_sorted 
                    samtools sort ${sam_file} -o  bam_files_sorted/${sam_file.baseName}.sorted.bam -O BAM
                    samtools index  bam_files_sorted/${sam_file.baseName}.sorted.bam
                    
                     """
}


// Checking the number of reads which mapped to the SARS-CoV-2 reference genome

process samtools_flagstat {
                                        conda '/home/bioinfo/anaconda3'
                                        publishDir "${params.outdir}" , mode: 'copy'
                                        

                                        input:
                                        tuple val(x), path(bam_file), path(bam_file_1)

                                        output:
                                        path "mapped_stats/${bam_file.baseName}.stats"

                                        script:
                                        """
                                        mkdir -p  mapped_stats
                                        samtools flagstat ${bam_file} > mapped_stats/${bam_file.baseName}.stats
                                        """
}




process samtools_index {
                                        conda '/home/bioinfo/anaconda3'
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
                                        conda '/home/bioinfo/anaconda3'
                                        publishDir "${params.outdir}" , mode: 'link'
                                        

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
                                        conda '/home/bioinfo/anaconda3'
                                        publishDir "${params.outdir}" , mode: 'link'
                                        

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
                                        conda '/home/bioinfo/anaconda3'
                                        publishDir "${params.outdir}" , mode: 'link'
                                        

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
                                        conda '/home/bioinfo/anaconda3'
                                        publishDir "${params.outdir}" , mode: 'copy'
                                        

                                        input: 
                                        tuple val(x), path(bam_file), path(metric_file),  path(reference)
                                       

                                        
                                        output:
                                        tuple val(x), path("bcftools_variant/${x}.vcf")

                                        script:
                                        """
                                        mkdir -p bcftools_variant 
                                        bcftools  mpileup -Ou -f ${reference} ${bam_file}|bcftools call -Ov -c  --ploidy 1 -o bcftools_variant/${x}.vcf
                                        """

}

process filter_bcf{
                                       conda '/home/bioinfo/anaconda3'
                                        publishDir "${params.outdir}" , mode: 'copy'
                                        

                                        input: 
                                        tuple val(x), path(vcf_file)
                                       
                                        output:
                                        tuple val(x), path("bcftools_filter/${x}.vcf.gz"), path("bcftools_filter/${x}.vcf.gz.csi")

                                        script:
                                        """
                                        mkdir -p bcftools_filter
                                        bcftools filter -Oz  ${vcf_file} -i 'TYPE="snp" && MIN(DP)>5 && QUAL>20' -o bcftools_filter/${x}.vcf.gz
                                        bcftools index bcftools_filter/${x}.vcf.gz

                                        """

}

process filter_bcf_indels{
                                        conda '/home/bioinfo/anaconda3'
                                        publishDir "${params.outdir}" , mode: 'copy'


                                        input:
                                        tuple val(x), path(vcf_file)

                                        output:
                                        tuple val(x), path("bcftools_filter_indels/${x}.vcf")

                                        script:
                                        """
                                        mkdir -p bcftools_filter_indels
                                        bcftools filter -O v  ${vcf_file} -i 'TYPE="indels" && MIN(DP)>5 && QUAL>20' -o bcftools_filter_indels/${x}.vcf

                                        """

}

// identify mutations against multiple different drugs used to treat TB 


process drug_snps {
                                        conda '/home/bioinfo/anaconda3'
                                        publishDir "${params.outdir}" , mode: 'copy'
                                        

                                        input: 
                                        tuple val(x), path(variant_file), path(variant_index)
                                       

                                        
                                        output:
                                        tuple val(x), 
                                        path("rpoB/${x}.vcf.gz"), 
                                        path("katG/${x}.vcf.gz"), 
                                        path("embB/${x}.vcf.gz"), 
                                        path("gyrA/${x}.vcf.gz"), 
                                        path("gyrB/${x}.vcf.gz"), 
                                        path("rrs/${x}.vcf.gz"), 
                                        path("fgd1/${x}.vcf.gz"), 
                                        path("fbiC/${x}.vcf.gz"), 
                                        path("pncA/${x}.vcf.gz")

                                        script:
                                        """
                                        mkdir -p rpoB katG embB gyrA gyrB rrs fgd1 fbiC pncA
                                        bcftools  view -Oz -r 'gi|448814763|ref|NC_000962.3|:759807-763325'  -o rpoB/${x}.vcf.gz ${variant_file}
                                        bcftools  view -Oz -r  'gi|448814763|ref|NC_000962.3|:2153889-2156111' -o katG/${x}.vcf.gz   ${variant_file}
                                        bcftools  view -Oz -r  'gi|448814763|ref|NC_000962.3|:4246514-4249810' -o embB/${x}.vcf.gz   ${variant_file}
                                        bcftools  view -Oz -r  'gi|448814763|ref|NC_000962.3|:7302-9818' -o gyrA/${x}.vcf.gz   ${variant_file}
                                        bcftools  view -Oz -r  'gi|448814763|ref|NC_000962.3|:5240-7267' -o gyrB/${x}.vcf.gz   ${variant_file}
                                        bcftools  view -Oz -r  'gi|448814763|ref|NC_000962.3|:1471846-1473382' -o rrs/${x}.vcf.gz  ${variant_file}
                                        bcftools  view -Oz -r  'gi|448814763|ref|NC_000962.3|:490783-491793' -o fgd1/${x}.vcf.gz  ${variant_file}
                                        bcftools  view -Oz -r  'gi|448814763|ref|NC_000962.3|:1302931-1305501' -o fbiC/${x}.vcf.gz  ${variant_file}
                                        bcftools  view -Oz -r  'gi|448814763|ref|NC_000962.3|:2288681-2289241' -o pncA/${x}.vcf.gz     ${variant_file}                            
                                        """

}




workflow {
my_reads=Channel.fromFilePairs("$params.reads")
my_reference=Channel.fromPath("$params.reference")
human_reads=Channel.fromPath("$params.human_reference")
//my_reads.view()
fastqc(my_reads)
//fastp_ch=fastp(my_reads)
//fastp_ch.view()
trimmomatic_ch=trimmomatic(my_reads)
//trimmomatic_ch.view()
//fastqc_trimmed(trimmomatic_ch)
bwa_human_ch=bwa_human(human_reads)
bwa_align_human_ch=bwa_human_ch.combine(trimmomatic_ch)
//bwa_align_human_ch.view()
bwa_human_sam=bwa_align_human(bwa_align_human_ch)
samtools_human_ch=samtools_human(bwa_human_sam)
unmapped_bam_ch =filter_unmapped_reads(samtools_human_ch)
split_file_ch=split_file(unmapped_bam_ch)
//split_file_ch.view()
bwa_ch=bwa_index(my_reference)
//bwa_ch.view()
bwa_combined_ch=bwa_ch.combine(split_file_ch)
//bwa_combined_ch.view()
bwa_aligned_ch=bwa_align(bwa_combined_ch)
//bwa_aligned_ch.view()
samtools_ch=samtools_TB(bwa_aligned_ch)
//samtools_ch.view()
samtools_flagstat(samtools_ch)
samtools_index(my_reference)
picard_sort_ch=picard_sort(samtools_ch)
picard_RG_ch=picard_read_groups(picard_sort_ch)
picard_duplicates_ch=picard_mark_duplicates(picard_RG_ch)
//lofreq_indels_ch=lofreq_indels(picard_duplicates_ch,my_reference)
//lofreq_indels_ch.view()
//lofreq_variant_ch=lofreq_call_variant(lofreq_indels_ch, my_reference)
//lofreq_filter_variant(lofreq_variant_ch)
picard_preprocess_ch=picard_duplicates_ch.combine(my_reference)
///picard_preprocess_ch.view()
bcf_variant_ch=variant_calling_bcf(picard_preprocess_ch)
//bcf_variant_ch.view()
filter_bcf_ch=filter_bcf(bcf_variant_ch)
//filter_bcf_indels(bcf_variant_ch)
drug_snps(filter_bcf_ch)
}
