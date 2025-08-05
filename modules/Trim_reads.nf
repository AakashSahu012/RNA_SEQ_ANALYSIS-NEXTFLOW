#!/usr/bin/env nextflow

process Trim_reads{
    tag "${sample_id}"
    publishDir "results/trimmed_reads", mode: 'copy'

    input:
    tuple val(sample_id) ,path (paired)
    
    output:
    tuple val(sample_id), path("${sample_id}_Trimmed_1.fastq.gz"), path("${sample_id}_Trimmed_2.fastq.gz")

    script:
    """
    java -jar /tools/trimmomatic/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 15 ${paired[0]} ${paired[1]} \\
    ${sample_id}_Trimmed_1.fastq.gz ${sample_id}_unpaired_1.fastq.gz \\
    ${sample_id}_Trimmed_2.fastq.gz ${sample_id}_unpaired_2.fastq.gz \\
    ILLUMINACLIP:/work/illuQia.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    
    """

}