#!/usr/bin/env nextflow

process Alignment{
    tag "${sample_id}"
    publishDir "${params.outdir}/alignments", mode: 'copy'

    input:
    tuple val (sample_id), path (reads)
    output:
    path "${sample_id}.bam" 

    script:
    """
    
    /tools/hisat2-2.2.1/hisat2 -x ${params.index} -1 ${reads[0]} -2 ${reads[1]} --threads 4 | \
    samtools view -Sb - | \
    samtools sort -@ 4 -o ${sample_id}.bam -

    """
}