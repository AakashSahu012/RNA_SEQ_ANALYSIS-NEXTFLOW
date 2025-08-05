#!/usr/bin/env nextflow

process DownloadSRA {
    tag "$sra_id"
    publishDir "/data/results", mode: 'copy'

    input:
    val sra_id

    output:
    tuple val(sra_id), path("${sra_id}_1.fastq.gz"), path("${sra_id}_2.fastq.gz")

    script:
    """
    fasterq-dump $sra_id --split-files --threads 15 --temp /data/tmp_sra
    gzip ${sra_id}_*.fastq
    rm -f ${sra_id}_1.fastq ${sra_id}_2.fastq


    """
}

