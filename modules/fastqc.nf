#!/usr/bin/env nextflow


process FastQC {
    tag "${sample_id}"

    publishDir "${params.outdir}/fastqc", mode: 'copy'

    input:
    tuple val(sample_id), path(paired)

    output:
    path "*.html"
    path "*.zip"

    script:
    """
    /tools/fastqc/fastqc -o . ${paired.join(' ')}
    """
}