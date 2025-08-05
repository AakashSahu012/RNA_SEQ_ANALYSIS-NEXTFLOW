#!/usr/bin/env nextflow
process FeatureCount {

    publishDir "${params.outdir}", mode: 'copy'

    input:
    file bam_list
    file gff

    output:
    file "combined_counts.txt"

    script:
    """
    /tools/subread-2.1.1-Linux-x86_64/bin/featureCounts -a ${gff} \\
                  -o combined_counts.txt \\
                  -T 2 \\
                  -t exon \\
                  -g Parent \\
                  -s 0 \\
                  -p \\
                  ${bam_list.join(' ')}
    """
}
