#!/usr/bin/env nextflow

include {DownloadSRA} from "./modules/DownloadSRA.nf"
include {FastQC} from "./modules/fastqc.nf"
include {Trim_reads} from "./modules/Trim_reads.nf"
include {Alignment} from "./modules/alignment.nf"
include {FeatureCount} from "./modules/FeatureCount.nf"

// // Channel to collect input files and group them by sample

params.sra_id = "sra_id.csv"
nextflow.enable.dsl=2

workflow {
    def sra_ch = Channel
        .fromPath('/data/RNA_seq/sra_id.csv')
        .splitCsv(header: true)
        .map { row ->
            println "Loaded SRA ID: ${row.sra_id}"  // Logging
            row.sra_id                             // Correct usage
        }
    sra_ch.view()
    download_sra_ch=DownloadSRA(sra_ch)
    // fastqc_ch =  Channel.fromFilePairs(download_sra_ch, flat: true).map().view()
    // fastqc_ch = Channel.fromFilePairs(download_sra_ch).view()
    // Trim_ch= Channel.fromFilePairs(download_sra_ch).view()
    // align_ch=Channel.fromFilePairs(Trim_ch).view()
    // bam_ch = Channel.fromPath("${params.bam}/*.bam").collect()
    fastqc_ch = FastQC(downloaded_ch)
    trimmed_ch   = Trim_reads(downloaded_ch)
    bam_ch=Alignment(trimmed_ch)
    bam_list=bam_ch.collect()
    gff_ch= Channel.value(file(params.gff))
    FeatureCount(bam_list,gff_ch)
}
    
