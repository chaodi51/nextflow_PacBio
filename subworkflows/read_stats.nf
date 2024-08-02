//
// Count total and mapped reads
//


// import modules
include { TOTAL_READS } from '../modules/total_reads.nf'
include { ALIGNED_READS } from '../modules/aligned_reads.nf'
// include { HOST_READS } from '../modules/host_reads.nf'
// include { STATS_TAB } from '../modules/stats_tab.nf'

process READ_STATS_FILE {
  publishDir params.outdir, mode:'copy'
  publishDir params.datalakedir, mode:'copy'
  label "ubuntu"

  input:
  val(total_reads)
  val(aligned_reads)

  output:
  path("readstats.tsv")

  script:
  """
  printf \"%s\\t%s\\t%s\\n" workflow_id total_reads aligned_reads > readstats.tsv
  printf \"%s\\t%s\\t%s\\n" \"$params.unique_id\" \"$total_reads\" \"$aligned_reads\" >> readstats.tsv
  """
}

workflow READ_STATS {
    take:
        // arbitrary names
        ch_sample // channel: fastq files
        ch_bamfile // channel: pbmm2 mapped bam files

    main:
        TOTAL_READS( ch_sample )
        ch_total_reads = TOTAL_READS.out.total_reads

        ALIGNED_READS( ch_bamfile )
        ch_aligned_reads = ALIGNED_READS.out.aligned_reads

        readstats_ch = READ_STATS_FILE(ch_total_reads, ch_aligned_reads)
        
    emit:
        read_stats = readstats_ch
    }
