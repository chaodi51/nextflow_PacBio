//
// coverage plot for multipl samples on each ref 
//


// import modules
include { GENOME_COV } from '../modules/genome_cov.nf'

workflow COV {
    take:
        pri_bam // channel: ch_align.pri_bam 

    main:
        GENOME_COV ( pri_bam )
        ch_gcov = GENOME_COV.out.gcov.collect()

    emit:
        gcov = ch_gcov
}
