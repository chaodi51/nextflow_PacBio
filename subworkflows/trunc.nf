//
// plot read starts and ends as potential truncation hotspots
//


// import modules
include { ST_END_COV } from '../modules/st_end_cov.nf'
//include { ANNO_REGION } from '../modules/anno_region.nf'

workflow TRUNC {
    take:
        pri_bam // channel: ch_align.pri_bam 

    main:
        ST_END_COV ( pri_bam )
        ch_st_end_cov = ST_END_COV.out.cov
    //    ch_bed_cov = ch_st_end_cov
    //        .map { tuple("${params.refdir}${it.simpleName.split('__')[1]}.bed", it) }

    //    ANNO_REGION( ch_bed_cov )
    //    ch_annotated_region = ANNO_REGION.out.anno_cov.collect()

    emit:
        st_end_cov = ch_st_end_cov
    //    bed_cov = ch_bed_cov
    //    annotated_region = ch_annotated_region
}
