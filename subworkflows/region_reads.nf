//
// Count reads mapped to specific regions on the reference
//


// import modules
include { COVERAGEBED } from '../modules/coveragebed.nf'
include { COVTAB } from '../modules/covtab.nf'

workflow REGION_READS {
    take:
        bam_beds // channel: [val(sample), .bam, .regions.bed, .ITR_ITR.bed., border.bam]

    main:
        COVERAGEBED( bam_beds )
        ch_region_counts = COVERAGEBED.out.counts.collect()

        COVTAB(ch_region_counts)
        ch_cov_tab = COVTAB.out.cov_tab

    emit:
        region_counts = ch_region_counts
        cov_tab = ch_cov_tab
    }
