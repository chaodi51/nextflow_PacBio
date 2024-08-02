//
// make ref indexes
//


// import modules
include { PBMM2_IDX } from '../modules/pbmm2_idx.nf'
include { FA_IDX } from '../modules/fa_idx.nf'

workflow REF {
    take:
        ch_ref_fa

    main:
        PBMM2_IDX( ch_ref_fa )
        ch_ref_mmi = PBMM2_IDX.out.mmi

        FA_IDX( ch_ref_fa )
        ch_ref_fai = FA_IDX.out.fai 

    emit:
        ref_mmi = ch_ref_mmi
        ref_fai = ch_ref_fai
    }
