//
// pull out read length and alignment length, 
// do lambda DNA normlization if needed
//


// import modules
include { READLEN } from '../modules/readlen.nf'
include { ALIGNLEN } from '../modules/alignlen.nf'

workflow LEN {
    take:
        ch_pri_bam // channel: primary bam files

    main:
        READLEN( ch_pri_bam )
        ch_pos_len = READLEN.out.pos_len_info
        ch_read_len = READLEN.out.read_len_count

        ALIGNLEN( ch_pri_bam )
        ch_align_len = ALIGNLEN.out.align_len_count

    emit:
        pos_len = ch_pos_len
        read_len = ch_read_len
        align_len = ch_align_len
    }
