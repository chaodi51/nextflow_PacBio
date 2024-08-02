//
// fastqc and multiqc
//


// import modules
include { FASTQC } from '../modules/fastqc.nf'
include { MULTIQC } from '../modules/multiqc.nf'
include { TIDYMULTIQC } from '../modules/tidymultiqc.nf'

workflow QC {
    take:
        ch_sample // channel: path to sample 

    main:
        FASTQC( ch_sample )
        ch_fastqc_html = FASTQC.out.fastqc_html.collect()
        ch_fastqc_zip = FASTQC.out.fastqc_zip.collect()

        MULTIQC( ch_fastqc_html, ch_fastqc_zip )
        ch_multiqc = MULTIQC.out.multiqc_json

        //Add tidymultiqc
        ch_tidymultiqc = TIDYMULTIQC(ch_multiqc)

     emit:
        multiqc_tsv = ch_tidymultiqc

    }
