//
// mapping and bam processing
//


// import modules
include { PBMM2_ALIGN } from '../modules/pbmm2_align.nf'

process READ_COUNT {
    label "samtools"

    input:
        path(bam)

    output:
        path("*.tsv")

    script:
    """
    samtools view -F3332 $bam | aligned_length.py > ${bam.baseName}.counts.tsv
    """
}

process REF_COUNTS {
    label "pandas"
  
    input:
    path(countfiles)
    val(primary)

    output:
    path("ref_counts.tsv"), emit: counts
    path("merged.tsv"), emit: merged

    script:
    """
    longest_alignment.py --outfile ref_counts.tsv --mergedfile merged.tsv --primaryref $primary --countfile ${countfiles.join(" ")}
    """
}

workflow ALIGN {
    take:
        sample
        reference
        primary_name 

    main:
        PBMM2_ALIGN( sample, reference )
        ch_bam = PBMM2_ALIGN.out.bam
        ch_bai = PBMM2_ALIGN.out.bai
        ch_countfile = READ_COUNT(PBMM2_ALIGN.out.bam)
        REF_COUNTS(ch_countfile.collect(), primary_name)
        ch_refcounts = REF_COUNTS.out.counts
        ch_merged = REF_COUNTS.out.merged
        
    emit:
        bam = ch_bam
        bai = ch_bai
        count = ch_refcounts
        merged = ch_merged
    }
