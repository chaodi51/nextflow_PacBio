//
// fastq creation
//

process PBINDEX {
    label 'pacbio'

    input:
    path(bam)

    output:
    path("*.pbi")

    script:
    """
    pbindex $bam
    """
}

process BAM2FASTQ {
    label 'pacbio'

    input:
    path(bam)
    path(bamindex)

    output:
    path("*.fastq.gz")

    script:
    """
    bam2fastq -o main $bam
    """
}

workflow MAKE_FASTQ {
    take:
        unalignedbam

    main:
        ch_index = PBINDEX(unalignedbam)
        ch_fastq = BAM2FASTQ(unalignedbam, ch_index)

    emit:
        fastq = ch_fastq
    }
