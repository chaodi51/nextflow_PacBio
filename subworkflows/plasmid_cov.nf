//
// information needed to make plasmid backbone coverage plot
//

process READ_NAMES {
    label "ubuntu"

    input:
        val(plasmid_name)
        path(long_align)

    output:
        path("*.txt")

    script:
    """
        grep $plasmid_name $long_align | cut -f 1 > readnames.txt
    """
}

process PLASMID_READS {
    label "pysam"

    input:
        path(readnames)
        path(inputbam)

    output:
        path("*.sam")

    script:
    """
       get_plasmid_reads.py --readnames $readnames --inputbam $inputbam --outputreads plasmid.sam
    """
}

process PLASMID_BAM {
    label "samtools"

    input:
        path(oldbam)
        path(samreads)

    output:
        path("plasmid.bam")

    script:
    """ 
        samtools view -H $oldbam > header.txt
        cat header.txt $samreads >> unsorted.sam
        samtools view -h unsorted.sam | samtools sort -o plasmid.bam
    """
}

process PCOVERAGE {
    label "samtools"
    label "hugemem"

    input:
        path(bam)

    output:
        path("plasmid_depth.tsv")

    script:
    """
        samtools depth -aa $bam > plasmid_depth.tsv
    """
}

workflow PLASMID_COV {
    take:
        longest_aligns
        plasmid_bam
        plasmid_name

    main:
        ch_readnames = READ_NAMES ( plasmid_name, longest_aligns )
        ch_plasmidreads = PLASMID_READS ( ch_readnames, plasmid_bam )
        ch_plasmidbam = PLASMID_BAM (plasmid_bam, ch_plasmidreads)
        ch_plasmidcov = PCOVERAGE ( ch_plasmidbam )
        
    emit:
        plasmid_cov = ch_plasmidcov
}
