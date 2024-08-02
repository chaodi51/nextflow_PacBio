process INDELQUAL {
    label "lofreq"
    label "himem"    
    
    input:
        path(bam)
        path(ref)

    output:
        path "${bam.simpleName}.indelqual.bam", emit: indelqual_bam
        path "${bam.simpleName}.indelqual.bam.bai", emit: indelqual_bai
        
    script:
    """ 
    lofreq indelqual --uniform 20,20 -f ${ref} -o ${bam.simpleName}.indelqual.bam ${bam}
    samtools index ${bam.simpleName}.indelqual.bam
    """
}

