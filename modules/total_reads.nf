process TOTAL_READS {
    label "ubuntu"
    
    input:
        path sample // channel: ch_sample

    output:
        stdout emit: total_reads

    script:
    """
    gunzip -c $sample | wc -l | tr -d "\\n"
    """
}
