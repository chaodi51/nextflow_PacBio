process ALIGNED_READS {
    label "samtools"

    input:
        path bamfile // channel: ch_bamfile

    output:
        stdout emit: aligned_reads

    script:
    """
    samtools flagstats $bamfile | grep "primary" | head -1 | cut -d " " -f1 | tr -d "\\n"
    """
}
