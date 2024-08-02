process HOST_READS {
    label "samtools"

    input:
        path bamfile // channel: ch_hostbamfile

    output:
        path "${bamfile.simpleName}.mapped.tsv", emit: host_reads

    script:
    """
    samtools flagstats $bamfile | grep "primary" | head -1|cut -d " " -f1 | \
    awk '{print "${bamfile.simpleName}\t"\$1 }'> ${bamfile.simpleName}.mapped.tsv
    """
}
