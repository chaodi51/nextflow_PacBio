process ALIGNLEN {
    publishDir params.outdir, mode:'copy'
    publishDir params.datalakedir, mode:'copy'
    label "pysam"
    
    input:
        path bamfile // channel: ch_pri_bam

    output:
        path "align_len_count.tsv", emit: align_len_count

    script:
    """
        ## use 'rlen' in samtools
        len.py ${bamfile} ${bamfile.simpleName}.align_len.txt
        cat ${bamfile.simpleName}.align_len.txt | sort -n | uniq -c | awk '{OFS=\"\\t\"} {print \"${params.unique_id}\", \$2, \$1}' | sed '1iworkflow_id\talignment_length\tcounts' > align_len_count.tsv
    """
}
