process READLEN {
    publishDir params.outdir, mode:'copy'
    publishDir params.datalakedir, mode:'copy'
    label "samtools"

    input:
        path bamfile // channel: ch_pri_bam

    output:
        path "read_len_counts.tsv", emit: read_len_count
        path "pos_len.tsv", emit: pos_len_info

    script:
    """
        samtools view ${bamfile} | cut -f 4,10 | awk '{print \$1, length(\$2)}' > read_len.txt
        cut -d " " -f 2 read_len.txt | sort -n | uniq -c | awk '{OFS=\"\\t\"} {print \"${params.unique_id}\", \$2, \$1}' | sed '1iworkflow_id\tread_length\tcounts' > read_len_counts.txt
        fill_missing.py --input read_len_counts.txt --output read_len_counts.tsv
        awk '{OFS=\"\\t\"} {print \"${params.unique_id}\", \$1, \$2}' read_len.txt | sed '1iworkflow_id\tposition\tread_length' > pos_len.tsv
    """
}
