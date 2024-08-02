process ST_END_COV {
    publishDir params.outdir, mode:'copy'
    publishDir params.datalakedir, mode:'copy'
    label "samtools"

    input:
        path bam // channel: pri_bam

    output:
        path "truncations.tsv", emit: cov

    script:
    """
    samtools view ${bam} | cut -f4,9 | awk '{print \$1"\\n"\$1+\$2-1}' | \
    sort | uniq -c | awk '{OFS=\"\\t\"} {print \"${params.unique_id}\", \$2, \$1}' | sort -nk2 | sed '1iworkflow_id\tposition\tcounts' > truncations.txt
    fill_missing.py --input truncations.txt --output truncations.tsv
    """
}
