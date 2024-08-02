process FASTQC {
    label "ncores"
    
    input:
        path(sample)

    output:
        path "fastqc/*.html", emit: fastqc_html
        path "fastqc/*.zip", emit: fastqc_zip

    script:
    """
        mkdir -p fastqc
        fastqc -t ${task.cpus} -q -o fastqc ${sample}
    """
}
