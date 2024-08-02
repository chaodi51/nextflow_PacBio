process MULTIQC {
    label "ncores"
    
    input:
        path(fastqc_html)
        path(fastqc_zip)

    output:
        path "*.html", emit: multiqci_html
        path "multiqc_data/*.json", emit: multiqc_json

    script:
    """
        multiqc . --filename multiqc.html
    """
}
