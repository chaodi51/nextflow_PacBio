process PLOT_TRUNC {
    publishDir "${params.outdir}/trunc/${params.run}", mode:'copy'
    label "tidyverse"
    
    input:
        path cov // channel: collection of '.anno.cov'
        path rmd

    output:
        path "${params.run}.trunc_hotspot.html", emit: trunc_html
        
    script:
    """ 
        cp -L ${rmd} plot.rmd
        Rscript -e 'rmarkdown::render("plot.rmd", output_format = "html_document", \
        params = list(files = "${cov}"), output_file = "${params.run}.trunc_hotspot.html")'
    """
}
