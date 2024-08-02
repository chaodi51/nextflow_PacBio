process PLOT_COV {
    publishDir "${params.outdir}/coverage/${params.run}", mode:'copy'
    label "tidyverse"
    
    input:
        path cov // channel: collection of '.cov'
        path rmd

    output:
        path "${params.run}.gcov.html", emit: gcov_html
        
    script:
    """ 
        cp -L ${rmd} plot.rmd
        Rscript -e 'rmarkdown::render("plot.rmd", output_format = "html_document", \
        params = list(covfiles = "${cov}"), output_file = "${params.run}.gcov.html")'
    """
}
