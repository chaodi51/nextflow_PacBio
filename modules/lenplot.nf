process LENPLOT {
    publishDir "${params.outdir}/length/${params.run}", mode:'copy'
    label "tidyverse"
    
    input:
        path align_len // channel: collection of align_len.txt
        path rmd

    output:
        path "${params.run}.len_hist.html", emit: hist
        
    script:
    """ 
        cp -L ${rmd} plot.rmd
        Rscript -e 'rmarkdown::render("plot.rmd", output_format = "html_document", \
        params = list(files = "${align_len}"), output_file = "${params.run}.len_hist.html")'
    """
}
