process CHIMERA_TRUNC {
    publishDir "${params.outdir}/chimera/alvis/${params.run}", mode:'copy'
    
    input:
        tuple path(chimeras), path(bed), path(rmd) // channel: {sample}__host_{aav}.chimeras.txt, rmd script

    output:
        path "${chimeras.simpleName}.trunc.html", emit: html
        stdout

    script:
    """ 
        echo ${chimeras}
        cp -L ${rmd} plot.rmd
        Rscript -e 'rmarkdown::render("plot.rmd", output_format = "html_document", \
        params = list(chimeras = "${chimeras}", annotation = "${bed}"), output_file = "${chimeras.simpleName}.trunc.html")'
    """
}