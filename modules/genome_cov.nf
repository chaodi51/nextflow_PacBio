process GENOME_COV {
    publishDir params.outdir, mode:'copy'
    publishDir params.datalakedir, mode:'copy'
    label "bedtools"

    input:
        path bam // channel: pri_bam

    output:
        path "${bam.simpleName}.cov", emit: gcov

    script:
    """
   
    bedtools genomecov -split -d -ibam ${bam} | cut -f2- - | \
        awk '{print \$0"\t${bam.simpleName.split('__')[0]}\t${bam.simpleName.split('__')[1]}"}' | \
        sed '1ipos\tdepth\tsample\tref' > ${bam.simpleName}.cov

    """
}
