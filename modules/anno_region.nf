process ANNO_REGION {
    publishDir "${params.outdir}/trunc/${params.run}", mode:'copy'
    label "biopython"

    input:
        tuple path(bed), path(cov) // channel: [ path(ref.bed), path(.cov) ]

    output:
        path "${cov.simpleName}.anno.cov", emit: anno_cov

    script:
    """
    anno_region.py ${bed} ${cov} ${cov.simpleName}.anno.cov
    """
}