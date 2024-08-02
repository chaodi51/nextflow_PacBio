process FIX_REGIONS {
    publishDir params.refdir, mode:'copy'
    label "pandas"

    input:
        path bed

    output:
        path "${bed.simpleName}.standardized.bed", emit: standardized_bed
        
    script:
    """
    fix_regions.py ${bed} ${bed.simpleName}.standardized.bed
    """
}