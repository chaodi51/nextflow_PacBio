process ALIGN_CHIMERA {
    label "minimap2"
    label "ncores"
    label "himem"

    input:
        path(ref)
        path(sample)

    output:
        path "${sample.simpleName}__${ref.baseName}.paf", emit: paf
        
    script:
    """
        minimap2 -x map-hifi -c -t ${task.cpus} ${ref} ${sample} > ${sample.simpleName}__${ref.baseName}.paf

    """
}
