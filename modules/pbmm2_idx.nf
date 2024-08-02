process PBMM2_IDX {
    label "pbmm2"

    input:
        path ref_fa

    output:
        path "${ref_fa.baseName}.mmi", emit: mmi

    script:
    """
        pbmm2 index $ref_fa "${ref_fa.baseName}.mmi"
    """
}
