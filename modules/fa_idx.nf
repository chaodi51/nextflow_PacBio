process FA_IDX {
    label "samtools"
    
    input:
        path ref_fa

    output:
        path "${ref_fa}.fai", emit: fai

    script:
    """
        samtools faidx $ref_fa
    """
}
