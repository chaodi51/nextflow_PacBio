process PBMM2_ALIGN {
    publishDir params.outdir, mode:'copy'
    publishDir params.datalakedir, mode:'copy'

    label "ncores"
    label "pbmm2"
    label "himem"

    input:
        path(sample)
        path(ref)

    output:
        path "*.bam", emit: bam
        path "*.bai", emit: bai

    script:
    """
        pbmm2 align ${ref} ${sample} ${ref.baseName}.bam \
        --sort --bam-index BAI -j ${task.cpus}
    """
}
