//
// methylation quantification
//

process PRIMROSE {
    publishDir params.outdir, mode:'copy'
    publishDir params.datalakedir, mode:'copy'
    label 'pacbio'

    input:
    path(hifibam)

    output:
    path("5mc.hifi_reads.bam")

    script:
    """
    primrose $hifibam 5mc.hifi_reads.bam
    """
}

process METHYL_ALIGN {
    publishDir params.outdir, mode:'copy'
    publishDir params.datalakedir, mode:'copy'
    label 'pacbio'
    label 'himem'

    input:
    path(methylbam)
    path(reference)

    output:
    path("aligned.5mc.hifi_reads.bam")

    script:
    """
    pbmm2 align $reference $methylbam aligned.5mc.hifi_reads.bam --sort -j 16 -J 16
    """
}

process METHYL_INDEX {
    label 'samtools'

    input:
    path(bam)

    output:
    path("*.bai")

    script:
    """
    samtools index $bam
    """
}

process CPG_SITES {
    label 'pacbio'

    input:
    path(mc_bam)
    path(mc_bai)

    output:
    path("*.bed")

    script:
    """
    /bin/pb-CpG-tools-v2.1.1-x86_64-unknown-linux-gnu/bin/aligned_bam_to_cpg_scores \
      --bam $mc_bam \
      --output-prefix methylation \
      --model /bin/pb-CpG-tools-v2.1.1-x86_64-unknown-linux-gnu/models/pileup_calling_model.v1.tflite \
      --threads 8
    """
}

process MAKE_TABLE {
    publishDir params.outdir, mode:'copy'
    publishDir params.datalakedir, mode:'copy'
    label 'pandas'

    input:
    path(table)

    output:
    path("methylation.tsv")

    script:
    """
    methylation_table.py --bed $table --uniqueid $params.unique_id
    """
}

workflow METHYLATION {
    take:
        unalignedbam
        reference

    main:
        hifibam_ch = unalignedbam
        methylbam_ch = PRIMROSE(hifibam_ch)
        aligned_bam_ch = METHYL_ALIGN(methylbam_ch, reference)
        aligned_bai_ch = METHYL_INDEX(aligned_bam_ch)
        methyl_bed_ch = CPG_SITES(aligned_bam_ch, aligned_bai_ch)
        methyl_tsv_ch = MAKE_TABLE(methyl_bed_ch)

    emit:
        tsv = methyl_tsv_ch
        primrose_bam = methylbam_ch
        aligned_primrose_bam = aligned_bam_ch
        aligned_primrose_bai = aligned_bai_ch
}
