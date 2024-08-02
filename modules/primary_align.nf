process PRIMARY_ALIGN {
    publishDir "${params.outdir}/pbmm2/${params.run}", mode:'copy'
    label "samtools"

    input:
        path bam

    output:
        path "${bam.baseName}.primary.bam", emit: pri_bam
        path "${bam.baseName}.primary.bam.bai", emit: pri_bai
        
    script:
    """
        samtools view -b -F 0x800 -F 0x100 $bam > ${bam.baseName}.primary.bam
        samtools index ${bam.baseName}.primary.bam
    """
}