process HOST_ALIGN {
    publishDir "${params.outdir}/pbmm2/${params.run}", mode:'copy'
    label "ncores"
    label "pbmm2"
    
    input:
        path(sample)
        path(host)

    output:
        path "${sample.simpleName}__host.bam", emit: host_bam
        path "${sample.simpleName}__host.bam.bai", emit: host_bai
        
    script:
    """
        pbmm2 align ${host} ${sample} ${sample.simpleName}__host.bam \
        --sort --bam-index BAI -j ${task.cpus}
    """
}
