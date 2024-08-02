
process PBMM2_ALIGN_UNMAPPED {

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
        --sort --bam-index BAI -j ${task.cpus} --unmapped
    """
}

process CAT_REFS {
    label "ubuntu"

    input:
        path(references)

    output:
        path("*.fasta")

    script:
    """
    cat ${references.join(" ")} >> combined_ref.fasta
    """
}

process GET_UNMAPPED {
    label "samtools"

    input:
        path(bam)

    output:
        path("*.tsv")

    script:
    """
    count=`samtools view -f4 $bam | wc -l | tr -d "\\n"`
    printf \"%s\\t%s\\n" unmapped \$count > unmapped.counts.tsv
    """
}

workflow UNALIGNED {
    take:
        sample
        references

    main:
        ch_comb_fa = CAT_REFS(references)
        PBMM2_ALIGN_UNMAPPED( sample, ch_comb_fa )
        ch_bam = PBMM2_ALIGN_UNMAPPED.out.bam
        ch_bai = PBMM2_ALIGN_UNMAPPED.out.bai
        unmapped_count_ch = GET_UNMAPPED(ch_bam)
        
    emit:
        count = unmapped_count_ch
    }
