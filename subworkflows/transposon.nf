
process ALIGN_TO_TRANSPOSONS {
    publishDir params.outdir, mode:'copy'
    publishDir params.datalakedir, mode:'copy'
    label "ncores"
    label "pbmm2"
    label "himem"

    input:
        path(sample_bam)
        path(transposon_fasta)

    output:
        path("*.bam"), emit: bam
        path("*.bai"), emit: bai
 
    script:
    """
        pbmm2 index $transposon_fasta ${transposon_fasta.baseName}.mmi
        pbmm2 align $transposon_fasta ${sample_bam} transposons.bam --sort --bam-index BAI -j ${task.cpus}
    """
}

process FLANKS_FASTQ {
    label "pysam"

    input:
        path(bam)

    output:
        path("*.fastq"), emit: fq
        path("*.tsv"), emit: tsv

    script:
    """
        transposon_flanks.py --bam $bam --fastq transposon_flanks.fastq --tsv tn_cov.tsv
    """
}

process ALIGN_TO_COMBINED {
    label "pbmm2"

    input:
        path(flanks_fastq)
        path(plasmid_fasta)

    output:
        path("*.bam")

    script:
    """
        pbmm2 index $plasmid_fasta ${plasmid_fasta.baseName}.mmi
        pbmm2 align $plasmid_fasta $flanks_fastq transposon_flanks.bam --sort --bam-index BAI
    """
}

process COMBINE_REF {
    label "ubuntu"

    input:
        path(fasta1)
        path(fasta2)

    output:
        path("*.fasta")

    script:
    """
        cat $fasta1 $fasta2 > combined.fasta
    """
}

process FLANK_TABLES {
    label "pysam"

    input:
       path(bam)

    output:
       path("*.tsv")

    script:
    """
        tn_flank_alignments.py --bam $bam --output flank
    """
}

process MERGE_TABLES {
    publishDir params.outdir, mode:'copy'
    publishDir params.datalakedir, mode:'copy'
    label "pandas"

    input:
        path(cov)
        path(flanks)

    output:
        path("*.tsv")

    script:
    """
        tn_merge_tables.py --cov $cov --uniqueid $params.unique_id --flanks ${flanks.join(" ")}
    """
}

workflow FIND_TRANSPOSONS {

     take:
         sample_bam
         transposon_fasta
         plasmid_fasta
         ecoli_fasta

     main:
         if(params.transposon_analysis) {
             ALIGN_TO_TRANSPOSONS ( sample_bam, transposon_fasta )
             ch_bam = ALIGN_TO_TRANSPOSONS.out.bam
             ch_bai = ALIGN_TO_TRANSPOSONS.out.bai
             FLANKS_FASTQ ( ch_bam )
             ch_flank_fq = FLANKS_FASTQ.out.fq
             ch_tn_cov = FLANKS_FASTQ.out.tsv
             ch_comb_ref = COMBINE_REF ( plasmid_fasta, ecoli_fasta ) 
             ch_flank_bam = ALIGN_TO_COMBINED ( ch_flank_fq, ch_comb_ref )
             ch_flank_tables = FLANK_TABLES ( ch_flank_bam ).collect()
             ch_merged_table = MERGE_TABLES ( ch_tn_cov, ch_flank_tables )
        } else {
              ch_bam = Channel.empty()
              ch_bai = Channel.empty()
              ch_merged_table = Channel.empty()
        }

     emit:
         bam = ch_bam
         bai = ch_bai
         metrics = ch_merged_table
}
