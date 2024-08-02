//
// calling SNV and indel using lofreq
//


// import modules
include { INDELQUAL } from '../modules/indelqual.nf' //Insert indel qualities
include { CALL_VAR } from '../modules/call_var.nf'
include { VAR_TSV } from '../modules/var_tsv.nf'

process DOWNSAMPLE {
    label "hugemem"
    label "samtools"

    input:
        path(bam)

    output:
        path("*.bam"), emit: bam
        path("*.bai"), emit: bai

    script:
    """
    samtools view -H $bam > header.txt
    samtools view -F2048 $bam | shuf -n $params.ndownsample > reads.sam
    cat header.txt reads.sam >> unsorted.sam
    samtools view -h unsorted.sam | samtools sort -o sorted.bam
    samtools index sorted.bam
    """
}

process CALL_SV {
    publishDir params.outdir, mode:'copy'
    publishDir params.datalakedir, mode:'copy'
    label "pacbio"
    label "hugemem"

    input:
        path(bam)
        path(fasta)

    output:
        path("struct_var.vcf")

    script:
    """
    pbsv discover $bam sample.svsig.gz
    pbsv call --ccs $fasta sample.svsig.gz struct_var.vcf
    """
}

workflow VAR {
    take:
        ch_bam
        ch_reference
        ch_refindex

    main:
        INDELQUAL( ch_bam, ch_reference )
        ch_indelqual_bam = INDELQUAL.out.indelqual_bam
        ch_indelqual_bai = INDELQUAL.out.indelqual_bai

        DOWNSAMPLE(ch_indelqual_bam)
        ch_subsamp_bam = DOWNSAMPLE.out.bam
        ch_subsamp_bai = DOWNSAMPLE.out.bai

        ch_sv_vcf = CALL_SV( ch_subsamp_bam, ch_reference)

        ch_var_vcf = CALL_VAR( ch_subsamp_bam, ch_subsamp_bai, ch_reference, ch_refindex )

        ch_var_tsv = VAR_TSV( ch_var_vcf )

    emit:
        snp_tsv = ch_var_tsv
        sv_vcf = ch_sv_vcf
}
