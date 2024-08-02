process MakePileup {
  label "samtools"
  label "himem"

  input:
  path(bam)
  path(fasta)

  output:
  path("*.dat")

  script:
  """
  samtools mpileup -d $params.variant_maxdepth -B -Q $params.min_consensus_baseq -aA -f $fasta $bam > pileup.dat
  """
}

process CallCns {
  label "varscan"

  input:
  path(pileup)

  output:
  path("*.vcf")

  script:
  """
  ( cat $pileup | varscan mpileup2cns --min-var-freq $params.min_variant_freq --strand-filter 0 --output-vcf > cns.vcf ) || ( echo -e "0\t0\t0\t0" | varscan mpileup2cns --strand-filter 0 --output-vcf > cns.vcf )
  """
}

process MakeVariationTable {
  publishDir params.outdir, mode:'copy'
  publishDir params.datalakedir, mode:'copy'
  label "pandas"

  input:
  path(pileup)
  val(thresholds)

  output:
  path("var_table.tsv")

  script:
  """
  pileup_to_consensus.py -p $pileup --workflow_id ${params.unique_id} -o var_table.tsv -t ${thresholds.join(" ")} -i $params.iupac_thresh
  """
}

process ConsensusFasta {
   publishDir params.outdir, mode:'copy'
   publishDir params.datalakedir, mode:'copy'
   label "pandas"

   input:
   path(consensus_tsv)

   output:
   path("*.fasta")

   script:
   """
   consensus_to_fasta.py --tsv $consensus_tsv --fasta consensus.fasta
   """
}

process PairwiseAlignment {
  publishDir params.outdir, mode:'copy'
  publishDir params.datalakedir, mode:'copy'

  input:
  path(ref_fa)
  path(consensus_fa)

  output:
  path("consensus.txt")

  script:
  """
  needle -awidth3 60 -auto  -asequence $ref_fa -bsequence $consensus_fa -outfile consensus.txt
  """
}

workflow CONSENSUS {

  take:
    ch_bam
    ch_fasta
    ch_threshlist

  main:
    pileup_ch = MakePileup(ch_bam, ch_fasta)
    final_cns_ch = MakeVariationTable(pileup_ch, ch_threshlist)
    cons_fa_ch = ConsensusFasta(final_cns_ch)
    pairwise_ch = PairwiseAlignment(ch_fasta, cons_fa_ch)

  emit:
    consensus = pairwise_ch
    var_table = final_cns_ch
}
