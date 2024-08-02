process AddMDString {
  label "samtools"

  input:
  path(bam)
  path(fasta)

  output:
  path("*.bam"), emit: bam
  path("*.bai"), emit: bai

  script:
  """
  samtools calmd -b $bam $fasta > with_md.bam
  samtools index with_md.bam
  """ 
}

process FormMetrics {

  input:
  path(bam)
  path(bamindex)
  path(fasta)

  output:
  path("*.tsv")

  script:
  """
  FormMetrics.py -f $fasta -b $bam -o seqqual_metrics.tsv
  """
}

process CatMetrics {
  label "ubuntu"

  input:
  path(metrics)

  output:
  path("*.tsv")

  script:
  """
  cat ${metrics.join(" ")} | sort -u -k1n -k2nr -k12nr >> tmp
  grep coverage tmp > header
  grep -v coverage tmp | awk '!a[\$1]++' > data
  cat header data >> combinedposmetrics.tsv
  """
}

process PositionalMetrics {
  publishDir params.outdir, mode:'copy'
  publishDir params.datalakedir, mode:'copy'
  label "pandas"

  input:
  path(metrics)

  output:
  path("positional_metrics.tsv")

  script:
  """
  #!/usr/bin/env python

  import pandas as pd
  x = pd.read_csv(\"$metrics\", sep=\"\\t\", index_col=0).reset_index().rename(columns={\"index\":\"position\"})
  x[\"workflow_id\"] = \"$params.unique_id\"
  x.to_csv(\"positional_metrics.tsv\", sep=\"\\t\", index=False)
  """
}

workflow POSITIONAL_METRICS {

  take:
    ch_bam
    ch_fasta

  main:
    AddMDString(ch_bam, ch_fasta)
    initial_metrics_ch = FormMetrics(AddMDString.out.bam, AddMDString.out.bai, ch_fasta)
    combined_metrics_ch = CatMetrics(initial_metrics_ch)
    final_pos_metrics_ch = PositionalMetrics(combined_metrics_ch)

  emit:
    pos_metrics = final_pos_metrics_ch
}
