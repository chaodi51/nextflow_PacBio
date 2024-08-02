process TIDYMULTIQC {
  publishDir params.outdir, mode:'copy'
  publishDir params.datalakedir, mode:'copy'
  label "pandas"

  input:
  path(json)

  output:
  path("multiqc_data.tsv")

  script:
  """
  tidymultiqc.py --multiqc-json $json --workflow_id ${params.unique_id} --output-tsv multiqc_data.tsv
  """
}
