process VAR_TSV {
    publishDir params.outdir, mode:'copy'
    publishDir params.datalakedir, mode:'copy'
    label "pandas"

    input:
        path(var_vcf)

    output:
        path("variants.tsv")

    script:
        """
        vcf_to_tsv.py --vcf $var_vcf --output variants.tsv --workflow_id $params.unique_id --frequency-called $params.freq_field --min-freq $params.min_variant_freq
        """
}
