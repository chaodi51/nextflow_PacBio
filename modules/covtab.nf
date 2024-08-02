process COVTAB {
    publishDir "${params.outdir}/regions_tab/${params.run}", mode:'copy'
    label "ubuntu"

    input:
        path region_counts // collection of count files

    output:
        path "${params.run}.region_counts.tsv", emit: cov_tab

    script:
    """
    for i in ${region_counts}; do
        sampleid=`echo \$i | sed 's/__.*//g'`

        awk '{print "'\${sampleid}'\t"\$1"\t"\$4"\t"\$2"\t"\$3"\t"\$5"\t"\$6}' \$i >> ${params.run}.region_counts.tsv
    done

    sed -i '1isampple\treference\tregion\tstart\tend\tread_counts\tpercent' ${params.run}.region_counts.tsv
    """
}
