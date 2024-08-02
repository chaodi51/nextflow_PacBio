process STATS_TAB {
    publishDir "${params.outdir}/stat_tab/${params.run}", mode:'copy'
    label "ubuntu"
    
    input:
        path total_reads // channel: ch_total_reads
        path all_aligned_reads // channel: ch_all_aligned_reads

    output:
        path "${params.run}.stats.tsv", emit: stats_tab

    script:
    """
    # total reads table
    for j in ${total_reads}; do 
        cat \$j >> total_reads.tsv
    done

    # mapped reads table
    for i in ${all_aligned_reads}; do
        sampleid=`cat \$i | cut -f1| awk 'split(\$0, a, "__") {print a[1]}'`
        ref=`cat \$i | cut -f1| awk 'split(\$0, a, "__") {print a[2]}'`
        tot_reads=`cat total_reads.tsv | grep \$sampleid |cut -f2`
        
        echo -en \${sampleid}'\\t' >> ${params.run}.stats.tsv
        echo -en \$tot_reads'\\t' >> ${params.run}.stats.tsv
        echo -en \$ref'\\t' >> ${params.run}.stats.tsv
        awk '{print \$2}' \$i >> ${params.run}.stats.tsv
    done

    sort -k1,1 ${params.run}.stats.tsv | awk '{print \$0"\t"\$4/\$2}' > foo
    sed 1i'Sample\tTotal_reads\tRef\tMapped_reads(primary)\tFraction' foo > ${params.run}.stats.tsv

    """
}


