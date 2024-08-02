process COVERAGEBED {
    label "bedtools"

    input:
        tuple val(sample),
              path(bam),
              path(region_bed),
              path(itr_bed),
              path(border_bed) // channel:  [val(sample), .bam, .regions.bed, .ITR_ITR.bed., border.bam]

    output:
        path "${bam.simpleName}.dissect.count", emit: counts

    script:
    """
    bedtools coverage -F 1.0 -a ${region_bed} -b ${bam} > tmp.txt
    bedtools coverage -f 1.0 -F 1.0 -a ${itr_bed} -b ${bam} >> tmp.txt  
    bedtools coverage -f 1.0 -a ${border_bed} -b ${bam} >> tmp.txt  
    cat tmp.txt | cut -f1-4,7 |  awk 'FNR==1 {tol=\$5} {print \$0"\t"\$5/tol}' > ${bam.simpleName}.dissect.count
    """
}
