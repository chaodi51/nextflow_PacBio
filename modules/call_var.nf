process CALL_VAR {
    label "lofreq"
    label "ncores"
    label "hugemem"
    
    input:
        path(bam) 
        path(bai) 
        path(ref) 
        path(ref_fai)

    output:
        path "*.vcf"
        
    script:
    """
    lofreq call-parallel --pp-threads 3 -B -N --call-indels -f ${ref} \
    -o variants.vcf ${bam}
    """
}
