process CHIMERA_REF {
    label "biopython"

    input:
        path(host) 
        path(ref)

    output:
        path "host_${ref.simpleName}.fasta", emit: comb_ref
        
    script:
    """
    mk_chimera_ref.py "host" ${ref.simpleName} ${host} ${ref} host_${ref.simpleName}.fasta
    """
}
