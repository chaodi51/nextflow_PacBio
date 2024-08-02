//
// chimeric reads detection and potential truncation sites plot
//


// import modules
include { CHIMERA_REF } from '../modules/chimera_ref.nf'
include { PBMM2_COMB_IDX } from '../modules/pbmm2_comb_idx.nf'
include { ALIGN_CHIMERA } from '../modules/align_chimera.nf'
include { ALVIS } from '../modules/alvis.nf'

process FORMAT {
    publishDir params.outdir, mode:'copy'
    publishDir params.datalakedir, mode:'copy'

    label "pandas"
    
    input:
        path(from_alvis)

    output:
        path("chimeras.tsv")

    script:
    """
        chimera_table.py --alvis $from_alvis --workflow-id $params.unique_id --output chimeras.tsv
    """
}

workflow CHIMERA {
    take:
        host 
        ref
        sample

    main:
        CHIMERA_REF( host, ref )
        ch_comb_ref = CHIMERA_REF.out.comb_ref

        PBMM2_COMB_IDX( ch_comb_ref )
        ch_comb_ref_mmi = PBMM2_COMB_IDX.out.mmi
        
        // sample_combref = Channel
        //    .fromPath ( "${projectDir}/conf/sample_combref.tsv", checkIfExists: true )
        //  .splitCsv ( header:true, sep: "\t" )
        //  .filter { it.run == params.run }
        //  .map { it -> tuple( "${it.Reference}", "${params.datadir}${it.Sample}.fastq.gz") }

        //sample_combref_pair = ch_comb_ref_mmi
        //    .filter { !it.name.contains('Lambda') }
        //    .map { it -> tuple("${it.simpleName}", it) }
        //    .join( sample_combref )
        //    .map { tuple(it[1], it[2]) }

        // sample_combref_pair = sample_ref_pair
        //     .filter { !it[1].contains('Lambda') }
        //     .map{ tuple(it[0], "${projectDir}/data/outputs/v${params.VERSION}/chimera/ref/${params.run}/host_${it[1].split('/')[-1]}") }

        ALIGN_CHIMERA( ch_comb_ref, sample )
        ch_paf = ALIGN_CHIMERA.out.paf

        ALVIS( ch_paf )
        ch_alvis_out = ALVIS.out.chimeras_txt

        ch_formatted_chimeras = FORMAT( ch_alvis_out )
        
    emit:
        alvis_chimeras_txt = ch_alvis_out
        alvis_chimeras_tsv = ch_formatted_chimeras
    }
