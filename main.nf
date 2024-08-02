#!/usr/bin/env nextflow 

/*
 * PacBio long reads sequencing analysis pipeline
 *
 * Authors:
 * - Chao Di <chao.di@sparktx.com>
 * - Nathan Johnson <nathan.johnson@sparktx.com>
 */

/* 
 * enables modules 
 */
nextflow.enable.dsl = 2

/*
 * Default pipeline parameters. They can be overriden on the command line eg.
 * given `params.foo` specify on the run command line `--foo some_value`.
 */

//params.VERSION = '1.0.0'
/* params for Tower */
//projectDir = "s3://sparktx-nextflow-pipeline/projects/nf_PacBio/"
params.run_id = ""
/* local params also work for Tower */
// this is the human genome, needed for chimeras
params.host = "s3://sparktx-nextflow-pipeline/projects/hd-short-read/ref/GCF_009914755.1_T2T-CHM13v2.0_genomic.fasta"
params.repcap_flavor = "pSpark100PK"
params.repcap = "s3://sparktx-nextflow-pipeline/projects/hd-short-read/ref/${params.repcap_flavor}.fasta"
params.helper = "s3://sparktx-nextflow-pipeline/projects/hd-short-read/ref/pAD2HPv2.fasta"
params.tn_fa = 's3://sparkds-nextflow-references/production/linear_contam/transposon_master.fasta'
params.ecoli_fa = 's3://sparkds-nextflow-references/production/linear_contam/dh10b.fasta'
params.iupac_thresh = 0.4
params.freq_field = "INFO_AF"
params.transposon_analysis = false

if (params.VERSION.startsWith("v")) {
    params.version = params.VERSION.substring(1)
} else {
    params.version = params.VERSION
}

if (params.run_id) {
  params.unique_id = params.run_id
} else {
  params.unique_id = workflow.sessionId
}

//params.outdir = "s3://sparktx-nextflow-pipeline/projects/nf_pacbio/outputs/${params.unique_id}/"
// params.outdir = "${projectDir}/debug"
params.outdir = "s3://sparktx-nextflow-pipeline/projects/nf_pacbio/outputs/v${params.version}/${params.unique_id}/"
//params.datalakedir = "s3://sparkds-datalake-groupdropin-bioinformatics/nf_pacbio/${params.unique_id}/"
params.datalakedir = "s3://sparkds-datalake-groupdropin-bioinformatics/nf_pacbio/v${params.version}/${params.unique_id}/"
// params.refpath = "${params.refdir}${params.reference}.fasta"
params.refpath = "${params.reference}"
params.plasmid_backbone = ""
params.rsconnect_url = "https://connect.sparkds.io/"
params.run_report = true

log.info """\
 P A C B I O - N F   P I P E L I N E
 ===================================
 run          : ${params.unique_id} 
 reads        : ${params.reads}
 outdir       : ${params.outdir}
 """


// import subworkflows
include { MAKE_FASTQ } from './subworkflows/make_fastq.nf'
include { METHYLATION } from './subworkflows/methylation.nf'
include { QC } from './subworkflows/qc.nf'
include { REF } from './subworkflows/prepare_reference.nf'
include { ALIGN } from './subworkflows/align.nf'
include { UNALIGNED } from './subworkflows/unaligned.nf'
include { READ_STATS } from './subworkflows/read_stats.nf'
include { REGION_READS } from './subworkflows/region_reads.nf'
include { LEN } from './subworkflows/length.nf'
include { COV } from './subworkflows/cov.nf'
include { TRUNC } from './subworkflows/trunc.nf'
include { VAR } from './subworkflows/variants.nf'
include { CHIMERA } from './subworkflows/chimera_alvis.nf'
include { POSITIONAL_METRICS } from './subworkflows/positional_metrics.nf'
include { CONSENSUS } from './subworkflows/consensus.nf'
include { PLASMID_COV } from './subworkflows/plasmid_cov.nf'
include { FIND_TRANSPOSONS } from './subworkflows/transposon.nf'

// import modules
// include { FIX_REGIONS } from './modules/fix_regions.nf'

include { FinalizeOutputs as FINALIZE_OUTPUTS } from './nf-finalize-outputs/modules/finalize-outputs/main.nf'
include { ReportRunner as REPORT_RUNNER } from './nf-report-runner/main'


process DIRNAME {
    label "ubuntu"

    input:
    val(ref)

    output:
    stdout

    script:
    """
    dirname $ref | sed 's/\\r\$//' | sed -z 's/\\n//g'
    """
}

process CONTIGNAME {
    label "ubuntu"

    input:
    path(ref)

    output:
    stdout

    script:
    """
    head -1 $ref | cut -c2- | sed 's/\\r\$//' | sed -z 's/\\n//g'
    """
}

process CAT_COUNTS {
    publishDir params.outdir, mode:'copy'
    publishDir params.datalakedir, mode:'copy'
    label "ubuntu"

    input:
    path(countfiles)
    
    output:
    path("reference_counts.tsv"), emit: reference_counts

    script:
    """
    printf \"%s\\t%s\\n\" "reference" "count" > reference_counts.tsv
    cat ${countfiles.join(" ")} >> reference_counts.tsv
    sed -i 's/GCF_009914755/human (T2T)/g' reference_counts.tsv 
    """
}

process SOFTCLIPPED {
  publishDir params.outdir, mode:'copy'
  publishDir params.datalakedir, mode:'copy'
  label 'samtoolspysam'

  input:
  path(bam)

  output:
  path("softclipped_reads.bam"), emit: sc_bam
  path("softclipped_reads.bam.bai"), emit: sc_bai

  script:
  """
  samtools index $bam
  samtools view -H $bam > header
  find_softclipped_reads.py --bam $bam > screads
  cat header screads | samtools view -h -o 'softclipped_reads.bam' -
  samtools index softclipped_reads.bam 
  """
}

process METADATA {
  publishDir params.outdir, mode:'copy'
  publishDir params.datalakedir, mode:'copy'
  label "ubuntu"

  input:
  val(sample)
  val(primary_ref)
  val(primary_contig)
  val(refdir)

  output:
  path("metadata.tsv"), emit: metadata_tsv

  script:
  """
  printf \"%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n\" workflow_id version sample_id group department program species lot requester email date project platform machine run fastq_path ref_path primary_ref primary_contig output_path >> metadata.tsv
  printf \"%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n\" \"$params.unique_id\" \"$params.VERSION\" \"${sample.baseName}\" \"$params.group\" \"$params.department\" \"$params.program\" \"$params.species\" \"$params.lot\" \"$params.requester\" \
    \"$params.email\" \"$params.date\" \"$params.project\" \"$params.platform\" \"$params.machine_id\" \"$params.run\" \"$params.reads\" \"$refdir\" \"${primary_ref.baseName}\" \"${primary_contig}\" \"$params.outdir\" >> metadata.tsv  
  """
}

process REPORT {
  secret 'RSTUDIO_CONNECT_API_USER'
  secret 'RSTUDIO_CONNECT_API_KEY'
  stageInMode 'copy'

  input:
  path(qc)
  path(read_stats)
  path(len)
  path(trunc)
  path(variant)
  path(positional)
  path(chimera)
  path(metadata)
  path(sample_report_template)
  path(toc_template)
  path(sample_report_script)
  path(toc_script)
  path(footer)
  path(methylation)
  path(countsfile)
  path(consensus)

  output:
  stdout

  script:
  """
  export API_USER=\$RSTUDIO_CONNECT_API_USER && export API_KEY=\$RSTUDIO_CONNECT_API_KEY
  deploy_nextflow_reports.py --version $params.version --sessionid $params.unique_id --debug --report-template $sample_report_template \
  --toc-template $toc_template --report-script $sample_report_script --toc-script $toc_script
  """
}

/* 
 * RUN MAIN WORKFLOW
 */
workflow {

//    plot_read_len = Channel.fromPath( "${baseDir}/bin/plot_read_len.rmd" ) 
//    plot_cov = Channel.fromPath( "${baseDir}/bin/plot_cov.rmd" ) 
//    plot_trunc = Channel.fromPath( "${baseDir}/bin/plot_trunc.rmd" ) 
//    plot_chimera_trunc = Channel.fromPath( "${baseDir}/bin/plot_chimera_trunc.rmd" )

    ch_ref_fa = Channel.fromPath( params.refpath, checkIfExists: true )
    ch_inputbam = Channel.fromPath( params.reads, checkIfExists: true )
    ch_transposon_fa = Channel.fromPath( params.tn_fa, checkIfExists: true )
    ch_ecoli_fa = Channel.fromPath( params.ecoli_fa, checkIfExists: true )
    samplename = "${file(params.reads).getSimpleName()}"
    refname = "${file(params.refpath).getBaseName()}"
    ch_refname = Channel.of(refname)
    ch_thresholds = Channel.of(params.consensus_thresh)

    refdir = DIRNAME(params.refpath)

    MAKE_FASTQ(ch_inputbam)
    ch_sample = MAKE_FASTQ.out.fastq

//    sample_ref_pair = Channel
//        .fromPath ( "${projectDir}/conf/master_mapping.tsv", checkIfExists: true )
//        .splitCsv ( header:true, sep: "\t" )
//        .filter { it.run == params.run }
//        .map { it -> tuple("${params.datadir}${it.Sample}.fastq.gz", "${params.refdir}${it.Reference}.mmi") }

    // subworkflows: fastqc and multiqc
    QC ( ch_sample )

    // subworkflow: make reference indexes 
    REF( ch_ref_fa )
    ch_ref_fai = REF.out.ref_fai
    primary_contig = CONTIGNAME(ch_ref_fa)

    // subworkflow: align reads
    ch_host = Channel.fromPath( params.host, checkIfExists: true )
    ch_repcap = Channel.fromPath( params.repcap, checkIfExists: true )
    ch_helper = Channel.fromPath( params.helper, checkIfExists: true )
//    ch_sample_host_pair = ch_sample.combine( ch_host )
    all_refs_ch = ch_ref_fa.mix(ch_host).mix(ch_repcap).mix(ch_helper)

    if(params.plasmid_backbone) {
        all_refs_ch = all_refs_ch.mix(Channel.fromPath( params.plasmid_backbone, checkIfExists: true ))
    }

    ALIGN( ch_sample.first(), all_refs_ch, ch_refname )
    UNALIGNED( ch_sample, all_refs_ch.collect())
    primary_bam_ch = ALIGN.out.bam.filter( {it.name.contains(refname)} )
    // need to remove unaligned etc per original ALIGN swf?
    // turns out only non-primary and supplementary alignments were previously removed

    if(params.plasmid_backbone) {
        plasmid_name = "${file(params.plasmid_backbone).getBaseName()}"
        ch_backbone = Channel.fromPath( params.plasmid_backbone, checkIfExists: true )
        ch_plasmid_bam = ALIGN.out.bam.filter( { ( (it.name.startsWith(plasmid_name)) & (!it.name.startsWith(refname) ) ) } )
        PLASMID_COV ( ALIGN.out.merged, ch_plasmid_bam, plasmid_name )
        plasmid_depth = PLASMID_COV.out.plasmid_cov
    }
    else {
        plasmid_depth = Channel.empty()
    }

    // subworkflow: stats table (total and mapped reads)
    READ_STATS( ch_sample, primary_bam_ch )
    ch_counts = CAT_COUNTS(ALIGN.out.count.mix(UNALIGNED.out.count).collect())

/*
    // subworkflow: count reads in specific regions, including ITR-ITR full-length
    // this is all of dubious usefulness because the underlying resources are never kept current
    ch_pri_bam = ch_align
        .pri_bam.filter { !it.name.contains("Lambda") }
        .map { tuple(it.simpleName.split('__')[1], it) }
    ch_ref_regions = Channel
        .fromPath( "${params.refdir}/elements/*.regions.bed", checkIfExists: true )
        .map { tuple(it.simpleName, it) }
    ch_ref_ITR_ITR = Channel
        .fromPath( "${params.refdir}/elements/*.ITR-ITR.bed", checkIfExists: true )
        .map { tuple(it.simpleName, it) }
    ch_ref_border = Channel
        .fromPath( "${params.refdir}/elements/*.border.bed", checkIfExists: true )
        .map { tuple(it.simpleName, it) }

    ch_bam_region_beds = ch_pri_bam
        .join(ch_ref_regions)
        .join(ch_ref_ITR_ITR)
        .join(ch_ref_border)

    REGION_READS( ch_bam_region_beds )
*/

    METHYLATION( ch_inputbam, ch_ref_fa )
    // subworkflow: read length, alignment length, lambda DNA normalization
    LEN( primary_bam_ch ) // a plot formerly got made here via bin/plot_read_len.rmd, recreate in report

    // subworkflow: coverage plots on each ref
    // Just replace with standard sequence quality metrics
    // COV( primary_bam_ch) // a plot formerly got made here via bin/plot_cov.rmd, recreate in report

    // subworkflow: plot read starts and ends as potential truncation hotspots
    TRUNC( primary_bam_ch ) // a plot formerly got made here via bin/plot_trunc.rmd, recreate in report

    // subworkflow: variants calling using lofreq
    // ch_bam_ref_pair = ch_align.pri_bam
    //    .map { tuple(it, "${params.refdir}/${it.simpleName.split('__')[1]}.fasta") }
    VAR( primary_bam_ch, ch_ref_fa, ch_ref_fai )

    // module: standardize the region annotations in the plasmids
    // this is also hopefully not truly necessary because of the same dependence on plasmid annotations that
    //    are not generally available
    // ch_bed_refs = Channel.fromPath( "${projectDir}/conf/liver_plasmid_regions.bed", checkIfExists: true )
    // FIX_REGIONS( ch_bed_refs )
    // ch_fixed_region_bed = FIX_REGIONS.out.standardized_bed

    // subworkflow: chimeric reads detection using alvis
    // ch_host_ref = ch_host.combine( ch_ref_fa.filter{ !it.name.contains("Lambda") } )
    CHIMERA( ch_host, ch_ref_fa, ch_sample ) // a plot formerly got made here via bin/plot_chimera_trunc.rmd, recreate in report

    CONSENSUS(primary_bam_ch, ch_ref_fa, ch_thresholds)

    FIND_TRANSPOSONS ( ch_inputbam, ch_transposon_fa, ch_ref_fa, ch_ecoli_fa )

    POSITIONAL_METRICS(primary_bam_ch, ch_ref_fa)

    SOFTCLIPPED(primary_bam_ch)

    ch_metadata = METADATA( ch_inputbam, ch_ref_fa, primary_contig, refdir ) 

    sample_report_template_ch = Channel.fromPath("${projectDir}/templates/template_nextflow_report.Rmd")
    toc_template_ch = Channel.fromPath("${projectDir}/templates/template_nextflow_aggregate.Rmd")
    sample_report_script_ch = Channel.fromPath("${projectDir}/Rscripts/deploy_nextflow_sample.R")
    toc_script_ch = Channel.fromPath("${projectDir}/Rscripts/deploy_nextflow_agg.R")
    footer_ch = Channel.fromPath("${projectDir}/templates/footer.html")

    static_files = Channel.empty().mix(
        QC.out.multiqc_tsv,
        READ_STATS.out.read_stats,
        LEN.out.read_len,
        TRUNC.out.st_end_cov,
        VAR.out.snp_tsv,
        POSITIONAL_METRICS.out.pos_metrics,
        CHIMERA.out.alvis_chimeras_tsv,
        METHYLATION.out.tsv,
        CONSENSUS.out.consensus,
        CAT_COUNTS.out.reference_counts,
        plasmid_depth,
        ch_ref_fa
    ).collect()
    rsconnect_tags = Channel.of('TDO', 'NGS/Proteomics', params.program, 'DNA-Sequencing', 'Adeno-associated Virus 2').collect()

    outputs = Channel.empty().mix(
        CAT_COUNTS.out.reference_counts,      // reference_counts.tsv
        METADATA.out.metadata_tsv,            // metadata.tsv
        LEN.out.align_len,                    // align_len_count.tsv
        LEN.out.read_len,                     // read_len_counts.tsv
        LEN.out.pos_len,                      // pos_len.tsv
        QC.out.multiqc_tsv,                   // multiqc_data.tsv
        READ_STATS.out.read_stats,            // readstats.tsv
        VAR.out.snp_tsv,                      // variants.tsv
        VAR.out.sv_vcf,                       // struct_var.vcf
        plasmid_depth,                        // plasmid coverage histogram if it exists, otherwise empty
        CHIMERA.out.alvis_chimeras_txt,       // chimeras.txt
        CHIMERA.out.alvis_chimeras_tsv,       // chimeras.tsv
        TRUNC.out.st_end_cov,                 // truncations.tsv
        CONSENSUS.out.consensus,              // consensus.txt
        CONSENSUS.out.var_table,              // var_table.tsv
        POSITIONAL_METRICS.out.pos_metrics,   // positional_metrics.tsv
        METHYLATION.out.tsv,                  // methylation.tsv
        METHYLATION.out.primrose_bam,         // 5mc.hifi_reads.bam
        METHYLATION.out.aligned_primrose_bam, // aligned.5mc.hifi_reads.bam
        METHYLATION.out.aligned_primrose_bai, // aligned.5mc.hifi_reads.bam.bai
        SOFTCLIPPED.out.sc_bam,               // softclipped_reads.bam 
        SOFTCLIPPED.out.sc_bai,               // softclipped_reads.bam.bai
        ALIGN.out.bam,                        // all other *.bam
        ALIGN.out.bai,                        // all other *.bai
        FIND_TRANSPOSONS.out.bam,             // transposons.bam
        FIND_TRANSPOSONS.out.bai,             // transposons.bam.bai
        FIND_TRANSPOSONS.out.metrics          // tn_metrics.tsv
    )

    if (params.run_report) {
        REPORT_RUNNER(params.datalakedir, params.unique_id, params.version, sample_report_template_ch, toc_template_ch, rsconnect_tags, params.rsconnect_url, ch_metadata, static_files, footer_ch)
        outputs = outputs.mix(REPORT_RUNNER.out.sample_report_url).collect() // https://connect.sparkds.io link to the report
        FINALIZE_OUTPUTS(REPORT_RUNNER.out.toc_url, params.run_id, params.outdir, outputs)
    }
    else {
        outputs = outputs.collect()
        FINALIZE_OUTPUTS(POSITIONAL_METRICS.out.pos_metrics, params.run_id, params.outdir, outputs)
    }
}
