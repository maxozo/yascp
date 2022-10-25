process VIREO_GT_FIX_HEADER
{
  tag "${pool_id}"
  publishDir  path: "${params.outdir}/infered_genotypes/${pool_id}/",
        mode: "${params.copy_mode}",
        overwrite: "true"
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      // println "container: /software/hgi/containers/wtsihgi-nf_genotype_match-1.0.sif\n"
      container "/software/hgi/containers/wtsihgi-nf_yascp_htstools-1.1.sif"
  } else {
      container "mercury/wtsihgi-nf_yascp_htstools-1.1"
  }
  //when: params.vireo.run_gtmatch_aposteriori
  // [W::hts_idx_load3] The index file is older than the data file: sorted_Expected_CRD_CMB13102395.bcf.gz.csi
  // [W::vcf_parse_info] INFO 'AD' is not defined in the header, assuming Type=String
  // [W::vcf_parse_info] INFO 'DP' is not defined in the header, assuming Type=String
  // [W::vcf_parse_info] INFO 'OTH' is not defined in the header, assuming Type=String
  // [W::vcf_parse_format] FORMAT 'AD' at 1:632317 is not defined in the header, assuming Type=String
  // [W::vcf_parse_format] FORMAT 'DP' at 1:632317 is not defined in the header, assuming Type=String
  // [W::vcf_parse_format] FORMAT 'PL' at 1:632317 is not defined in the header, assuming Type=String
  // Error: The INFO field is not defined in the header: AD


  // ##INFO=<ID=AD,Number=A,Type=Integer,Description="alternative allele  (variant-by-cell) of reads">
  // ##INFO=<ID=DP,Number=1,Type=Integer,Description="depth UMIs for each variant in each cell">
  //
  // 
  // 
  // 
  // 


  label 'process_low'

  input:
    tuple val(pool_id), path(vireo_gt_vcf)

  output:
    tuple val(pool_id), path("${vireo_fixed_vcf}"), path("${vireo_fixed_vcf}.tbi"), emit: gt_pool

  script:
  sorted_vcf = "${pool_id}_vireo_srt.vcf.gz"
  vireo_fixed_vcf = "${pool_id}_headfix_vireo.vcf.gz"
  """
    # fix header of vireo VCF
    bcftools view -h ${vireo_gt_vcf} > init_head.txt
    sed -i '/^##fileformat=VCFv.*/a ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">' init_head.txt
    head -n -1 init_head.txt > header.txt
    echo '##INFO=<ID=AD,Number=A,Type=Integer,Description="alternative allele  (variant-by-cell) of reads">' >> header.txt
    echo '##INFO=<ID=DP,Number=1,Type=Integer,Description="depth UMIs for each variant in each cell">' >> header.txt
    echo '##INFO=<ID=PL,Number=1,Type=Integer,Description="depth UMIs for each variant in each cell">' >> header.txt
    echo '##INFO=<ID=OTH,Number=1,Type=Integer,Description="????">' >> header.txt
    echo '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="???">' >> header.txt
    echo '##FORMAT=<ID=AD,Number=G,Type=Integer,Description="????n">' >> header.txt
    echo '##FORMAT=<ID=DP,Number=G,Type=Integer,Description="????n">' >> header.txt
    tail -n1 init_head.txt >> header.txt

    # sort VCF file (bcftools sort bails out with an error)
    bcftools view ${vireo_gt_vcf} | \
    awk '\$1 ~ /^#/ {print \$0;next} {print \$0 | "sort -k1,1V -k2,2n"}' | \
    bcftools view -Oz -o ${sorted_vcf} -

    bcftools reheader -h header.txt ${sorted_vcf} | \
    bcftools view -Oz -o pre_${vireo_fixed_vcf}
    tabix -p vcf pre_${vireo_fixed_vcf}
    bcftools +fixref pre_${vireo_fixed_vcf} -Oz -o ${vireo_fixed_vcf} -- -d -f ${params.reference_assembly_fasta_dir}/genome.fa -m flip
    tabix -p vcf ${vireo_fixed_vcf}
  """
}
process REPLACE_GT_DONOR_ID2{

    publishDir  path: "${params.outdir}/deconvolution/vireo_gt_fix/${samplename}/",
          pattern: "GT_replace_*",
          mode: "${params.copy_mode}",
          overwrite: "true"
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "/software/hgi/containers/mercury_scrna_deconvolution_62bd56a-2021-12-15-4d1ec9312485.sif"
        //// container "/software/hgi/containers/mercury_scrna_deconvolution_latest.img"
    } else {
        container "mercury/scrna_deconvolution:62bd56a"
    }

  label 'process_medium'

  input:
    tuple val(samplename), path(gt_donors), path(vireo_sample_summary),path(vireo___exp_sample_summary),path(vireo__donor_ids),path(vcf_file),path(donor_gt_csi)
    
  output:
    tuple val(samplename), path("GT_replace_donor_ids.tsv"), emit: sample_donor_ids
    tuple val(samplename), path("GT_replace_GT_donors.vireo.vcf.gz"), path(vcf_file),path(donor_gt_csi), emit: sample_donor_vcf
    tuple val(samplename), path("GT_replace_GT_donors.vireo.vcf.gz"), emit: infered_vcf
    path("GT_replace_${samplename}.sample_summary.txt"), emit: sample_summary_tsv
    path("GT_replace_${samplename}__exp.sample_summary.txt"), emit: sample__exp_summary_tsv
    path("GT_replace_${samplename}_assignments.tsv"), emit: assignments

  script:

    in=""

    """
      bcftools query -l GT_donors.vireo.vcf.gz > donors_in_vcf.tsv
      replace_donors.py -id ${samplename} ${in} --input_file ${params.input_data_table}
      bcftools reheader --samples replacement_assignments.tsv -o GT_replace_GT_donors.vireo.vcf.gz GT_donors.vireo.vcf.gz 
    """
}

process REPLACE_GT_DONOR_ID{


    publishDir  path: "${params.outdir}/deconvolution/vireo_gt_fix/${samplename}/",
          pattern: "GT_replace_*",
          mode: "${params.copy_mode}",
          overwrite: "true"
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "/software/hgi/containers/mercury_scrna_deconvolution_62bd56a-2021-12-15-4d1ec9312485.sif"
        //// container "/software/hgi/containers/mercury_scrna_deconvolution_latest.img"
    } else {
        container "mercury/scrna_deconvolution:62bd56a"
    }

  label 'process_medium'

  input:
    tuple val(samplename), path(gt_donors), path(vireo_sample_summary),path(vireo___exp_sample_summary),path(vireo__donor_ids),path(vcf_file),path(donor_gt_csi)
    path(gt_match_results)
  output:
    path("test.out", emit: replacements)
    tuple val(samplename), path("GT_replace_donor_ids.tsv"), emit: sample_donor_ids
    tuple val(samplename), path("GT_replace_GT_donors.vireo.vcf.gz"), path(vcf_file),path(donor_gt_csi), emit: sample_donor_vcf
    path("GT_replace_${samplename}.sample_summary.txt"), emit: sample_summary_tsv
    path("GT_replace_${samplename}__exp.sample_summary.txt"), emit: sample__exp_summary_tsv
    path("GT_replace_${samplename}_assignments.tsv"), emit: assignments

  script:
    if(params.genotype_phenotype_mapping_file==''){
      in=""
    }else if (params.use_phenotype_ids_for_gt_match){
      in="--genotype_phenotype_mapping ${params.genotype_phenotype_mapping_file}"
      // in=""
    }else{
      in=""
    }

    """
      echo ${samplename} > test.out
      gunzip -k -d --force GT_donors.vireo.vcf.gz
      replace_donors.py -id ${samplename} ${in} --input_file ${params.input_data_table}
      bgzip GT_replace_GT_donors.vireo.vcf
    """
}

process GT_MATCH_POOL_IBD
{
  tag "${pool_id}_ibd"

  publishDir  path: "${params.outdir}/gtmatch/${pool_id}",
          mode: "${params.copy_mode}",
          overwrite: "true"

  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "/software/hgi/containers/wtsihgi-nf_yascp_plink1-1.0.img"
  } else {
      container "mercury/wtsihgi-nf_yascp_plink1-1.0"
  }

  label 'process_tiny'

  input:
    tuple val(pool_id), path(vireo_gt_vcf)
    val(mode)
    val(mode2)

  output:
    tuple val(pool_id),path("*_${pool_id}.genome*"), emit:plink_ibd

  script:
    """
      plink --vcf ${vireo_gt_vcf} --genome unbounded --const-fid dummy --out ${mode2}_${mode}_${pool_id}
    """
}

process GT_MATCH_POOL_AGAINST_PANEL
{
  tag "${pool_id}_vs_${panel_id}"

  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      // println "container: /software/hgi/containers/wtsihgi-nf_genotype_match-1.0.sif\n"
      container "/software/hgi/containers/wtsihgi-nf_yascp_htstools-1.1.sif"
  } else {
      container "mercury/wtsihgi-nf_yascp_htstools-1.1"
  }

  label 'process_long'
  //when: params.vireo.run_gtmatch_aposteriori

  input:
    tuple val(pool_id), path(vireo_gt_vcf), path(vireo_gt_tbi), val(panel_id), path(ref_gt_vcf), path(ref_gt_csi)

  output:
    tuple val(pool_panel_id), path("${gt_check_output_txt}"), emit:gtcheck_results

  script:
  pool_panel_id = "pool_${pool_id}_panel_${panel_id}"
  panel_filnam = "${ref_gt_vcf}" - (~/\.[bv]cf(\.gz)?$/)
  gt_check_output_txt = "${pool_id}_gtcheck_${panel_filnam}.txt"
  """
    bcftools gtcheck --no-HWE-prob -g ${ref_gt_vcf} ${vireo_gt_vcf} > ${gt_check_output_txt}
  """
}

process ASSIGN_DONOR_FROM_PANEL
{
  // sum gtcheck discrepancy scores from multiple ouputput files of the same panel
  tag "${pool_panel_id}"
  publishDir  path: "${params.outdir}/gtmatch/${pool_id}",
          pattern: "*.csv",
          mode: "${params.copy_mode}",
          overwrite: "true"
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      // println "container: /software/hgi/containers/wtsihgi-nf_genotype_match-1.0.sif\n"
      container "/software/hgi/containers/wtsihgi-nf_genotype_match-1.0.sif"
  } else {
      container "mercury/wtsihgi-nf_genotype_match-1.0"
  }

  input:
    tuple val(pool_panel_id), path(gtcheck_output_files)

  output:
    tuple val(pool_id), path("${assignment_table_out}"), emit: gtcheck_assignments
    path("${score_table_out}", emit: gtcheck_scores)

  label 'process_low'

  script:
  (_, pool_id) = ("${pool_panel_id}" =~ /^pool_(\S+)_panel_/)[0]
  score_table_out = "${pool_panel_id}_gtcheck_score_table.csv"
  assignment_table_out = "${pool_panel_id}_gtcheck_donor_assignments.csv"

  """
    gtcheck_assign.py ${pool_panel_id} ${gtcheck_output_files}
  """
}

process ASSIGN_DONOR_OVERALL
{
  // decide final donor assignment across different panels from per-panel donor assignments
  tag "${pool_panel_id}"

  publishDir  path: "${params.outdir}/gtmatch/${pool_id}",
          pattern: "*.csv",
          mode: "${params.copy_mode}",
          overwrite: "true"

  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      // println "container: /software/hgi/containers/wtsihgi-nf_genotype_match-1.0.sif\n"
      container "/software/hgi/containers/mercury_scrna_deconvolution_62bd56a-2021-12-15-4d1ec9312485.sif"
  } else {
      container "mercury/wtsihgi-nf_genotype_match-1.0"
  }

  input:
    tuple val(pool_id), path(gtcheck_assign_files)

  output:
    tuple val(pool_id), path("${donor_assignment_file}"), emit: donor_assignments
    path(stats_assignment_table_out), emit: donor_match_table
    tuple val(pool_id),path(stats_assignment_table_out), emit: donor_match_table_with_pool_id
    path("*.csv")

  label 'process_tiny'

  script:
  donor_assignment_file = "${pool_id}_gt_donor_assignments.csv"
  stats_assignment_table_out = "stats_${pool_id}_gt_donor_assignments.csv"
  """
    gtcheck_assign_summary.py ${donor_assignment_file} ${gtcheck_assign_files}
  """
}




process REPLACE_GT_ASSIGNMENTS_WITH_PHENOTYPE{
  label 'process_low'
  publishDir  path: "${params.outdir}/gtmatch/",
          pattern: "*_assignments.csv",
          mode: "${params.copy_mode}",
          overwrite: "true"

  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "/software/hgi/containers/mercury_scrna_deconvolution_62bd56a-2021-12-15-4d1ec9312485.sif"
      //// container "/software/hgi/containers/mercury_scrna_deconvolution_latest.img"
  } else {
      container "mercury/scrna_deconvolution:62bd56a"
  }

  input:
    path(gt_match_results)

  output:
    path(gt_match_results, emit: donor_match_table)

  script:
    """
      perform_replacement.py --genotype_phenotype_mapping ${params.genotype_phenotype_mapping_file} --assignemts ${gt_match_results}

    """

}

process ENHANCE_STATS_FILE{

  tag "${pool_id}"

  publishDir  path: "${params.outdir}/gtmatch/${pool_id}",
        mode: "${params.copy_mode}",
        overwrite: "true"


  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "/software/hgi/containers/mercury_scrna_deconvolution_62bd56a-2021-12-15-4d1ec9312485.sif"
      //// container "/software/hgi/containers/mercury_scrna_deconvolution_latest.img"
  } else {
      container "mercury/scrna_deconvolution:62bd56a"
  }

  label 'process_small'


  input:
    tuple val(pool_id),path(ibd_table),path(stats_table),val(expected_ids),path(withinn_pool_ibd)
    val(condition)


  output:
    
    tuple val(pool_id),path("PiHAT_Stats_File*"), emit: stats_table_PiHat_enhanced optional true
    path ('Max_PiHAT_For_Expected*'), emit: max_PiHAT_For_Expected optional true
    path('Done.tmp'), emit: done_validation
  script:
    if (params.extra_sample_metadata==''){
      md_inp = ""
    }else{
      md_inp = "-md ${params.extra_sample_metadata}"
    }
    """
      add_PiHat_to_GT_match.py -mt ${stats_table} -ph ${ibd_table} -c ${condition} -e ${expected_ids} -m ${params.genotype_phenotype_mapping_file} -id ${pool_id} -wpi ${withinn_pool_ibd} ${md_inp} || echo 'we dont have expected samples in this cohort'
      echo 'Done' > Done.tmp
    """

}

process ENHANCE_VIREO_METADATA_WITH_DONOR{
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "/software/hgi/containers/mercury_scrna_deconvolution_62bd56a-2021-12-15-4d1ec9312485.sif"
        //// container "/software/hgi/containers/mercury_scrna_deconvolution_latest.img"
    } else {
        container "mercury/scrna_deconvolution:62bd56a"
    }
  label 'process_small'


  input:
    path(extra_sample_metadata)
    path(donor_n_cells)
    path(out_gt)

  output:
    path('replaced_vireo_exp__donor_n_cells_out.tsv'), emit: replaced_vireo_exp__donor_n_cells_out

  script:
    """
      enhance_vireo_with_metadata.py --Extra_Metadata_Donors ${extra_sample_metadata} --vireo_data ${donor_n_cells}
    """
}



process COMBINE_MATCHES_IN_EXPECTED_FORMAT{
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "/software/hgi/containers/mercury_scrna_deconvolution_62bd56a-2021-12-15-4d1ec9312485.sif"
        //// container "/software/hgi/containers/mercury_scrna_deconvolution_latest.img"
    } else {
        container "mercury/scrna_deconvolution:62bd56a"
    }
  label 'process_small'

  input:
    path(stats_files)

  output:
    path('All_Infered_Expected.csv'), emit: all_Infered_Expected

  script:
    if (params.cohorts_to_drop_from_GT_Relatednes_check==''){
      md_in = ""
    }else{
      md_in = "-dr ${params.cohorts_to_drop_from_GT_Relatednes_check}"
    }
    """
      combine_all_GTmatched_in_expected_format.py --files "${stats_files}" ${md_in}
    """
}


workflow MATCH_GT_VIREO {
  take:
    ch_pool_id_vireo_vcf
    ch_ref_vcf

  main:
    // ch_ref_vcf.subscribe { println "match_genotypes: ch_ref_vcf = ${it}" }

    // VIREO header causes problems downstream
    VIREO_GT_FIX_HEADER(ch_pool_id_vireo_vcf)
    VIREO_GT_FIX_HEADER.out.gt_pool
      .combine(ch_ref_vcf)
      .set { ch_gt_pool_ref_vcf }
    // ch_gt_pool_ref_vcf.subscribe { println "match_genotypes: ch_gt_pool_ref_vcf = ${it}\n" }

    // now match genotypes against a panels
    GT_MATCH_POOL_AGAINST_PANEL(ch_gt_pool_ref_vcf)

    // group by panel id
    GT_MATCH_POOL_AGAINST_PANEL.out.gtcheck_results
      .groupTuple()
      .set { gt_check_by_panel }
    

    ASSIGN_DONOR_FROM_PANEL(gt_check_by_panel)
    ASSIGN_DONOR_FROM_PANEL.out.gtcheck_assignments
      .groupTuple()
      .set{ ch_donor_assign_panel }
    // ch_donor_assign_panel.subscribe {println "ASSIGN_DONOR_OVERALL: ch_donor_assign_panel = ${it}\n"}

    ASSIGN_DONOR_OVERALL(ch_donor_assign_panel)

  emit:
    pool_id_donor_assignments_csv = ASSIGN_DONOR_OVERALL.out.donor_assignments
    donor_match_table = ASSIGN_DONOR_OVERALL.out.donor_match_table
    donor_match_table_with_pool_id = ASSIGN_DONOR_OVERALL.out.donor_match_table_with_pool_id
    gt_pool = VIREO_GT_FIX_HEADER.out.gt_pool
}
