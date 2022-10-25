// match deconvoluted donors by genotype to a reference panel

include { MATCH_GT_VIREO; GT_MATCH_POOL_IBD } from '../modules/nf-core/modules/genotypes/main'
include {GT_IBD_MATCH_RELATIONSHIPS} from '../modules/nf-core/modules/ibd_relationsips/main'
include {COMBINE_MATCHES_IN_EXPECTED_FORMAT} from '../modules/nf-core/modules/genotypes/main'
include {Relationships_Between_Infered_Expected; Relationships_Between_Infered_Expected as Relationships_Between_Infered_GT_Matched} from '../modules/nf-core/modules/infered_expected_relationship/main'
include {SUBSET_WORKF} from "$projectDir/modules/nf-core/modules/subset_genotype/main"

workflow match_genotypes {
  take:
    ch_pool_id_vireo_vcf
    merged_expected_genotypes

  main:
    Channel.fromPath(
      params.genotype_input.tsv_donor_panel_vcfs,
      followLinks: true,
      checkIfExists: true
    )
    .splitCsv(header: true, sep: '\t')
    .map { row -> tuple(row.label, file(row.vcf_file_path), file("${row.vcf_file_path}.csi")) }
    .set { ch_ref_vcf }

    MATCH_GT_VIREO(ch_pool_id_vireo_vcf, ch_ref_vcf)
    
    // «««««««««
    // This channel creates an input that contains the GT matched results as per input
    // «««««««««
    COMBINE_MATCHES_IN_EXPECTED_FORMAT(MATCH_GT_VIREO.out.donor_match_table.collect())
    COMBINE_MATCHES_IN_EXPECTED_FORMAT.out.all_Infered_Expected.splitCsv(header: true, sep: '\t').map { row -> tuple(row.experiment_id, row.donor_vcf_ids) }
    .set { gt_matched_samples }
    
    // «««««««««
    // This channel creates an input that contains the expected vcfs as per input
    // «««««««««
    Channel.fromPath(params.input_data_table,      
      followLinks: true,
      checkIfExists: true
    ).splitCsv(header: true, sep: '\t').map { row -> tuple(row.experiment_id, row.donor_vcf_ids) }
    .set { donors_in_pools }
    
    MATCH_GT_VIREO.out.gt_pool.set{vireo_GT_Genotypes}
    // «««««««««
    // compare genotypes within a pool (identity by descent)
    GT_MATCH_POOL_IBD(ch_pool_id_vireo_vcf,'Withing_pool','InferedOnly')
    GT_MATCH_POOL_IBD.out.plink_ibd.set{idb_pool}
    // compare genotypes with the expected/matched genotypes to estimate the relationship between matches.
    // #### Note: this code takes in the subset genotypes (either determines as an expected inputs or determined as final matches [dependant on flag] and also the final GT match results to later extract the PI_HAT value)

    Channel.fromPath(
      params.genotype_input.tsv_donor_panel_vcfs,
      followLinks: true,
      checkIfExists: true
    ).splitCsv(header: true, sep: '\t')
    .map { row -> tuple(row.label, file(row.vcf_file_path), file("${row.vcf_file_path}.csi")) }
    .set { ch_ref_vcf }

    SUBSET_WORKF(ch_ref_vcf,gt_matched_samples)
    merged_GT_Matched_genotypes = SUBSET_WORKF.out.merged_expected_genotypes

    Relationships_Between_Infered_Expected(donors_in_pools,merged_expected_genotypes,vireo_GT_Genotypes,'Expected',MATCH_GT_VIREO.out.donor_match_table_with_pool_id,idb_pool)
    Relationships_Between_Infered_GT_Matched(gt_matched_samples,merged_GT_Matched_genotypes,vireo_GT_Genotypes,'GT_Matched', Relationships_Between_Infered_Expected.out.donor_match_table,idb_pool)
    
    Relationships_Between_Infered_Expected.out.done_validation.set{ou1}
    Relationships_Between_Infered_GT_Matched.out.done_validation.set{ou2}
    ou1.mix(ou2).set{ou3}
    // Now based on these two files we will enhance the stats file with PiHat values.
    // 


  emit:
    pool_id_donor_assignments_csv = MATCH_GT_VIREO.out.pool_id_donor_assignments_csv
    donor_match_table = MATCH_GT_VIREO.out.donor_match_table
    out_finish_val = ou3
}
