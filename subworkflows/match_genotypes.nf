// match deconvoluted donors by genotype to a reference panel

include { MATCH_GT_VIREO; GT_MATCH_POOL_IBD } from "$projectDir/modules/nf-core/modules/genotypes/main"
include {COMBINE_MATCHES_IN_EXPECTED_FORMAT} from "$projectDir/modules/nf-core/modules/genotypes/main"
include {Relationships_Between_Infered_Expected; Relationships_Between_Infered_Expected as Relationships_Between_Infered_GT_Matched} from '../modules/nf-core/modules/infered_expected_relationship/main'
include {SUBSET_WORKF} from "$projectDir/modules/nf-core/modules/subset_genotype/main"
include {CONCORDANCE_CALCLULATIONS} from "$projectDir/modules/nf-core/modules/concordance/main"

workflow match_genotypes {
  take:
    vireo_out_sample_donor_vcf
    merged_expected_genotypes
    gt_pool
    gt_math_pool_against_panel_input
    genome
    ch_ref_vcf
    cellsnp_cell_vcfs2
    cell_assignments
  main:

       
    MATCH_GT_VIREO(gt_math_pool_against_panel_input)
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
    
    // «««««««««
    // compare genotypes within a pool (identity by descent)
    GT_MATCH_POOL_IBD(vireo_out_sample_donor_vcf,'Withing_pool','InferedOnly')
    GT_MATCH_POOL_IBD.out.plink_ibd.set{idb_pool}
    // compare genotypes with the expected/matched genotypes to estimate the relationship between matches.
    // #### Note: this code takes in the subset genotypes (either determines as an expected inputs or determined as final matches [dependant on flag] and also the final GT match results to later extract the PI_HAT value)

    SUBSET_WORKF(ch_ref_vcf,gt_matched_samples,'GTMatchedSubset',genome)
    merged_GT_Matched_genotypes = SUBSET_WORKF.out.merged_expected_genotypes

    // Here we need to account for the fact that no GT is expected.
    // donors_in_pools.subscribe { println "donors_in_pools: $it" } //Needs an if statement which determines if pool has no donors expected.
    // Now based on these two files we will enhance the stats file with PiHat values and
    Relationships_Between_Infered_Expected(donors_in_pools,merged_expected_genotypes,gt_pool,'InferedExpected',MATCH_GT_VIREO.out.donor_match_table_with_pool_id,idb_pool)
    outfile_for_final_gt = Relationships_Between_Infered_Expected.out.donor_match_table
    
    Relationships_Between_Infered_GT_Matched(gt_matched_samples,merged_GT_Matched_genotypes,gt_pool,'InferedGTMatched', Relationships_Between_Infered_Expected.out.donor_match_table,idb_pool)
    outfile_for_final_gt = Relationships_Between_Infered_GT_Matched.out.donor_match_table
    Relationships_Between_Infered_Expected.out.done_validation.set{ou1}
    Relationships_Between_Infered_GT_Matched.out.done_validation.set{ou2}
    ou1.mix(ou2).set{ou3}

    // Now we calculate the concordance/discordance scores of each of the cells against the donors.
    // merged_GT_Matched_genotypes.subscribe { println "merged_GT_Matched_genotypes: $it" }
    // merged_expected_genotypes.subscribe { println "merged_expected_genotypes: $it" }
    input3 = merged_GT_Matched_genotypes.join(merged_expected_genotypes, remainder: true)
    input3.map { row -> 
      if(row[1]==null){
        tuple(row[0], "$projectDir/assets/fake_file.fq","$projectDir/assets/fake_file2.fq",row[2],row[[3]]) 
      }else{
         tuple(row[0], row[1],row[2],row[3],row[[4]]) 
      }
    }.set{input32}

    // input32.subscribe { println "input32: $it" }
    input4 = input32.combine(cellsnp_cell_vcfs2, by: 0)
    // input4.subscribe { println "input4: $it" }
    input5 = input4.combine(MATCH_GT_VIREO.out.donor_match_table_with_pool_id, by:0)
    // input5.subscribe { println "input5: $it" }
    input6 = input5.combine(cell_assignments, by:0)
    // input6.subscribe { println "input6: $it" }
    CONCORDANCE_CALCLULATIONS(input6)


  emit:
    pool_id_donor_assignments_csv = MATCH_GT_VIREO.out.pool_id_donor_assignments_csv
    donor_match_table = MATCH_GT_VIREO.out.donor_match_table
    out_finish_val = ou3
    donor_match_table_enhanced = outfile_for_final_gt
}
