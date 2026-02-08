
include { SUBSET_GENOTYPE } from '../subset_genotype/main'
include {JOIN_CHROMOSOMES;JOIN_STUDIES_MERGE} from '../subset_genotype/main'
include {JOIN_STUDIES_MERGE as JOIN_INFERED_EXPECTED_MERGE} from '../subset_genotype/main'
include { GT_MATCH_POOL_IBD as GT_MATCH_INFERED_EXPECTED; ENHANCE_STATS_FILE } from '../genotypes/main'

workflow RELATIONSHIPS_BETWEEN_INFERED_EXPECTED {
    take:
      donors_in_pools
      merged_expected_genotypes
      vireo_GT_Genotypes
      mode
      donor_match_table
      idb_pool
    main:
      
      vireo_GT_Genotypes.map { row -> tuple("${row[0]}", row[1],row[2]) }.set { vireo_GT_Genotypes }
      // vireo_GT_Genotypes.view()
      merged_expected_genotypes.map { row -> tuple("${row[0]}", row[1],row[2]) }.set { merged_expected_genotypes }
      idb_pool.map{ row -> tuple("${row[0]}", row[1]) }.set { idb_pool }
      // merged_expected_genotypes.view()
      vireo_GT_Genotypes.mix(merged_expected_genotypes).set{all_Expected_and_infeared}
      
      all_Expected_and_infeared.unique().groupTuple(by: 0).set{grouped}
      // grouped.view()
      JOIN_INFERED_EXPECTED_MERGE(grouped,'InferedMerge',mode)
      // have to do this for each of the pools too. 
      JOIN_INFERED_EXPECTED_MERGE.out.merged_expected_genotypes.map { row -> tuple(row[0], row[1]) }
      .set { sample_name_vcf_no_csi }

      
      GT_MATCH_INFERED_EXPECTED(sample_name_vcf_no_csi,'Expected_Infered',mode)
    
      GT_MATCH_INFERED_EXPECTED.out.plink_ibd.combine(donor_match_table, by: 0).set{ibd_genome_mix}
      ibd_genome_mix.combine(donors_in_pools, by: 0).set{ibd_genome_expected_mix}
      ibd_genome_expected_mix.combine(idb_pool, by: 0).set{ibd_genome_expected_mix2}

      
      ENHANCE_STATS_FILE(ibd_genome_expected_mix2,mode)
     
    emit:
        donor_match_table = ENHANCE_STATS_FILE.out.stats_table_PiHat_enhanced
        done_validation = ENHANCE_STATS_FILE.out.done_validation
}
