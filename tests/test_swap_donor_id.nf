nextflow.enable.dsl=2

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

include { SWAP_DONOR_IDS; SUBSET_GENOTYPE} from '../modules/nf-core/modules/subset_genotype/main'

workflow TEST_SWAP_DONOR_IDS {
    Channel.fromPath(params.genotype_input.sample_id_swap_table)
      .set { ch_id_swap_table }

    Channel.fromPath(params.input_data_table, followLinks: true, checkIfExists: true)
    .set {ch_input_table}

    ch_input_table.splitCsv(header: true, sep: params.input_tables_column_delimiter)
      .map{row->tuple(row.experiment_id, params.genotype_input.full_vcf_file, row.donor_vcf_ids)}
      .set{ ch_exp_vcf_donorslist }

    ch_exp_vcf_donorslist.subscribe onNext: { println "test_swap_donor_id: $it" }, onComplete: { println 'Done' }
    SWAP_DONOR_IDS(ch_exp_vcf_donorslist, ch_id_swap_table)

    SWAP_DONOR_IDS.out
    .set {ch_exp_vcf_donorids}

    ch_exp_vcf_donorids.subscribe onNext: { println "ch_exp_vcf_donorids: $it" }, onComplete: { println 'Done' }
    SUBSET_GENOTYPE(ch_exp_vcf_donorids)
}
