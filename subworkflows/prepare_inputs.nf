nextflow.enable.dsl=2

include {PREP_COLLECTMETADATA; MERGE_METADATA} from "$projectDir/modules/local/merge_metadata/main"
include { YASCP_INPUTS } from "$projectDir/modules/local/prepere_yascp_inputs/main"

workflow DECONV_INPUTS{
    // This is a function that prpeares the inputs for deconvolution that is utilised twice.
    take:
        cellbender_path
        PREPARE_INPUTS
    main:
        cellbender_path
            .map{row->tuple(row[0], file("${row[1]}".replaceFirst(/.*results/,"${params.outdir}")))}
            .set{ch_experiment_filth5} // this channel is used for task 'split_donor_h5ad'
            
        PREPARE_INPUTS.out.ch_experiment_bam_bai_barcodes.map { experiment, bam, bai, barcodes -> tuple(experiment,
                            bam,
                            bai)}.set{pre_ch_experiment_bam_bai_barcodes}

        cellbender_path.map{row->tuple(row[0], file("${row[1]}".replaceFirst(/.*results/,"${params.outdir}")+'/barcodes.tsv.gz'))}.set{barcodes}

        channel__file_paths_10x= cellbender_path.map{row->tuple(row[0],
                                                    file("${row[1]}".replaceFirst(/.*results/,"${params.outdir}")+'/barcodes.tsv.gz'),
                                                    file("${row[1]}".replaceFirst(/.*results/,"${params.outdir}")+'/features.tsv.gz'),
                                                    file("${row[1]}".replaceFirst(/.*results/,"${params.outdir}")+'/matrix.mtx.gz'))}

        channel__file_paths_10x_single = cellbender_path.map{row->tuple(row[0],
                                                    file("${row[1]}".replaceFirst(/.*results/,"${params.outdir}"))
                                                    )}
        pre_ch_experiment_bam_bai_barcodes.combine(barcodes, by: 0).set{ch_experiment_bam_bai_barcodes}
    emit:
        channel__file_paths_10x
        ch_experiment_bam_bai_barcodes
        ch_experiment_filth5
        channel__file_paths_10x_single
}

workflow PREPARE_INPUTS {
	// this workflow processes the outputs from cellbender to perform the data preparation
    take: channel_input_data_table_pre
    main:

        log.info "... Prepearing inputs based on the 10x folder required for downstream analysis..."

		YASCP_INPUTS(channel_input_data_table_pre)
		channel_input_data_table = YASCP_INPUTS.out.input_file_corectly_formatted
        

        channel_input_data_table
            .splitCsv(header: true, sep: params.input_tables_column_delimiter)
            .map{row->tuple(row.experiment_id, row.n_pooled)}
            .set{ch_experiment_npooled}

        chanel_cr_outs = channel_input_data_table
            .splitCsv(header: true, sep: params.input_tables_column_delimiter)
            .map{row -> file(row.data_path_10x_format)}

        channel_input_data_table
            .splitCsv(header: true, sep: params.input_tables_column_delimiter)
            .map{ row ->
                def bam_file = file("${row.data_path_10x_format}/possorted_genome_bam.bam")
                def bam_file_atac = file("${row.data_path_10x_format}/gex_possorted_bam.bam")
                def barcodes_file_gz = file("${row.data_path_10x_format}/filtered_feature_bc_matrix/barcodes.tsv.gz")
                def barcodes_file = file("${row.data_path_10x_format}/filtered_feature_bc_matrix/barcodes.tsv")
                def selected_bam = bam_file.exists() ? bam_file : bam_file_atac
                
                // Check if barcodes.tsv.gz exists, if not, fall back to barcodes.tsv
                def selected_barcodes_file = barcodes_file_gz.exists() ? barcodes_file_gz : barcodes_file

                return tuple(row.experiment_id, selected_bam, selected_barcodes_file.toString())
            }
            .set{pre_ch_experiment_bam_barcodes}

        channel__file_paths_10x =  channel_input_data_table
            .splitCsv(header: true, sep: params.input_tables_column_delimiter)
            .map{row -> tuple(
            row.experiment_id,
            file("${row.data_path_10x_format}/filtered_feature_bc_matrix/barcodes.tsv.gz"),
            file("${row.data_path_10x_format}/filtered_feature_bc_matrix/features.tsv.gz"),
            file("${row.data_path_10x_format}/filtered_feature_bc_matrix/matrix.mtx.gz")
        )}

        channel__metadata = channel_input_data_table
            .splitCsv(header: true, sep: params.input_tables_column_delimiter)
            .map { row -> 
                def metrics_csv = file("${row.data_path_10x_format}/metrics_summary.csv")
                def fallback_tsv = file("${row.data_path_10x_format}/summary.csv")
                
                tuple(
                    row.experiment_id,
                    metrics_csv.exists() ? metrics_csv : fallback_tsv
                )
            }

        channel_dsb = channel_input_data_table
            .splitCsv(header: true, sep: params.input_tables_column_delimiter)
            .map{row -> tuple(row.experiment_id, file(row.data_path_10x_format+'/raw_feature_bc_matrix'),file(row.data_path_10x_format+'/filtered_feature_bc_matrix'))}

        PREP_COLLECTMETADATA(channel__metadata)
        channel__metadata=MERGE_METADATA(PREP_COLLECTMETADATA.out.metadata.collect())

        channel_input_data_table
            .splitCsv(header: true, sep: params.input_tables_column_delimiter)
            .map{row -> tuple(row.experiment_id, file(row.data_path_10x_format+'/raw_feature_bc_matrix'))}.set{ch_experimentid_paths10x_raw}

        channel_input_data_table
            .splitCsv(header: true, sep: params.input_tables_column_delimiter)
            .map{row -> tuple(row.experiment_id, file(row.data_path_10x_format+'/filtered_feature_bc_matrix'))}.set{ch_experimentid_paths10x_filtered}


        if (params.split_h5ad_per_donor.run) {
            log.info "params.split_h5ad_per_donor.run=true: create pre_ch_experiment_filth5 from params.input_tables_column_delimiter"
            channel_input_data_table
                .splitCsv(header: true, sep: params.input_tables_column_delimiter)
            .map{row->tuple(row.experiment_id, "${row.data_path_10x_format}/filtered_feature_bc_matrix")}
            .set{pre_ch_experiment_filth5} // this channel is used for task 'split_donor_h5ad'
        }
        else {
            // create dummy channel
            log.info "using dummy channel for pre_ch_experiment_filth5, because unused."
            pre_ch_experiment_filth5 = Channel.from("foo").map { foo -> tuple("foo1","foo2") }
        }

        if (params.genotype_input.run_with_genotype_input) {
            log.info "You have selected params.vireo.run_with_genotype_input=true -> will run Vireo with genotype input. Input VCF and list of donors per experiment_id gathered from params.input_n_pooled_table)"
            if (params.genotype_input.subset_genotypes){
                log.info("----We will subset genotypes to the donors listed in the donor_vcf_ids for use in the Deconvolution----")
                // If we subset the genotypes, th ids in the donor_vcf_ids will be used to generate an individual vcf per pool.
                channel_input_data_table
                    .splitCsv(header: true, sep: params.input_tables_column_delimiter)
                    .map{row->tuple(row.experiment_id, params.genotype_input.full_vcf_file,params.genotype_input.full_vcf_file+'.csi', row.donor_vcf_ids)}
                    .set{pre_ch_experiment_donorsvcf_donorslist}
            }else{
                log.info('----We are using full VCF for each donor without subsetting----')
                // When we do not subset the genotypes, a full vcf will be fed in Vireo and soupocell assignments, but the number of pooled individuals will be used to detect the correct number of individuals in pool.
                // in this case we have to be aware that the last number is a number of donors pooled instead of IDs as per above
                channel_input_data_table
                    .splitCsv(header: true, sep: params.input_tables_column_delimiter)
                    .map{row->tuple(row.experiment_id, params.genotype_input.full_vcf_file,params.genotype_input.full_vcf_file+'.csi', row.n_pooled)}
                    .set{pre_ch_experiment_donorsvcf_donorslist}
            }

        } else {
            // Create dummy channel since we are not using VCFs, hence no genotype input
            log.info "We are not using genotypes in this pipeline run"
            pre_ch_experiment_donorsvcf_donorslist = Channel.from("foo").map { foo -> tuple("foo1","foo2","foo3",'foo4') }
        }

        pre_ch_experiment_filth5
            .set{ch_experiment_filth5}

        pre_ch_experiment_bam_barcodes
            .map { a,b,c -> tuple(a, file(b), file("${b}.bai"), file(c))}
            .set {ch_experiment_bam_bai_barcodes}

        pre_ch_experiment_donorsvcf_donorslist
            .set{ch_experiment_donorsvcf_donorslist}

		log.info " ----Prepearing the inputs ---"

    emit:
		ch_experiment_bam_bai_barcodes
		ch_experiment_npooled
		ch_experiment_filth5
		ch_experiment_donorsvcf_donorslist
        ch_experimentid_paths10x_raw
        ch_experimentid_paths10x_filtered
        channel__file_paths_10x
        channel__metadata
        channel_input_data_table
        channel_dsb
        chanel_cr_outs

}
