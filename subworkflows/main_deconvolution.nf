nextflow.enable.dsl=2

// main deconvolution modules, common to all input modes:

include { CELLSNP } from "$projectDir/modules/nf-core/modules/cellsnp/main"
include { SUBSET_GENOTYPE } from "$projectDir/modules/nf-core/modules/subset_genotype/main"
include { VIREO } from "$projectDir/modules/nf-core/modules/vireo/main"
include { GUZIP_VCF } from "$projectDir/modules/nf-core/modules/guzip_vcf/main"
include { SOUPORCELL } from "$projectDir/modules/nf-core/modules/souporcell/main"
include { SPLIT_DONOR_H5AD } from "$projectDir/modules/nf-core/modules/split_donor_h5ad/main"
include { PLOT_DONOR_CELLS } from "$projectDir/modules/nf-core/modules/plot_donor_cells/main"
include {MULTIPLET} from "$projectDir/modules/nf-core/modules/multiplet/main"
include {SOUPORCELL_VS_VIREO} from "$projectDir/modules/nf-core/modules/plot_souporcell_vs_vireo/main"
include { REPLACE_GT_ASSIGNMENTS_WITH_PHENOTYPE; ENHANCE_VIREO_METADATA_WITH_DONOR } from "$projectDir/modules/nf-core/modules/genotypes/main"
include { match_genotypes } from './match_genotypes'
include {REPLACE_GT_DONOR_ID } from "$projectDir/modules/nf-core/modules/genotypes/main"
include {SUBSET_WORKF} from "$projectDir/modules/nf-core/modules/subset_genotype/main"
include {REPLACE_GT_DONOR_ID2 } from "$projectDir/modules/nf-core/modules/genotypes/main"


workflow  main_deconvolution {

    take:
		ch_experiment_bam_bai_barcodes
		ch_experiment_npooled
		ch_experiment_filth5
		ch_experiment_donorsvcf_donorslist
        channel__file_paths_10x

    main:
		log.info "#### running DECONVOLUTION workflow #####"
        channel_input_data_table = Channel.fromPath(params.input_data_table, followLinks: true, checkIfExists: true)
        
        if (params.run_with_genotype_input) {
            if (params.genotype_input.subset_genotypes){
                log.info "---We are subsetting genotypes----"
                // We have to produce a single vcf file for each individual pool.
                // Therefore we create 2 channels:
                // 1) All the expected vcf ids listed in the donor table
                Channel.fromPath(params.input_data_table,      
                followLinks: true,
                checkIfExists: true
                ).splitCsv(header: true, sep: '\t').map { row -> tuple(row.experiment_id, row.donor_vcf_ids) }
                .set { donors_in_pools }
                
                // 2) All the vcfs provided to us. 
                Channel.fromPath(
                params.genotype_input.tsv_donor_panel_vcfs,
                followLinks: true,
                checkIfExists: true
                ).splitCsv(header: true, sep: '\t')
                .map { row -> tuple(row.label, file(row.vcf_file_path), file("${row.vcf_file_path}.csi")) }
                .set { ch_ref_vcf }

                // This will subsequently result in a joint vcf file for all the cohorts listed for each of the pools that can be used in VIREO and/or GT matching algorythm.
                SUBSET_WORKF(donors_in_pools,ch_ref_vcf)
                merged_expected_genotypes = SUBSET_WORKF.out.merged_expected_genotypes
                // SUBSET_GENOTYPE(ch_experiment_donorsvcf_donorslist.map { experiment, donorsvcf,donortbi, donorslist -> tuple(experiment,
                //                 file(donorsvcf),file(donortbi),
                //                 donorslist)})
            }
        }


        
        ch_poolid_donor_assignment = Channel.empty()
        ch_experiment_donorsvcf_donorslist.map { experiment, donorsvcf, donorstbi,donorslist -> tuple(experiment, donorslist.replaceAll(/,/, " ").replaceAll(/"/, ""))}.set{donors_in_lane}

        if (params.existing_cellsnp != ''){
            expl1 = Channel.fromPath( "${params.existing_cellsnp}/*/*.vcf.gz")
            expl1.map{row -> tuple("${row[-2]}".replaceAll('cellsnp_',''), "${row}".replaceAll('/cellSNP.base.vcf.gz',''))}.set{cellsnp_output_dir}
            
        }else{
            log.info('Running CELLSNP')
            channel_input_data_table
                .splitCsv(header: true, sep: params.input_tables_column_delimiter)
                .map{row->tuple(row.experiment_id, "${row.data_path_10x_format}/possorted_genome_bam.bam" ,row.data_path_10x_format+'/filtered_feature_bc_matrix/barcodes.tsv.gz')}
                .set{pre_ch_experiment_bam_barcodes}
            pre_ch_experiment_bam_barcodes
                .map { a,b,c -> tuple(a, file(b), file("${b}.bai"), file(c))}
                .set {ch_experiment_bam_bai_barcodes}

            CELLSNP(ch_experiment_bam_bai_barcodes,
                Channel.fromPath(params.cellsnp.vcf_candidate_snps).collect())
            cellsnp_output_dir = CELLSNP.out.cellsnp_output_dir
        }

        MULTIPLET(
            params.output_dir,
            channel__file_paths_10x,
            params.sample_qc.cell_filters.filter_multiplets.expected_multiplet_rate,
            params.sample_qc.cell_filters.filter_multiplets.n_simulated_multiplet,
            params.sample_qc.cell_filters.filter_multiplets.multiplet_threshold_method,
            params.sample_qc.cell_filters.filter_multiplets.scale_log10
        )


        // Here we run Vireo software to perform the donor deconvolution. Note that we have coded the pipeline to be capable in using
        // the full genotypes as an input and also subset to the individuals provided as an input in the donor_vcf_ids column. The
        // VIREO:
        if (params.run_with_genotype_input) {
            log.info "---running Vireo with genotype input----"
            // for each experiment_id to deconvolute, subset donors vcf to its donors and subset genomic regions.
            // Here we subset the genotypes. This happens if the input.nf contains subset_genotypes = true
            log.info "---We are using subset genotypes running Vireo----"
            cellsnp_output_dir.combine(ch_experiment_npooled, by: 0)
                .combine(merged_expected_genotypes, by: 0).set{full_vcf}
        }
        // Vireo without genotype input:
        else {
            // Vireo can also be run without the genotypes, and the performance is equally good then running with.
            // here we run it without the genotypes and the donors are labeled as donor0, donor1 etc, dependant on the number of donors set in the input file n_pooled column
            log.info "-----running Vireo without genotype input----"
            cellsnp_output_dir.combine(ch_experiment_npooled, by: 0).set{full_vcf}
            full_vcf.map {experiment, cellsnp, npooled -> tuple(experiment, cellsnp, npooled,[],[])}.set{full_vcf}
        }

        // When all the channels are prpeared accordingly we exacute the vireo with the prpeared channel.
        full_vcf.filter { experiment, cellsnp, npooled, t,ti -> npooled != '1' }.set{full_vcf2}
        full_vcf.filter { experiment, cellsnp, npooled, t,ti -> npooled == '1' }.set{not_deconvoluted}


        VIREO(full_vcf2)
        REPLACE_GT_DONOR_ID2(VIREO.out.all_required_data)
        vireo_out_sample_donor_vcf = REPLACE_GT_DONOR_ID2.out.infered_vcf
        vireo_out_sample_summary_tsv = REPLACE_GT_DONOR_ID2.out.sample_summary_tsv
        vireo_out_sample__exp_summary_tsv = REPLACE_GT_DONOR_ID2.out.sample__exp_summary_tsv
        vireo_out_sample_donor_ids = REPLACE_GT_DONOR_ID2.out.sample_donor_ids

        if (params.run_with_genotype_input) {
            match_genotypes(vireo_out_sample_donor_vcf,merged_expected_genotypes)
        }


        if (params.souporcell.run){
            // YASCP pipeline is also capable in running SOUPORCELL instead of VIREO. If activated SOUPORCELL will be used.
            // yascp currently doesnt have an option to take souporcell assignments as downstream instead of vireo but this will be added shortly.
            // This runs the Souporcell Preprocessing

            if (params.souporcell.use_raw_barcodes) {
                // read raw cellranger barcodes per pool for souporcell
                log.info """Here use the raw barcodes"""
            } else if (! params.souporcell.use_raw_barcodes) {
                ch_experiment_bam_bai_barcodes
                    .combine(ch_experiment_npooled, by: 0)
                    .set{ch_experiment_bam_bai_barcodes_npooled}
            }

            // This runs the Souporcell
            // Similarly to the VIREO Soupocell can be run with and without genotypes and folowing prpeares the inputs accordingly to each option.
            // Soupocell cant digest a gz file gence we extract the data.
            if (params.run_with_genotype_input) {

                if (params.genotype_input.subset_genotypes){
                    // this will run the soupocell with the subset genotypes. This happens if the input.nf contains subset_genotypes = true
                    log.info "---We are using subset genotypes running Suporcell----"
                    GUZIP_VCF(SUBSET_GENOTYPE.out.samplename_subsetvcf)
                    // ch_experiment_bam_bai_barcodes_npooled.combine(SUBSET_GENOTYPE.out.samplename_subsetvcf, by: 0).set{full_vcf}
                    ch_experiment_bam_bai_barcodes_npooled.combine(GUZIP_VCF.out.souporcell_vcf, by: 0).set{full_vcf}

                }else{
                    log.info "---We are using a full genotype input for Suporcell----"
                    // this however currently doesnt work and the individuals have to be provided.
                    // here just add the full vcf path to each of the ch_experiment_bam_bai_barcodes_npooled
                    GUZIP_VCF(tuple('full_vcf', file(params.genotype_input.full_vcf_file)))
                    GUZIP_VCF.out.souporcell_vcf.map { sample, vcf -> vcf }.set{vcf_file}
                    ch_experiment_bam_bai_barcodes_npooled.map {
                        samplename, bam_file, bai_file, barcodes_tsv_gz, souporcell_n_clusters ->
                        tuple(samplename, bam_file, bai_file, barcodes_tsv_gz, souporcell_n_clusters)
                        }.set{full_vcf}
                    full_vcf.combine(vcf_file).set{full_vcf}

                }
                full_vcf.combine(donors_in_lane, by: 0).set{full_vcf}
            }
            else{
                // Soupocell can also run without the genotypes. This will prpeare channels to run it withoutt.
                log.info "-----running Suporcell without genotype input----"
                // here make add an empty entry [] to the ch_experiment_bam_bai_barcodes_npooled
                ch_experiment_bam_bai_barcodes_npooled.map {
                    samplename, bam_file, bai_file, barcodes_tsv_gz, souporcell_n_clusters ->
                    tuple(samplename, bam_file, bai_file, barcodes_tsv_gz, souporcell_n_clusters,[],[])
                    }.set{full_vcf}
            }

            // When all the channels are prepeared then we can run the soupocell accordingly.
            SOUPORCELL(full_vcf,
                Channel.fromPath(params.souporcell.reference_fasta).collect())
            // Regardless if the Soupocell is run with or without genotypes we still need to match the donor ids with the cluster ids since this does not happen automatically is Soupocell.
        }

        ch_experiment_bam_bai_barcodes
            .map { samplename, bam_file, bai_file, barcodes_tsv_gz -> tuple(samplename, file(bam_file)) }
            .combine(vireo_out_sample_donor_ids, by: 0 )
            .set { ch_experiment_bam_vireo_donor_ids }

        // If sample is not deconvoluted we will use scrublet to detect the doublets and remove them.
        not_deconvoluted.map{ experiment, donorsvcf, npooled,t,t2 -> tuple(experiment, 'None')}.set{not_deconvoluted2}
        file_cellmetadata = MULTIPLET.out.file__cellmetadata
        scrublet_paths = MULTIPLET.out.scrublet_paths
        split_channel = vireo_out_sample_donor_ids.combine(ch_experiment_filth5, by: 0)
        split_channel2 = not_deconvoluted2.combine(ch_experiment_filth5, by: 0)
        // combining these 2 channels in one
        split_channel3 = split_channel.mix(split_channel2)
        // adding the scrublet paths to the channel.
        split_channel4 = split_channel3.combine(scrublet_paths, by: 0)

        split_channel5 = split_channel4.map{
            val_sample, val_donor_ids_tsv, val_filtered_matrix_h5, path_scrublet ->
            [  val_sample,file(val_donor_ids_tsv),file(val_filtered_matrix_h5),path_scrublet,params.outdir]
        }

        SPLIT_DONOR_H5AD(split_channel5)

        // collect file paths to h5ad files in tsv tables:
        SPLIT_DONOR_H5AD.out.donors_h5ad_tsv
        .collectFile(name: "donors_h5ad.tsv",
                newLine: false, sort: true,
                seed: "experiment_id\tdonor\th5ad_filepath\n",
                storeDir:params.outdir+'/deconvolution/filepaths')

        // paste experiment_id and donor ID columns with __ separator
        SPLIT_DONOR_H5AD.out.exp__donors_h5ad_tsv
        .collectFile(name: "exp__donors_h5ad.tsv",
                newLine: false, sort: true,
                seed: "experiment_id\th5ad_filepath\n",
                storeDir:params.outdir+'/deconvolution/filepaths')

        SPLIT_DONOR_H5AD.out.donors_h5ad_assigned_tsv
        .collectFile(name: "donors_h5ad_assigned.tsv",
                newLine: false, sort: true,
                seed: "experiment_id\tdonor\th5ad_filepath\n",
                storeDir:params.outdir+'/deconvolution/filepaths')

        // paste experiment_id and donor ID columns with __ separator
        out_h5ad =SPLIT_DONOR_H5AD.out.exp__donors_h5ad_assigned_tsv
        .collectFile(name: "exp__donors_h5ad_assigned.tsv",
                newLine: false, sort: true,
                seed: "experiment_id\th5ad_filepath\n",
                storeDir:params.outdir+'/deconvolution/filepaths')

        SPLIT_DONOR_H5AD.out.h5ad_tsv
        .collectFile(name: "cellranger_as_h5ad.tsv",
                newLine: true, sort: true, // only one line in each file to collate, without ending new line character, so add it here.
                seed: "experiment_id\th5ad_filepath", // don't need \n here since newLine: true
                storeDir:params.outdir+'/deconvolution/filepaths')

        // all vireo() outputs collected -> plot_donor_ncells():
        vireo_out_sample_summary_tsv
        .collectFile(name: "vireo_donor_n_cells.tsv",
                newLine: false, sort: true,
                seed: "experiment_id\tdonor\tn_cells\n",
                storeDir:params.outdir+'/deconvolution/filepaths')
        .set{ch_vireo_donor_n_cells_tsv} // donor column: donor0, .., donorx, doublet, unassigned

        out_split = SPLIT_DONOR_H5AD.out.donor_n_cells

        // paste experiment_id and donor ID columns with __ separator
        vireo_out_sample__exp_summary_tsv = out_split
            .collectFile(name: "vireo_exp__donor_n_cells.tsv",
                    newLine: false, sort: true,
                    seed: "experiment_id\tn_cells\n",
                    storeDir:params.outdir+'/deconvolution/filepaths')


        if (params.run_with_genotype_input & params.genotype_input.posterior_assignment) {
            if (params.extra_sample_metadata!='' & params.add_donor_metadata){
                // Here we have sample level metadata but we have chosen to keep the donor ids.
                // in this scenario we enhance the donor level vireo metadata file to add the donor metadata to the h5ads and eventually to the donor and teanche report
                ENHANCE_VIREO_METADATA_WITH_DONOR(params.extra_sample_metadata,vireo_out_sample__exp_summary_tsv,REPLACE_GT_DONOR_ID.out.assignments.collect())
                vireo_out_sample__exp_summary_tsv = ENHANCE_VIREO_METADATA_WITH_DONOR.out.replaced_vireo_exp__donor_n_cells_out
            }
        }


        PLOT_DONOR_CELLS(ch_vireo_donor_n_cells_tsv)

        if (params.souporcell.run && params.vireo.run) {
            // This is an orringinal comparison of Vireo and soupocell, this needs to be cleaned up.
            SOUPORCELL_VS_VIREO(
                vireo_out_sample_donor_ids
                .combine(SOUPORCELL.out.souporcell_output_files.map {a,b,c,d -> tuple(a,b)},
                    by: 0))
        }

    emit:
        out_h5ad
        vireo_out_sample__exp_summary_tsv
        vireo_out_sample_donor_vcf
        ch_poolid_csv_donor_assignments = ch_poolid_donor_assignment
        sample_possorted_bam_vireo_donor_ids = ch_experiment_bam_vireo_donor_ids
}
