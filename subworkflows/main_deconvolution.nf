nextflow.enable.dsl=2

// main deconvolution modules, common to all input modes:

include { CELLSNP } from "$projectDir/modules/nf-core/modules/cellsnp/main"
include { SUBSET_GENOTYPE; SUBSET_GENOTYPE_SWAP_DONOR_ID } from '../modules/nf-core/modules/subset_genotype/main'
include { VIREO } from '../modules/nf-core/modules/vireo/main'
include { GUZIP_VCF } from '../modules/nf-core/modules/guzip_vcf/main'
include { SOUPORCELL } from '../modules/nf-core/modules/souporcell/main'
include { SPLIT_DONOR_H5AD } from '../modules/nf-core/modules/split_donor_h5ad/main'
include { PLOT_DONOR_CELLS } from '../modules/nf-core/modules/plot_donor_cells/main'
include {MULTIPLET} from "../modules/nf-core/modules/multiplet/main"


workflow  main_deconvolution {

    take:
		ch_experiment_bam_bai_barcodes
		ch_experiment_npooled
		ch_experiment_filth5
		ch_experiment_donorsvcf_donorslist
    ch_donorid_swap_table
    ch_vcf_dir
    channel__file_paths_10x

    main:
		log.info "#### running workflow main_deconvolution() ..."


        if (params.run_with_genotype_input) {
            if (params.genotype_input.subset_genotypes){
                log.info "---We are subsetting genotypes----"
                ch_experiment_donorsvcf_donorslist
                .map {
                  experiment, donorsvcf, donorslist -> tuple(
                  experiment, file(donorsvcf), donorslist
                  )
                }
                .set { ch_exp_vcf_donors }
                if (params.genotype_input.swap_donor_ids) {
                  SUBSET_GENOTYPE_SWAP_DONOR_ID(ch_exp_vcf_donors, ch_vcf_dir, ch_donorid_swap_table)
                } else {
                  SUBSET_GENOTYPE(ch_exp_vcf_donors)
                }
            }

        }

        ch_experiment_donorsvcf_donorslist.map { experiment, donorsvcf, donorslist -> tuple(experiment, donorslist.replaceAll(/,/, " ").replaceAll(/"/, ""))}.set{donors_in_lane}

        CELLSNP(ch_experiment_bam_bai_barcodes,
            Channel.fromPath(params.cellsnp.vcf_candidate_snps).collect())

        MULTIPLET(
            params.output_dir,
            channel__file_paths_10x,
            params.sample_qc.cell_filters.filter_multiplets.expected_multiplet_rate,
            params.sample_qc.cell_filters.filter_multiplets.n_simulated_multiplet,
            params.sample_qc.cell_filters.filter_multiplets.multiplet_threshold_method,
            params.sample_qc.cell_filters.filter_multiplets.scale_log10
        )

        // cellsnp() outputs -> vireo():
        if (params.vireo.run){

            // cellsnp() from pipeline provided inputs:


            // Vireo:
            if (params.run_with_genotype_input) {
                log.info "---running Vireo with genotype input----"
                // for each experiment_id to deconvolute, subset donors vcf to its donors and subset genomic regions.
                if (params.genotype_input.subset_genotypes){
                    log.info "---We are using subset genotypes running Vireo----"
                    CELLSNP.out.cellsnp_output_dir.combine(ch_experiment_npooled, by: 0)
                        .combine(SUBSET_GENOTYPE.out.samplename_subsetvcf, by: 0).set{full_vcf}

                }else{
                    log.info "---We are using a full genotype input for Vireo----"
                    CELLSNP.out.cellsnp_output_dir.combine(ch_experiment_npooled, by: 0).set{full_vcf}
                    full_vcf.map { experiment, cellsnpvcf, npooled -> tuple(experiment,cellsnpvcf,npooled,file(params.genotype_input.full_vcf_file))}.set {full_vcf}
                }
            }
            // Vireo without genotype input:
            else {
                log.info "-----running Vireo without genotype input----"
                CELLSNP.out.cellsnp_output_dir.combine(ch_experiment_npooled, by: 0).set{full_vcf}
                full_vcf.map {experiment, cellsnp, npooled -> tuple(experiment, cellsnp, npooled,[])}.set{full_vcf}

            }
            // TODO - filter out the full_vcf to only samples with more than 1 pooled id.
            // full_vcf.view()
            full_vcf.filter { experiment, cellsnp, npooled, t -> npooled != '1' }.set{full_vcf2}
            full_vcf.filter { experiment, cellsnp, npooled, t -> npooled == '1' }.set{not_deconvoluted}

            // full_vcf2.view()
            VIREO(full_vcf2)
            vireo_out_sample_summary_tsv = VIREO.out.sample_summary_tsv
            vireo_out_sample__exp_summary_tsv = VIREO.out.sample__exp_summary_tsv
            vireo_out_sample_donor_ids = VIREO.out.sample_donor_ids
        }


        if (params.souporcell.run){
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
            if (params.run_with_genotype_input) {

                if (params.genotype_input.subset_genotypes){
                    log.info "---We are using subset genotypes running Suporcell----"
                    // here combine the ch_experiment_bam_bai_barcodes_npooled with the output of subset
                    GUZIP_VCF(SUBSET_GENOTYPE.out.samplename_subsetvcf)
                    // ch_experiment_bam_bai_barcodes_npooled.combine(SUBSET_GENOTYPE.out.samplename_subsetvcf, by: 0).set{full_vcf}
                    ch_experiment_bam_bai_barcodes_npooled.combine(GUZIP_VCF.out.souporcell_vcf, by: 0).set{full_vcf}


                }else{
                    log.info "---We are using a full genotype input for Suporcell----"
                    // this however currently doesnt work and the individuals have to be provided.

                    // here just add the full vcf path to each of the ch_experiment_bam_bai_barcodes_npooled
                    GUZIP_VCF(tuple('full_vcf', file(params.genotype_input.full_vcf_file)))
                    GUZIP_VCF.out.souporcell_vcf.map { sample, vcf -> vcf }.set{vcf_file}
                    // we would now flatten this and take the the file of the tuple to be used next.
                    ch_experiment_bam_bai_barcodes_npooled.map {
                        samplename, bam_file, bai_file, barcodes_tsv_gz, souporcell_n_clusters ->
                        tuple(samplename, bam_file, bai_file, barcodes_tsv_gz, souporcell_n_clusters)
                        }.set{full_vcf}
                    full_vcf.combine(vcf_file).set{full_vcf}

                }
                full_vcf.combine(donors_in_lane, by: 0).set{full_vcf}
            }
            else{
                log.info "-----running Suporcell without genotype input----"
                // here make add an empty entry [] to the ch_experiment_bam_bai_barcodes_npooled
                ch_experiment_bam_bai_barcodes_npooled.map {
                    samplename, bam_file, bai_file, barcodes_tsv_gz, souporcell_n_clusters ->
                    tuple(samplename, bam_file, bai_file, barcodes_tsv_gz, souporcell_n_clusters,[],[])
                    }.set{full_vcf}
            }


            // Now that channel is created run suporcell
            full_vcf.filter { samplename, bam_file, bai_file, barcodes_tsv_gz, souporcell_n_clusters, t1, t2 -> souporcell_n_clusters != '1' }.set{full_vcf2}
            full_vcf.filter { samplename, bam_file, bai_file, barcodes_tsv_gz, souporcell_n_clusters, t1, t2 -> souporcell_n_clusters == '1' }.set{not_deconvoluted}

            SOUPORCELL(full_vcf,
                Channel.fromPath(params.souporcell.reference_fasta).collect())

        }

        if (params.vireo.run){

                // If sample is not deconvoluted we will use scrublet to detect the doublets and remove them.
                not_deconvoluted.map{ experiment, donorsvcf, npooled,t -> tuple(experiment, 'None')}.set{not_deconvoluted}
                file_cellmetadata = MULTIPLET.out.file__cellmetadata
                scrublet_paths = MULTIPLET.out.scrublet_paths
                split_channel = vireo_out_sample_donor_ids.combine(ch_experiment_filth5, by: 0)
                split_channel2 = not_deconvoluted.combine(ch_experiment_filth5, by: 0)
                split_channel = split_channel.mix(split_channel2)
                split_channel = split_channel.combine(scrublet_paths, by: 0)

                // run a multiplet detection and when splitting the donor specific h5ad remove these from non deconvoluted samples

                SPLIT_DONOR_H5AD(split_channel)

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

                // paste experiment_id and donor ID columns with __ separator
                vireo_out_sample__exp_summary_tsv = SPLIT_DONOR_H5AD.out.donor_n_cells
                .collectFile(name: "vireo_exp__donor_n_cells.tsv",
                        newLine: false, sort: true,
                        seed: "experiment_id\tn_cells\n",
                        storeDir:params.outdir+'/deconvolution/filepaths')

                vireo_out_sample_summary_tsv.view()
                vireo_out_sample__exp_summary_tsv.view()

                PLOT_DONOR_CELLS(ch_vireo_donor_n_cells_tsv)


        }else{
            out_h5ad =Channel.fromPath(params.cellsnp.vcf_candidate_snps).collect()
            vireo_out_sample__exp_summary_tsv=Channel.fromPath(params.cellsnp.vcf_candidate_snps).collect()
        }

    emit:
        out_h5ad
        vireo_out_sample__exp_summary_tsv

}
