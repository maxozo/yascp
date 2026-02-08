nextflow.enable.dsl=2

// main deconvolution modules, common to all input modes:

include { CELLSNP;
          CAPTURE_CELLSNP_FILES;
          DYNAMIC_DONOR_EXCLUSIVE_SNP_SELECTION; MPILEUP;
          ASSESS_CALL_RATE } from "$projectDir/modules/local/cellsnp/main"
include { VIREO;
          POSTPROCESS_SUMMARY;
          REMOVE_DUPLICATED_DONORS_FROM_GT;
          VIREO_SUBSAMPLING;
          VIREO_SUBSAMPLING_PROCESSING;
          GENOTYPE_MATCHER;
          CAPTURE_VIREO } from "$projectDir/modules/local/vireo/main"
include { SUBSET_BAM_PER_BARCODES_AND_VARIANTS } from "$projectDir/modules/local/subset_bam_per_barcodes_and_variants/main"
include { FREEBAYES } from "$projectDir/modules/local/freebayes/main"
include { PREPROCESS_GENOTYPES } from "$projectDir/modules/local/genotypes/main"
include { SPLIT_DONOR_H5AD; PREP_ASSIGNMENTS_FILE } from "$projectDir/modules/local/split_donor_h5ad/main"
include { PLOT_DONOR_CELLS } from "$projectDir/modules/local/plot_donor_cells/main"
include { ENHANCE_VIREO_METADATA_WITH_DONOR } from "$projectDir/modules/local/genotypes/main"
include { MATCH_GENOTYPES } from './match_genotypes'
include { ENHANCE_STATS_GT_MATCH } from "$projectDir/modules/local/genotypes/main"
include { SUBSET_WORKF } from "$projectDir/modules/local/subset_genotype/main"
include { REPLACE_GT_DONOR_ID2 } from "$projectDir/modules/local/genotypes/main"
include { STAGE_FILE } from "$projectDir/modules/local/retrieve_resources/retrieve_resources"
include { GT_MATCH_POOL_IBD } from "$projectDir/modules/local/genotypes/main"
include {VIREO_GT_FIX_HEADER;
         VIREO_ADD_SAMPLE_PREFIX; 
         MERGE_GENOTYPES_IN_ONE_VCF_IDX_PAN; 
         MERGE_GENOTYPES_IN_ONE_VCF_FREEBAYES} from "$projectDir/modules/local/genotypes/main"

include {COLLECT_FILE as COLLECT_FILE1;
        COLLECT_FILE as COLLECT_FILE2;
        COLLECT_FILE as COLLECT_FILE3;
        COLLECT_FILE as COLLECT_FILE4;
        COLLECT_FILE as COLLECT_FILE5;
        COLLECT_FILE as COLLECT_FILE6} from "$projectDir/modules/local/collect_file/main"

workflow  MAIN_DECONVOLUTION {


    take:
		ch_experiment_bam_bai_barcodes
		ch_experiment_npooled
		ch_experiment_filth5
		ch_experiment_donorsvcf_donorslist
        scrublet_paths
        ch_ref_vcf2
        genome

    main:
		log.info "#### running DECONVOLUTION workflow #####"
        Channel.empty().set { ch_versions }


        vcf_candidate_snps = STAGE_FILE(params.cellsnp.vcf_candidate_snps)

        if (params.genotype_input.run_with_genotype_input) {
            log.info """---Running with genotype input ---"""
            // We have to produce a single vcf file for each individual pool.
            // Therefore we create 2 channels:
            // 1) All the expected vcf ids listed in the donor table
            log.info "# selecting SNPs and subsetting GT #"
            Channel.fromPath(params.input_data_table,      
                followLinks: true,
                checkIfExists: true
                ).splitCsv(header: true, sep: '\t').map { row -> tuple(row.experiment_id, row.donor_vcf_ids) }
                .set { donors_in_pools }
            
            // This will subsequently result in a joint vcf file per cohort per donors listed for each of the pools that can be used in VIREO and/or GT matching algorythm.
            PREPROCESS_GENOTYPES(ch_ref_vcf2)
            ch_ref_vcf = PREPROCESS_GENOTYPES.out.vcf_tuple
            ch_versions = ch_versions.mix(PREPROCESS_GENOTYPES.out.versions)
            SUBSET_WORKF(ch_ref_vcf,donors_in_pools,'AllExpectedGT',genome)
            ch_versions = ch_versions.mix(SUBSET_WORKF.out.subset_workf_versions)
            merged_expected_genotypes = SUBSET_WORKF.out.merged_expected_genotypes
            merged_expected_genotypes2 = merged_expected_genotypes.combine(vcf_candidate_snps)
            GT_MATCH_POOL_IBD(SUBSET_WORKF.out.samplename_subsetvcf_ibd,'Withing_expected','Expected')
            ch_versions = ch_versions.mix(GT_MATCH_POOL_IBD.out.versions)

            if (params.use_bam_derived_cellsnp_panel){
                // Here should add an option to derive panel from actual bam files as different technologies has different coverages and ceirtain panels may not work. 
                MPILEUP(ch_experiment_bam_bai_barcodes,params.reference_assembly_fasta_dir)
                mpileup_out_chanel = MPILEUP.out.pileup

            }


            // This takes the subset genotypes expected in the pool and select informative SNPs between them and adds them to the panel. 
            DYNAMIC_DONOR_EXCLUSIVE_SNP_SELECTION(params.add_snps_to_pile_up_based_on_genotypes_provided,merged_expected_genotypes2)
            ch_versions = ch_versions.mix(DYNAMIC_DONOR_EXCLUSIVE_SNP_SELECTION.out.versions)
            cellsnp_panels = DYNAMIC_DONOR_EXCLUSIVE_SNP_SELECTION.out.cellsnp_pool_panel

            // merged_expected_genotypes.subscribe { println "merged_expected_genotypes: $it" }
            informative_uninformative_sites = DYNAMIC_DONOR_EXCLUSIVE_SNP_SELECTION.out.informative_uninformative_sites

        }else{
            ch_ref_vcf = Channel.of()
        }

        ch_poolid_donor_assignment = Channel.empty()
        ch_experiment_donorsvcf_donorslist.map { experiment, donorsvcf, donorstbi,donorslist -> tuple(experiment, donorslist.replaceAll(/,/, " ").replaceAll(/"/, ""))}.set{donors_in_lane}
        ch_experiment_bam_bai_barcodes.combine(ch_experiment_npooled, by: 0).set{cellsnp_with_npooled}

        if (params.genotype_input.run_with_genotype_input) {
            // Here we provide sites to pile up the snps within the pool. 
            // As a starting point instead of the default cellsnp panel users are encouraged to:
            if (params.provide_within_pool_donor_specific_sites_for_pilup){
                log.info """---Running with provide_within_pool_donor_specific_sites_for_pilup ---"""
                cellsnp_with_npooled_pre = cellsnp_with_npooled.join(cellsnp_panels, remainder: true)
                // 
            }else{
                log.info """---Running WITHOUT provide_within_pool_donor_specific_sites_for_pilup ---"""
                cellsnp_with_npooled_pre = cellsnp_with_npooled.combine(vcf_candidate_snps)
            }
            
            cellsnp_with_npooled_pre2 = cellsnp_with_npooled_pre.combine(vcf_candidate_snps)
            if (params.use_bam_derived_cellsnp_panel){
                cellsnp_with_npooled_pre = cellsnp_with_npooled.combine(mpileup_out_chanel, by: 0)
                cellsnp_with_npooled_pre2 = cellsnp_with_npooled_pre.combine(vcf_candidate_snps)
                cellsnp_with_npooled_pre2.subscribe { println "cellsnp_with_npooled_pre2: $it" }
            }
            
            cellsnp_with_npooled2 = cellsnp_with_npooled_pre2.map { experiment, bam, bai, barcodes, nrPool, panel, default_cellsnp ->
                                                                    [
                                                                        experiment,
                                                                        bam,
                                                                        bai,
                                                                        barcodes,
                                                                        nrPool,
                                                                        panel == null ? default_cellsnp : panel
                                                                    ]
                                                                }
            cellsnp_with_npooled2.subscribe { println "cellsnp_with_npooled2: $it" }
        }else{
            log.info """---Running without genotype input ---"""
            cellsnp_with_npooled2 = cellsnp_with_npooled.combine(vcf_candidate_snps)
        }

        log.info('Capturing some of the existing CELLSNP files')
        CAPTURE_CELLSNP_FILES(params.existing_cellsnp)
        CAPTURE_CELLSNP_FILES.out.cellsnp_loc.splitCsv(header: false, sep: ' ')
            .map{row->tuple(row[0], "${row[1]}")}
            .set{cellsnp_output_dir1}
        cellsnp_output_dir1.join(cellsnp_with_npooled2, remainder: true).set{filter_channel}
        // filter_channel.subscribe { println "cellsnp_output_dir1 filter_channel: $it" }
        filter_channel.filter{ it[1] == null }.map{row -> tuple(row[0], row[2],row[3],row[4],row[5],row[6])}.set{cellsnp_with_npooled}
        
        CAPTURE_CELLSNP_FILES.out.cellsnp_loc.splitCsv(header: false, sep: ' ')
            .map{row->tuple(row[0], file("${row[1]}/cellSNP.cells.vcf.gz"))}
            .set{cellsnp_cell_vcfs}           
        CAPTURE_CELLSNP_FILES.out.cellsnp_loc.splitCsv(header: false, sep: ' ')
            .map{row->tuple(row[0])}
            .set{cellsnp_cap_ids} 
        // cellsnp_with_npooled.subscribe { println "cellsnp_with_npooled: $it" }
        CELLSNP(cellsnp_with_npooled)
        ch_versions = ch_versions.mix(CELLSNP.out.versions)
        if (params.genotype_input.run_with_genotype_input) {
            // Here we assess how many informative sites has been called on. 
            CELLSNP.out.cell_vcfs.combine(DYNAMIC_DONOR_EXCLUSIVE_SNP_SELECTION.out.informative_uninformative_sites, by: 0).set{assess_call_rate_input}
            ASSESS_CALL_RATE(assess_call_rate_input)
            COLLECT_FILE4(ASSESS_CALL_RATE.out.variants_description.collect(),"all_variants_description.tsv",params.outdir+'/deconvolution/concordances',1,'')
        }
        for_bam_pileups = CELLSNP.out.for_bam_pileups
        
        cellsnp_cap_ids.combine(ch_experiment_bam_bai_barcodes, by: 0).combine(vcf_candidate_snps).set{cellsnp_cap_ids2}
        cellsnp_cap_ids2.map{row -> tuple(row[0], row[4],row[1])}.set{cellsnp_cap_ids2}
        for_bam_pileups2 = cellsnp_cap_ids2.mix(for_bam_pileups)

        cellsnp_output_dir2 = CELLSNP.out.cellsnp_output_dir
        cellsnp_output_dir=cellsnp_output_dir1.concat(cellsnp_output_dir2)
        cellsnp_cell_vcfs2=cellsnp_cell_vcfs.concat(CELLSNP.out.cell_vcfs)

        cellsnp_output_dir.combine(ch_experiment_npooled, by: 0).set{full_vcf2}
        full_vcf2.map {experiment, cellsnp, npooled -> tuple(experiment, cellsnp, npooled,[],[])}.set{full_vcf2}
         // Here we run Vireo software to perform the donor deconvolution. Note that we have coded the pipeline to be capable in using
         // the full genotypes as an input and also subset to the individuals provided as an input in the donor_vcf_ids column.
        if (params.genotype_input.vireo_with_gt && params.genotype_input.run_with_genotype_input) {
            log.info "---running Vireo with genotype input----"
            // for each experiment_id to deconvolute, subset donors vcf to its donors and subset genomic regions.
            // Here we subset the genotypes. This happens if the input.nf contains subset_genotypes = true
            // log.info "---We are using subset genotypes running Vireo----"
            // We need to make sure that the expected genotypes dont contain repeated genotypes - donors sequenced twice.

            if (params.genotype_phenotype_mapping_file!=''){
                REMOVE_DUPLICATED_DONORS_FROM_GT(merged_expected_genotypes,params.genotype_phenotype_mapping_file,params.input_data_table)
                ch_versions = ch_versions.mix(REMOVE_DUPLICATED_DONORS_FROM_GT.out.versions)
                merged_expected_genotypes2= REMOVE_DUPLICATED_DONORS_FROM_GT.out.merged_expected_genotypes
            }else{
                merged_expected_genotypes2 = merged_expected_genotypes
            }

            cellsnp_output_dir.combine(ch_experiment_npooled, by: 0)
                .combine(merged_expected_genotypes2, by: 0).set{full_vcf_pre}

            full_vcf_pre.join(full_vcf2, remainder: true).set{filter_channel_vcf}
            // just take the null chanel, and add that to full vcf 
            filter_channel_vcf.filter{ it[1] == null }.map{row -> tuple(row[0], row[2],row[3],row[4],row[5])}.set{full_vcf4}
            full_vcf_pre.mix(full_vcf4).set{full_vcf}
        } else {
            // Vireo can also be run without the genotypes, and the performance is equally good then running with, however have to be careful in what sites file is provided.
            // Here we run it without the genotypes and the donors are labeled as donor0, donor1 etc, dependant on the number of donors set in the input file n_pooled column
            log.info "-----running Vireo without genotype input----"
            full_vcf = full_vcf2
        }


        CAPTURE_VIREO(params.existing_vireo)
        CAPTURE_VIREO.out.output_dir.flatten().map{row->tuple("${row}".replaceFirst(/.*vireo_/,""), "${row}/donor_ids.tsv")}
            .set{vireo_out_sample_donor_ids_cap}  
        CAPTURE_VIREO.out.output_dir.flatten().map{row->tuple("${row}")}
            .set{vireo_paths_0} 
        CAPTURE_VIREO.out.output_dir.flatten().map{row->tuple("${row}".replaceFirst(/.*vireo_/,""), "${row}/GT_donors.vireo.vcf.gz")}
            .set{vireo_out_sample_donor_vcf_cap}  
        CAPTURE_VIREO.out.output_dir.flatten().map{row->tuple("${row}".replaceFirst(/.*vireo_/,""), "${row}/summary.tsv")}
            .set{vireo_out_sample_summary_tsv_cap} 
        CAPTURE_VIREO.out.output_dir.flatten().map{row->tuple("${row}".replaceFirst(/.*vireo_/,""), "${row}/GT_donors.vireo.vcf.gz", "${row}/donor_ids.tsv")}
            .set{all_required_data_cap_pre} 

        // When all the channels are prpeared accordingly we exacute the vireo with the prpeared channel.
        full_vcf.filter { experiment, cellsnp, npooled, t,ti -> npooled != '1' }.set{full_vcf2}
        full_vcf2.join(vireo_out_sample_donor_ids_cap, remainder: true).set{filter_channel_vireo}

        all_required_data_cap_pre.join(full_vcf2, remainder: false).set{all_required_data_cap_pre}
        all_required_data_cap_pre.map{row -> tuple(row[0], row[1],row[2],row[5],row[6])}.set{all_required_data_cap}

        filter_channel_vireo.filter{ it[5] == null }.map{row -> tuple(row[0], row[1],row[2],row[3],row[4])}.set{full_vcf3}
        full_vcf.filter { experiment, cellsnp, npooled, t,ti -> npooled == '1' }.set{not_deconvoluted}
        
        not_deconvoluted.join(for_bam_pileups2, remainder: false).set{for_bam_pileups_pre}
        for_bam_pileups_pre.map{row -> tuple(row[0], row[5],row[6],row[7])}.set{for_bam_pileups_3}

        Channel.of(1..params.vireo.subsample_times).set{itterations}
        full_vcf2.combine(itterations).set{vireo_extra_repeats}
        // full_vcf3.subscribe { println "1:: vireo_paths_map $it" }
        VIREO(full_vcf3)
        ch_versions = ch_versions.mix(VIREO.out.versions)
        VIREO.out.summary.mix(vireo_out_sample_summary_tsv_cap).set{summary_files}
        POSTPROCESS_SUMMARY(summary_files)
        
        //But to make it consistent we still randomly assign donor IDs per pool for each of the names.
        vireo_with_gt = Channel.of(params.genotype_input.vireo_with_gt)
        VIREO.out.all_required_data.set{replacement_input3}
        replacement_input3.mix(all_required_data_cap).set{replacement_input2}

        replacement_input2.combine(POSTPROCESS_SUMMARY.out.summary_tsvs, by: 0).set{replacement_input}

        replacement_input.combine(vireo_with_gt).set{vir_repl_input}

        REPLACE_GT_DONOR_ID2(vir_repl_input)
        // combine the outputs of the capture_vireo and Vireo
        vireo_paths_0.concat(REPLACE_GT_DONOR_ID2.out.output_dir).set{vireo_paths}

        VIREO_GT_FIX_HEADER(REPLACE_GT_DONOR_ID2.out.infered_vcf,genome)
        VIREO_ADD_SAMPLE_PREFIX(VIREO_GT_FIX_HEADER.out.infered_vcf)
        
        vireo_out_sample_donor_vcf = VIREO_GT_FIX_HEADER.out.infered_vcf
        vireo_out_sample_summary_tsv = REPLACE_GT_DONOR_ID2.out.sample_summary_tsv
        vireo_out_sample__exp_summary_tsv = REPLACE_GT_DONOR_ID2.out.sample__exp_summary_tsv
        vireo_out_sample_donor_ids = REPLACE_GT_DONOR_ID2.out.sample_donor_ids

        ch_experiment_bam_bai_barcodes
            .map { samplename, bam_file, bai_file, barcodes_tsv_gz -> tuple(samplename, file(bam_file)) }
            .combine(vireo_out_sample_donor_ids, by: 0 )
            .set { ch_experiment_bam_vireo_donor_ids }

        // If sample is not deconvoluted we will use scrublet to detect the doublets and remove them.
        not_deconvoluted.map{ experiment, donorsvcf, npooled,t,t2 -> tuple(experiment, 'None')}.set{not_deconvoluted2}
        split_channel = vireo_out_sample_donor_ids.combine(ch_experiment_filth5, by: 0)
        split_channel2 = not_deconvoluted2.combine(ch_experiment_filth5, by: 0)

        // combining these 2 channels in one
        split_channel3 = split_channel.mix(split_channel2)

        // adding the scrublet paths to the channel.
        split_channel4 = split_channel3

        split_channel5 = split_channel4.map{
            val_sample, val_donor_ids_tsv, val_filtered_matrix_h5 ->
            [  val_sample,file(val_donor_ids_tsv),file(val_filtered_matrix_h5),params.outdir]
        }

        SPLIT_DONOR_H5AD(split_channel5)
        PREP_ASSIGNMENTS_FILE(split_channel5)
        cell_assignments = PREP_ASSIGNMENTS_FILE.out.cell_assignments
        if (params.bam_pileup_per_donor){
            bam_subset_chanel = SPLIT_DONOR_H5AD.out.sample_donor_level_barcodes.transpose().combine(for_bam_pileups_3, by: 0)
            SUBSET_BAM_PER_BARCODES_AND_VARIANTS(bam_subset_chanel)  // This process subsets the ba
            ch_versions = ch_versions.mix(SUBSET_BAM_PER_BARCODES_AND_VARIANTS.out.versions)

            chromosomes =  Channel.of(1,2,3,4,5,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21)
            freebayes_pre = SUBSET_BAM_PER_BARCODES_AND_VARIANTS.out.freebayes_input

            freebayes_in = freebayes_pre.combine(chromosomes)

            FREEBAYES(freebayes_in,genome)
            ch_versions = ch_versions.mix(FREEBAYES.out.versions)
            vireo_out_sample_donor_vcf = FREEBAYES.out.freebayes_vcf

            FREEBAYES.out.gt_pool.groupTuple(by:0).set{fbb1}
            fbb1.map{row->tuple(row[0], row[1], row[2], row[3][0])}.set{fbb1}
            MERGE_GENOTYPES_IN_ONE_VCF_IDX_PAN(fbb1,'freebayes')
            ch_versions = ch_versions.mix(MERGE_GENOTYPES_IN_ONE_VCF_IDX_PAN.out.versions)
            MERGE_GENOTYPES_IN_ONE_VCF_IDX_PAN.out.gt_pool.groupTuple(by:0).set{fbb2}
            MERGE_GENOTYPES_IN_ONE_VCF_FREEBAYES(fbb2,'freebayes')
            ch_versions = ch_versions.mix(MERGE_GENOTYPES_IN_ONE_VCF_FREEBAYES.out.versions)
            gt_math_pool_against_panel_input2 = MERGE_GENOTYPES_IN_ONE_VCF_FREEBAYES.out.gt_pool.combine(ch_ref_vcf)
            vir_inp =  MERGE_GENOTYPES_IN_ONE_VCF_FREEBAYES.out.vir_input

        }else{
            gt_math_pool_against_panel_input2 = Channel.of()
            vir_inp =  Channel.of()
        }

        vir_inp2 = vir_inp.collect()
        vir_inp3 = vireo_paths.collect()
        vireo_paths2 = gt_matcher_inp = vir_inp3.mix(vir_inp2)
              
        GENOTYPE_MATCHER(gt_matcher_inp)

        matched_donors = GENOTYPE_MATCHER.out.matched_donors

        if (params.genotype_input.run_with_genotype_input) {
            if (params.vireo.do_vireo_subsampling){
                VIREO_SUBSAMPLING(vireo_extra_repeats)
                ch_versions = ch_versions.mix(VIREO_SUBSAMPLING.out.versions)
                VIREO_SUBSAMPLING.out.output_dir.concat(VIREO.out.output_dir_subsampling).set{tuple_1}
                tuple_1.groupTuple(by:0).set{vspp0}
                VIREO_SUBSAMPLING_PROCESSING(vspp0)
                subsampling_donor_swap = VIREO_SUBSAMPLING_PROCESSING.out.subsampling_donor_swap
            }else{
                subsampling_donor_swap = Channel.from("$projectDir/assets/fake_file.fq")
            }

            VIREO_GT_FIX_HEADER.out.gt_pool
                .combine(ch_ref_vcf)
                .set { gt_math_pool_against_panel_input3 }
            gt_math_pool_against_panel_input = gt_math_pool_against_panel_input3.mix(gt_math_pool_against_panel_input2)

            MATCH_GENOTYPES(vireo_out_sample_donor_vcf,merged_expected_genotypes,VIREO_GT_FIX_HEADER.out.gt_pool,gt_math_pool_against_panel_input,genome,ch_ref_vcf,cellsnp_cell_vcfs2,cell_assignments,subsampling_donor_swap,informative_uninformative_sites)
            gt_matches = MATCH_GENOTYPES.out.donor_match_table.collect()


            ENHANCE_STATS_GT_MATCH(MATCH_GENOTYPES.out.donor_match_table_enhanced,params.input_data_table)
            COLLECT_FILE5(ENHANCE_STATS_GT_MATCH.out.assignments.collect(),"assignments_all_pools.tsv",params.outdir+'/deconvolution/vireo_processed',1,'')
            COLLECT_FILE6(ENHANCE_STATS_GT_MATCH.out.assignments.collect(),"assignments_all_pools.tsv",params.outdir+'/deconvolution/gtmatch',1,'')
            assignments_all_pools = COLLECT_FILE5.out.output_collection
            gt_matches = Channel.from("$projectDir/assets/fake_file.fq")
            assignments_all_pools = Channel.from("$projectDir/assets/fake_file.fq")
        }else{
            gt_matches = Channel.from("$projectDir/assets/fake_file.fq")
            assignments_all_pools = Channel.from("$projectDir/assets/fake_file.fq")
        }

        header_seed="experiment_id\th5ad_filepath"
        out_h5ad = COLLECT_FILE1(SPLIT_DONOR_H5AD.out.exp__donors_h5ad_assigned_tsv.collect(),"exp__donors_h5ad_assigned.tsv",0,0,header_seed)

        header_seed="experiment_id\tdonor\tn_cells"
        ch_vireo_donor_n_cells_tsv = COLLECT_FILE2(REPLACE_GT_DONOR_ID2.out.sample_summary_tsv.collect(),"vireo_donor_n_cells.tsv",0,0,header_seed)

        // paste experiment_id and donor ID columns with __ separator
        header_seed="experiment_id\tn_cells"
        vireo_out_sample__exp_summary_tsv = COLLECT_FILE3(SPLIT_DONOR_H5AD.out.donor_n_cells.collect(),"vireo_exp__donor_n_cells.tsv",0,0,header_seed)

        if (params.genotype_input.run_with_genotype_input & params.genotype_input.posterior_assignment) {
            if (params.extra_sample_metadata!=''){
                // Here we have sample level metadata but we have chosen to keep the donor ids.
                // in this scenario we enhance the donor level vireo metadata file to add the donor metadata to the h5ads and eventually to the donor and teanche report
                ENHANCE_VIREO_METADATA_WITH_DONOR(params.extra_sample_metadata,vireo_out_sample__exp_summary_tsv,REPLACE_GT_DONOR_ID.out.assignments.collect())
                vireo_out_sample__exp_summary_tsv = ENHANCE_VIREO_METADATA_WITH_DONOR.out.replaced_vireo_exp__donor_n_cells_out
            }
        }
        PLOT_DONOR_CELLS(ch_vireo_donor_n_cells_tsv)

    emit:
        out_h5ad
        vireo_out_sample__exp_summary_tsv
        vireo_out_sample_donor_vcf
        ch_poolid_csv_donor_assignments = ch_poolid_donor_assignment
        sample_possorted_bam_vireo_donor_ids = ch_experiment_bam_vireo_donor_ids
        assignments_all_pools
        vireo_paths2
        matched_donors
        deconvolution_versions = ch_versions
}
