
// Load base.config by default for all pipelines - typically included in the nextflow config.
include { main_deconvolution } from './main_deconvolution.nf'
include { CELLTYPIST } from '../modules/nf-core/modules/celltypist/main'

workflow deconvolution {
    take:
        ch_experiment_bam_bai_barcodes
        ch_experiment_npooled
        ch_experiment_filth5
        ch_experiment_donorsvcf_donorslist
        ch_donorid_swap_table
        ch_vcf_dir
        channel__file_paths_10x
    main:
        // // run main deconvolution pipeline on prepared input channels:
        main_deconvolution(ch_experiment_bam_bai_barcodes,
                ch_experiment_npooled,
                ch_experiment_filth5,
                ch_experiment_donorsvcf_donorslist,
                ch_donorid_swap_table,
                ch_vcf_dir,
                channel__file_paths_10x)
        out_h5ad = main_deconvolution.out.out_h5ad
        vireo_out_sample__exp_summary_tsv=main_deconvolution.out.vireo_out_sample__exp_summary_tsv
    emit:
        out_h5ad
        vireo_out_sample__exp_summary_tsv

}
