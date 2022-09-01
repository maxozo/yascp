include { GATHER_DATA;  SPLIT_DATA_BY_STUDY} from '../modules/nf-core/modules/gather_data/main'
include { ENCRYPT_DIR; ENCRYPT_TARGET } from '../modules/local/encrypt'
include { TRANSFER;SUMMARY_STATISTICS_PLOTS } from '../modules/nf-core/modules/summary_statistics_plots/main'
include { split_bam_by_donor } from "../modules/local/cellranger_bam_per_donor"

workflow data_handover{
    take:
        outdir
        qc_input
        ch_poolid_csv_donor_assignments
        sample_possorted_bam_vireo_donor_ids

    main:
        log.info 'running data handover'
        cram_encrypted_dirs = Channel.empty()
        GATHER_DATA(outdir,qc_input.collect(),params.input_data_table)
        if (params.encrypt){
            ENCRYPT_DIR(GATHER_DATA.out.outfiles_dataset)
        }

        if (params.split_bam){
            split_bam_by_donor(sample_possorted_bam_vireo_donor_ids, params.reference_assembly_fasta_dir)
            cram_encrypted_dirs = ENCRYPT_TARGET(split_bam_by_donor.out.possorted_cram_files)
        }

        SPLIT_DATA_BY_STUDY(
          outdir,
          ENCRYPT_DIR.out.encrypted_dir,
          cram_encrypted_dirs.collect(),
          ch_poolid_csv_donor_assignments.collect()
          )

        SUMMARY_STATISTICS_PLOTS(outdir,GATHER_DATA.out.outfiles_dataset,params.input_data_table)

        // We also generate a report.
        // If we run it in sanger we transfer the data to the local website.
        TRANSFER(SUMMARY_STATISTICS_PLOTS.out.summary_plots)

}
