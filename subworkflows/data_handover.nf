include { GATHER_DATA} from "$projectDir/modules/local/gather_data/main"
include { TRANSFER;SUMMARY_STATISTICS_PLOTS } from "$projectDir/modules/local/summary_statistics_plots/main"
include { SUBSET_BAM_PER_BARCODES } from "$projectDir/modules/local/subset_bam_per_barcodes_and_variants/main"

workflow DATA_HANDOVER{
    take:
        outdir
        input_channel
        qc_input
        ch_poolid_csv_donor_assignments
        sample_possorted_bam_vireo_donor_ids
        genome
    main:
        log.info 'running data handover'

        GATHER_DATA(outdir,qc_input.collect(),input_channel)
        gh_out  = GATHER_DATA.out.outfiles_dataset

        if (params.split_bam){
            GATHER_DATA.out.barcodes_files.flatten().map{sample -> tuple("${sample}".replaceFirst(/.*\//,"").replaceFirst(/\..*/,""),"${sample}".replaceFirst(/.*\//,"").replaceFirst(/\.tsv.*/,""),sample)}.set{barcodes}
            barcodes.combine(sample_possorted_bam_vireo_donor_ids, by: 0).set{full_split_chanel_input}

            GATHER_DATA.out.barcodes_files.flatten().map{sample -> tuple("${sample}".replaceFirst(/.*\//,"").replaceFirst(/\..*/,""),"${sample}".replaceFirst(/.*\//,"").replaceFirst(/\.tsv.*/,""),sample)}.set{barcodes}
            barcodes.combine(sample_possorted_bam_vireo_donor_ids, by: 0).set{full_split_chanel_input}
            SUBSET_BAM_PER_BARCODES(full_split_chanel_input,genome)
        }

        SUMMARY_STATISTICS_PLOTS(outdir,gh_out,params.input_data_table)
        TRANSFER(SUMMARY_STATISTICS_PLOTS.out.summary_plots,params.rsync_to_web_file,outdir)

}
