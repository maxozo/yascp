
// Load base.config by default for all pipelines - typically included in the nextflow config.
include { CELLBENDER } from '../modules/local/cellbender/main'
include { SPLIT_CITESEQ_GEX } from '../modules/local/citeseq/main'

include {CAPTURE_CELLBENDER_FILES} from "$projectDir/modules/local/cellbender/functions"

workflow AMBIENT_RNA {
    take:
        ch_experimentid_paths10x_raw
		ch_experimentid_paths10x_filtered
        channel__metadata
        
    main:
        Channel.empty().set { ch_versions }
        log.info params.input_data_table
        log.info """---Running Cellbender pipeline ---"""
        
        // LOGIC TO capture CELLBENDER THAT HAS ALREADY BEEN PERFORMED AND ONLY PROCESS NEW SAMPLES
        // NEEDED SINCE CELLBENDER IS A LONG PROCESS AND QUITE OFTEN THIS IS REDUNDANT IF 50/100 SAMPLES HAVE ALREADY PERFORMED THIS PROCESS.
        // Capture previously processed CellBender outputs and raw input paths that were already used

        CAPTURE_CELLBENDER_FILES(params.cellbender_location,"${params.outdir}/preprocessing",params.input_data_table)
        
        CAPTURE_CELLBENDER_FILES.out.alt_input_unfiltered.flatten()
            .map{ sample -> tuple("${sample}".replaceFirst(/.*\/captured\/unfiltered\//,"").replaceFirst(/\/.*/,""), sample) }
            .set{ ch_previously_processed_raw_inputs }
        
        CAPTURE_CELLBENDER_FILES.out.alt_input_filtered.flatten()
            .map{ sample -> tuple("${sample}".replaceFirst(/.*\/captured\/filtered\//,"").replaceFirst(/\/.*/,""), sample) }
            .set{ ch_previously_filtered_outputs }

        // CAPTURE_CELLBENDER_FILES.out.cb_to_use_downstream.flatten()
        //     .map{ sample -> tuple("${sample}".replaceFirst(/.*\/cellbender\//,"").replaceFirst(/\/.*/,""), sample) }
        //     .set{ ch_previously_filtered_outputs }

        ch_experimentid_paths10x_raw
            .join(ch_previously_processed_raw_inputs, by: [0], remainder: true)
            .set{ ch_joined_raw_with_processed }

        ch_experimentid_paths10x_filtered
            .join(ch_previously_filtered_outputs, by: [0], remainder: true)
            .set{ ch_joined_filtered_with_processed }

        ch_joined_raw_with_processed
            .filter{ it[1] != null }
            .filter{ it[2] == null }
            .map{ row -> tuple(row[0], row[1]) }
            .set{ ch_new_raw_inputs }

        ch_joined_filtered_with_processed
            .filter{ it[1] != null }
            .filter{ it[2] == null }
            .map{ row -> tuple(row[0], row[1]) }
            .set{ ch_new_filtered_inputs }


        ch_joined_raw_with_processed
            .filter{ it[1] != null }
            .filter{ it[2] != null }
            .map{ row -> tuple(row[0], row[2]) }
            .set{ ch_previously_processed_raw_inputs_matched }


        ch_joined_filtered_with_processed
            .filter{ it[1] != null }
            .filter{ it[2] != null }
            .map{ row -> tuple(row[0], row[2]) }
            .set{ ch_previously_processed_filtered_outputs }

        CELLBENDER(ch_new_raw_inputs,ch_new_filtered_inputs,channel__metadata)

        ch_versions = ch_versions.mix(CELLBENDER.out.cellbender_versions)

        // AMBIENTNESS_QUANTIFICATION()

        cellbender_path_processed = CELLBENDER.out.cellbender_path
        cellbender_path_raw_processed = CELLBENDER.out.cellbender_path_raw
        
        cellbender_path = cellbender_path_processed.concat(ch_previously_processed_filtered_outputs)
        cellbender_path_raw = cellbender_path_raw_processed.concat(ch_previously_processed_raw_inputs_matched)

        ch_previously_processed_raw_inputs_matched.subscribe { println "✅ SKIP raw: $it" }
        ch_previously_processed_filtered_outputs.subscribe { println "✅ SKIP filtered: $it" }

    emit:
        cellbender_path
        cellbender_path_raw
        cellbender_versions = ch_versions
}
