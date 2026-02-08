include {
    UMAP_CALCULATE;
    UMAP_GATHER;
    UMAP_PLOT_SWARM;
    UMAP_CALCULATE_AND_PLOT;GENERATE_FINAL_UMAPS;
} from "./functions.nf"

workflow UMAP {
    take:
        outdir
        anndata
        metadata
        pcs
        reduced_dims
        use_pcs_as_reduced_dims
        n_neighbors
        umap_init
        umap_min_dist
        umap_spread
        colors_quantitative
        colors_categorical
        method
    main:
        UMAP_CALCULATE(
            outdir,
            anndata,
            metadata,
            pcs,
            reduced_dims,
            use_pcs_as_reduced_dims,
            n_neighbors,
            umap_init,
            umap_min_dist,
            umap_spread,
            method
        )
        UMAP_CALCULATE.out.adata_out.collect().set{all_umaps}
        UMAP_CALCULATE.out.outdir_anndata.groupTuple()
            .reduce([:]) { map, tuple ->  // 'map' is used to collect values;
                                          // 'tuple' is the record
                def file_id = tuple[0]    // the first item is the 'iter_id'
                def group = map[file_id]  // the aggregation for cur 'iter_id'
                if( !group ) group = [ file_id ]  // if new, create a new entry
                group[1] = tuple[1][0]  // cannot use uniq here because nf
                group[2] = tuple[2][0]  // stages the same file names in
                group[3] = tuple[3][0]  // different working dirs...
                group[4] = tuple[4][0]  // the solution is to use first item
                group[5] = tuple[5][0]
                group[6] = tuple[6]     // list of umaps to merge
                map[file_id] = group    // set back into the map
                return map // return it so it will be used in the next iteration
            }.flatMap { it.values() }.set{umap_gather_input} // tricky part: get the list of values of
                                     // in the map, each value is the
                                     // aggregation build above
                                     // the 'flatMap' emits each of these
                                     // aggregation list as a single item
   
        // Gather step.
        // Gather by tuple ... if we just to a collect, then will get all
        // umap_calculate calls, not split by reduced_dims. See link below:
        // http://nextflow-io.github.io/patterns/index.html#_process_outputs_into_groups
        UMAP_GATHER(
            umap_gather_input
        )

        if (params.run_celltype_assignment){
            GENERATE_FINAL_UMAPS(UMAP_GATHER.out.anndata,params.outdir)
        }

        // Make plots
        UMAP_PLOT_SWARM(
            UMAP_GATHER.out.outdir,
            UMAP_GATHER.out.anndata,
            // cluster.out.reduced_dims,
            colors_quantitative,
            colors_categorical,
            '20'
        )
    emit:
        // Return merged input data file.
        outdir = UMAP_GATHER.out.outdir
        anndata = UMAP_GATHER.out.anndata
        metadata = UMAP_GATHER.out.metadata
        pcs = UMAP_GATHER.out.pcs
        reduced_dims = UMAP_GATHER.out.reduced_dims

}