
include {
    SCCAF_ASSESS_CLUSTERING;
    SCCAF_OPTIMIZE_CLUSTERING;
} from "./functions.nf"

workflow SCCAF {
  take:
    outdir
    anndata
    external_clustering
    min_accuracy

  main:
   
    SCCAF_ASSESS_CLUSTERING(
        outdir,
        anndata,
        external_clustering
    )

    SCCAF_OPTIMIZE_CLUSTERING(
        outdir,
        anndata,
        external_clustering,
        min_accuracy
    )
}
