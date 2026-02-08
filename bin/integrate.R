#!/usr/bin/env Rscript
# 1. integrate all donors by SCT to get one integrated SCT object
# 2a. integrate all donors by CITE to get one integrated CITE object
# 2b. integrate all donors by CITE_bgRemoved to get one integrated CITE_bgRemoved object
# 3a. integrate SCT+CITE
# 3b. integrate SCT+CITE_bgRemoved

#### libraries ####
library("Matrix")
library(Seurat)
library(ggplot2)
library("RColorBrewer")
library(future)
# options(future.globals.maxSize= 300 * 1024^3)
# if (future::supportsMulticore()) {
#   future::plan(future::multicore)
# } else {
#   future::plan(future::multisession)
# }
future.seed=TRUE
args = commandArgs(trailingOnly=TRUE)

# args=vector(mode='list', length=6); args[[1]]='pct_counts_gene_group__mito_transcript'; args[[2]]=5; args[[3]]=30; args[[4]]=12;args[[5]]= 10;args[[6]]= 9
if (args[1]=='NONE') {
  vars_to_regress = c()
  reg_name = 'regress__NONE'
}else{
  vars_to_regress = strsplit(as.character(args[1]), split=";")
  vars_to_regress = unlist(vars_to_regress)
  reg_name = paste0('regress__',args[1])
}
print(vars_to_regress)
# integrate sct
#  NONE 5 30 12 10 9
#  k.anchor= strtoi(5)
#  dims= strtoi(30)
#  ndim_sct= strtoi(12)
#  ndim_citeBgRemoved= strtoi(10)
#  ndim_cite_integrated = strtoi(9)
k.anchor=strtoi(args[2])
print(k.anchor)
dims=strtoi(args[3])

# Number of dimensions to use for SCT integrated data (check elbow plot to see if it's correct)
ndim_sct = strtoi(args[4])
# Number of dimensions to use for cite_bgRemoved integrated data (check elbow plot to see if it's correct)
ndim_citeBgRemoved = strtoi(args[5])
# Number of dimensions to use for cite integrated data (check elbow plot to see if it's correct)
ndim_cite_integrated = strtoi(args[6])
data_dir <- getwd()
outdir <- paste0('.')
figdir <- paste0(outdir,'/figures','__',reg_name,'/')
dir.create(figdir,showWarnings = F)
tmp_rds_file <- paste0(reg_name,'__','all_samples_integrated.RDS')
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
# Select the sc eqtlgen rds files
wg2_dir <- getwd()
cite_files <- list.files(pattern='.RDS',
                         path=Sys.glob('.'),recursive=T, full.names=T)
#####

sobj_list <- list()
sample_names <- c()
i <- 1
i2 <- 1
for(f in cite_files){
  i2 <- i2+1

  sample_id <- gsub('.withADT|.RDS','',basename(f))
  sample_name <- sample_id
  print(paste0(sample_id,': sample'))
  sobj_per_donor <- readRDS(f)  #// Here we actually dont have the doublets anymore since they were not merged back in the files.
  # print(sobj_per_donor)
  for(donor in names(sobj_per_donor)){
    print(donor)
    don = sobj_per_donor[[donor]]
    # Then add this pool-donor seurat object to list so that sobj_list
    # will contain seurat objects for all donors of all pools
    # if ("CITE" %in% Assays(don)){
      sobj_list[[i]] <- don
      sample_names <- c(sample_names,  paste0(sample_name,'-',donor))     
      i <- i+1
    # }else{
    #   print("CITE assay is not present")
    # }
  }
  # if (i2>5){
  #     break
  # }
}


print("---- Integrate ------")
integrate_sct <- function(slist, reference_samples, k.anchor=5, dims=30){
  # Integrate a list of seurat objects together using the SCT assay using rpca
  # slist: list of seurat objects
  # reference_samples: list of names to use as reference (reference parameter in FindIntegrationAnchors)
  # k.anchor: Number of anchors to use in FindIntegrationAnchors
  # dims: Number of dimeonsions to use in FindIntegrationAnchors
  # See https://satijalab.org/seurat/articles/integration_introduction.html and
  #     https://satijalab.org/seurat/articles/integration_rpca.html for more info
  
  # In case SCT is not default assay, change it to SCT
  for(i in 1:length(slist)){
    DefaultAssay(slist[[i]]) <- 'SCT'
  }

  print('Using reference samples:')
  print(reference_samples)
  
  print('Select integration features')
  features <- SelectIntegrationFeatures(object.list = slist, nfeatures = 3000)
  
  print('Prep SCT integration')
  slist <- PrepSCTIntegration(object.list = slist, anchor.features = features)
  
  # print('Run PCA')
  # slist <- lapply(X = slist, FUN = RunPCA, features = features,verbose=F)
  
  # reference parameter of FindIntegrationAnchors expects indexes, not names, so
  # look up the index of the samples given to reference_samples
  reference_sample_index <- which(names(slist) %in% reference_samples)
  print('Find anchors')
  # Use rpca because it is faster for large dataset, see https://satijalab.org/seurat/articles/integration_rpca.html
  Csparse_validate = "CsparseMatrix_validate"
  anchors <- FindIntegrationAnchors(object.list = slist, normalization.method = "SCT",
                                    anchor.features = features, dims = 1:dims, 
                                    k.anchor = k.anchor, reduction = 'rpca',
                                    reference=reference_sample_index)
  print('Integrate data')
  integrated.sct <- IntegrateData(anchorset = anchors, normalization.method = "SCT", dims = 1:dims)
  rm(anchors, features)
  gc() 
  return(integrated.sct)
}



integrate_cite <- function(slist,reference_samples, normalize, assay, k.anchor=5, dims=20){


  # slist =sobj_list;random_donors_for_reference=reference_samples;normalize=F; assay='CITE_bgRemoved';k.anchor=5; dims=20
                                                     
  # Integrate a list of seurat objects together using assay (e.g. CITE or CITE_bgRemoved)
  # slist: list of seurat objects
  # reference_samples: list of names to use as reference (reference parameter in FindIntegrationAnchors)
  # normalize: Whether to normalize CITE data or not
  # assay: Which assay to use
  # k.anchor: Number of anchors to use in FindIntegrationAnchors
  # dims: Number of dimeonsions to use in FindIntegrationAnchors. Lower than default because less proteins
  
  print('Using reference samples:')
  print(reference_samples)
  
  for(i in 1:length(slist)){

    DefaultAssay(slist[[i]]) <- assay
  }
  
  # Normalize date if paramter is set T

  slist_test <- lapply(X = slist, FUN = function(x) {
    if(normalize==T){
      x <- NormalizeData(x, verbose = FALSE, normalization.method='CLR', margin=2) 
    }
  })
  
  
  print('Select integration features')
  seurat.features.adt <- rownames(slist[[1]])
  
  print(paste('# of selected features before isotype removal:', length(seurat.features.adt)))
  # remove isotype controls from features as these should not be used for integration
  seurat.features.adt <- seurat.features.adt[!grepl('isotype',seurat.features.adt,
                                                    ignore.case = T)]
  print(paste('# of selected features after isotype removal:', length(seurat.features.adt)))
  print('Scale data and run PCA')
  slist <- lapply(X = slist, FUN = function(x) {
    x <- ScaleData(x, features = seurat.features.adt, verbose = FALSE)
    x <- RunPCA(x, features = seurat.features.adt, verbose = FALSE)
  })
  
  print('Find anchors')
  anchors.adt <- FindIntegrationAnchors(object.list = slist, 
                                        reduction = "rpca", 
                                        dims = 1:dims)
  
  print('Integrate data')
  integrated.adt <- IntegrateData(anchorset = anchors.adt, dims = 1:dims)
  
  # NOTE: With rpca integration for ADT we do scale after integration
  # unlike SCT integration for expression above
  print('Scale data and run PCA and umap')
  integrated.adt <- ScaleData(integrated.adt)
  return(integrated.adt)
}

names(sobj_list) <- sample_names
#####
Csparse_validate = "CsparseMatrix_validate"
pools_to_integrate <- unique(gsub('-donor[0-9]+','',names(sobj_list)))
# from each pool randomly select one donor to use as reference with one donor from all other pools
random_donors_for_reference <- c()
pools_to_use <- c()
# Set a seed so that when rerunning the code the same samples get sampled
set.seed(94235)
# For each pool, randomly sample 1 donor
for(pool in pools_to_integrate){
  random_donor <- sample(names(sobj_list)[grepl(pool, names(sobj_list))],1)
  random_donors_for_reference <- c(random_donors_for_reference, random_donor)
}

# integrate sct
print('---- Integrating SCT---')
dir.create(paste0(figdir,'/2.elbow_plots/'),showWarnings = F)

# Integrate SCT
# slist=sobj_list
# reference_samples=random_donors_for_reference
random_donor_integration_sct <- integrate_sct(sobj_list, random_donors_for_reference)
random_donor_integration_sct <- RunPCA(random_donor_integration_sct, verbose = FALSE)
# Number of dimensions to use for SCT integrated data (check elbow plot to see if it's correct)
ElbowPlot(random_donor_integration_sct, ndims=50)+geom_vline(xintercept = ndim_sct+0.5,
                                                             lty=2,colour='red')
ggsave(paste0(figdir,'/2.elbow_plots/SCT.pdf'), width=5, height=5)
random_donor_integration_sct <- RunUMAP(random_donor_integration_sct, reduction = "pca", dims = 1:ndim_sct)

print('---- Integrating CITE---')
# integrate CITE
random_donor_integration_cite <- integrate_cite(sobj_list, 
                                              random_donors_for_reference,
                                              normalize=F,
                                              assay='CITE')
# Number of dimensions to use for cite integrated data (check elbow plot to see if it's correct)
random_donor_integration_cite <- RunPCA(random_donor_integration_cite, verbose = FALSE)
ElbowPlot(random_donor_integration_cite, ndims=50)+geom_vline(xintercept = ndim_cite_integrated+0.5,
                                                              lty=2,colour='red')
random_donor_integration_cite <- RunUMAP(random_donor_integration_cite, reduction = "pca", 
                                         dims = 1:ndim_cite_integrated)


print('---- Integrating cite_bgRemoved---')
# integrate cite_bgRemoved
random_donor_integration_citeBgRemoved <- integrate_cite(sobj_list, 
                                                         random_donors_for_reference,
                                                         normalize=F,
                                                         assay='CITE_bgRemoved')
# Number of dimensions to use for cite_bgRemoved integrated data (check elbow plot to see if it's correct)
random_donor_integration_citeBgRemoved <- RunPCA(random_donor_integration_citeBgRemoved, verbose = FALSE)
ElbowPlot(random_donor_integration_citeBgRemoved, ndims=50)+geom_vline(xintercept = ndim_citeBgRemoved+0.5,
                                                                       lty=2,colour='red')
ggsave(paste0(figdir,'/2.elbow_plots/CITE_bgRemoved.pdf'), width=14, height=4)
random_donor_integration_citeBgRemoved <- RunUMAP(random_donor_integration_citeBgRemoved, reduction = "pca", 
                                                  dims = 1:ndim_citeBgRemoved)
                                              
#####
print("---- Integration done -----")
#### PCA and UMAP ####







# Seperate Code - adding datasets together
#####
dir.create(paste0(figdir,'/2.initial_umaps/'),showWarnings = F)
##### plot each umap #####
n_pools <- length(unique(random_donor_integration_sct$Pool))
n_celltype <- length(unique(random_donor_integration_sct$Azimuth.predicted.celltype.l2))

p1 <- DimPlot(random_donor_integration_sct, group.by='Pool',  cols=myPalette(n_pools), shuffle=T)
p2 <- DimPlot(random_donor_integration_sct, group.by='Azimuth.predicted.celltype.l2',  cols=myPalette(n_celltype))
p1+p2       
ggsave(paste0(figdir,'/2.initial_umaps/all_samples_sct.pdf'), width=14, height=4)


p1 <- DimPlot(random_donor_integration_citeBgRemoved, group.by='Pool',  cols=myPalette(n_pools), shuffle=T)
p2 <- DimPlot(random_donor_integration_citeBgRemoved, group.by='Azimuth.predicted.celltype.l2',  cols=myPalette(n_celltype))
p1+p2       
ggsave(paste0(figdir,'/2.initial_umaps/all_samples_citeBgRemoved.pdf'), width=14, height=4)


p1 <- DimPlot(random_donor_integration_cite, group.by='Pool',  cols=myPalette(n_pools), shuffle=T)
p2 <- DimPlot(random_donor_integration_cite, group.by='Azimuth.predicted.celltype.l2',  cols=myPalette(n_celltype))
p1+p2       
ggsave(paste0(figdir,'/2.initial_umaps/all_samples_citeBgRemoved.pdf'), width=14, height=4)
#####


##### add the assays together #####
random_donor_integration_sct[["CITE.integrated"]] <- random_donor_integration_cite[["integrated"]]
random_donor_integration_sct[["CITE_bgRemoved.integrated"]] <- random_donor_integration_citeBgRemoved[["integrated"]]
random_donor_integration_sct[["pca.cite.integrated"]] <- random_donor_integration_cite[["pca"]]
random_donor_integration_sct[["pca.cite_bgRemoved.integrated"]] <- random_donor_integration_citeBgRemoved[["pca"]]
random_donor_integration_sct[["umap.cite.integrated"]] <- random_donor_integration_cite[["umap"]]
random_donor_integration_sct[["umap.cite_bgRemoved.integrated"]] <- random_donor_integration_citeBgRemoved[["umap"]]

  # Rename SCT assays
random_donor_integration_sct <- RenameAssays(object = random_donor_integration_sct, integrated = 'SCT.integrated')
names(random_donor_integration_sct@reductions)[[which(names(random_donor_integration_sct@reductions)=='umap')]] <- 'umap.SCT.integrated'
names(random_donor_integration_sct@reductions)[[which(names(random_donor_integration_sct@reductions)=='pca')]] <- 'pca.SCT.integrated'

# remove predicted_ADT because we have our own CITE measurements
random_donor_integration_sct@assays$predicted_ADT <- NULL
#####
print('---- Done, lets save the Object ----')
saveRDS(random_donor_integration_sct, file=tmp_rds_file)

