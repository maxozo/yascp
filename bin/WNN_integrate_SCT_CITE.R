#!/usr/bin/env Rscript
# Use weighted nearest neighbours to integrate SCT with CITE and SCT with CITE(background removed)

#### libraries ####
library(ggbeeswarm)
library(Seurat)
library(ggplot2)
library("RColorBrewer")
library(tidyverse)
library(patchwork)
library(future)
# options(future.globals.maxSize= 128849018880)
if (future::supportsMulticore()) {
  future::plan(future::multicore)
} else {
  future::plan(future::multisession)
}
library("RColorBrewer")
#####
# .libPaths('/lustre/scratch123/hgi/teams/hgi/mo11/tmp_projects/jaguar_yascp/nieks_pipeline/yascp_run/work/6e/dd58d21d32b7a0eefad878415f35e6/pack')
# remotes::install_version("Matrix", "1.6.1")
# library(Matrix)
#### set up directories, colors paths ####

# outdir <- paste0('./out')
# dir.create(outdir,showWarnings = F)
# figdir <- paste0(outdir,'/figures/')
# dir.create(figdir,showWarnings = F)
# tmp_rds_dir <- paste0(outdir,'/tmp_rds_files/')
# dir.create(tmp_rds_dir,showWarnings = F)
# tmp_rds_file <- paste0(tmp_rds_dir,'all_samples_integrated.RDS')
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))


args = commandArgs(trailingOnly=TRUE)
tmp_rds_file = args[1]
n2 = strsplit(as.character(tmp_rds_file), split="__all_samples")[[1]][1]
outname = paste0(n2,'__all_samples.wnn.integrated.RDS')

integrated_data <- readRDS(tmp_rds_file)
myPallette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
sct_ndim <- 12
cite_ndim <- 9
cite_bgRemoved_ndim <- 10
#####

#### Double check that correct dimensions are chose ####
ElbowPlot(integrated_data, ndims=50, reduction='pca.SCT.integrated')+geom_vline(xintercept = sct_ndim+0.5,
                                                             lty=2,colour='red')

ElbowPlot(integrated_data, ndims=50, reduction='pca.cite.integrated')+geom_vline(xintercept = cite_ndim+0.5,
                                                                                lty=2,colour='red')
ElbowPlot(integrated_data, ndims=50, reduction='pca.cite_bgRemoved.integrated')+geom_vline(xintercept = cite_bgRemoved_ndim+0.5,
                                                                                 lty=2,colour='red')
#####

# Identify multimodal neighbors. These will be stored in the neighbors slot, 
# and can be accessed using integrated_data[['weighted.nn.sct_cite']]
# The WNN graph can be accessed at integrated_data[["wknn.sct_cite"]], 
# and the SNN graph used for clustering at integrated_data[["wsnn.sct_cite"]]
# Cell-specific modality weights can be accessed at integrated_data$RNA.weight
Csparse_validate = "CsparseMatrix_validate"
integrated_data <- FindMultiModalNeighbors(
  integrated_data, reduction.list = list("pca.SCT.integrated", "pca.cite.integrated"), 
  dims.list = list(1:sct_ndim, 1:cite_ndim), 
  modality.weight.name = "weight.SCT_cite",
  weighted.nn.name = 'weighted.nn.sct_cite',
  knn.graph.name = 'wknn.sct_cite',
  snn.graph.name = 'wsnn.sct_cite'
)


# Do the same also for the CITE bg removed data (using the relevant reduction)
integrated_data <- FindMultiModalNeighbors(
  integrated_data, reduction.list = list("pca.SCT.integrated", "pca.cite_bgRemoved.integrated"), 
  dims.list = list(1:sct_ndim, 1:cite_bgRemoved_ndim), modality.weight.name = "weight.SCT_citeBgRemoved",
  weighted.nn.name = 'weighted.nn.sct_citeBgRemoved',
  knn.graph.name = 'wknn.sct_citeBgRemoved',
  snn.graph.name = 'wsnn.sct_citeBgRemoved'
)

#### calculate umaps ####
integrated_data <- RunUMAP(integrated_data, nn.name = "weighted.nn.sct_cite", 
                           reduction.name = "wnn.umap.SCT_cite", reduction.key = "wnnUMAPsctCITE_")
integrated_data <- RunUMAP(integrated_data, nn.name = "weighted.nn.sct_citeBgRemoved", 
                           reduction.name = "wnn.umap.SCT_citeBgRemoved", reduction.key = "wnnUMAPsctCITEbgRemoved_")
####

#### plot all different umaps ####
# theme that removes the x-axis ticks and labels
blank_theme <-  theme(axis.title=element_blank(),
                     axis.text=element_blank(),
                     axis.ticks=element_blank())

p1 <- DimPlot(integrated_data, reduction = 'umap.SCT.integrated', group.by = 'Azimuth.predicted.celltype.l1', cols=myPalette(8)) + NoLegend()+
  blank_theme+ggtitle('SCT umap')
p2 <- DimPlot(integrated_data, reduction = 'umap.cite.integrated', group.by = 'Azimuth.predicted.celltype.l1',  cols=myPalette(8)) + NoLegend()+
  blank_theme+ggtitle('CITE umap')
p3 <- DimPlot(integrated_data, reduction = 'umap.cite_bgRemoved.integrated', group.by = 'Azimuth.predicted.celltype.l1',  cols=myPalette(8))+
  blank_theme+ggtitle('CITE(background removed) umap')
p4 <- DimPlot(integrated_data, reduction = 'wnn.umap.SCT_cite', group.by = 'Azimuth.predicted.celltype.l1', cols=myPalette(8)) + NoLegend()+
  blank_theme+ggtitle('SCT + CITE umap')
p5 <- DimPlot(integrated_data, reduction = 'wnn.umap.SCT_citeBgRemoved', group.by = 'Azimuth.predicted.celltype.l1',  cols=myPalette(8)) + NoLegend()+
  blank_theme+ggtitle('SCT + CITE(background removed) umap')

p1+p2+p3+p4+p5
#####



# Identify clusters using low resolution
# FindNeighbors has only been run for SCT+CITE and SCT+cite(bg remove),
# so also run FindNeighbors for the other 
DefaultAssay(integrated_data) <- 'SCT'
ElbowPlot(integrated_data, reduction='pca.SCT.integrated', ndims=50)+geom_vline(xintercept=12.5, colour='red',lty=2)
integrated_data <- FindNeighbors(integrated_data, reduction='pca.SCT.integrated', dims=1:12)
DefaultAssay(integrated_data) <- 'CITE.integrated'
integrated_data <- FindNeighbors(integrated_data, reduction='pca.cite.integrated',
                             graph=c('CITE.integrated_nn', 'CITE.integrated_snn'))
DefaultAssay(integrated_data) <- 'CITE_bgRemoved.integrated'
integrated_data <- FindNeighbors(integrated_data, reduction='pca.cite_bgRemoved.integrated',
                             graph=c('CITE_bgRemoved.integrated_nn', 
                                     'CITE_bgRemoved.integrated_snn'))

for(graph in names(integrated_data@graphs)[grepl('snn',names(integrated_data@graphs))]){
  print(graph)
  integrated_data <- FindClusters(object = integrated_data, resolution = 0.5, 
                              algorithm = 3, 
                              graph.name = graph, 
                              verbose = F)
  integrated_data@meta.data[[paste0(graph,'_res.0.5')]] <- integrated_data$seurat_clusters
  integrated_data$seurat_clusters <- NULL
  print(paste0('New column added: ',colnames(integrated_data@meta.data)[ncol(integrated_data@meta.data)]))
}



saveRDS(integrated_data, file=outname)




# sobj = integrated_data
# azimuth_dir <- paste0(outdir,'/azimuth/')
# dir.create(azimuth_dir,showWarnings = F)

# unique_donors <- unique(paste(sobj$Pool, sobj$donor.vireo))
# donor_rename <- data.frame('pool_donor'=unique_donors, 
#                             'new_name'=paste0('donor_',1:length(unique_donors)))
# sobj$new_donor_name <- donor_rename[match(paste(sobj$Pool, sobj$donor.vireo),
#                                                   donor_rename$pool_donor),]$new_name

# count_l1 <- all_samples@meta.data %>% group_by(new_donor_name) %>% count(Azimuth.predicted.celltype.l1)
# count_l2 <- all_samples@meta.data %>% group_by(new_donor_name) %>% count(Azimuth.predicted.celltype.l2)
# count_l3 <- all_samples@meta.data %>% group_by(new_donor_name) %>% count(Azimuth.predicted.celltype.l3)

# make_wide_table <- function(long_table){
#   wide_table <- data.frame(spread(long_table, key = colnames(long_table)[2], value = n))
#   wide_table[is.na(wide_table)] <- 0
#   rownames(wide_table) <- wide_table$new_donor_name
#   wide_table$new_donor_name <- NULL 
#   wide_table$total <- rowSums(wide_table)
#   return(wide_table)
# }

# count_l1_wide <- make_wide_table(count_l1)
# count_l2_wide <- make_wide_table(count_l2)
# count_l3_wide <- make_wide_table(count_l3)
# write.table(count_l1_wide, paste0(azimuth_dir,'/azimuth.l1.countPerDonor.tsv'),
#             sep='\t',col.names=NA, quote=F)
# write.table(count_l2_wide, paste0(azimuth_dir,'/azimuth.l2.countPerDonor.tsv'),
#             sep='\t',col.names=NA, quote=F)
# write.table(count_l3_wide, paste0(azimuth_dir,'/azimuth.l3.countPerDonor.tsv'),
#             sep='\t',col.names=NA, quote=F)


# dir.create(azimuth_dir,'figures/',showWarnings = F)
# count_l1$total <- count_l1_wide[match(count_l1$new_donor_name, rownames(count_l1_wide)),]$total
# ggplot(count_l1, aes(Azimuth.predicted.celltype.l1, n, fill=total))+
#   geom_quasirandom(shape=21)+
#   geom_hline(yintercept=150, colour='red',lty=2)+
#   scale_fill_binned(type="viridis")+
#   theme_pubr()+
#   theme(legend.position='right')
# ggsave(paste0(azimuth_dir,'figures/azimuth_count_l1.png'), width=8, height=5)

# count_l2$total <- count_l2_wide[match(count_l2$new_donor_name, rownames(count_l2_wide)),]$total
# ggplot(count_l2, aes(Azimuth.predicted.celltype.l2, n, fill=total))+
#   geom_quasirandom(shape=21)+
#   geom_hline(yintercept=150, colour='red',lty=2)+
#   scale_fill_binned(type="viridis")+
#   theme_pubr()+
#   theme(legend.position='right',
#         axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# ggsave(paste0(azimuth_dir,'figures/azimuth_count_l2.png'), width=24, height=5)

# count_l3$total <- count_l3_wide[match(count_l3$new_donor_name, rownames(count_l3_wide)),]$total
# ggplot(count_l3, aes(Azimuth.predicted.celltype.l3, n, fill=total))+
#   geom_quasirandom(shape=21)+
#   geom_hline(yintercept=150, colour='red',lty=2)+
#   scale_fill_binned(type="viridis")+
#   theme_pubr()+
#   theme(legend.position='right',
#         axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# ggsave(paste0(azimuth_dir,'figures/azimuth_count_l3.png'), width=36, height=5)
