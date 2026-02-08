#!/usr/bin/env Rscript
# Add VDJ information to the seurat object

#### libraries ####
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(ggpubr)
library(immunarch)
library("RColorBrewer")
library('Seurat')
library(future)
# options(future.globals.maxSize = 60 * 1024^3)
if (future::supportsMulticore()) {
  future::plan(future::multicore)
} else {
  future::plan(future::multisession)
}

cat('Libraries are loaded')
myPallette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
cols <- myPallette(10)
####



#### set up directories, colors paths, read data ####
data_dir <- './'
# args = list() 
# args[1]='figures__regress__pct_counts_gene_group__mito_transcript'
# args[2]='out'
# args[3] ='tmp_rds_files'
# args[4]='regress__pct_counts_gene_group__mito_transcript__all_samples.wnn.integrated.RDS'
args = commandArgs(trailingOnly=TRUE)

wnn_integrated_file = args[1]


myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))


all_samples <- readRDS(wnn_integrated_file)
vdj_files <- unique(list.files('./',
                               pattern = "filtered_contig_annotations.csv", recursive = TRUE,
                               full.names = T))

# add a sample_id column to metadata
all_samples$sample_id <- gsub('.*_cellranger','cellranger',rownames(all_samples@meta.data))
#####


#### Add VDJ info to seurat object ####
tcr_names <- c()
bcr_names <- c()
immdata_TCR <- list('data'=list(), 'meta'=data.frame())
immdata_BCR <- list('data'=list(), 'meta'=data.frame())
tcr_i <- 1
bcr_i <- 1
total  <- 0
all_samples$clonotype <- ''
all_samples$vdj_type <- ''

for (file_path in vdj_files) {
  vdj_type <- ''
  # Get wether it is TCR or BCR based on file path name
  if(grepl('vdj_t',file_path) | grepl('vdj_TCR',file_path)){
    vdj_type <- 'TCR'
  }else if(grepl('vdj_b',file_path) | grepl('vdj_BCR',file_path)){
    vdj_type <- 'BCR'
  } else{stop('No VDJ type in file path')}
  
  # 'cellranger multi' has different file structure than 'cellranger count', 
  # so extract sample_id in different ways
  if(grepl('multi', file_path)){
    sample_id <- basename(dirname(dirname(file_path)))
  }else{
    sample_id <- basename(dirname(file_path))
  }

  if(file.info(file_path)$size==0){
    print(paste(file_path,'is empty, skip'))
    next
  }
  print(file_path)
  
  # Use immunearch to load vdj data, as they do a lot of VDJ preprocessing.
  # See https://immunarch.com/ for more info
  immunearch_data <- suppressWarnings(suppressMessages(repLoad(file_path)))

  immunearch_data$data$filtered_contig_annotations$Barcode_single <- paste0(gsub('-1','',
                                                                              gsub(';.*','',immunearch_data$data$filtered_contig_annotations$Barcode)),
                                                                            '_',sample_id)
  
  # Filter the cells with VDJ info to only include cells that passed QC from SCT/CITE processing
  immunearch_data$data$filtered_contig_annotations <- immunearch_data$data$filtered_contig_annotations[immunearch_data$data$filtered_contig_annotations$Barcode_single %in% all_samples@meta.data$Barcode,]
  immunearch_data$data$filtered_contig_annotations$donor <- all_samples@meta.data[match(immunearch_data$data$filtered_contig_annotations$Barcode_single,
                                                                                        all_samples@meta.data$Barcode),]$donor.vireo
  # add clonotype ID to the seurat metadata of each cell (if the cell has a clonotype)
  # this is quite slow, should rewrite this at some point
  if(nrow(immunearch_data$data$filtered_contig_annotations) != 0){
    for(i in 1:nrow(immunearch_data$data$filtered_contig_annotations)){
      row <- immunearch_data$data$filtered_contig_annotations[i,]
      for(bcode in strsplit(row$Barcode,';')[[1]]){
        bcode <- paste0(gsub('-1',paste0('_',sample_id),bcode))
        if(bcode %in% all_samples$Barcode){
          all_samples@meta.data[all_samples$Barcode==bcode,]$clonotype <- row$raw_clonotype_id
          all_samples@meta.data[all_samples$Barcode==bcode,]$vdj_type <- vdj_type
        }
      }
    }  
  }
  
  
  if(!sample_id %in% all_samples$sample_id){
    next
  }
  this_sample <- all_samples@meta.data[all_samples$sample_id==sample_id,]
  this_sample <- this_sample[1,c('Pool'),drop=F]
  rownames(this_sample) <- NULL
  
  # Make a separate metadata for each of the donors
  for(donor in unique(immunearch_data$data$filtered_contig_annotations$donor)){
    immunarch_this_donor <- immunearch_data$data$filtered_contig_annotations[which(immunearch_data$data$filtered_contig_annotations$donor==donor),]
    this_sample_with_donor <- this_sample
    
    this_sample_with_donor$vdj_type <- vdj_type
    this_sample_with_donor$sample_name <- donor
    
    colnames(this_sample_with_donor) <- c('Sample','VDJ_type','donor')
    this_sample_with_donor$Sample <- paste0(this_sample_with_donor$Sample,'_',this_sample_with_donor$donor)
    if(vdj_type=='TCR'){
      immdata_TCR$data[[tcr_i]] <- immunarch_this_donor
      immdata_TCR$meta <- rbind(immdata_TCR$meta, this_sample_with_donor)
      rownames(immdata_TCR$meta) <- immdata_TCR$Sample
      tcr_i <- tcr_i + 1
      tcr_names <- c(tcr_names, this_sample_with_donor$Sample)
      
    }else if(vdj_type=='BCR'){
      immdata_BCR$data[[bcr_i]] <- immunarch_this_donor
      immdata_BCR$meta <- rbind(immdata_BCR$meta, this_sample_with_donor)
      bcr_i <- bcr_i+1
      bcr_names <- c(bcr_names, this_sample_with_donor$Sample)
      rownames(immdata_BCR$meta) <- immdata_BCR$Sample
      
    }else{
      stop('ERROR: vdj type missing')
    }
  }
  
  total <- total+1
}

vdj_dir <- paste0('./vdj/')
dir.create(vdj_dir,showWarnings = F)
names(immdata_TCR$data) <- tcr_names
names(immdata_BCR$data) <- bcr_names

n2 = strsplit(as.character(wnn_integrated_file), split="__all_samples")[[1]][1]
outname = paste0(n2,'__all_samples_integrated.vdj.RDS')
outname_BCR = paste0(n2,'__all_samples_integrated.BCR.RDS')
outname_TCR = paste0(n2,'__all_samples_integrated.TCR.RDS')

saveRDS(immdata_BCR, file=outname_BCR)
saveRDS(immdata_TCR, file=outname_TCR)

saveRDS(all_samples, file=outname)

