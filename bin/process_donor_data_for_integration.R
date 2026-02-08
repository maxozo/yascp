#!/usr/bin/env Rscript
library("Matrix")
library(Seurat)
library(ggplot2)
library("RColorBrewer")
library(future)
options(future.globals.maxSize = 60 * 1024^3)
args = commandArgs(trailingOnly=TRUE)
# STAT3_A1_BM vireo_STAT3_A1_BM matched_donors.txt STAT3_A1_BM.withADT.RDS NONE
# args=vector(mode='list', length=6); args[1]='FOXP3_G1_BM'; args[[2]]='vireo_FOXP3_G1_BM'; args[[3]]='matched_donors.txt'; args[[4]]='FOXP3_G1_BM.withADT.RDS';args[[5]]= 'NONE'
# FOXP3_G1_BM vireo_FOXP3_G1_BM matched_donors.txt FOXP3_G1_BM.withADT.RDS NONE
sample_name <- 'STAT3_A1_BM'
vireo_input = 'vireo_STAT3_A1_BM'
matched_donor_files = 'matched_donors.txt'
rds_file = 'STAT3_A1_BM.withADT.RDS'
vars_regress ='NONE'

sample_name = args[1]
vireo_input = args[2]
matched_donor_files = args[3]
rds_file = args[4]
vars_regress = args[5]
print(sample_name)
drop_doublets = TRUE
print('vars_regress')
print(vars_regress)
print('vireo_input')
print(vireo_input)
# LDP13 vireo_LDP13  matched_donors.txt LDP13.withADT.RDS
# matched_donor_files <- 'matched_donors.txt'
matched_donors <- read.table(matched_donor_files,header=T,sep='\t')

# vireo_input = '/lustre/scratch123/hgi/teams/hgi/mo11/tmp_projects/sle_project/sle_full_run/work/77/a20c457346abcd381106c81624951b/vireo_LDP1'
# rds_file = '/lustre/scratch123/hgi/teams/hgi/mo11/tmp_projects/sle_project/sle_full_run/work/77/a20c457346abcd381106c81624951b/LDP1.withADT.RDS'

# sample_name <- gsub('.withADT|.RDS','',basename(rds_file))
sobj <- readRDS(rds_file)
sobj@assays$RNA <- sobj[!grepl("^MT-", rownames(sobj)), ]@assays$RNA
sobj@assays$RNA <- sobj[!grepl("^RP[SL]", rownames(sobj)), ]@assays$RNA
sobj$Pool <- sample_name

if (vars_regress=='NONE' ||  is.na(vars_regress)) {
  vars_to_regress = c()
  reg_name = 'regress__NONE'
}else{
  vars_to_regress = strsplit(as.character(vars_regress), split=";")
  vars_to_regress = unlist(vars_to_regress)
  reg_name = paste0('regress__',vars_regress)
}

# Add the vireo donor assignment to the seurat object
df <- read.table(paste0(vireo_input,'/donor_ids.tsv'), header=T, sep='\t')

df$sample_name <- sample_name
df$sample <- sample_name

# This is still some old code from when there were donors in multiple batches but when we didn't
# have genotypes, at some point his will be updated. Get the matched genotypes

df$doublet_logLikRatio <- NULL
if (drop_doublets){
    print('--- Dropping Doublets ---')
    df <- df[!grepl('doublet|unassigned', df$donor_id),]
}else{
    print('--- Not Dropping Doublets ---')
    df <- df[!grepl('unassigned', df$donor_id),]
    matched_donors[nrow(matched_donors) + 1,] = c(paste0("doublet",count),sample_name,"doublet")
}

df$matched.donor <- matched_donors[match(paste0(df$sample_name, 
                                                df$donor_id), 
                                        paste0(matched_donors$sample, 
                                                matched_donors$old_donor)),]$new_donor
this_donor_cells <- df

sobj@meta.data$donor.vireo <- this_donor_cells[match(paste0(rownames(sobj@meta.data),'_',sample_name),
                                                    paste0(this_donor_cells$cell,
                                                            '_',this_donor_cells$sample)),]$matched.donor

# Subset the Seurat object to exclude donors with NA in donor.vireo
sobj <- sobj[, rownames(sobj@meta.data[!is.na(sobj$donor.vireo), ])]

# Split the Seurat object by donor
sobj_per_donor <- SplitObject(sobj, split.by = "donor.vireo")
donors_to_drop <- c()

for (donor in names(sobj_per_donor)) {
    print(donor)

    # SCTransform with glmGamPoi method
    sobj_per_donor[[donor]] <- SCTransform(sobj_per_donor[[donor]], 
                                           vars.to.regress = vars_to_regress, 
                                           method = "glmGamPoi",
                                           verbose = FALSE, return.only.var.genes = FALSE)

    # Normalize and calculate cell cycle scores
    sobj_per_donor[[donor]] <- CellCycleScoring(
        NormalizeData(sobj_per_donor[[donor]]),
        s.features = cc.genes$s.genes,
        g2m.features = cc.genes$g2m.genes,
        assay = 'SCT',
        set.ident = TRUE
    )

    # Run PCA on cell cycle genes
    sobj_per_donor[[donor]] <- RunPCA(sobj_per_donor[[donor]], features = c(cc.genes$s.genes, cc.genes$g2m.genes), verbose = FALSE)
    p1 <- DimPlot(sobj_per_donor[[donor]]) + ggtitle('Before cell cycle correction (cell cycle genes only)')

    # Check number of PCs
    dims2 <- dim(sobj_per_donor[[donor]]@assays$RNA$counts)
    num_pcs <- min(20, dims2[2] - 1)
    if (num_pcs < 20) {
        donors_to_drop <- append(donors_to_drop, donor)
    }

    # Run PCA on all variable genes
    sobj_per_donor[[donor]] <- RunPCA(sobj_per_donor[[donor]], dims = 1:10, npcs = 20, verbose = FALSE)
    p2 <- DimPlot(sobj_per_donor[[donor]]) + ggtitle('Before cell cycle correction (all genes)')
    sobj_per_donor[[donor]]@project.name <- paste0(sample_name, '-', donor)
}

# Drop the specified donors
sobj_per_donor <- sobj_per_donor[!names(sobj_per_donor) %in% donors_to_drop]

# Save the result
saveRDS(sobj_per_donor, file = paste0('normalised__', rds_file))
