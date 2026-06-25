library(Seurat)
library(sctransform)
library(scDblFinder)
library(dplyr)
library(ggplot2)
library(Matrix)


create_seurat_obj = function(file_h5) {
  
  data.raw <- Read10X_h5(file_h5)
  data <- CreateSeuratObject(data.raw,project = "Data",names.field=2,names.delim="-")
  
  data[["RNA"]] <- as(data[["RNA"]], "Assay")
  
  data[["percent.ribo.raw"]] <- PercentageFeatureSet(data, pattern = "^rp[sl][[:digit:]]", assay="RNA")  
  data[["percent.mt.raw"]] <-PercentageFeatureSet(data, pattern = "^mt-", assay="RNA")
  
  data_filter <- subset(data, subset = nFeature_RNA > 500 & nFeature_RNA < 7000  & nCount_RNA > 500 & 
                          nCount_RNA < 25000 & percent.mt.raw < 15)
  
  return(data_filter)
  
}

apply_transform_NEW <- function(data) {
    
  data<-SCTransform(data,return.only.var.genes = F,verbose=F)
  data<-RunPCA(data)
  
  data<-FindNeighbors(object = data, dims = 1:15, reduction = "pca")
  data<-FindClusters(object=data, resolution = 0.3)
  data<-RunUMAP(object=data,dims = 1:15)

  return(data)
}


#### Zhao et al.d dataset ####

z_4dpf = create_seurat_obj("matrix_cellbender_filtered.h5" )


sce <- as.SingleCellExperiment(z_4dpf)
dbl <- scDblFinder(sce)
z_4dpf$DoubletScore <- dbl$scDblFinder.score
z_4dpf$DoubletClass <- dbl$scDblFinder.class
z_4dpf_woDoublet =  subset(z_4dpf, subset = DoubletClass == "singlet")

z_4dpf_woDoublet <- SCTransform(z_4dpf_woDoublet,return.only.var.genes = F,verbose=F)
z_4dpf_woDoublet <- RunPCA(z_4dpf_woDoublet)

z_4dpf_woDoublet_process <- FindNeighbors(object = z_4dpf_woDoublet, dims = 1:30, reduction = "pca", k.param = 30)
z_4dpf_woDoublet_process <- FindClusters(object=z_4dpf_woDoublet_process, resolution = 0.3)

z_4dpf_woDoublet_process <- RunUMAP(object=z_4dpf_woDoublet_process ,dims = 1:30, n.neighbors = 30)


### tgfb dataset ###

WT_tgfb = create_seurat_mito_filtering("output_WT_tgfb_filtered_SEURAT.h5")
mut_tgfb = create_seurat_mito_filtering("output_mut_tgfb_filtered_SEURAT.h5")

WT_tgfb_apply = apply_transform_(WT_tgfb)
mut_tgfb_apply = apply_transform(mut_tgfb)

list_tgfb_raw = list(WT_tgfb, mut_tgfb_apply)
list_tgfb_apply = list(WT_tgfb_apply, mut_tgfb_apply)

list_tgfb_woDoublets <- lapply(seq_along(list_tgfb_apply), function(i) {
  
  sce <- as.SingleCellExperiment(list_tgfb_apply[[i]])
  dbl <- scDblFinder(sce)
  
  raw <- list_tgfb_raw[[i]]
  raw$DoubletScore <- dbl$scDblFinder.score
  raw$DoubletClass <- dbl$scDblFinder.class
  
  subset(raw, subset = DoubletClass == "singlet")
})

list_tgfb_woDoublets <- lapply(seq_along(list_tgfb_woDoublets), function(i) {
  SCTransform(list_tgfb_woDoublets[[i]], variable.features.n = 3000, return.only.var.genes = FALSE, verbose = FALSE)
})

list_tgfb_woDoublets[[1]]$Cond = "WT"
list_tgfb_woDoublets[[2]]$Cond = "Mut"

var.features <- SelectIntegrationFeatures(object.list = list_tgfb_woDoublets, nfeatures = 3000)
merged_tgfb <- merge(list_tgfb_woDoublets[[1]], y= list_tgfb_woDoublets[-1], add.cell.ids = c("WT","mut"), project = "merged_rep", merge.data = TRUE)
VariableFeatures(merged_tgfb) <- var.features

merged_tgfb <- RunPCA(merged_tgfb, verbose = FALSE)

pdsred_values <- as.matrix(GetAssayData(merged_tgfb, slot = "scale.data")["pDsRed-express-1", ])
pca_embeddings <- Embeddings(merged_tgfb, "pca")[, 1:14]
weight_factor <- 7

combined <- cbind(pdsred_values * weight_factor, pca_embeddings)
colnames(combined) <- paste0("Combined_", 1:ncol(combined))
merged_tgfb[["combined"]] <- CreateDimReducObject(embeddings = combined, key = "Combined_", assay = DefaultAssay(merged_tgfb))

merged_tgfb_process <- RunUMAP(merged_tgfb, reduction = "combined", dims = 1:14)
merged_tgfb_process  <- FindNeighbors(merged_tgfb_process, reduction = "combined", dims = 1:14) %>% FindClusters(resolution = 0.3)
