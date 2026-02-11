library(Seurat)
library(sctransform)
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