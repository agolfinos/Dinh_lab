
gene = "LAMP3"
p1 = SpatialFeaturePlot(data2, features = c(gene))
p2 = SpatialFeaturePlot(data3, features = c(gene))
gene = "CXCL9"

p3 = SpatialFeaturePlot(data2, features = c(gene))
p4 = SpatialFeaturePlot(data3, features = c(gene))

cowplot::plot_grid(p1, p3, p2, p4, nrow = 2)

DefaultAssay(data1) <- 'SCT'
SpatialFeaturePlot(data1, features = c("KRT17", "CXCL9", 'LAMP3', 'CD8A'))

DefaultAssay(data4) <- 'SCT'
SpatialFeaturePlot(data4, features = c("LAMP3"))

DefaultAssay(data3) <- 'SCT'
SpatialFeaturePlot(data3, features = c("KRT17", "CXCL9", 'LAMP3', 'CD8A'))





library(Seurat)
library(tidyverse)
library(scales)
library(pals)
library(cowplot)
library(fgsea)
library(msigdbr)
library(data.table)
library(mvhspatialplots)

setwd("/Volumes/hqdinh2-1/My_Projects/HNC_SPORE/")

data1 <- readRDS("/Volumes/hqdinh2//Projects/RawData_FromUWBiotech/HNCVisium_2022-05-05/CK17-27_analysis/ck17_27_seurat.RDS")
data2 <- readRDS("/Volumes/hqdinh2/Projects/RawData_FromUWBiotech/HNCVisium_2022-05-05/CK17-5_analysis/ck17_5_seurat.RDS")
data3 <- readRDS("/Volumes/hqdinh2/Projects/RawData_FromUWBiotech/HNCVisium_2022-05-05/CK17-7_analysis/ck17_7_seurat.RDS")
data4 <- readRDS("/Volumes/hqdinh2/Projects/RawData_FromUWBiotech/HNCVisium_2022-05-05/CK17-21_analysis/ck17_21_seurat.RDS")

DefaultAssay(data1) <- 'SCT'
DefaultAssay(data2) <- 'SCT'
DefaultAssay(data3) <- 'SCT'
DefaultAssay(data4) <- 'SCT'

# Run SCTransform integration workflow: https://satijalab.org/seurat/archive/v3.0/integration.html

data1 <- FindVariableFeatures(data1)
data2 <- FindVariableFeatures(data2)
data3 <- FindVariableFeatures(data3)
data4 <- FindVariableFeatures(data4)

ABC.list <- list(data1, data2, data3, data4)

genes.common <- Reduce(intersect, list(rownames(data1), rownames(data2), rownames(data3), rownames(data4)))
list.features <- SelectIntegrationFeatures(object.list = ABC.list,
                                           nfeatures = 2000,
                                           assay = c("SCT", "SCT", "SCT", "SCT"))

options(future.globals.maxSize = 24000 * 1024^2)

ABC.list <- PrepSCTIntegration(object.list = ABC.list,
                               anchor.features = list.features,
                               assay = "SCT",
                               verbose = F)
ABC.anchors <- FindIntegrationAnchors(object.list = ABC.list,
                                      normalization.method = "SCT",
                                      anchor.features = list.features,
                                      verbose = F)                       
save(ABC.anchors, file = "tmp.rda")
ABC <- IntegrateData(anchorset = ABC.anchors,
                     features.to.integrate = genes.common,
                     normalization.method = "SCT", 
                     verbose = F)
ABC <- RunPCA(ABC, verbose = F)
ABC <- FindNeighbors(ABC, reduction = "pca", dims = 1:30)
ABC <- FindClusters(ABC, resolution = 0.8, verbose = F)
ABC <- RunUMAP(ABC, reduction = "pca", dims = 1:30)

  