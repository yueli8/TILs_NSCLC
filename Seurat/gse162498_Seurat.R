library(Seurat)
library(ggplot2)
library(cowplot)
library(scater)
library(scran)
library(BiocParallel)
library(BiocNeighbors)
library(data.table)
library(dplyr)
library(Matrix)
library(clustree)
setwd("~/gse162498")

a1 <- Read10X(data.dir = "~/gse162498/P34_Tumor_raw_feature_bc_matrix")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
#pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
#Determine the ‘dimensionality’ of the dataset
#pbmc <- JackStraw(pbmc, num.replicate = 100,dims = 40)
#pbmc <- ScoreJackStraw(pbmc, dims = 1:40)
#JackStrawPlot(pbmc, dims = 1:40)
#ElbowPlot(pbmc,ndims = 40)
#確定下面的dims
pbmc <- RunUMAP(pbmc, dims = 1:20)
saveRDS(pbmc, file = "P34_T.rds")

a1 <- Read10X(data.dir = "~/gse162498/P35_Tumor_raw_feature_bc_matrix")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
#pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
#Determine the ‘dimensionality’ of the dataset
#pbmc <- JackStraw(pbmc, num.replicate = 100,dims = 40)
#pbmc <- ScoreJackStraw(pbmc, dims = 1:40)
#JackStrawPlot(pbmc, dims = 1:40)
#ElbowPlot(pbmc,ndims = 40)
#確定下面的dims
pbmc <- RunUMAP(pbmc, dims = 1:20)
saveRDS(pbmc, file = "P35_T.rds")

a1 <- Read10X(data.dir = "~/gse162498/P42_Tumor_raw_feature_bc_matrix")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
#pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
#Determine the ‘dimensionality’ of the dataset
#pbmc <- JackStraw(pbmc, num.replicate = 100,dims = 40)
#pbmc <- ScoreJackStraw(pbmc, dims = 1:40)
#JackStrawPlot(pbmc, dims = 1:40)
#ElbowPlot(pbmc,ndims = 40)
#確定下面的dims
pbmc <- RunUMAP(pbmc, dims = 1:20)
saveRDS(pbmc, file = "P42_T.rds")

a1 <- Read10X(data.dir = "~/gse162498/P43_Tumor_raw_feature_bc_matrix")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
#pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
#Determine the ‘dimensionality’ of the dataset
#pbmc <- JackStraw(pbmc, num.replicate = 100,dims = 40)
#pbmc <- ScoreJackStraw(pbmc, dims = 1:40)
#JackStrawPlot(pbmc, dims = 1:40)
#ElbowPlot(pbmc,ndims = 40)
#確定下面的dims
pbmc <- RunUMAP(pbmc, dims = 1:20)
saveRDS(pbmc, file = "P43_T.rds")

a1 <- Read10X(data.dir = "~/gse162498/P46_Tumor_raw_feature_bc_matrix")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
#pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
#Determine the ‘dimensionality’ of the dataset
#pbmc <- JackStraw(pbmc, num.replicate = 100,dims = 40)
#pbmc <- ScoreJackStraw(pbmc, dims = 1:40)
#JackStrawPlot(pbmc, dims = 1:40)
#ElbowPlot(pbmc,ndims = 40)
#確定下面的dims
pbmc <- RunUMAP(pbmc, dims = 1:20)
saveRDS(pbmc, file = "P46_T.rds")

a1 <- Read10X(data.dir = "~/gse162498/P47_Tumor_raw_feature_bc_matrix")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
#pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
#Determine the ‘dimensionality’ of the dataset
#pbmc <- JackStraw(pbmc, num.replicate = 100,dims = 40)
#pbmc <- ScoreJackStraw(pbmc, dims = 1:40)
#JackStrawPlot(pbmc, dims = 1:40)
#ElbowPlot(pbmc,ndims = 40)
#確定下面的dims
pbmc <- RunUMAP(pbmc, dims = 1:20)
saveRDS(pbmc, file = "P47_T.rds")

a1 <- Read10X(data.dir = "~/gse162498/P55_Tumor_raw_feature_bc_matrix")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
#pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
#Determine the ‘dimensionality’ of the dataset
#pbmc <- JackStraw(pbmc, num.replicate = 100,dims = 40)
#pbmc <- ScoreJackStraw(pbmc, dims = 1:40)
#JackStrawPlot(pbmc, dims = 1:40)
#ElbowPlot(pbmc,ndims = 40)
#確定下面的dims
pbmc <- RunUMAP(pbmc, dims = 1:20)
saveRDS(pbmc, file = "P55_T.rds")

a1 <- Read10X(data.dir = "~/gse162498/P57_Tumor_raw_feature_bc_matrix")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
#pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
#Determine the ‘dimensionality’ of the dataset
#pbmc <- JackStraw(pbmc, num.replicate = 100,dims = 40)
#pbmc <- ScoreJackStraw(pbmc, dims = 1:40)
#JackStrawPlot(pbmc, dims = 1:40)
#ElbowPlot(pbmc,ndims = 40)
#確定下面的dims
pbmc <- RunUMAP(pbmc, dims = 1:20)
saveRDS(pbmc, file = "P57_T.rds")

a1 <- Read10X(data.dir = "~/gse162498/P57_Blood_raw_feature_bc_matrix")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
#pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
#Determine the ‘dimensionality’ of the dataset
#pbmc <- JackStraw(pbmc, num.replicate = 100,dims = 40)
#pbmc <- ScoreJackStraw(pbmc, dims = 1:40)
#JackStrawPlot(pbmc, dims = 1:40)
#ElbowPlot(pbmc,ndims = 40)
#確定下面的dims
pbmc <- RunUMAP(pbmc, dims = 1:20)
saveRDS(pbmc, file = "P57_B.rds")

a1 <- Read10X(data.dir = "~/gse162498/P58_Tumor_raw_feature_bc_matrix")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
#pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
#Determine the ‘dimensionality’ of the dataset
#pbmc <- JackStraw(pbmc, num.replicate = 100,dims = 40)
#pbmc <- ScoreJackStraw(pbmc, dims = 1:40)
#JackStrawPlot(pbmc, dims = 1:40)
#ElbowPlot(pbmc,ndims = 40)
#確定下面的dims
pbmc <- RunUMAP(pbmc, dims = 1:20)
saveRDS(pbmc, file = "P58_T.rds")

a1 <- Read10X(data.dir = "~/gse162498/P58_Blood_raw_feature_bc_matrix")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
#pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
#Determine the ‘dimensionality’ of the dataset
#pbmc <- JackStraw(pbmc, num.replicate = 100,dims = 40)
#pbmc <- ScoreJackStraw(pbmc, dims = 1:40)
#JackStrawPlot(pbmc, dims = 1:40)
#ElbowPlot(pbmc,ndims = 40)
#確定下面的dims
pbmc <- RunUMAP(pbmc, dims = 1:20)
saveRDS(pbmc, file = "P58_B.rds")

a1 <- Read10X(data.dir = "~/gse162498/P60_Blood_raw_feature_bc_matrix")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
#pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
#Determine the ‘dimensionality’ of the dataset
#pbmc <- JackStraw(pbmc, num.replicate = 100,dims = 40)
#pbmc <- ScoreJackStraw(pbmc, dims = 1:40)
#JackStrawPlot(pbmc, dims = 1:40)
#ElbowPlot(pbmc,ndims = 40)
#確定下面的dims
pbmc <- RunUMAP(pbmc, dims = 1:20)
saveRDS(pbmc, file = "P60_B.rds")

a1 <- Read10X(data.dir = "~/gse162498/P60_Tumor_raw_feature_bc_matrix")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
#pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
#Determine the ‘dimensionality’ of the dataset
#pbmc <- JackStraw(pbmc, num.replicate = 100,dims = 40)
#pbmc <- ScoreJackStraw(pbmc, dims = 1:40)
#JackStrawPlot(pbmc, dims = 1:40)
#ElbowPlot(pbmc,ndims = 40)
#確定下面的dims
pbmc <- RunUMAP(pbmc, dims = 1:20)
saveRDS(pbmc, file = "P60_T.rds")

a1 <- Read10X(data.dir = "~/gse162498/P60_Juxta_raw_feature_bc_matrix")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
#pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
#Determine the ‘dimensionality’ of the dataset
#pbmc <- JackStraw(pbmc, num.replicate = 100,dims = 40)
#pbmc <- ScoreJackStraw(pbmc, dims = 1:40)
#JackStrawPlot(pbmc, dims = 1:40)
#ElbowPlot(pbmc,ndims = 40)
#確定下面的dims
pbmc <- RunUMAP(pbmc, dims = 1:20)
saveRDS(pbmc, file = "P60_J.rds")

a1 <- Read10X(data.dir = "~/gse162498/P61_Blood_raw_feature_bc_matrix")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
#pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
#Determine the ‘dimensionality’ of the dataset
#pbmc <- JackStraw(pbmc, num.replicate = 100,dims = 40)
#pbmc <- ScoreJackStraw(pbmc, dims = 1:40)
#JackStrawPlot(pbmc, dims = 1:40)
#ElbowPlot(pbmc,ndims = 40)
#確定下面的dims
pbmc <- RunUMAP(pbmc, dims = 1:20)
saveRDS(pbmc, file = "P61_B.rds")

a1 <- Read10X(data.dir = "~/gse162498/P61_Tumor_raw_feature_bc_matrix")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
#pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
#Determine the ‘dimensionality’ of the dataset
#pbmc <- JackStraw(pbmc, num.replicate = 100,dims = 40)
#pbmc <- ScoreJackStraw(pbmc, dims = 1:40)
#JackStrawPlot(pbmc, dims = 1:40)
#ElbowPlot(pbmc,ndims = 40)
#確定下面的dims
pbmc <- RunUMAP(pbmc, dims = 1:20)
saveRDS(pbmc, file = "P61_T.rds")

a1 <- Read10X(data.dir = "~/gse162498/P61_Juxta_raw_feature_bc_matrix")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
#pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
#Determine the ‘dimensionality’ of the dataset
#pbmc <- JackStraw(pbmc, num.replicate = 100,dims = 40)
#pbmc <- ScoreJackStraw(pbmc, dims = 1:40)
#JackStrawPlot(pbmc, dims = 1:40)
#ElbowPlot(pbmc,ndims = 40)
#確定下面的dims
pbmc <- RunUMAP(pbmc, dims = 1:20)
saveRDS(pbmc, file = "P61_J.rds")

P34_T<-readRDS(file="P34_T.rds")
P35_T<-readRDS(file="P35_T.rds")
P42_T<-readRDS(file="P42_T.rds")
P43_T<-readRDS(file="P43_T.rds")
P46_T<-readRDS(file="P46_T.rds")
P47_T<-readRDS(file="P47_T.rds")
P55_T<-readRDS(file="P55_T.rds")
P57_B<-readRDS(file="P57_B.rds")
P57_T<-readRDS(file="P57_T.rds")
P58_B<-readRDS(file="P58_B.rds")
P58_T<-readRDS(file="P58_T.rds")
P60_B<-readRDS(file="P60_B.rds")
P60_J<-readRDS(file="P60_J.rds")
P60_T<-readRDS(file="P60_T.rds")
P61_B<-readRDS(file="P61_B.rds")
P61_J<-readRDS(file="P61_J.rds")
P61_T<-readRDS(file="P61_T.rds")

P34_T<-RenameCells(P34_T,add.cell.id="P34_T",for.merge=T)
P34_T@meta.data$tech<-"Tumor"
P34_T@meta.data$celltype<-"Tumor_P34"

P35_T<-RenameCells(P35_T,add.cell.id="P35_T",for.merge=T)
P35_T@meta.data$tech<-"Tumor"
P35_T@meta.data$celltype<-"Tumor_P35"

P42_T<-RenameCells(P42_T,add.cell.id="P42_T",for.merge=T)
P42_T@meta.data$tech<-"Tumor"
P42_T@meta.data$celltype<-"Tumor_P42"

P43_T<-RenameCells(P43_T,add.cell.id="P43_T",for.merge=T)
P43_T@meta.data$tech<-"Tumor"
P43_T@meta.data$celltype<-"Tumor_P43"

P46_T<-RenameCells(P46_T,add.cell.id="P46_T",for.merge=T)
P46_T@meta.data$tech<-"Tumor"
P46_T@meta.data$celltype<-"Tumor_P46"

P47_T<-RenameCells(P47_T,add.cell.id="P47_T",for.merge=T)
P47_T@meta.data$tech<-"Tumor"
P47_T@meta.data$celltype<-"Tumor_P47"

P55_T<-RenameCells(P55_T,add.cell.id="P55_T",for.merge=T)
P55_T@meta.data$tech<-"Tumor"
P55_T@meta.data$celltype<-"Tumor_P55"

P57_B<-RenameCells(P57_B,add.cell.id="P57_B",for.merge=T)
P57_B@meta.data$tech<-"Blood"
P57_B@meta.data$celltype<-"Blood_P57"

P57_T<-RenameCells(P57_T,add.cell.id="P57_T",for.merge=T)
P57_T@meta.data$tech<-"Tumor"
P57_T@meta.data$celltype<-"Tumor_P57"

P58_B<-RenameCells(P58_B,add.cell.id="P58_B",for.merge=T)
P58_B@meta.data$tech<-"Blood"
P58_B@meta.data$celltype<-"Blood_P58"

P58_T<-RenameCells(P58_T,add.cell.id="P58_T",for.merge=T)
P58_T@meta.data$tech<-"Tumor"
P58_T@meta.data$celltype<-"Tumor_P58"

P60_B<-RenameCells(P60_B,add.cell.id="P60_B",for.merge=T)
P60_B@meta.data$tech<-"Blood"
P60_B@meta.data$celltype<-"Blood_P60"

P60_J<-RenameCells(P60_J,add.cell.id="P60_J",for.merge=T)
P60_J@meta.data$tech<-"Juxta"
P60_J@meta.data$celltype<-"Juxta_P60"

P60_T<-RenameCells(P60_T,add.cell.id="P60_T",for.merge=T)
P60_T@meta.data$tech<-"Tumor"
P60_T@meta.data$celltype<-"Tumor_P60"

P61_B<-RenameCells(P61_B,add.cell.id="P61_B",for.merge=T)
P61_B@meta.data$tech<-"Blood"
P61_B@meta.data$celltype<-"Blood_P61"

P61_J<-RenameCells(P61_J,add.cell.id="P61_J",for.merge=T)
P61_J@meta.data$tech<-"Juxta"
P61_J@meta.data$celltype<-"Juxta_P61"

P61_T<-RenameCells(P61_T,add.cell.id="P61_T",for.merge=T)
P61_T@meta.data$tech<-"Tumor"
P61_T@meta.data$celltype<-"Tumor_P61"

P5758_T<-merge(P57_T,P58_T)
P6061_T<-merge(P60_T,P61_T)
P57586061_T<-merge(P5758_T,P6061_T)

P5758_B<-merge(P57_B,P58_B)
P6061_B<-merge(P60_B,P61_B)
P57586061_B<-merge(P5758_B,P6061_B)

P57586061<-merge(P57586061_T,P57586061_B)

library(Seurat)
library(ggplot2)

saveRDS(P57586061, file="P57586061_before_integrate.rds")
hms<-P57586061

hms<-readRDS(file="P57586061_before_integrate.rds")
#before integrate
hms[["percent.mt"]] <- PercentageFeatureSet(hms, pattern = "^Mt-")
VlnPlot(hms, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
pancreas <- NormalizeData(object = hms, normalization.method = "LogNormalize", scale.factor = 1e4)
pancreas <- FindVariableFeatures(pancreas, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
pancreas <- ScaleData(pancreas, verbose = FALSE)
pancreas <- RunPCA(pancreas, npcs = 30, verbose = FALSE)
pancreas <- RunUMAP(pancreas, reduction = "pca", dims = 1:30)
p1 <- DimPlot(pancreas, reduction = "umap", group.by = "tech")
p2 <- DimPlot(pancreas, reduction = "umap", group.by = "celltype", label = TRUE, repel = TRUE) + 
  NoLegend()
plot_grid(p1,p2)
#integrate
pancreas.list <- SplitObject(pancreas, split.by = "celltype")
for (i in 1: length(pancreas.list)) {
  pancreas.list[[i]] <- NormalizeData(pancreas.list[[i]], verbose = FALSE)
  pancreas.list[[i]] <- FindVariableFeatures(pancreas.list[[i]], selection.method = "vst", nfeatures = 2000, 
                                             verbose = FALSE)
}
#reference.list <- pancreas.list[c("Blood_P57","Tumor_P57","Blood_P58","Tumor_P58", 
#"Blood_P60","Tumor_P60", "Blood_P61", "Tumor_P61","Juxta_P60","Juxta_P61")]
reference.list <- pancreas.list[c("Tumor_P57","Tumor_P58","Tumor_P60","Tumor_P61",
                                  "Blood_P57","Blood_P58", "Blood_P60", "Blood_P61" )]
pancreas.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, dims = 1:30)
DefaultAssay(pancreas.integrated) <- "integrated"
pancreas.integrated <- ScaleData(pancreas.integrated, verbose = FALSE)
pancreas.integrated <- RunPCA(pancreas.integrated, npcs = 30, verbose = FALSE)
pancreas.integrated <- RunUMAP(pancreas.integrated, reduction = "pca", dims = 1:30)
p1 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "tech")
p2 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "celltype")
plot_grid(p1,p2)
saveRDS(pancreas.integrated, file = "hms_after_integrated.rds")
saveRDS(pancreas.integrated, file = "P57586061_after_integrated.rds")

hms_individual_integrated<-readRDS(file="hms_after_integrated.rds")
p1 <- DimPlot(hms_individual_integrated, reduction = "umap", group.by = "celltype")
p1

#pbmc <- JackStraw(hms_individual_integrated, num.replicate = 100,dims = 40)
#pbmc <- ScoreJackStraw(pbmc, dims = 1:30)
#JackStrawPlot(pbmc, dims = 1:30)
#ElbowPlot(pbmc,ndims = 30)
#find how many 15cluster
#ElbowPlot(hms_individual_integrated)
hms_neighbor<- FindNeighbors(hms_individual_integrated, dims = 1:30)
obj <- FindClusters(hms_neighbor, resolution = seq(0.5,1.2,by=0.1))
#resolution設置在0.4-1.2之間,越大clusters越多,查看拐點
clustree(obj)
hms_cluster <- FindClusters( hms_neighbor, resolution = 1.2)
head(Idents(hms_cluster), 5)
hms_cluster<- RunUMAP(hms_cluster, dims = 1:30)
DimPlot(hms_cluster, reduction = "umap")
saveRDS(hms_cluster, file = "hms_cluster_test_1.2.rds")

scRNA.markers <- FindAllMarkers(hms_cluster, 
                                only.pos = TRUE,  #特异性高表达marker
                                min.pct = 0.05, 
                                logfc.threshold = 0.05
)
write.table(scRNA.markers,file="cellMarkers.txt",sep="\t",row.names=F,quote=F)

#挑选每个细胞亚群中特意高表达的20个基因
top20 <- scRNA.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC) 
write.csv(file="top20_cell_markers.csv",top20)
#整理成表格，只显示基因名字
top20_table=unstack(top20, gene ~ cluster)
names(top20_table)=gsub("X","cluster",names(top20_table))
write.csv(file="top20_marker_genes_1.2.csv",top20_table,row.names=F)

#细胞及细胞中基因与RNA数量
slotNames(hms_cluster)
#assay
hms_cluster@assays
dim(hms_cluster@meta.data)
View(hms_cluster@meta.data)

hms_cluster<-readRDS(file="hms_cluster_test_1.2.rds")
DimPlot(hms_cluster, reduction = "umap")

new.cluster.ids <- c("T", "Epithelial","T", "T","T", "T","T", "T","B",
                     "T", "T","T", "T","T", "T","T", "T","Epithelial",
                     "Epithelial","Macrophage","T", "T","Macrophage","Epithelial",
                     "T","Monocyte","T", "T","Natural_Killer","B","B","Progenitor") 

names(new.cluster.ids) <- levels(hms_cluster)
hms_cluster_id<- RenameIdents(hms_cluster, new.cluster.ids)
DimPlot(hms_cluster_id, reduction = "umap", label = TRUE, pt.size = 0.5) 
DimPlot(hms_cluster_id, reduction = "umap", label = FALSE, pt.size = 0.5) 
saveRDS(hms_cluster_id, file = "hms_cluster_id_test.rds")

hms_cluster_id<-readRDS("hms_cluster_id_test.rds")

T<-subset(hms_cluster_id, idents=c('T'))
DimPlot(T, reduction = "umap")
saveRDS(T, file="T.rds")

Epithelial<-subset(hms_cluster_id, idents=c('Epithelial'))
DimPlot(Epithelial, reduction = "umap")
saveRDS(Epithelial, file="Epithelial.rds")

B<-subset(hms_cluster_id, idents=c('B'))
DimPlot(B, reduction = "umap")
saveRDS(B, file="B.rds")

Secretory<-subset(hms_cluster_id, idents=c('Secretory'))
DimPlot(Secretory, reduction = "umap")
saveRDS(Secretory, file="Secretory.rds")

Macrophage<-subset(hms_cluster_id, idents=c('Macrophage'))
DimPlot(Macrophage, reduction = "umap")
saveRDS(Macrophage, file="Macrophage.rds")

Monocyte<-subset(hms_cluster_id, idents=c('Monocyte'))
DimPlot(Monocyte, reduction = "umap")
saveRDS(Monocyte, file="Monocyte.rds")

Natural_Killer<-subset(hms_cluster_id, idents=c('Natural_Killer'))
DimPlot(Natural_Killer, reduction = "umap")
saveRDS(Natural_Killer, file="Natural_Killer.rds")

Progenitor<-subset(hms_cluster_id, idents=c('Progenitor'))
DimPlot(Progenitor, reduction = "umap")
saveRDS(Progenitor, file="Progenitor.rds")

#input each cluster
T<-readRDS("T.rds")
Epithelial<-readRDS("Epithelial.rds")
B_cell<-readRDS("B.rds")
Secretory<-readRDS("Secretory.rds")
Macrophage<-readRDS("Macrophage.rds")
Monocyte<-readRDS("Monocyte.rds")
Natural_Killer<-readRDS("Natural_Killer.rds")
Progenitor<-readRDS("Progenitor.rds")

#DimPlot
DimPlot(T, reduction = "umap", split.by = "tech")
DimPlot(Epithelial, reduction = "umap", split.by = "tech")
DimPlot(B, reduction = "umap", split.by = "tech")
DimPlot(Secretory, reduction = "umap", split.by = "tech")
DimPlot(Macrophage, reduction = "umap", split.by = "tech")
DimPlot(Monocyte, reduction = "umap", split.by = "tech")
DimPlot(Natural_Killer, reduction = "umap", split.by = "tech")
DimPlot(Progenitor, reduction = "umap", split.by = "tech")

#regroup T cell
T<-readRDS("T.rds")
#pbmc <- JackStraw(T, num.replicate = 100,dims = 40)
#pbmc <- ScoreJackStraw(pbmc, dims = 1:30)
#JackStrawPlot(pbmc, dims = 1:30)
#ElbowPlot(pbmc,ndims = 30)
#find how many 15cluster
hms_neighbor<- FindNeighbors(pbmc, dims = 1:30)
obj <- FindClusters(hms_neighbor, resolution = seq(0.5,1.2,by=0.1))
#resolution設置在0.4-1.2之間,越大clusters越多,查看拐點
clustree(obj)
hms_cluster <- FindClusters( hms_neighbor, resolution = 1.2)
head(Idents(hms_cluster), 5)
hms_cluster<- RunUMAP(hms_cluster, dims = 1:30)
DimPlot(hms_cluster, reduction = "umap")
saveRDS(hms_cluster, file = "T_test_1.2.rds")

scRNA.markers <- FindAllMarkers(hms_cluster, 
                                only.pos = TRUE,  #特异性高表达marker
                                min.pct = 0.05, 
                                logfc.threshold = 0.05)
write.table(scRNA.markers,file="TcellMarkers.txt",sep="\t",row.names=F,quote=F)

#挑选每个细胞亚群中特意高表达的20个基因
top20 <- scRNA.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC) 
write.csv(file="top20_cell_markers.csv",top20)
#整理成表格，只显示基因名字
top20_table=unstack(top20, gene ~ cluster)
names(top20_table)=gsub("X","cluster",names(top20_table))
write.csv(file="Ttop20_marker_genes_1.2.csv",top20_table,row.names=F)

hms_cluster<-readRDS("T_test_1.2.rds")
DimPlot(hms_cluster, reduction = "umap")

#细胞及细胞中基因与RNA数量
slotNames(hms_cluster)
#assay
hms_cluster@assays
dim(hms_cluster@meta.data)
View(hms_cluster@meta.data)

new.cluster.ids <- c("Naive_CD4_T", "Natural_Killer_T","Natural_Killer_T","T_Helper",
                     "Effector_Memory_CD4_T","Resident_Memory_CD8_T","Effector_Memory_CD4_T", 
                     "Terminally_Exhausted_CD8_T","Regulatory_T","Effector_Memory_CD4_T",
                     "Regulatory_T","Terminally_Exhausted_CD8_T","B","Macrophage","Regulatory_T",
                     "Cytotoxic_CD8_T","Effector_Memory_CD8_T", "Macrophage",
                     "Cycling_Natural_Killer", "Terminally_Exhausted_CD8_T",
                     "Pre_Exhausted_CD8_T","Epithelial")
names(new.cluster.ids) <- levels(hms_cluster)
hms_cluster_id<- RenameIdents(hms_cluster, new.cluster.ids)
DimPlot(hms_cluster_id, reduction = "umap", label = TRUE, pt.size = 0.5) 
DimPlot(hms_cluster_id, reduction = "umap", label = FALSE, pt.size = 0.5) 
saveRDS(hms_cluster_id, file = "T_cluster_id_test.rds")

T_cluster_id<-readRDS("T_cluster_id_test.rds")

#only T cells
Only_T<-subset(hms_cluster_id, idents=c("Naive_CD4_T", "Natural_Killer_T","T_Helper",
                                        "Effector_Memory_CD4_T","Resident_Memory_CD8_T",
                                        "Terminally_Exhausted_CD8_T","Regulatory_T",
                                        "Cytotoxic_CD8_T","Effector_Memory_CD8_T",
                                        "Pre_Exhausted_CD8_T"))
DimPlot(Only_T, reduction = "umap", label = TRUE, pt.size = 0.5) 
DimPlot(Only_T, reduction = "umap", label = FALSE, pt.size = 0.5) 
saveRDS(Only_T, file = "Only_T_cluster_id_test.rds")

#Plot
RidgePlot(CD8_T, features = c("CD8A","CD8B","CD38","CD69","ENTPD1","GZMA","GZMH","MYO1F","SYNE1","TSC22D3","XCL2"),
          cols = c("green3","orangered"), group.by="tech", ncol = 4) + theme(axis.title.y = element_blank())

markers.to.plot<-c("CD8A","CD8B","CD38","CD69","ENTPD1","GZMA","GZMH","MYO1F","SYNE1","TSC22D3","XCL2")
DoHeatmap(subset(hms_cluster_id,downsample=50000),features = markers.to.plot,size=5)

genes11_heatmap<-DotPlot(hms_cluster_id,features = c("CD8A","CD8B","CD38","CD69","ENTPD1",
                                                     "GZMA","GZMH","MYO1F","SYNE1","TSC22D3","XCL2"))+RotatedAxis()
genes11_heatmap<-genes11_heatmap$data
write.csv(genes11_heatmap,"genes11_heatmap")

#expression level in each patients
CD8_T<-readRDS("CD8_T.rds")
a<-DoHeatmap(subset(CD8_T,downsample=50000),features = markers.to.plot,size=5,group.by="celltype")
a1<-a$data
write.table(a1,"a1")
               
Only_T <- readRDS(file="Only_T_cluster_id_test.rds") 

DimPlot(Only_T, reduction = "umap", label = TRUE, pt.size = 0.5) 
#DimPlot(Only_T, reduction = "umap", label = FALSE, pt.size = 0.5) 

#doheatmap
markers.to.plot<-c("CD8A","CD8B","CD38","CD69","ENTPD1","GZMA","GZMH","MYO1F","SYNE1","TSC22D3","XCL2")
p1<-DoHeatmap(subset(Only_T,downsample=50000),features = markers.to.plot,size=5)
p1
a<-p1$data
write.table(a,"a")

markers.to.plot<-c("CD8A","CD8B","CD38","CD69","ENTPD1","GZMA","GZMH","MYO1F","SYNE1","TSC22D3","XCL2")
p1<-DoHeatmap(subset(Only_T,downsample=50000),features = markers.to.plot,size=5,group.by="tech")
p1
a<-p1$data
write.table(a,"a")


genes11_heatmap<-DotPlot(Only_T,features = c("CD8A","CD8B","CD38","CD69","ENTPD1",
                                                     "GZMA","GZMH","MYO1F","SYNE1","TSC22D3","XCL2"))+RotatedAxis()

#Plot
RidgePlot(Only_T, features = c("CD8A","CD8B","CD38","CD69","ENTPD1","GZMA","GZMH","MYO1F","SYNE1","TSC22D3","XCL2"),
          cols = c("green3","orangered"), group.by="tech", ncol = 4) + theme(axis.title.y = element_blank())

markers.to.plot<-c("CD8A","CD8B","CD38","CD69","ENTPD1","GZMA","GZMH","MYO1F","SYNE1","TSC22D3","XCL2")
DoHeatmap(subset(Only_T,downsample=50000),features = markers.to.plot,size=5)

genes11_heatmap<-DotPlot(hms_cluster_id,features = c("CD8A","CD8B","CD38","CD69","ENTPD1",
                                                     "GZMA","GZMH","MYO1F","SYNE1","TSC22D3","XCL2"))+RotatedAxis()
genes11_heatmap<-genes11_heatmap$data
write.csv(genes11_heatmap,"genes11_heatmap")

#expression level in each patients
CD8_T<-readRDS("CD8_T.rds")
a<-DoHeatmap(subset(CD8_T,downsample=50000),features = markers.to.plot,size=5,group.by="celltype")
a1<-a$data
write.table(a1,"a1")

#Dotplot
setwd("~/gse162498/data_finally")
hms_cluster<-readRDS("hms_cluster_id_test.rds")
DimPlot(hms_cluster, reduction = "umap", label = TRUE, pt.size = 0.5) 
DimPlot(hms_cluster, reduction = "umap", label = FALSE, pt.size = 0.5) 
features = c("NNMT","ANXA1","IL7R","KLRB1","NKG7","AIF1","S100A4","STMN1",
             "FTL", "APOE","APOC1","JCHAIN", "CD79A", 
             "SCGB1A1","SFTPC", "BATF","CXCL13",
             "FOXP3", "CCR7","IL32","GZMB")
DotPlot(hms_cluster,features=features)+RotatedAxis()

hms_cluster<-readRDS("Only_T_cluster_id_test.rds")
DimPlot(hms_cluster, reduction = "umap", label = TRUE, pt.size = 0.5) 
DimPlot(hms_cluster, reduction = "umap", label = FALSE, pt.size = 0.5) 
features = c("CCR7","SELL","NFKBIZ","SLC2A3","CXCL13","KLRB1","IL7R","ANXA1",
             "GZMK", "CCL4", "GZMB","CD8A","LAG3","BATF",
             "FOXP3","NKG7","PRF1","GNLY","CD8B","CD79A",
            "IFI6","ISG15","IFITM1")
DotPlot(hms_cluster,features=features)+RotatedAxis()

#feature plot
Only_T <- readRDS(file="Only_T_cluster_id_test.rds") 
DimPlot(Only_T, reduction = "umap", label = TRUE, pt.size = 0.5) 
markers.to.plot<-c("CD8A","CD8B","CD38","CD69","ENTPD1","GZMA","GZMH","MYO1F","SYNE1","TSC22D3","XCL2")

#第二次重新分羣,發表的數據
setwd("~/gse162498")
hms_cluster<-readRDS(file="hms_cluster_test_1.2.rds")
DimPlot(hms_cluster, reduction = "umap")
DimPlot(hms_cluster, reduction = "umap",label = TRUE)
new.cluster.ids <- c("naive cd4 T", "Epithelial","treg", "naive cd4 T","naive cd4 T", 
                     "nkt","nkt", "effect cd8","B","chronic activation cd4",
                     "naive cd4 T", "nkt","naive cd4 T", "Effect memory cd4 T","treg",
                     "naive memory cd8 T","nkt", "Epithelial",
                     "Epithelial","Macrophage","nkt", "cycling nk","b","Epithelial",
                     "pre_exhausted cd8 T","dc","microglial", "endothelial",
                     "NK","B","B","Progenitor") 
names(new.cluster.ids) <- levels(hms_cluster)
hms_cluster_id<- RenameIdents(hms_cluster, new.cluster.ids)
DimPlot(hms_cluster_id, reduction = "umap", label = TRUE, pt.size = 0.5) 
DimPlot(hms_cluster_id, reduction = "umap", label = FALSE, pt.size = 0.5) 
saveRDS(hms_cluster_id, file = "cluster01_id_test.rds")

hms_cluster<-readRDS(file="hms_cluster_test_1.2.rds")
DimPlot(hms_cluster, reduction = "umap")
DimPlot(hms_cluster, reduction = "umap",label = TRUE)
new.cluster.ids <- c("T", "Epithelial","T", "T","T", 
                     "T","T", "T","B","T",
                     "T", "T","T", "T","T",
                     "T","T", "Epithelial",
                     "Epithelial","Macrophage","T", "nk","B","Epithelial",
                     "T","dc","microglial", "endothelial",
                     "nk","B","B","Progenitor") 
names(new.cluster.ids) <- levels(hms_cluster)
hms_cluster_id<- RenameIdents(hms_cluster, new.cluster.ids)
DimPlot(hms_cluster_id, reduction = "umap", label = TRUE, pt.size = 0.5) 
DimPlot(hms_cluster_id, reduction = "umap", label = FALSE, pt.size = 0.5) 
saveRDS(hms_cluster_id, file = "cluster02_id_test.rds")
DotPlot(hms_cluster_id,features = c("NNMT","IGFBP7","MGP","VEGFA","SOX4","IL7R","KLF2","S100A4","STMN1","TUBB",
"NKG7","CXCL13", "FTL","APOC1","APOE","JCHAIN","CD79A","SCGB1A1","GZMK","CST7",
"BATF","FOXP3","IL32")) +RotatedAxis()                          
                                               
                                               