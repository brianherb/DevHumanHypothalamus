
setwd('/local/projects-t3/idea/bherb/Hypothalamus/PubRes')

source('./Code/HypoPub_LoadPackages.R')
source('./Code/HypoPub_LoadFunction.R')
source('./Code/HypoPub_LoadReferenceFiles.R')


###################
#### Tri1 data ####
###################

### Tri1 samples - CS12 might be a problem - very low counts 

Tri1Samples = c("CS12","CS13","CS14","CS15","CS20","CS22_1","CS22_2")

CS22_1_data <- Read10X(data.dir = "/local/projects-t3/idea/bherb/Hypothalamus/CS22_Hypothalamus")
CS22_2_data <- Read10X(data.dir = "/local/projects-t3/idea/bherb/Hypothalamus/CS22_2_hypothalamus")
CS20_data <- Read10X(data.dir = "/local/projects-t3/idea/bherb/Hypothalamus/CS20_hypothalamus")
CS15_data <- Read10X(data.dir = "/local/projects-t3/idea/bherb/Hypothalamus/CS15_forebrain")
CS14_data <- Read10X(data.dir = "/local/projects-t3/idea/bherb/Hypothalamus/CS14cortex_GRCh38")
CS13_data <- Read10X(data.dir = "/local/projects-t3/idea/bherb/Hypothalamus/CS13prosencephalon_GRCh38")
CS12_data <- Read10X(data.dir = "/local/projects-t3/idea/bherb/Hypothalamus/CS12Diencephalon_GRCh38")

Tri1RawDat = list(CS12_data,CS13_data,CS14_data,CS15_data,CS20_data,CS22_1_data,CS22_2_data)
names(Tri1RawDat) = Tri1Samples

## cell number tracking 

Tri1CellTracking = data.frame(sample=Tri1Samples,Raw=0,AfterCutoffs=0,AfterDoublets=0)
rownames(Tri1CellTracking) = Tri1Samples

Tri1AllDat = vector(mode='list',length=length(Tri1Samples))
names(Tri1AllDat) = Tri1Samples
Tri1NormDat = vector(mode='list',length=length(Tri1Samples))
names(Tri1NormDat) = Tri1Samples


for(i in Tri1Samples){

Tri1AllDat[[i]] <- CreateSeuratObject(counts = Tri1RawDat[[i]])
Tri1CellTracking[i,'Raw'] = ncol(Tri1AllDat[[i]])
Tri1AllDat[[i]]@meta.data$sample=i
Tri1AllDat[[i]] = RenameCells(Tri1AllDat[[i]],add.cell.id=i)
Tri1AllDat[[i]] <- PercentageFeatureSet(Tri1AllDat[[i]], pattern = "^MT-", col.name = "percent.mt")
## QC plots
pdf(file=paste("./TestPlots/",i,"_all_QC.pdf",sep=''),width=12,height=8)
print(VlnPlot(Tri1AllDat[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
)
print(FeatureScatter(Tri1AllDat[[i]], feature1 = "nCount_RNA", feature2 = "percent.mt"))
print(FeatureScatter(Tri1AllDat[[i]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA"))
dev.off()
Tri1NormDat[[i]] = subset(Tri1AllDat[[i]], subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & nCount_RNA > 1000& nCount_RNA < 15000 & percent.mt < 20)
Tri1CellTracking[i,'AfterCutoffs'] = ncol(Tri1NormDat[[i]])
Tri1NormDat[[i]] <- SCTransform(Tri1NormDat[[i]], vars.to.regress = "percent.mt", verbose = FALSE)
Tri1NormDat[[i]] <- RunPCA(Tri1NormDat[[i]], verbose = FALSE)
Tri1NormDat[[i]] <- RunUMAP(Tri1NormDat[[i]], dims = 1:50, verbose = FALSE)
Tri1NormDat[[i]] <- FindNeighbors(Tri1NormDat[[i]], dims = 1:50, verbose = FALSE)
Tri1NormDat[[i]] <- FindClusters(Tri1NormDat[[i]], verbose = FALSE, resolution=0.5)

DefaultAssay(Tri1NormDat[[i]]) = 'RNA'

pdf(file=paste("./TestPlots/",i,"_Check_Cluster_PreDupRm.pdf",sep=''),width=12,height=8)
print(DimPlot(Tri1NormDat[[i]], reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE))
print(FeaturePlot(Tri1NormDat[[i]], features = c("AQP4","AGT"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(Tri1NormDat[[i]], features = c("FREM2","FZD2"), pt.size = 0.2, ncol = 2) )
print(FeaturePlot(Tri1NormDat[[i]], features = c("ITM2A","IGFBP7"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(Tri1NormDat[[i]], features = c("MATN4","SCRG1"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(Tri1NormDat[[i]], features = c("SMAGP","HPGD"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(Tri1NormDat[[i]], features = c("ERMN","MOBP"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(Tri1NormDat[[i]], features = c("SLC44A1","GPR17"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(Tri1NormDat[[i]], features = c("C1QB","HEXB"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(Tri1NormDat[[i]], features = c("SLC17A6","GAD1"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(Tri1NormDat[[i]], features = c("AQP4","RAX"), pt.size = 0.2, ncol = 2))
dev.off()

## doublet id and removal -  
DefaultAssay(Tri1NormDat[[i]]) = 'SCT'
TMPsce = Seurat::as.SingleCellExperiment(Tri1NormDat[[i]])
TMPsce <- scDblFinder(TMPsce,clust.method='overcluster')
Tri1NormDat[[i]]@meta.data$scDblFinder = colData(TMPsce)$scDblFinder.class
Tri1NormDat[[i]] = subset(Tri1AllDat[[i]],cells = colnames(Tri1NormDat[[i]])[which(Tri1NormDat[[i]]@meta.data$scDblFinder=='singlet')]) #7894
Tri1CellTracking[i,'AfterDoublets'] = ncol(Tri1NormDat[[i]])
Tri1NormDat[[i]] <- SCTransform(Tri1NormDat[[i]], vars.to.regress = "percent.mt", verbose = FALSE)
Tri1NormDat[[i]] <- RunPCA(Tri1NormDat[[i]], verbose = FALSE)
Tri1NormDat[[i]] <- RunUMAP(Tri1NormDat[[i]], dims = 1:50, verbose = FALSE)
Tri1NormDat[[i]] <- FindNeighbors(Tri1NormDat[[i]], dims = 1:50, verbose = FALSE)
Tri1NormDat[[i]] <- FindClusters(Tri1NormDat[[i]], verbose = FALSE, resolution=0.5)

pdf(file=paste("./TestPlots/",i,"_Check_Cluster_PostDupRm.pdf",sep=''),width=12,height=8)
print(DimPlot(Tri1NormDat[[i]], reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE))
print(FeaturePlot(Tri1NormDat[[i]], features = c("AQP4","AGT"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(Tri1NormDat[[i]], features = c("FREM2","FZD2"), pt.size = 0.2, ncol = 2) )
print(FeaturePlot(Tri1NormDat[[i]], features = c("ITM2A","IGFBP7"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(Tri1NormDat[[i]], features = c("MATN4","SCRG1"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(Tri1NormDat[[i]], features = c("SMAGP","HPGD"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(Tri1NormDat[[i]], features = c("ERMN","MOBP"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(Tri1NormDat[[i]], features = c("SLC44A1","GPR17"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(Tri1NormDat[[i]], features = c("C1QB","HEXB"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(Tri1NormDat[[i]], features = c("SLC17A6","GAD1"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(Tri1NormDat[[i]], features = c("AQP4","RAX"), pt.size = 0.2, ncol = 2))
dev.off()

cat(paste('\n\n',i,' done normalization','\n\n',sep=''))

}


save(Tri1AllDat,file='./SeuratObj/Tri1Samples_AllCells.rda',compress=TRUE)

save(Tri1NormDat,file='./SeuratObj/Tri1Samples_PostSCT.rda',compress=TRUE)

write.csv(Tri1CellTracking,file='./SeuratObj/Tri1Samples_CellCounts.csv')



###################
#### Tri2 data ####
###################

load("/local/projects-t3/idea/bherb/Hypothalamus/AllRef_Tri2_integrated.rda") ## Tri2 samples with Reference samples - previously computed, use just as starting point 

Tri2Samples = c("GW18","GW19","GW20","GW22","GW25")

GW25_data <- Read10X(data.dir = "/local/projects-t3/idea/bherb/Hypothalamus/GW25_3V_hypo") #8437
GW22_data <- Read10X(data.dir = "/local/projects-t3/idea/bherb/Hypothalamus/GW22T_hypo") #3085
GW20_data <- Read10X(data.dir = "/local/projects-t3/idea/bherb/Hypothalamus/GW20_hypothalamus") #8795
GW19_data <- Read10X(data.dir = "/local/projects-t3/idea/bherb/Hypothalamus/GW19_hypothalamus") #3283
GW18_data <- Read10X(data.dir = "/local/projects-t3/idea/bherb/Hypothalamus/GW18_hypothalamus")

Tri2RawDat = list(GW18_data,GW19_data,GW20_data,GW22_data,GW25_data)
names(Tri2RawDat) = Tri2Samples

## cell number tracking 

Tri2CellTracking = data.frame(sample=Tri2Samples,Raw=0,AfterCutoffs=0,AfterDoublets=0)
rownames(Tri2CellTracking) = Tri2Samples

Tri2AllDat = vector(mode='list',length=length(Tri2Samples))
names(Tri2AllDat) = Tri2Samples
Tri2NormDat = vector(mode='list',length=length(Tri2Samples))
names(Tri2NormDat) = Tri2Samples


for(i in Tri2Samples){

Tri2AllDat[[i]] <- CreateSeuratObject(counts = Tri2RawDat[[i]])
Tri2CellTracking[i,'Raw'] = ncol(Tri2AllDat[[i]])
Tri2AllDat[[i]]@meta.data$sample=i
Tri2AllDat[[i]] = RenameCells(Tri2AllDat[[i]],add.cell.id=i)
Tri2AllDat[[i]] <- PercentageFeatureSet(Tri2AllDat[[i]], pattern = "^MT-", col.name = "percent.mt")
## QC plots
pdf(file=paste("./TestPlots/",i,"_all_QC.pdf",sep=''),width=12,height=8)
print(VlnPlot(Tri2AllDat[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
)
print(FeatureScatter(Tri2AllDat[[i]], feature1 = "nCount_RNA", feature2 = "percent.mt"))
print(FeatureScatter(Tri2AllDat[[i]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA"))
dev.off()
## just for reference - AllRefTri2.integrated previously made
Tri2AllDat[[i]]@meta.data$Maj_Clust_From_Ref = AllRefTri2.integrated@meta.data[colnames(Tri2AllDat[[i]]),'Maj_Clust_From_Ref']
Tri2NormDat[[i]] = subset(Tri2AllDat[[i]], subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & nCount_RNA > 1000& nCount_RNA < 15000 & percent.mt < 20)
Tri2CellTracking[i,'AfterCutoffs'] = ncol(Tri2NormDat[[i]])
Tri2NormDat[[i]] <- SCTransform(Tri2NormDat[[i]], vars.to.regress = "percent.mt", verbose = FALSE)
Tri2NormDat[[i]] <- RunPCA(Tri2NormDat[[i]], verbose = FALSE)
Tri2NormDat[[i]] <- RunUMAP(Tri2NormDat[[i]], dims = 1:50, verbose = FALSE)
Tri2NormDat[[i]] <- FindNeighbors(Tri2NormDat[[i]], dims = 1:50, verbose = FALSE)
Tri2NormDat[[i]] <- FindClusters(Tri2NormDat[[i]], verbose = FALSE, resolution=0.5)

DefaultAssay(Tri2NormDat[[i]]) = 'RNA'

pdf(file=paste("./TestPlots/",i,"_Check_Cluster_PreDupRm.pdf",sep=''),width=12,height=8)
DimPlot(Tri2NormDat[[i]], label = TRUE) + NoLegend()
DimPlot(Tri2NormDat[[i]], reduction = "umap", group.by = "Maj_Clust_From_Ref", label = TRUE, repel = TRUE)
print(FeaturePlot(Tri2NormDat[[i]], features = c("AQP4","AGT"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(Tri2NormDat[[i]], features = c("FREM2","FZD2"), pt.size = 0.2, ncol = 2) )
print(FeaturePlot(Tri2NormDat[[i]], features = c("ITM2A","IGFBP7"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(Tri2NormDat[[i]], features = c("MATN4","SCRG1"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(Tri2NormDat[[i]], features = c("SMAGP","HPGD"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(Tri2NormDat[[i]], features = c("ERMN","MOBP"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(Tri2NormDat[[i]], features = c("SLC44A1","GPR17"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(Tri2NormDat[[i]], features = c("C1QB","HEXB"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(Tri2NormDat[[i]], features = c("SLC17A6","GAD1"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(Tri2NormDat[[i]], features = c("AQP4","RAX"), pt.size = 0.2, ncol = 2))
dev.off()

## doublet id and removal -  
DefaultAssay(Tri2NormDat[[i]]) = 'SCT'
TMPsce = Seurat::as.SingleCellExperiment(Tri2NormDat[[i]])
TMPsce <- scDblFinder(TMPsce,clust.method='overcluster')
Tri2NormDat[[i]]@meta.data$scDblFinder = colData(TMPsce)$scDblFinder.class
Tri2NormDat[[i]] = subset(Tri2AllDat[[i]],cells = colnames(Tri2NormDat[[i]])[which(Tri2NormDat[[i]]@meta.data$scDblFinder=='singlet')]) #7894
Tri2CellTracking[i,'AfterDoublets'] = ncol(Tri2NormDat[[i]])
Tri2NormDat[[i]]@meta.data$Maj_Clust_From_Ref = AllRefTri2.integrated@meta.data[colnames(Tri2NormDat[[i]]),'Maj_Clust_From_Ref']
Tri2NormDat[[i]] <- SCTransform(Tri2NormDat[[i]], vars.to.regress = "percent.mt", verbose = FALSE)
Tri2NormDat[[i]] <- RunPCA(Tri2NormDat[[i]], verbose = FALSE)
Tri2NormDat[[i]] <- RunUMAP(Tri2NormDat[[i]], dims = 1:50, verbose = FALSE)
Tri2NormDat[[i]] <- FindNeighbors(Tri2NormDat[[i]], dims = 1:50, verbose = FALSE)
Tri2NormDat[[i]] <- FindClusters(Tri2NormDat[[i]], verbose = FALSE, resolution=0.5)

pdf(file=paste("./TestPlots/",i,"_Check_Cluster_PostDupRm.pdf",sep=''),width=12,height=8)
DimPlot(Tri2NormDat[[i]], label = TRUE) + NoLegend()
DimPlot(Tri2NormDat[[i]], reduction = "umap", group.by = "Maj_Clust_From_Ref", label = TRUE, repel = TRUE)
print(FeaturePlot(Tri2NormDat[[i]], features = c("AQP4","AGT"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(Tri2NormDat[[i]], features = c("FREM2","FZD2"), pt.size = 0.2, ncol = 2) )
print(FeaturePlot(Tri2NormDat[[i]], features = c("ITM2A","IGFBP7"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(Tri2NormDat[[i]], features = c("MATN4","SCRG1"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(Tri2NormDat[[i]], features = c("SMAGP","HPGD"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(Tri2NormDat[[i]], features = c("ERMN","MOBP"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(Tri2NormDat[[i]], features = c("SLC44A1","GPR17"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(Tri2NormDat[[i]], features = c("C1QB","HEXB"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(Tri2NormDat[[i]], features = c("SLC17A6","GAD1"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(Tri2NormDat[[i]], features = c("AQP4","RAX"), pt.size = 0.2, ncol = 2))
dev.off()

cat(paste('\n\n',i,' done normalization','\n\n',sep=''))

}

save(Tri2AllDat,file='./SeuratObj/Tri2Samples_AllCells.rda',compress=TRUE)

save(Tri2NormDat,file='./SeuratObj/Tri2Samples_PostSCT.rda',compress=TRUE)

## final cell counts - GW18=5201, GW19=3066, GW20=6118, GW22=2925, GW25=769  

write.csv(Tri2CellTracking,file='./SeuratObj/Tri2Samples_CellCounts.csv')


