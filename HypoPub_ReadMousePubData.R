setwd('/local/projects-t3/idea/bherb/Hypothalamus/PubRes')

source('./Code/HypoPub_LoadPackages.R')
source('./Code/HypoPub_LoadFunction.R')
source('./Code/HypoPub_LoadReferenceFiles.R')

## Mouse 

#### Blackshaw developmental samples 
BlackMeta = read.csv('/local/projects-t3/idea/bherb/Hypothalamus/Blackshaw/data/GSE132355_E10-P45_umap_data.csv',row.names=1)

BlackSamples = c('E10','E11','E12','E13','E14','E15_1','E15_2','E16_1','E16_2','E18','P4','P8','P14','P45_1','P45_2')

BlackRawDat = vector(mode='list',length=length(BlackSamples))
names(BlackRawDat) = BlackSamples

for(i in BlackSamples){
BlackRawDat[[i]] = Read10X(data.dir = paste("/local/projects-t3/idea/bherb/Hypothalamus/Blackshaw/data/",i,sep=''))
cat(i)
}

table(flexsplit(rownames(BlackMeta),'_')[,1])

 E10   E11   E12   E13   E14   E15  E15L   E16 E16V1   E18   P14    P4   P45  P45X    P8 
 3927  7171  9310  9677 16543 10867  2698 16792  1134  8865  6944  4521  8084  8941 13677 

test= BlackRawDat[['E15_1']]
testSum = colSums(test)

summary(testSum[flexsplit(rownames(BlackMeta)[grep('E15_',rownames(BlackMeta))],'_')[,2]])
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  223.0   367.0   547.0   710.7   867.0 16397.0 

test= BlackRawDat[['E15_2']]
testSum = colSums(test)

summary(testSum[flexsplit(rownames(BlackMeta)[grep('E15_',rownames(BlackMeta))],'_')[,2]])
    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
   0.000    0.000    0.000    5.472    0.000 5015.000 

test= BlackRawDat[['E15_1']]
testSum = colSums(test)

summary(testSum[flexsplit(rownames(BlackMeta)[grep('E15L_',rownames(BlackMeta))],'_')[,2]])
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   0.00    0.00    1.00   16.76    1.00 3969.00 

test= BlackRawDat[['E15_2']]
testSum = colSums(test)

summary(testSum[flexsplit(rownames(BlackMeta)[grep('E15L_',rownames(BlackMeta))],'_')[,2]])
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  221.0   280.0   329.0   419.2   426.0 10539.0 


### so E15_1 is E15_ and E15_2 is E15L_  - rename 
test= BlackRawDat[['E16_1']]
testSum = colSums(test)

summary(testSum[flexsplit(rownames(BlackMeta)[grep('E16_',rownames(BlackMeta))],'_')[,2]])
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    226     737    1041    1273    1508   28437

test= BlackRawDat[['E16_2']]
testSum = colSums(test)

summary(testSum[flexsplit(rownames(BlackMeta)[grep('E16_',rownames(BlackMeta))],'_')[,2]])
   #Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
     NA      NA      NA     NaN      NA      NA   16792 

test= BlackRawDat[['E16_1']]
testSum = colSums(test)

summary(testSum[flexsplit(rownames(BlackMeta)[grep('E16V1_',rownames(BlackMeta))],'_')[,2]])
 #   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
     NA      NA      NA     NaN      NA      NA    1134 

test= BlackRawDat[['E16_2']]
testSum = colSums(test)

summary(testSum[flexsplit(rownames(BlackMeta)[grep('E16V1_',rownames(BlackMeta))],'_')[,2]])
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  229.0   404.2   783.0  1272.2  1908.0 10273.0 


### so E16_1 is E16_ and E16_2 is E16V1_  - rename 
test= BlackRawDat[['P45_1']]
testSum = colSums(test)

summary(testSum[flexsplit(rownames(BlackMeta)[grep('P45_',rownames(BlackMeta))],'_')[,2]])
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    237     530     939    1598    1948   57791 

test= BlackRawDat[['P45_2']]
testSum = colSums(test)

summary(testSum[flexsplit(rownames(BlackMeta)[grep('P45_',rownames(BlackMeta))],'_')[,2]])
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   0.00    0.00    1.00   26.65    2.00 6282.00 

test= BlackRawDat[['P45_1']]
testSum = colSums(test)

summary(testSum[flexsplit(rownames(BlackMeta)[grep('P45X_',rownames(BlackMeta))],'_')[,2]])
    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
    0.00     0.00     0.00    31.71     1.00 23992.00 

test= BlackRawDat[['P45_2']]
testSum = colSums(test)

summary(testSum[flexsplit(rownames(BlackMeta)[grep('P45X_',rownames(BlackMeta))],'_')[,2]])
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    224     523     973    1686    2098   30753 

### so P45_1 is P45_ and P45_2 is P45X_  - rename 

BlackSamples = c('E10','E11','E12','E13','E14','E15','E15L','E16','E16V1','E18','P4','P8','P14','P45','P45X')

names(BlackRawDat) = BlackSamples

### loop through metadata file here, and reduce data frames here 

for(i in BlackSamples){
tmpCells = flexsplit(rownames(BlackMeta)[grep(paste(i,'_',sep=''),rownames(BlackMeta))],'_')[,2]
BlackRawDat[[i]] = BlackRawDat[[i]][,tmpCells]
cat(i)
}

### convert to human genes  - for now, just grab first gene - about 10k na's, about 1000 dup

for(i in BlackSamples){

tmpGene = dimnames(BlackRawDat[[i]])[[1]]
tmpGene2 = MMtoHS(tmpGene)
naInd = which(is.na(tmpGene2))
dupMatch = match(tmpGene2,tmpGene2)
dupInd = which(dupMatch!=c(1:length(tmpGene2)))
dimnames(BlackRawDat[[i]])[[1]] = tmpGene2
BlackRawDat[[i]] = BlackRawDat[[i]][-dupInd,]
cat(i)
}

## cell number tracking 

BlackAllDat = vector(mode='list',length=length(BlackSamples))
names(BlackAllDat) = BlackSamples
BlackNormDat = vector(mode='list',length=length(BlackSamples))
names(BlackNormDat) = BlackSamples

for(i in BlackSamples){

BlackAllDat[[i]] <- CreateSeuratObject(counts = BlackRawDat[[i]][-2,]) #-2 for NA gene
BlackAllDat[[i]]@meta.data$sample=i
BlackAllDat[[i]] = RenameCells(BlackAllDat[[i]],add.cell.id=i)
BlackAllDat[[i]] <- PercentageFeatureSet(BlackAllDat[[i]], pattern = "^MT-", col.name = "percent.mt")
BlackAllDat[[i]]@meta.data = cbind(BlackAllDat[[i]]@meta.data, BlackMeta[colnames(BlackAllDat[[i]]),])
## QC plots
pdf(file=paste("./TestPlots/Blackshaw_",i,"_QC.pdf",sep=''),width=12,height=8)
print(VlnPlot(BlackAllDat[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
)
print(FeatureScatter(BlackAllDat[[i]], feature1 = "nCount_RNA", feature2 = "percent.mt"))
print(FeatureScatter(BlackAllDat[[i]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA"))
dev.off()
## just for Blackerence - AllBlackBlack.integrated previously made
#BlackNormDat[[i]] = subset(BlackAllDat[[i]], subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & nCount_RNA > 1000& nCount_RNA < 15000 & percent.mt < 20)
BlackNormDat[[i]] = BlackAllDat[[i]]
if(sum(BlackAllDat[[i]]@meta.data$percent.mt)==0){
BlackNormDat[[i]] <- SCTransform(BlackNormDat[[i]], verbose = FALSE)
  } else {
BlackNormDat[[i]] <- SCTransform(BlackNormDat[[i]], vars.to.regress = "percent.mt", verbose = FALSE)
}
BlackNormDat[[i]] <- RunPCA(BlackNormDat[[i]], verbose = FALSE)
BlackNormDat[[i]] <- RunUMAP(BlackNormDat[[i]], dims = 1:50, verbose = FALSE)
BlackNormDat[[i]] <- FindNeighbors(BlackNormDat[[i]], dims = 1:50, verbose = FALSE)
BlackNormDat[[i]] <- FindClusters(BlackNormDat[[i]], verbose = FALSE, resolution=0.5)

DefaultAssay(BlackNormDat[[i]]) = 'RNA'

pdf(file=paste("./TestPlots/Blackshaw_",i,"_Check_Cluster.pdf",sep=''),width=12,height=8)
print(DimPlot(BlackNormDat[[i]], label = TRUE) + NoLegend())
print(DimPlot(BlackNormDat[[i]], reduction = "umap", group.by = "Cluster", label = TRUE, repel = TRUE))
print(DimPlot(BlackNormDat[[i]], reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE))
dev.off()
cat(i)
}

save(BlackAllDat,file='./SeuratObj/BlackSamples_AllCells.rda',compress=TRUE)
save(BlackNormDat,file='./SeuratObj/BlackSamples_PostSCT.rda',compress=TRUE)












#### Reference samples - do same preproc 


RefSamples = c('Mof','Chen','Mick.male','Mick.female','Kim_FC1','Kim_MC1','Kim_MC2','Kim_MC3','Wen','Cam')

Mof_data <- Read10X(data.dir = "/local/projects-t3/idea/bherb/Hypothalamus/Moffitt") #27998 genes x  31299 cells 
Mof_meta = read.csv('/local/projects-t3/idea/bherb/Hypothalamus/Moffitt/Moffitt_metadata.csv')
# cell names as rows 
Mof_meta$Major_Cell_Type = Mof_meta$Cell_class
Mof_meta$Sub_Cell_Type = as.character(Mof_meta$Non_neuronal_cluster)
Mof_meta$Sub_Cell_Type[which(Mof_meta$Non_neuronal_cluster=='')] = as.character(Mof_meta$Neuronal_cluster)[which(Mof_meta$Non_neuronal_cluster=='')]
Mof_meta$Sub_Cell_Type[which(Mof_meta$Sub_Cell_Type=='')] = NA
rownames(Mof_meta) = paste('Mof_',Mof_meta$Cell,sep='')

ChenCount = read.table('/local/projects-t3/idea/bherb/Hypothalamus/Chen/GSE87544_Merged_17samples_14437cells_count.txt',header=TRUE) #  23284 genes x 14437 cells (first col is gene name )
rownames(ChenCount) = toupper(ChenCount$Gene)
Chen_data = ChenCount[,-1]
Chen_meta = read.csv('/local/projects-t3/idea/bherb/Hypothalamus/Chen/GSE87544_1443737Cells.SVM.cluster.identity.renamed.csv',row.names=1)
Chen_meta$Major_Cell_Type = Chen_meta$SVM_clusterID
Chen_meta$Sub_Cell_Type = Chen_meta$SVM_clusterID
rownames(Chen_meta) = paste('Chen_',rownames(Chen_meta),sep='')


Mick.male_data <- Read10X(data.dir = "/local/projects-t3/idea/bherb/Hypothalamus/Mickelsen/GSM3562050") #28692 genes x 3439 cells
Mick.male_meta = read.csv("/local/projects-t3/idea/bherb/Hypothalamus/Mickelsen/Mickelsen_Male_Meta.csv",row.names=1)
rownames(Mick.male_meta) = paste('Mick.male_',rownames(Mick.male_meta),sep='')

Mick.female_data <- Read10X(data.dir = "/local/projects-t3/idea/bherb/Hypothalamus/Mickelsen/GSM3562051") #28692 genes x 3793 cells
Mick.female_meta=read.csv("/local/projects-t3/idea/bherb/Hypothalamus/Mickelsen/Mickelsen_Female_Meta.csv",row.names=1)
rownames(Mick.female_meta) = paste('Mick.female_',rownames(Mick.female_meta),sep='')


Kim.FC1_data <- Read10X(data.dir = "/local/projects-t3/idea/bherb/Hypothalamus/Kim/2018_0428") # 28000 genes x 10030 cells
Kim.FC1_meta = read.csv("/local/projects-t3/idea/bherb/Hypothalamus/Kim/Kim_FC1_Meta.csv",row.names=1)
rownames(Kim.FC1_meta) = paste('Kim_FC1_',rownames(Kim.FC1_meta),sep='')

Kim.MC1_data <- Read10X(data.dir = "/local/projects-t3/idea/bherb/Hypothalamus/Kim/2018_0420") # 28000 genes x 5270 cells
Kim.MC1_meta = read.csv("/local/projects-t3/idea/bherb/Hypothalamus/Kim/Kim_MC1_Meta.csv",row.names=1)
rownames(Kim.MC1_meta) = paste('Kim_MC1_',rownames(Kim.MC1_meta),sep='')

Kim.MC2_data <- Read10X(data.dir = "/local/projects-t3/idea/bherb/Hypothalamus/Kim/2018_0429") # 28000 genes x 7926 cells
Kim.MC2_meta = read.csv("/local/projects-t3/idea/bherb/Hypothalamus/Kim/Kim_MC2_Meta.csv",row.names=1)
rownames(Kim.MC2_meta) = paste('Kim_MC2_',rownames(Kim.MC2_meta),sep='')

Kim.MC3_data <- Read10X(data.dir = "/local/projects-t3/idea/bherb/Hypothalamus/Kim/2018_0520") # 28000 genes x 7344 cells
Kim.MC3_meta = read.csv("/local/projects-t3/idea/bherb/Hypothalamus/Kim/Kim_MC3_Meta.csv",row.names=1)
rownames(Kim.MC3_meta) = paste('Kim_MC3_',rownames(Kim.MC3_meta),sep='')

Wen_data <- Read10X(data.dir = "/local/projects-t3/idea/bherb/Hypothalamus/Wen/GSM3880020") #31053 genes x  11437 cells
Wen_meta = read.csv("/local/projects-t3/idea/bherb/Hypothalamus/Wen/Wen_Meta.csv",row.names=1)
rownames(Wen_meta) = paste('Wen_',rownames(Wen_meta),sep='')

Cam_Norm_data = read.table('/local/projects-t3/idea/bherb/Hypothalamus/Campbell/GSE93374_Merged_all_020816_BatchCorrected_LNtransformed_doubletsremoved_Data.txt')

Cam_data = read.table('/local/projects-t3/idea/bherb/Hypothalamus/Campbell/GSE93374_Merged_all_020816_DGE.txt')
camgenes = toupper(rownames(Cam_data))
dupind= which(camgenes=='PMS2') ## this is the only gene duplicated
Cam_data = Cam_data[-dupind[1],]
rownames(Cam_data) = camgenes[-dupind[1]]

### limit to cells present in Cam_Norm
Cam_data = Cam_data[,colnames(Cam_Norm_data)] # 26773 genes x  20921 cells


Cam_meta = read.table('/local/projects-t3/idea/bherb/Hypothalamus/Campbell/GSE93374_cell_metadata.txt',sep='\t',header=TRUE)
Cam.meta.neu = flexsplit(Cam_meta$Neuron.Subclusters,'.') #30
Cam.meta.sub = flexsplit(Cam_meta$All.Cell.Subclusters,'.') #36
Cam.meta.all = flexsplit(Cam_meta$All.Cell.Clusters,'.') #20
## add col to meta, full names 
Cam_meta = Cam_meta[,c(1:11)]
colnames(Cam_meta) = flexsplit(colnames(Cam_meta),'.')[,2]
rownames(Cam_meta) = Cam_meta$ID
Cam_meta$clust_all_name = Cam.meta.all[match(Cam_meta$clust_all,Cam.meta.all[1:36,1]),2]
Cam_meta$clust_sub_name = Cam.meta.sub[match(Cam_meta$clust_all_micro,Cam.meta.sub[1:20,1]),2]
Cam_meta$clust_neuro_name = Cam.meta.neu[match(Cam_meta$clust_neurons,Cam.meta.neu[1:30,1]),2]
Cam_meta$Major_Cell_Type = Cam_meta$clust_all_name
Cam_meta$Sub_Cell_Type = as.character(Cam_meta$clust_all_name)
Cam_meta$Sub_Cell_Type[which(!is.na(Cam_meta$clust_neuro_name))]=as.character(Cam_meta$clust_neuro_name[which(!is.na(Cam_meta$clust_neuro_name))])
rownames(Cam_meta) = paste('Cam_',rownames(Cam_meta),sep='')


### watch for mt, meta 

RefSamples = c('Mof','Chen','Mick.male','Mick.female','Kim_FC1','Kim_MC1','Kim_MC2','Kim_MC3','Wen','Cam')
RefRegions = c('POA','ALL','LHA','LHA','VMHvl','VMHvl','VMHvl','VMHvl','SCN','ARC')

RefRawDat = list(Mof_data,Chen_data,Mick.male_data,Mick.female_data,Kim.FC1_data,Kim.MC1_data,Kim.MC2_data,Kim.MC3_data,Wen_data,Cam_data)
names(RefRawDat) = RefSamples

RefMeta = list(Mof_meta,Chen_meta,Mick.male_meta,Mick.female_meta,Kim.FC1_meta,Kim.MC1_meta,Kim.MC2_meta,Kim.MC3_meta,Wen_meta,Cam_meta)
names(RefMeta) = RefSamples

## cell number tracking 

RefCellTracking = data.frame(sample=RefSamples,Raw=0,AfterCutoffs=0,AfterDoublets=0)
rownames(RefCellTracking) = RefSamples

RefAllDat = vector(mode='list',length=length(RefSamples))
names(RefAllDat) = RefSamples
RefNormDat = vector(mode='list',length=length(RefSamples))
names(RefNormDat) = RefSamples


for(i in RefSamples){

RefAllDat[[i]] <- CreateSeuratObject(counts = RefRawDat[[i]])
RefCellTracking[i,'Raw'] = ncol(RefAllDat[[i]])
RefAllDat[[i]]@meta.data$sample=i
RefAllDat[[i]]@meta.data$region=RefRegions[which(RefSamples==i)]
RefAllDat[[i]] = RenameCells(RefAllDat[[i]],add.cell.id=i)
RefAllDat[[i]] <- PercentageFeatureSet(RefAllDat[[i]], pattern = "^MT-", col.name = "percent.mt")
RefAllDat[[i]]@meta.data = cbind(RefAllDat[[i]]@meta.data, RefMeta[[i]][colnames(RefAllDat[[i]]),])
## QC plots
pdf(file=paste("./TestPlots/",i,"_all_QC.pdf",sep=''),width=12,height=8)
print(VlnPlot(RefAllDat[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
)
print(FeatureScatter(RefAllDat[[i]], feature1 = "nCount_RNA", feature2 = "percent.mt"))
print(FeatureScatter(RefAllDat[[i]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA"))
dev.off()
## just for reference - AllRefRef.integrated previously made
RefNormDat[[i]] = subset(RefAllDat[[i]], subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & nCount_RNA > 1000& nCount_RNA < 15000 & percent.mt < 20)
RefCellTracking[i,'AfterCutoffs'] = ncol(RefNormDat[[i]])
if(sum(RefAllDat[[i]]@meta.data$percent.mt)==0){
RefNormDat[[i]] <- SCTransform(RefNormDat[[i]], verbose = FALSE)
  } else {
RefNormDat[[i]] <- SCTransform(RefNormDat[[i]], vars.to.regress = "percent.mt", verbose = FALSE)
}
RefNormDat[[i]] <- RunPCA(RefNormDat[[i]], verbose = FALSE)
RefNormDat[[i]] <- RunUMAP(RefNormDat[[i]], dims = 1:50, verbose = FALSE)
RefNormDat[[i]] <- FindNeighbors(RefNormDat[[i]], dims = 1:50, verbose = FALSE)
RefNormDat[[i]] <- FindClusters(RefNormDat[[i]], verbose = FALSE, resolution=0.5)

DefaultAssay(RefNormDat[[i]]) = 'RNA'

pdf(file=paste("./TestPlots/",i,"_Check_Cluster_PreDupRm.pdf",sep=''),width=12,height=8)
print(DimPlot(RefNormDat[[i]], label = TRUE) + NoLegend())
print(DimPlot(RefNormDat[[i]], reduction = "umap", group.by = "Major_Cell_Type", label = TRUE, repel = TRUE))
print(DimPlot(RefNormDat[[i]], reduction = "umap", group.by = "Sub_Cell_Type", label = TRUE, repel = TRUE))
print(FeaturePlot(RefNormDat[[i]], features = c("AQP4","AGT"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(RefNormDat[[i]], features = c("FREM2","FZD2"), pt.size = 0.2, ncol = 2) )
print(FeaturePlot(RefNormDat[[i]], features = c("ITM2A","IGFBP7"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(RefNormDat[[i]], features = c("MATN4","SCRG1"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(RefNormDat[[i]], features = c("SMAGP","HPGD"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(RefNormDat[[i]], features = c("ERMN","MOBP"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(RefNormDat[[i]], features = c("SLC44A1","GPR17"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(RefNormDat[[i]], features = c("C1QB","HEXB"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(RefNormDat[[i]], features = c("SLC17A6","GAD1"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(RefNormDat[[i]], features = c("AQP4","RAX"), pt.size = 0.2, ncol = 2))
dev.off()

## doublet id and removal -  
DefaultAssay(RefNormDat[[i]]) = 'SCT'
TMPsce = Seurat::as.SingleCellExperiment(RefNormDat[[i]])
TMPsce <- scDblFinder(TMPsce,clust.method='overcluster')
RefNormDat[[i]]@meta.data$scDblFinder = colData(TMPsce)$scDblFinder.class
RefNormDat[[i]] = subset(RefAllDat[[i]],cells = colnames(RefNormDat[[i]])[which(RefNormDat[[i]]@meta.data$scDblFinder=='singlet')]) #7894
RefCellTracking[i,'AfterDoublets'] = ncol(RefNormDat[[i]])
RefNormDat[[i]] <- SCTransform(RefNormDat[[i]], vars.to.regress = "percent.mt", verbose = FALSE)
RefNormDat[[i]] <- RunPCA(RefNormDat[[i]], verbose = FALSE)
RefNormDat[[i]] <- RunUMAP(RefNormDat[[i]], dims = 1:50, verbose = FALSE)
RefNormDat[[i]] <- FindNeighbors(RefNormDat[[i]], dims = 1:50, verbose = FALSE)
RefNormDat[[i]] <- FindClusters(RefNormDat[[i]], verbose = FALSE, resolution=0.5)

pdf(file=paste("./TestPlots/",i,"_Check_Cluster_PostDupRm.pdf",sep=''),width=12,height=8)
print(DimPlot(RefNormDat[[i]], label = TRUE) + NoLegend())
print(DimPlot(RefNormDat[[i]], reduction = "umap", group.by = "Major_Cell_Type", label = TRUE, repel = TRUE))
print(DimPlot(RefNormDat[[i]], reduction = "umap", group.by = "Sub_Cell_Type", label = TRUE, repel = TRUE))
print(FeaturePlot(RefNormDat[[i]], features = c("AQP4","AGT"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(RefNormDat[[i]], features = c("FREM2","FZD2"), pt.size = 0.2, ncol = 2) )
print(FeaturePlot(RefNormDat[[i]], features = c("ITM2A","IGFBP7"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(RefNormDat[[i]], features = c("MATN4","SCRG1"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(RefNormDat[[i]], features = c("SMAGP","HPGD"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(RefNormDat[[i]], features = c("ERMN","MOBP"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(RefNormDat[[i]], features = c("SLC44A1","GPR17"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(RefNormDat[[i]], features = c("C1QB","HEXB"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(RefNormDat[[i]], features = c("SLC17A6","GAD1"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(RefNormDat[[i]], features = c("AQP4","RAX"), pt.size = 0.2, ncol = 2))
dev.off()

cat(paste('\n\n',i,' done normalization','\n\n',sep=''))

}

save(RefAllDat,file='./SeuratObj/RefSamples_AllCells.rda',compress=TRUE)

save(RefNormDat,file='./SeuratObj/RefSamples_PostSCT.rda',compress=TRUE)

write.csv(RefCellTracking,file='./SeuratObj/RefSamples_CellCounts.csv')


## merging samples  - Do  Ref, Tri2, Tri1, Tri2+Tri1, Ref+Tri2, Ref+All, Ref+Tri1?, Tri2Neurons, Tri2Neurons+CS22 neurons, RefNeuron + Tri2 Neuron, Ref Oligo + Tri2 Oligo 

### create reference dataset 

RefCombDat = RefNormDat

RefCombDat[['Mof']] = subset(RefCombDat[['Mof']],cells = colnames(RefCombDat[['Mof']])[sample(c(1:ncol(RefCombDat[['Mof']])),2000)])
RefCombDat[['Chen']] = subset(RefCombDat[['Chen']],cells = colnames(RefCombDat[['Chen']])[sample(c(1:ncol(RefCombDat[['Chen']])),2000)])
RefCombDat[['Mick.male']] = subset(RefCombDat[['Mick.male']],cells = colnames(RefCombDat[['Mick.male']])[sample(c(1:ncol(RefCombDat[['Mick.male']])),1000)])
RefCombDat[['Mick.female']] = subset(RefCombDat[['Mick.female']],cells = colnames(RefCombDat[['Mick.female']])[sample(c(1:ncol(RefCombDat[['Mick.female']])),1000)])
RefCombDat[['Kim_FC1']] = subset(RefCombDat[['Kim_FC1']],cells = colnames(RefCombDat[['Kim_FC1']])[sample(c(1:ncol(RefCombDat[['Kim_FC1']])),500)])
RefCombDat[['Kim_MC1']] = subset(RefCombDat[['Kim_MC1']],cells = colnames(RefCombDat[['Kim_MC1']])[sample(c(1:ncol(RefCombDat[['Kim_MC1']])),500)])
RefCombDat[['Kim_MC2']] = subset(RefCombDat[['Kim_MC2']],cells = colnames(RefCombDat[['Kim_MC2']])[sample(c(1:ncol(RefCombDat[['Kim_MC2']])),500)])
RefCombDat[['Kim_MC3']] = subset(RefCombDat[['Kim_MC3']],cells = colnames(RefCombDat[['Kim_MC3']])[sample(c(1:ncol(RefCombDat[['Kim_MC3']])),500)])
RefCombDat[['Wen']] = subset(RefCombDat[['Wen']],cells = colnames(RefCombDat[['Wen']])[sample(c(1:ncol(RefCombDat[['Wen']])),2000)])
RefCombDat[['Cam']] = subset(RefCombDat[['Cam']],cells = colnames(RefCombDat[['Cam']])[sample(c(1:ncol(RefCombDat[['Cam']])),2000)])


for (i in 1:length(RefCombDat)) {
    RefCombDat[[i]] <- NormalizeData(RefCombDat[[i]], verbose = FALSE)
    RefCombDat[[i]] <- FindVariableFeatures(RefCombDat[[i]], selection.method = "vst", nfeatures = 2000, verbose = FALSE)
}

options(future.globals.maxSize = 4000 * 1024^2)

AllRef.anchors <- FindIntegrationAnchors(object.list = RefCombDat, dims = 1:30)
AllRef.integrated <- IntegrateData(anchorset = AllRef.anchors, dims = 1:30)
DefaultAssay(AllRef.integrated) <- "integrated"
AllRef.integrated <- ScaleData(AllRef.integrated, verbose = FALSE)
AllRef.integrated <- RunPCA(AllRef.integrated, npcs = 30, verbose = FALSE)
AllRef.integrated <- RunUMAP(AllRef.integrated, reduction = "pca", dims = 1:30)

AllRef.integrated@meta.data$BinaryNeu = NA

AllRef.integrated@meta.data$BinaryNeu[c(grep('excit',AllRef.integrated@meta.data$Major_Cell_Type,ignore.case =TRUE),grep('glu',AllRef.integrated@meta.data$Major_Cell_Type,ignore.case =TRUE))] = 'Excitatory'
AllRef.integrated@meta.data$BinaryNeu[c(grep('inhib',AllRef.integrated@meta.data$Major_Cell_Type,ignore.case =TRUE),grep('gaba',AllRef.integrated@meta.data$Major_Cell_Type,ignore.case =TRUE))] = 'Inhibitory'

AllRef.integrated@meta.data$Major_Cell_Type[which(is.na(AllRef.integrated@meta.data$Major_Cell_Type))] = "Unk"

AllRef.integrated@meta.data$Sub_Cell_Type[which(is.na(AllRef.integrated@meta.data$Sub_Cell_Type))] = "Unk"


pdf('./TestPlots/AllRef_Major_Sub_clusters.pdf',width=12,height=8)
print(DimPlot(AllRef.integrated, label = TRUE) + NoLegend())
for(i in c("sample","region","BinaryNeu","SVM_clusterID","clust_all_name")){
print(DimPlot(AllRef.integrated, reduction = "umap", group.by = i, label = TRUE, repel = TRUE))
}
print(FeaturePlot(AllRef.integrated, features = c("SLC17A6","GAD1"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(AllRef.integrated, features = c("RAX","CHGB"), pt.size = 0.2, ncol = 2))
dev.off()

### clustering to identify cell types 



AllRef.integrated <- FindNeighbors(AllRef.integrated, dims = 1:30, verbose = FALSE)

AllRef.integrated = BuildClusterTree(object=AllRef.integrated)

for(i in seq(from=0.5, to=1, by=0.1)){
AllRef.integrated <- FindClusters(AllRef.integrated, verbose = FALSE, resolution=i)
}

pdf('./TestPlots/AllRef_Check_clusters.pdf',width=12,height=8)
print(DimPlot(AllRef.integrated, label = TRUE) + NoLegend())
for(i in c("sample","region","BinaryNeu","SVM_clusterID","clust_all_name")){
print(DimPlot(AllRef.integrated, reduction = "umap", group.by = i, label = TRUE, repel = TRUE))
}
for (k in seq(from=0.5,to=1,by=0.1)){
print(DimPlot(AllRef.integrated, reduction = "umap", group.by = paste("integrated_snn_res.",k,sep=""), label = TRUE, repel = TRUE))
}
dev.off()


pdf('./TestPlots/AllRef_check_ThalamusGenes.pdf',width=12,height=8)
print(DimPlot(AllRef.integrated, label = TRUE) + NoLegend())
for(i in c("sample","region","BinaryNeu","SVM_clusterID","clust_all_name")){
print(DimPlot(AllRef.integrated, reduction = "umap", group.by = i, label = TRUE, repel = TRUE))
}
for(k in ThalGenes){
print(FeaturePlot(AllRef.integrated, features = k, pt.size = 0.4, order=TRUE))
}
dev.off()



### choose res = 0.6 for cell type designations 


Microglia = 7
Tanycyte_A1B1 = 4
Tanycyte_B2 = 18
Mat_Oligo_1 = 13
Mat_Oligo_2 = 2
Mat_Oligo_3 = 11
Maturing_Oligo = 15
Immature_Oligo = 3
Ependymal = 17
Astrocyte = 5
Neuron = 0,1,6,8,9,12
Fibroblast = 19
Epithelial = 10, 14, 16

Remove? Single sample small clusters - 20, 21, 22 

AllRef.integrated@meta.data$Maj_Clust_Int=NA

AllRef.integrated@meta.data$Maj_Clust_Int[which(!is.na(match(AllRef.integrated@meta.data$integrated_snn_res.0.6,c("0","1","6","8","9","12"))))] = "Neuron"

AllRef.integrated@meta.data$Maj_Clust_Int[which(AllRef.integrated@meta.data$integrated_snn_res.0.6==4)]="Tanycyte_A12B1"

AllRef.integrated@meta.data$Maj_Clust_Int[which(AllRef.integrated@meta.data$integrated_snn_res.0.6==18)]="Tanycyte_B2"

AllRef.integrated@meta.data$Maj_Clust_Int[which(AllRef.integrated@meta.data$integrated_snn_res.0.6==5)] = "Astrocyte"

AllRef.integrated@meta.data$Maj_Clust_Int[which(!is.na(match(AllRef.integrated@meta.data$integrated_snn_res.0.6,c("10","14","16"))))] = "Epithelial"

AllRef.integrated@meta.data$Maj_Clust_Int[which(AllRef.integrated@meta.data$integrated_snn_res.0.6==13)]="Mat_Oligo_1"

AllRef.integrated@meta.data$Maj_Clust_Int[which(AllRef.integrated@meta.data$integrated_snn_res.0.6==2)]="Mat_Oligo_2"

AllRef.integrated@meta.data$Maj_Clust_Int[which(AllRef.integrated@meta.data$integrated_snn_res.0.6==11)]="Mat_Oligo_3"

AllRef.integrated@meta.data$Maj_Clust_Int[which(AllRef.integrated@meta.data$integrated_snn_res.0.6==15)]="Maturing_Oligo"

AllRef.integrated@meta.data$Maj_Clust_Int[which(AllRef.integrated@meta.data$integrated_snn_res.0.6==3)]="Immature_Oligo"

AllRef.integrated@meta.data$Maj_Clust_Int[which(AllRef.integrated@meta.data$integrated_snn_res.0.6==17)]="Ependymal"

AllRef.integrated@meta.data$Maj_Clust_Int[which(AllRef.integrated@meta.data$integrated_snn_res.0.6==19)]="Fibroblast"

AllRef.integrated@meta.data$Maj_Clust_Int[which(AllRef.integrated@meta.data$integrated_snn_res.0.6==7)]="Microglia"

AllRef.integrated@meta.data$Maj_Clust_Int[which(!is.na(match(AllRef.integrated@meta.data$integrated_snn_res.0.6,c("20","21","22"))))] = "Unknown"


pdf('./TestPlots/AllRef_Check_Maj_Cell_clustering.pdf',width=12,height=8)
print(DimPlot(AllRef.integrated, reduction = "umap", group.by = "integrated_snn_res.0.6", label = TRUE, repel = TRUE))
print(DimPlot(AllRef.integrated, reduction = "umap", group.by = "sample", label = TRUE, repel = TRUE))
print(DimPlot(AllRef.integrated, reduction = "umap", group.by = "Maj_Clust_Int", label = TRUE, repel = TRUE))
print(DimPlot(AllRef.integrated, reduction = "umap", group.by = "SVM_clusterID", label = TRUE, repel = TRUE))
print(DimPlot(AllRef.integrated, reduction = "umap", group.by = "clust_all_name", label = TRUE, repel = TRUE))
dev.off()

## get rid of Unknown cells 
AllRef.integrated = subset(AllRef.integrated,cells=colnames(AllRef.integrated)[which(AllRef.integrated@meta.data$Maj_Clust_Int!='Unknown')]) #11899

for(i in RefSamples){
tmpCell = intersect(colnames(AllRef.integrated),colnames(RefCombDat[[i]]))
RefCombDat[[i]] = subset(RefCombDat[[i]],cells = tmpCell)
RefCombDat[[i]]@meta.data$Maj_Clust_Int = AllRef.integrated@meta.data[tmpCell,'Maj_Clust_Int']
}

save(RefCombDat,file="./SeuratObj/Mus_Ref_12000cells_Indv_Objects.rda",compress=TRUE)

save(AllRef.integrated,file="./SeuratObj/Mus_Ref_12000cells_seurat_obj.rda",compress=TRUE)

### Merge Ref with Tri2 

load('./SeuratObj/Tri2Samples_PostSCT.rda')
load("./SeuratObj/Mus_Ref_12000cells_Indv_Objects.rda")

RefTri2NormDat = c(RefCombDat,Tri2NormDat)

for (i in 1:length(RefTri2NormDat)) {
    RefTri2NormDat[[i]] <- NormalizeData(RefTri2NormDat[[i]], verbose = FALSE)
    RefTri2NormDat[[i]] <- FindVariableFeatures(RefTri2NormDat[[i]], selection.method = "vst", nfeatures = 2000, verbose = FALSE)
}

RefTri2.anchors <- FindIntegrationAnchors(object.list = RefTri2NormDat, dims = 1:30)
RefTri2.integrated <- IntegrateData(anchorset = RefTri2.anchors, dims = 1:30)
DefaultAssay(RefTri2.integrated) <- "integrated"
RefTri2.integrated <- ScaleData(RefTri2.integrated, verbose = FALSE)
RefTri2.integrated <- RunPCA(RefTri2.integrated, npcs = 30, verbose = FALSE)
RefTri2.integrated <- RunUMAP(RefTri2.integrated, reduction = "pca", dims = 1:30)



pdf('./TestPlots/RefTri2_Major_Sub_clusters.pdf',width=12,height=8)
print(DimPlot(RefTri2.integrated, label = TRUE) + NoLegend())
for(i in c("sample","region","Maj_Clust_Int","SVM_clusterID","clust_all_name")){
print(DimPlot(RefTri2.integrated, reduction = "umap", group.by = i, label = TRUE, repel = TRUE))
}
print(FeaturePlot(RefTri2.integrated, features = c("SLC17A6","GAD1"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(RefTri2.integrated, features = c("RAX","CHGB"), pt.size = 0.2, ncol = 2))
dev.off()


RefTri2.integrated <- RunPCA(RefTri2.integrated, npcs = 50, verbose = FALSE)

pdf('./TestPlots/RefTri2_Test_PCs.pdf',width=12,height=8)
for (i in seq(from=10,to=50,by=10)) {
RefTri2.integrated <- RunUMAP(RefTri2.integrated, reduction = "pca", dims = 1:i)
print(DimPlot(RefTri2.integrated, reduction = "umap", group.by = 'sample', label = TRUE, repel = TRUE))
print(DimPlot(RefTri2.integrated, reduction = "umap", group.by = 'Maj_Clust_Int', label = TRUE, repel = TRUE))
}
dev.off()


RefTri2.integrated <- FindNeighbors(RefTri2.integrated, dims = 1:50, verbose = FALSE)

RefTri2.integrated = BuildClusterTree(object=RefTri2.integrated)

for(i in seq(from=0.1, to=2.0, by=0.1)){
RefTri2.integrated <- FindClusters(RefTri2.integrated, verbose = FALSE, resolution=i)
}

RefTri2.integrated  <- CellCycleScoring(RefTri2.integrated, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

pdf('./TestPlots/RefTri2_Check_clusters.pdf',width=12,height=8)
print(DimPlot(RefTri2.integrated, label = TRUE) + NoLegend())
for(i in c("sample","region","Phase","Maj_Clust_From_Ref","SVM_clusterID","clust_all_name")){
print(DimPlot(RefTri2.integrated, reduction = "umap", group.by = i, label = TRUE, repel = TRUE))
}
for (k in seq(from=0.1,to=2,by=0.1)){
print(DimPlot(RefTri2.integrated, reduction = "umap", group.by = paste("integrated_snn_res.",k,sep=""), label = TRUE, repel = TRUE))
}
dev.off()

## res = 0.3 captures major cell types 

subplots(RefTri2.integrated,'integrated_snn_res.0.3',file='./TestPlots/RefTri2_Res03_Indiv_Clusters.pdf')

neurons - 0,12,13,14,17
neuro progenitor - 6
oligo progenitor / dividing  - 9
immature oligo - 3 
maturing oligo - 10
mature oligo - 5 
glia progenitor - 4 
tanycyte - 2 
Astrocyte - 1 
Ependymal - 11
Epithelial - 8 
Fibroblast - 15
Microglia - 7
UNK - 16


RefTri2.integrated@meta.data$Cell_Type_From_Ref

RefTri2.integrated@meta.data[which(!is.na(match(RefTri2.integrated@meta.data$integrated_snn_res.0.3,c("0","12","13","14","17")))),'Cell_Type_From_Ref'] = "Neuron"
RefTri2.integrated@meta.data[which(RefTri2.integrated@meta.data$integrated_snn_res.0.3==8),'Cell_Type_From_Ref']="Epithelial"
RefTri2.integrated@meta.data[which(RefTri2.integrated@meta.data$integrated_snn_res.0.3==15),'Cell_Type_From_Ref']="Fibroblast"
RefTri2.integrated@meta.data[which(RefTri2.integrated@meta.data$integrated_snn_res.0.3==9),'Cell_Type_From_Ref']="DividingProgenitor"
RefTri2.integrated@meta.data[which(RefTri2.integrated@meta.data$integrated_snn_res.0.3==3),'Cell_Type_From_Ref']="ImmatureOligo"
RefTri2.integrated@meta.data[which(RefTri2.integrated@meta.data$integrated_snn_res.0.3==5),'Cell_Type_From_Ref']="MatureOligo"
RefTri2.integrated@meta.data[which(RefTri2.integrated@meta.data$integrated_snn_res.0.3==10),'Cell_Type_From_Ref']="MaturingOligo"
RefTri2.integrated@meta.data[which(RefTri2.integrated@meta.data$integrated_snn_res.0.3==1),'Cell_Type_From_Ref']="Astrocyte"
RefTri2.integrated@meta.data[which(RefTri2.integrated@meta.data$integrated_snn_res.0.3==2),'Cell_Type_From_Ref']="Tanycyte"
RefTri2.integrated@meta.data[which(RefTri2.integrated@meta.data$integrated_snn_res.0.3==11),'Cell_Type_From_Ref']="Ependymal"
RefTri2.integrated@meta.data[which(RefTri2.integrated@meta.data$integrated_snn_res.0.3==7),'Cell_Type_From_Ref']="Microglia"
RefTri2.integrated@meta.data[which(RefTri2.integrated@meta.data$integrated_snn_res.0.3==16),'Cell_Type_From_Ref']="UKN"
RefTri2.integrated@meta.data[which(RefTri2.integrated@meta.data$integrated_snn_res.0.3==6),'Cell_Type_From_Ref']="NeuroProgenitor"
RefTri2.integrated@meta.data[which(RefTri2.integrated@meta.data$integrated_snn_res.0.3==4),'Cell_Type_From_Ref']="GlialProgenitor"

pdf('./TestPlots/RefTri2_Transfered_Cell_Types.pdf',width=12,height=8)
print(DimPlot(RefTri2.integrated, label = TRUE) + NoLegend())
for(i in c("sample","Cell_Type_From_Ref","region","Maj_Clust_Int","SVM_clusterID","clust_all_name")){
print(DimPlot(RefTri2.integrated, reduction = "umap", group.by = i, label = TRUE, repel = TRUE))
}
dev.off()


save(RefTri2.integrated,file="./SeuratObj/Tri2_with_Ref_seurat_obj.rda",compress=TRUE)

