
library(Seurat)

setwd('/local/projects-t3/idea/bherb/Hypothalamus/PubRes')



###########################
#### Trimester 1 only 

### running Kmeans and GENIE3 on Trimester 1 samples (3/1/21)

load('./SeuratObj/Tri1Samples_PostSCT.rda')

Tri1CombDat = Tri1NormDat[2:7] ## don't include CS12 - fewer genes included in SCT than other samples. 

## large var gene list 
for (i in 1:length(Tri1CombDat)) {
    Tri1CombDat[[i]] <- NormalizeData(Tri1CombDat[[i]], verbose = FALSE)
    Tri1CombDat[[i]] <- FindVariableFeatures(Tri1CombDat[[i]], selection.method = "vst", nfeatures = 10000, verbose = FALSE)
}

for(i in 1:length(Tri1CombDat)) {
  if(i==1){
geneTot = VariableFeatures(Tri1CombDat[[i]])
  } else {
geneTot = c(geneTot,VariableFeatures(Tri1CombDat[[i]]))
  }
}

geneTot = names(table(geneTot))[which(table(geneTot)>3)]

for(i in 1:length(Tri1CombDat)){
geneTot = intersect(geneTot,as.character(rownames(Tri1CombDat[[i]])))
} #7003

for (i in 1:length(Tri1CombDat)) {
    Tri1CombDat[[i]] <- NormalizeData(Tri1CombDat[[i]], verbose = FALSE)
    #Tri1CombDat[[i]] <- FindVariableFeatures(Tri1CombDat[[i]], selection.method = "vst", nfeatures = 2000, verbose = FALSE)
}

options(future.globals.maxSize = 4000 * 1024^2)

AllTri1.anchors7003 <- FindIntegrationAnchors(object.list = Tri1CombDat, anchor.features = geneTot,dims = 1:30)
save(AllTri1.anchors7003,file='./Analysis/AllTri1.anchors_7003genes.rda',compress=TRUE)
AllTri1.integrated7003 <- IntegrateData(anchorset = AllTri1.anchors7003, dims = 1:30,features = geneTot)
DefaultAssay(AllTri1.integrated7003) <- "integrated"
AllTri1.integrated7003 <- ScaleData(AllTri1.integrated7003, verbose = FALSE)
AllTri1.integrated7003 <- RunPCA(AllTri1.integrated7003, npcs = 50, verbose = FALSE)
AllTri1.integrated7003 <- RunUMAP(AllTri1.integrated7003, reduction = "pca", dims = 1:50)
AllTri1.integrated7003 <- RunTSNE(AllTri1.integrated7003, reduction = "pca", dims = 1:50)

save(AllTri1.integrated7003,file='./SeuratObj/Tri1Samples_integrated7003.rda',compress=TRUE)

pdf('./TestPlots/AllTri1_integrated7003.pdf',width=12,height=8)
print(DimPlot(AllTri1.integrated7003, label = TRUE) + NoLegend())
for(i in c("sample")){
print(DimPlot(AllTri1.integrated7003, reduction = "tsne", group.by = i, label = TRUE, repel = TRUE))
}
dev.off()

dat7003 = as.matrix(AllTri1.integrated7003@assays$integrated[1:7003,1:20574])

orig7003 = as.matrix(AllTri1.integrated7003@assays$RNA[rownames(dat7003),1:20574])

## replace all zero counts and zero any negative numbers 
for(i in 1:ncol(dat7003)){
zeroInd = which(orig7003[,i]==0 | dat7003[,i]<0)
dat7003[zeroInd,i] = 0
if(i%%100==0) cat(paste(i,', ',sep=''))
}

dist7003 = as.matrix(dist(AllTri1.integrated7003@reductions[['pca']]@cell.embeddings[,1:50],diag=TRUE))

write.csv(t(dat7003),'/local/projects-t3/idea/bherb/Hypothalamus/PubRes/Networks/Tri1_7003_ExpressionMatrix.csv')

dat7003knn5 = smoother_aggregate_nearest_nb(mat = dat7003, D = dist7003, k = 5)

## round? 14 digits... need ~3? 

write.csv(t(dat7003knn5),'/local/projects-t3/idea/bherb/Hypothalamus/PubRes/Networks/Tri1_7003knn5_ExpressionMatrix.csv') 


######################
#### Trimester 2 only 

### gene list derived from eariler var gene call for individual sample genie3 calculations

geneTot = unique(c(as.character(HSTF$Name),as.character(genieGW25[,2]),as.character(genieGW25[,3]),as.character(genieGW22[,2]),as.character(genieGW22[,3]),as.character(genieGW20[,2]),as.character(genieGW20[,3]),as.character(genieGW19[,2]),as.character(genieGW19[,3]),as.character(genieGW18[,2]),as.character(genieGW18[,3])))


load('./SeuratObj/Tri2Samples_PostSCT.rda')

Tri2CombDat = Tri2NormDat

for(i in 1:5){
geneTot = intersect(geneTot,as.character(rownames(Tri2NormDat[[i]])))
} #7425


for (i in 1:length(Tri2CombDat)) {
    Tri2CombDat[[i]] <- NormalizeData(Tri2CombDat[[i]], verbose = FALSE)
    #Tri2CombDat[[i]] <- FindVariableFeatures(Tri2CombDat[[i]], selection.method = "vst", nfeatures = 2000, verbose = FALSE)
}

options(future.globals.maxSize = 4000 * 1024^2)

AllTri2.anchors7425 <- FindIntegrationAnchors(object.list = Tri2CombDat, anchor.features = geneTot,dims = 1:30)
save(AllTri2.anchors7425,file='./Analysis/AllTri2.anchors_7425genes.rda',compress=TRUE)
AllTri2.integrated7425 <- IntegrateData(anchorset = AllTri2.anchors, dims = 1:30,features = geneTot)
DefaultAssay(AllTri2.integrated7425) <- "integrated"
AllTri2.integrated7425 <- ScaleData(AllTri2.integrated7425, verbose = FALSE)
AllTri2.integrated7425 <- RunPCA(AllTri2.integrated7425, npcs = 50, verbose = FALSE)
AllTri2.integrated7425 <- RunUMAP(AllTri2.integrated7425, reduction = "pca", dims = 1:50)

save(AllTri2.integrated7425,file='./SeuratObj/Tri2Samples_integrated7425.rda',compress=TRUE)

dat7425 = as.matrix(AllTri2.integrated7425@assays$integrated[1:7425,1:25000])

dist7425 = as.matrix(dist(AllTri2.integrated7425@reductions[['pca']]@cell.embeddings[,1:50],diag=TRUE))

write.csv(t(dat7425),'/local/projects-t3/idea/bherb/Hypothalamus/PubRes/Networks/Tri2_7425_ExpressionMatrix.csv')

dat7425knn5 = smoother_aggregate_nearest_nb(mat = dat7425, D = dist7425, k = 5)

write.csv(t(dat7425knn5),'/local/projects-t3/idea/bherb/Hypothalamus/PubRes/Networks/Tri2_7425knn5_ExpressionMatrix.csv')


######################
### All samples 


totVar = unique(c(rownames(AllTri1.integrated7003),rownames(AllTri2.integrated7425))) #10441

rm(AllTri1.integrated7003,AllTri2.integrated7425)

load('./SeuratObj/AllSamples_integrated.rda') #All.integrated

obj = All.integrated

obj = subset(obj,cells=colnames(obj)[which(obj@meta.data$Clean_Set=='Yes')],features=totVar)

AllCombDat <- SplitObject(obj, split.by = "sample")

features <- totVar

for (i in 1:length(AllCombDat)) {
    AllCombDat[[i]] <- NormalizeData(AllCombDat[[i]], verbose = FALSE)
    AllCombDat[[i]] <- ScaleData(AllCombDat[[i]], features = features, verbose = FALSE)
    AllCombDat[[i]] <- RunPCA(AllCombDat[[i]], features = features, verbose = FALSE)
}

options(future.globals.maxSize = 80000 * 1024^2)

All.anchors <- FindIntegrationAnchors(object.list = AllCombDat, anchor.features = features, reduction = "rpca")

save(All.anchors,file='./SeuratObj/AllSamples10441_Anchors.rda',compress=TRUE)

## here - next step requires >100 Gb memory!??

All.combined <- IntegrateData(anchorset = All.anchors)

DefaultAssay(All.combined) <- "integrated"

save(All.combined,file='./SeuratObj/AllSamples10441_Integrated.rda',compress=TRUE)

# Run the standard workflow for visualization and clustering
All.combined <- ScaleData(All.combined, verbose = FALSE)
All.combined <- RunPCA(All.combined, npcs = 30, verbose = FALSE)
All.combined <- RunUMAP(All.combined, reduction = "pca", dims = 1:30)
All.combined <- FindNeighbors(All.combined, reduction = "pca", dims = 1:30)
All.combined <- FindClusters(All.combined, resolution = 0.5)




