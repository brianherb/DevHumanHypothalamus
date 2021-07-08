
library(Seurat)
library(SingleCellExperiment)
library(slingshot)
library(RColorBrewer)
library(gam)
library(ggplot2)

setwd('/local/projects-t3/idea/bherb/Hypothalamus/PubRes')

load('./SeuratObj/AllSamples_integrated.rda') #All.integrated

## all data?  - starts with 'Dividing_OligDominant_2'

All.integrated@active.ident=as.factor(All.integrated$SubClass)

tmpIdent = unique(All.integrated@active.ident)[grep('Olig',unique(All.integrated@active.ident),fixed=TRUE)]

tmpIdent = tmpIdent[-grep('Astro',tmpIdent)]

Oligo_Tmp = subset(All.integrated,idents=tmpIdent)

DefaultAssay(Oligo_Tmp) = 'integrated' ## 4258 cells

Oligo_PCA = Oligo_Tmp@reductions[["pca"]]@cell.embeddings[,1:30]

Oligo_Clust = Oligo_Tmp$SubClass

#Oligo_SCE = Seurat::as.SingleCellExperiment(Oligo_Tmp) ## just 2000 var genes - try for now? - transfered PCA
Oligo_SCE = SingleCellExperiment(Oligo_Tmp[['integrated']],colData=Oligo_Tmp@meta.data)

reducedDims(Oligo_SCE) = list(PCA=Oligo_PCA)


Oligo_SCE <- slingshot(Oligo_SCE, clusterLabels = 'SubClass', reducedDim = 'PCA')


save(Oligo_SCE,file='./Analysis/Oligo_oligo_SCE_lineage.rda')

colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(Oligo_SCE$slingPseudotime_1, breaks=100)]
plotcol2 <- colors[as.numeric(as.factor(Oligo_SCE$SubClass))]

pdf(file='./TestPlots/Oligo_OligoLineage_PC1_2_SubclassHannah.pdf',width=12,height=8)
plot(reducedDims(Oligo_SCE)$PCA, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(Oligo_SCE), lwd=2, col='black')
plot(reducedDims(Oligo_SCE)$PCA, col = brewer.pal(9,'Set1')[Oligo_SCE$seurat_clusters], pch=16, asp = 1)
lines(SlingshotDataSet(Oligo_SCE), lwd=2, type = 'lineages', col = 'black')
plot(reducedDims(Oligo_SCE)$PCA, col = plotcol2, pch=16, asp = 1)
lines(SlingshotDataSet(Oligo_SCE), lwd=2, col='black')
dev.off()


lin1 <- getLineages(Oligo_PCA, Oligo_Clust, start.clus = 'Dividing_OligDominant_2')


Oligo.t <- Oligo_SCE$slingPseudotime_1

# for time, only look at the 100 most variable genes
#Y <- log1p(assays(GW18AllGE_Ctx_SCE)$norm)
#var100 <- names(sort(apply(Y,1,var),decreasing = TRUE))[1:100]
#Y <- Y[var100,]
Oligo.Y = as.matrix(Oligo_SCE@assays@data[[1]][1:2000,1:3519])

# 


# fit a GAM with a loess term for pseudotime
Oligo.pval <- apply(Oligo.Y,1,function(z){
    d <- data.frame(z=z, t=Oligo.t)
    suppressWarnings({
      tmp <- suppressWarnings(gam(z ~ lo(t), data=d))
    })
    p <- summary(tmp)[3][[1]][2,3]
    p
})

##FDR? 

Oligo.fdr = p.adjust(Oligo.pval,method='fdr')

Oligo.cor = NA

for(i in 1:nrow(Oligo.Y)){
Oligo.cor[i] = cor(Oligo.t,Oligo.Y[i,],use="complete.obs")
if(i%%100==0) cat(paste(i,', ',sep=''))
}

rownames(Oligo.Y)[which(Oligo.cor>0.65 & Oligo.fdr<0.05)]
rownames(Oligo.Y)[which(Oligo.cor<(-0.85) & Oligo.fdr<0.05)]


## skip for now
Y2 = Oligo_Tmp@assays$SCT[1:20309,1:2370]

# fit a GAM with a loess term for pseudotime
Oligo.pval2 <- apply(Y2,1,function(z){
    d <- data.frame(z=z, t=t)
    suppressWarnings({
      tmp <- suppressWarnings(gam(z ~ lo(t), data=d))
    })
    p <- summary(tmp)[3][[1]][2,3]
    p
})

##FDR? 

Oligo.fdr2 = p.adjust(Oligo.pval2,method='fdr')

Oligo.cor2 = NA

for(i in 1:nrow(Y2)){
Oligo.cor2[i] = cor(t,Y2[i,])
if(i%%100==0) cat(paste(i,', ',sep=''))
}

rownames(Y2)[which(Oligo.cor2>0.7 & Oligo.fdr2<0.05)]
rownames(Y2)[which(Oligo.cor2<(-0.5) & Oligo.fdr2<0.05)]


## also calculate cor 

## check pseudotime

All.integrated@meta.data$Oligo_pseudotime = 0

All.integrated@meta.data[colnames(Oligo_SCE),'Oligo_pseudotime'] = Oligo_SCE$slingPseudotime_1


pdf('./TestPlots/Oligo_Oligo_Lineage_SubclassHannah.pdf',width=12,height=8)
print(DimPlot(All.integrated, label = TRUE) + NoLegend())
for(i in c("sample","MajorClass","SubClass")){
print(DimPlot(All.integrated, reduction = "umap", group.by = i, label = TRUE, repel = TRUE))
}
print(FeaturePlot(All.integrated, features = 'Oligo_pseudotime', pt.size = 0.4, order=TRUE))
for(k in rownames(Y)[which(Oligo.cor>0.65 & Oligo.fdr<0.05)]){
print(FeaturePlot(All.integrated, features = k, pt.size = 0.4, order=TRUE))
}
for(k in rownames(Y)[which(Oligo.cor<(-0.85) & Oligo.fdr<0.05)]){
print(FeaturePlot(All.integrated, features = k, pt.size = 0.4, order=TRUE))
}
dev.off()

### line plots 


OligLingenes = unique(c(rownames(Oligo.Y)[which(Oligo.cor>0.4 & Oligo.fdr<0.05)],rownames(Oligo.Y)[which(Oligo.cor<(-0.5) & Oligo.fdr<0.05)]))

OligLingenes = intersect(OligLingenes,TF_NP$Transcription_factors)

OligLindata=cbind(Oligo.t,t(Oligo.Y[OligLingenes,]))
OligLindata=as.data.frame(OligLindata)
colnames(OligLindata) = c('pseudotime',OligLingenes)

OligLindata$cluster = All.integrated@meta.data[colnames(Oligo_SCE),'SubClass']

OligLindataM = reshape2::melt(data=OligLindata, id.vars = "pseudotime", measure.vars = OligLingenes)

pdf('./TestPlots/Tri12_Oligo_TestPseudotime_SmoothLines.pdf',width=12,height=8)
print(ggplot(OligLindataM, aes(x = pseudotime, y = value, colour = variable)) + geom_smooth())
dev.off()


### genes marking early oligo lineage

All.integrated@meta.data$EarlyOligo = 'NO'


All.integrated@meta.data[which(All.integrated@meta.data$MajorClass=='Dividing_Progenitor' & All.integrated@meta.data$ColorByMod=='Oligodendrocyte'),'EarlyOligo'] = 'YES'

All.integrated@active.ident=factor(All.integrated$EarlyOligo)

DefaultAssay(All.integrated) = 'RNA'
EarlyOligoDEG = FindAllMarkers(All.integrated)


pdf('./TestPlots/Oligo_Early_Lineage_Genes.pdf',width=12,height=8)
print(DimPlot(All.integrated, label = TRUE) + NoLegend())
for(i in c("sample","MajorClass","SubClass")){
print(DimPlot(All.integrated, reduction = "umap", group.by = i, label = TRUE, repel = TRUE))
}
print(FeaturePlot(All.integrated, features = 'Oligo_pseudotime', pt.size = 0.4, order=TRUE))
for(k in EarlyOligoDEG$gene[which(EarlyOligoDEG$avg_log2FC<0)[1:20]]){
print(FeaturePlot(All.integrated, features = k, pt.size = 0.4, order=TRUE))
}

dev.off()


pdf('./TestPlots/Oligo_Early_Lineage_Genes_Set2.pdf',width=12,height=8)
print(DimPlot(All.integrated, label = TRUE) + NoLegend())
for(i in c("sample","MajorClass","SubClass")){
print(DimPlot(All.integrated, reduction = "umap", group.by = i, label = TRUE, repel = TRUE))
}
print(FeaturePlot(All.integrated, features = 'Oligo_pseudotime', pt.size = 0.4, order=TRUE))
for(k in EarlyOligoDEG$gene[which(EarlyOligoDEG$avg_log2FC<0)[21:40]]){
print(FeaturePlot(All.integrated, features = k, pt.size = 0.4, order=TRUE))
}

dev.off()



### early lineage DEG's? 


All.integrated.Early = subset(All.integrated,cells=colnames(All.integrated)[which(All.integrated@meta.data$MajorClass=='Dividing_Progenitor')])

All.integrated.Early = All.integrated.Early[,-which(is.na(All.integrated.Early@meta.data$ColorByMod))]

All.integrated.Early@active.ident=factor(All.integrated.Early$ColorByMod)

DefaultAssay(All.integrated.Early) = 'RNA'
EarlyDEG = FindAllMarkers(All.integrated.Early)

EarlyDEGTF = FindAllMarkers(All.integrated.Early,features=intersect(HSTF$Name,rownames(All.integrated.Early)))

write.csv(EarlyDEG,file='./Analysis/EarlyLineages_DEG.csv')


write.csv(EarlyDEGTF,file='./Analysis/EarlyLineages_DEG_TF.csv')


AstEarly = EarlyDEGTF$gene[which( EarlyDEGTF$avg_log2FC>0 & EarlyDEGTF$cluster=='Astrocyte')][1:10]

#EarlyDEGTF$pct.2<0.1 &

EpendEarly = EarlyDEGTF$gene[which( EarlyDEGTF$avg_log2FC>0 & EarlyDEGTF$cluster=='Ependymal')][1:10]


ImNeuroEarly = EarlyDEGTF$gene[which( EarlyDEGTF$avg_log2FC>0 & EarlyDEGTF$cluster=='ImmatureNeuron')][1:10]

MatNeuroEarly = EarlyDEGTF$gene[which( EarlyDEGTF$avg_log2FC>0 & EarlyDEGTF$cluster=='MatureNeuron')][1:10]

OligoEarly = EarlyDEGTF$gene[which( EarlyDEGTF$avg_log2FC>0 & EarlyDEGTF$cluster=='Oligodendrocyte')][1:10]

EarlyTF = unique(c(AstEarly,EpendEarly,ImNeuroEarly,MatNeuroEarly,OligoEarly))


ClustOrder = rev(sort(unique(as.character(All.integrated.Early@active.ident))))


All.integrated.Early@active.ident=factor(All.integrated.Early@active.ident,levels=ClustOrder)


pdf('./TestPlots/All_EarlyTFs_DotPlot.pdf',width=20,height=9)
print(DotPlot(All.integrated.Early, features =  EarlyTF,idents=c('Astrocyte','Ependymal','ImmatureNeuron','MatureNeuron','Oligodendrocyte')) + RotatedAxis())
dev.off()


ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

color_list <- ggplotColours(n=5)

TestLin = sort(unique(All.integrated.Early@meta.data$ColorByMod))[-5]

pdf('./TestPlots/Color_By_Modules_IndvPlots2_justEarly.pdf',width=12,height=8)
for(j in 1:length(TestLin)){
All.integrated.Early@meta.data$TmpCol = NA
All.integrated.Early@meta.data$TmpCol[which(All.integrated.Early@meta.data$ColorByMod==TestLin[j])] = TestLin[j]
#print(FeaturePlot(All.integrated, features = 'TmpCol', pt.size = 0.4,cols=c(color_list[j],"lightgrey"),order=TRUE))
print(DimPlot(All.integrated.Early, reduction = "umap", group.by = 'TmpCol', cols=c(color_list[j],"gray90"), label = TRUE, repel = TRUE))
}
dev.off()











### modules or GRN's? 

## Top ME's -  








#### check correlation with gene modules 

MEs = readRDS('./Analysis/AllSamp_mergeColors.k100.MEs_CoreCellGene.rds') 

Eigen = MEs[['eigengenes']]

OligoEigen = Eigen[colnames(All.integrated)[which(All.integrated$Oligo_pseudotime!=0)],]

OligoCor = rep(0,ncol(OligoEigen))

names(OligoCor) = colnames(OligoEigen)

OligoPseudo = All.integrated@meta.data$Oligo_pseudotime[which(All.integrated$Oligo_pseudotime!=0)]

naind = which(is.na(OligoEigen[,1]))

OligoPseudo = OligoPseudo[-naind]

OligoEigen = OligoEigen[-naind,]


for(i in names(OligoCor)){
OligoCor[i] = cor(as.vector(OligoEigen[,i]),OligoPseudo)

}

## missing things that increase in middle? 

## line plots? 


OligLinMEs = c(names(sort(OligoCor))[1:5],names(sort(OligoCor,decreasing=TRUE))[1:5])

OligLinMEdata=cbind(All.integrated@meta.data[rownames(OligoEigen),'Oligo_pseudotime'],OligoEigen[,OligLinMEs])
OligLinMEdata=as.data.frame(OligLinMEdata)
colnames(OligLinMEdata) = c('pseudotime',OligLinMEs)

OligLinMEdata$cluster = All.integrated@meta.data[rownames(OligLinMEdata),'SubClass']

OligLinMEdataM = reshape2::melt(data=OligLinMEdata, id.vars = "pseudotime", measure.vars = OligLinMEs)

pdf('./TestPlots/Tri12_Oligo_TestPseudotime_MEs_SmoothLines.pdf',width=12,height=8)
print(ggplot(OligLinMEdataM, aes(x = pseudotime, y = value, colour = variable)) + geom_smooth())
dev.off()


## GRNs 















### dev lineage? 

## Dividing_1 first


All.integrated@active.ident=as.factor(All.integrated$SubClass)

devIdent = unique(All.integrated@active.ident)[c(grep('Dividing',unique(All.integrated@active.ident),fixed=TRUE),grep('Radial',unique(All.integrated@active.ident),fixed=TRUE))]

devIdent = devIdent[-grep('OligDom',devIdent)]

devIdent = devIdent[-grep('Micro',devIdent)]

devIdent = devIdent[-grep('VLMC',devIdent)]

All_dev = subset(All.integrated,idents=devIdent)

DefaultAssay(All_dev) = 'integrated' ## 4258 cells

All_dev_PCA = All_dev@reductions[["pca"]]@cell.embeddings[,1:30]

All_dev_Clust = All_dev$SubClass

#All_SCE = Seurat::as.SingleCellExperiment(All_dev) ## just 2000 var genes - try for now? - transfered PCA
All_dev_SCE = SingleCellExperiment(All_dev[['integrated']],colData=All_dev@meta.data)

reducedDims(All_dev_SCE) = list(PCA=All_dev_PCA)


All_dev_SCE <- slingshot(All_dev_SCE, clusterLabels = 'SubClass', reducedDim = 'PCA')

save(All_dev_SCE,file='./Analysis/All_dev_SCE_lineage.rda')

colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(All_dev_SCE$slingPseudotime_1, breaks=100)]
plotcol2 <- colors[as.numeric(as.factor(All_dev_SCE$SubClass))]

pdf(file='./TestPlots/All_DevLineage_PC1_2_SubclassHannah.pdf',width=12,height=8)
plot(reducedDims(All_dev_SCE)$PCA, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(All_dev_SCE), lwd=2, col='black')
plot(reducedDims(All_dev_SCE)$PCA, col = brewer.pal(9,'Set1')[All_dev_SCE$seurat_clusters], pch=16, asp = 1)
lines(SlingshotDataSet(All_dev_SCE), lwd=2, type = 'lineages', col = 'black')
plot(reducedDims(All_dev_SCE)$PCA, col = plotcol2, pch=16, asp = 1)
lines(SlingshotDataSet(All_dev_SCE), lwd=2, col='black')
dev.off()


lin1 <- getLineages(All_dev_PCA, All_dev_Clust, start.clus = 'Dividing_1')


t <- All_dev_SCE$slingPseudotime_1

# for time, only look at the 100 most variable genes
#Y <- log1p(assays(GW18AllGE_Ctx_SCE)$norm)
#var100 <- names(sort(apply(Y,1,var),decreasing = TRUE))[1:100]
#Y <- Y[var100,]
Y = as.matrix(All_dev_SCE@assays@data[[1]][1:2000,1:11318])

# 


# fit a GAM with a loess term for pseudotime
Dev.pval <- apply(Y,1,function(z){
    d <- data.frame(z=z, t=t)
    suppressWarnings({
      dev <- suppressWarnings(gam(z ~ lo(t), data=d))
    })
    p <- summary(dev)[3][[1]][2,3]
    p
})

##FDR? 

Dev.fdr = p.adjust(Dev.pval,method='fdr')


Dev.cor = NA

for(i in 1:nrow(Y)){
Dev.cor[i] = cor(t,Y[i,],use="complete.obs")
if(i%%100==0) cat(paste(i,', ',sep=''))
}

rownames(Y)[which(Dev.cor>0.65 & Dev.fdr<0.05)]
rownames(Y)[which(Dev.cor<(-0.85) & Dev.fdr<0.05)]



## check pseudotime

All.integrated@meta.data$Dev_pseudotime1 = 0

All.integrated@meta.data[colnames(All_dev_SCE),'Dev_pseudotime1'] = All_dev_SCE$slingPseudotime_1

All.integrated@meta.data$Dev_pseudotime2 = 0

All.integrated@meta.data[colnames(All_dev_SCE),'Dev_pseudotime2'] = All_dev_SCE$slingPseudotime_2

All.integrated@meta.data$Dev_pseudotime3 = 0

All.integrated@meta.data[colnames(All_dev_SCE),'Dev_pseudotime3'] = All_dev_SCE$slingPseudotime_3

All.integrated@meta.data$Dev_pseudotime4 = 0

All.integrated@meta.data[colnames(All_dev_SCE),'Dev_pseudotime4'] = All_dev_SCE$slingPseudotime_4

All.integrated@meta.data$Dev_pseudotimeMean = 0

All.integrated@meta.data[colnames(All_dev_SCE),'Dev_pseudotimeMean'] = rowMeans(as.data.frame(colData(All_dev_SCE)[,c('slingPseudotime_1','slingPseudotime_2','slingPseudotime_3','slingPseudotime_4')]),na.rm=TRUE)




pdf('./TestPlots/All_Dev_Lineage_SubclassHannah.pdf',width=12,height=8)
print(DimPlot(All.integrated, label = TRUE) + NoLegend())
for(i in c("sample","MajorClass","SubClass")){
print(DimPlot(All.integrated, reduction = "umap", group.by = i, label = TRUE, repel = TRUE))
}
for(i in c("Dev_pseudotime1","Dev_pseudotime2","Dev_pseudotime3","Dev_pseudotime4","Dev_pseudotimeMean")){
print(FeaturePlot(All.integrated, features = i, pt.size = 0.4, order=TRUE))
}
for(k in rownames(Y)[which(Dev.cor>0.65 & Dev.fdr<0.05)]){
print(FeaturePlot(All.integrated, features = k, pt.size = 0.4, order=TRUE))
}
for(k in rownames(Y)[which(Dev.cor<(-0.85) & Dev.fdr<0.05)]){
print(FeaturePlot(All.integrated, features = k, pt.size = 0.4, order=TRUE))
}
dev.off()





## Try slingshot on Eigengenes 

## start - Dividing_1


All.integrated@active.ident=as.factor(All.integrated$SubClass)

DefaultAssay(All.integrated) = 'integrated' ## 4258 cells

EigenAll_PCA = All.integrated@reductions[["pca"]]@cell.embeddings[rownames(Eigen),1:30]

EigenAll_Clust = All.integrated$SubClass[rownames(Eigen)]

#Oligo_SCE = Seurat::as.SingleCellExperiment(Oligo_Tmp) ## just 2000 var genes - try for now? - transfered PCA
EigenAll_SCE = SingleCellExperiment(t(Eigen),colData=All.integrated@meta.data[rownames(Eigen),c('sample','MajorClass','SubClass','Neuron','Region')])

reducedDims(EigenAll_SCE) = list(PCA=EigenAll_PCA)


EigenAll_SCE <- slingshot(EigenAll_SCE, clusterLabels = 'SubClass', reducedDim = 'PCA')

save(EigenAll_SCE,file='./Analysis/EigenAll_SCE_lineage.rda')

colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(EigenAll_SCE$slingPseudotime_1, breaks=100)]
plotcol2 <- colors[as.numeric(as.factor(EigenAll_SCE$SubClass))]

pdf(file='./TestPlots/EigenAll_Lineage_SubclassHannah.pdf',width=12,height=8)
plot(reducedDims(EigenAll_SCE)$PCA, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(EigenAll_SCE), lwd=2, col='black')
plot(reducedDims(EigenAll_SCE)$PCA, col = brewer.pal(9,'Set1')[EigenAll_SCE$seurat_clusters], pch=16, asp = 1)
lines(SlingshotDataSet(EigenAll_SCE), lwd=2, type = 'lineages', col = 'black')
plot(reducedDims(EigenAll_SCE)$PCA, col = plotcol2, pch=16, asp = 1)
lines(SlingshotDataSet(EigenAll_SCE), lwd=2, col='black')
dev.off()


lin1 <- getLineages(EigenAll_PCA, EigenAll_Clust, start.clus = 'Dividing_OligDominant_2')


EigenAll.t <- EigenAll_SCE$slingPseudotime_1





All.integrated@meta.data$Dev_pseudotime1 = 0

All.integrated@meta.data[colnames(All_dev_SCE),'Dev_pseudotime1'] = All_dev_SCE$slingPseudotime_1

All.integrated@meta.data$Dev_pseudotime2 = 0

All.integrated@meta.data[colnames(All_dev_SCE),'Dev_pseudotime2'] = All_dev_SCE$slingPseudotime_2

All.integrated@meta.data$Dev_pseudotime3 = 0

All.integrated@meta.data[colnames(All_dev_SCE),'Dev_pseudotime3'] = All_dev_SCE$slingPseudotime_3

All.integrated@meta.data$Dev_pseudotime4 = 0

All.integrated@meta.data[colnames(All_dev_SCE),'Dev_pseudotime4'] = All_dev_SCE$slingPseudotime_4

All.integrated@meta.data$Dev_pseudotimeMean = 0

All.integrated@meta.data[colnames(All_dev_SCE),'Dev_pseudotimeMean'] = rowMeans(as.data.frame(colData(All_dev_SCE)[,c('slingPseudotime_1','slingPseudotime_2','slingPseudotime_3','slingPseudotime_4')]),na.rm=TRUE)



slingTimes = colnames(colData(EigenAll_SCE))[grep('slingPseud',colnames(colData(EigenAll_SCE)))]

##104 col


pdf('./TestPlots/All_Eigen_Lineage_SubclassHannah.pdf',width=12,height=8)
print(DimPlot(All.integrated, label = TRUE) + NoLegend())
for(i in c("sample","MajorClass","SubClass")){
print(DimPlot(All.integrated, reduction = "umap", group.by = i, label = TRUE, repel = TRUE))
}
for(i in slingTimes){

All.integrated@meta.data[,i] = 0

All.integrated@meta.data[colnames(EigenAll_SCE),i] = colData(EigenAll_SCE)[,i]

print(FeaturePlot(All.integrated, features = i, pt.size = 0.4, order=TRUE))
}
dev.off()











### Junque 


#### Trimester 2 Oligo lineage - I did this with monocle at some point, but might be simplier with slingshot since I'm expecting a single lineage

AllTri2.integrated@active.ident=as.factor(AllTri2.integrated$SubClass)

AllTri2_Tmp = subset(AllTri2.integrated,idents=unique(AllTri2.integrated@active.ident)[grep('Olig',unique(AllTri2.integrated@active.ident),fixed=TRUE)])

DefaultAssay(AllTri2_Tmp) = 'integrated' #3001 cells 

#AllTri2_Tmp@assays$RNA = AllTri2_Tmp@assays$RNA[rownames(AllTri2_Tmp@assays$integrated),]

#AllTri2_Tmp@assays$SCT = AllTri2_Tmp@assays$SCT[rownames(AllTri2_Tmp@assays$integrated),]

#AllTri2_Tmp[['SCT']] <- NULL
#AllTri2_Tmp[['RNA']] <- NULL


AllTri2_PCA = AllTri2_Tmp@reductions[["pca"]]@cell.embeddings[,1:30]

AllTri2_Clust = AllTri2_Tmp$Maj_Clust_From_Ref

#AllTri2_SCE = Seurat::as.SingleCellExperiment(AllTri2_Tmp) ## just 2000 var genes - try for now? - transfered PCA
AllTri2_SCE = SingleCellExperiment(AllTri2_Tmp[['integrated']],colData=AllTri2_Tmp@meta.data)

reducedDims(AllTri2_SCE) = list(PCA=AllTri2_PCA)


AllTri2_SCE <- slingshot(AllTri2_SCE, clusterLabels = 'Maj_Clust_From_Ref', reducedDim = 'PCA')



colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(AllTri2_SCE$slingPseudotime_1, breaks=100)]

pdf(file='./TestPlots/AllTri2_OligoLineage_PC1_2.pdf',width=12,height=8)
plot(reducedDims(AllTri2_SCE)$PCA, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(AllTri2_SCE), lwd=2, col='black')
plot(reducedDims(AllTri2_SCE)$PCA, col = brewer.pal(9,'Set1')[AllTri2_SCE$seurat_clusters], pch=16, asp = 1)
lines(SlingshotDataSet(AllTri2_SCE), lwd=2, type = 'lineages', col = 'black')
dev.off()


lin1 <- getLineages(AllTri2_PCA, AllTri2_Clust, start.clus = 'ImmatureOligo')


t <- AllTri2_SCE$slingPseudotime_1

# for time, only look at the 100 most variable genes
#Y <- log1p(assays(GW18AllGE_Ctx_SCE)$norm)
#var100 <- names(sort(apply(Y,1,var),decreasing = TRUE))[1:100]
#Y <- Y[var100,]
Y = as.matrix(AllTri2_SCE@assays@data[[1]][1:2000,1:2370])

# fit a GAM with a loess term for pseudotime
Oligo.pval <- apply(Y,1,function(z){
    d <- data.frame(z=z, t=t)
    suppressWarnings({
      tmp <- suppressWarnings(gam(z ~ lo(t), data=d))
    })
    p <- summary(tmp)[3][[1]][2,3]
    p
})

##FDR? 

Oligo.fdr = p.adjust(Oligo.pval,method='fdr')

Oligo.cor = NA

for(i in 1:nrow(Y)){
Oligo.cor[i] = cor(t,Y[i,])
if(i%%100==0) cat(paste(i,', ',sep=''))
}

rownames(Y)[which(Oligo.cor>0.7 & Oligo.fdr<0.05)]
rownames(Y)[which(Oligo.cor<(-0.5) & Oligo.fdr<0.05)]


Y2 = AllTri2_Tmp@assays$SCT[1:20309,1:2370]

# fit a GAM with a loess term for pseudotime
Oligo.pval2 <- apply(Y2,1,function(z){
    d <- data.frame(z=z, t=t)
    suppressWarnings({
      tmp <- suppressWarnings(gam(z ~ lo(t), data=d))
    })
    p <- summary(tmp)[3][[1]][2,3]
    p
})

##FDR? 

Oligo.fdr2 = p.adjust(Oligo.pval2,method='fdr')

Oligo.cor2 = NA

for(i in 1:nrow(Y2)){
Oligo.cor2[i] = cor(t,Y2[i,])
if(i%%100==0) cat(paste(i,', ',sep=''))
}

rownames(Y2)[which(Oligo.cor2>0.7 & Oligo.fdr2<0.05)]
rownames(Y2)[which(Oligo.cor2<(-0.5) & Oligo.fdr2<0.05)]


## also calculate cor 

## check pseudotime

AllTri2.integrated@meta.data$Oligo_pseudotime = 0

AllTri2.integrated@meta.data[colnames(AllTri2_SCE),'Oligo_pseudotime'] = AllTri2_SCE$slingPseudotime_1


pdf('./TestPlots/AllTri2_Oligo_Lineage.pdf',width=12,height=8)
print(DimPlot(AllTri2.integrated, label = TRUE) + NoLegend())
for(i in c("sample","Maj_Clust_From_Ref","neuron_name")){
print(DimPlot(AllTri2.integrated, reduction = "umap", group.by = i, label = TRUE, repel = TRUE))
}
print(FeaturePlot(AllTri2.integrated, features = 'Oligo_pseudotime', pt.size = 0.4, order=TRUE))
for(k in rownames(Y)[which(Oligo.cor>0.7 & Oligo.fdr<0.05)]){
print(FeaturePlot(AllTri2.integrated, features = k, pt.size = 0.4, order=TRUE))
}
for(k in rownames(Y)[which(Oligo.cor<(-0.5) & Oligo.fdr<0.05)]){
print(FeaturePlot(AllTri2.integrated, features = k, pt.size = 0.4, order=TRUE))
}
dev.off()



### looks very good - captures genes that are turned on, fewer are specific to oligo progenitors  - look for genes that correlate with PDGFRA

pdf('./TestPlots/AllTri2_Oligo_Lineage_SCT.pdf',width=12,height=8)
print(DimPlot(AllTri2.integrated, label = TRUE) + NoLegend())
for(i in c("sample","Maj_Clust_From_Ref","neuron_name")){
print(DimPlot(AllTri2.integrated, reduction = "umap", group.by = i, label = TRUE, repel = TRUE))
}
print(FeaturePlot(AllTri2.integrated, features = 'Oligo_pseudotime', pt.size = 0.4, order=TRUE))
for(k in rownames(Y2)[which(Oligo.cor2>0.7 & Oligo.fdr2<0.05)]){
print(FeaturePlot(AllTri2.integrated, features = k, pt.size = 0.4, order=TRUE))
}
for(k in rownames(Y2)[which(Oligo.cor2<(-0.5) & Oligo.fdr2<0.05)]){
print(FeaturePlot(AllTri2.integrated, features = k, pt.size = 0.4, order=TRUE))
}
dev.off()


PDGFRA.cor = NA

for(i in 1:nrow(Y)){
PDGFRA.cor[i] = cor(Y[74,],Y[i,],method='spearman')
if(i%%100==0) cat(paste(i,', ',sep=''))
} ### top 30 aren't great except XRCC2, and that is just DNA repair 


pdf('./TestPlots/AllTri2_Oligo_Lineage_PDGFRAlike.pdf',width=12,height=8)
print(DimPlot(AllTri2.integrated, label = TRUE) + NoLegend())
for(i in c("sample","Maj_Clust_From_Ref","neuron_name")){
print(DimPlot(AllTri2.integrated, reduction = "umap", group.by = i, label = TRUE, repel = TRUE))
}
print(FeaturePlot(AllTri2.integrated, features = 'Oligo_pseudotime', pt.size = 0.4, order=TRUE))
for(k in rownames(Y)[order(PDGFRA.cor,decreasing=TRUE)][1:30]){
print(FeaturePlot(AllTri2.integrated, features = k, pt.size = 0.4, order=TRUE))
}
dev.off()


###edit below - for making smooth lines across pseudotime 

Oligo.t <- AllTri2_SCE$slingPseudotime_1

Oligo.Y = AllTri2_SCE@assays@data[[1]]

Tri2OligLingenes = c(unique(c(rownames(Y2)[which(Oligo.cor2>0.7 & Oligo.fdr2<0.05)],rownames(Y)[which(Oligo.cor>0.7 & Oligo.fdr<0.05)])),unique(c(rownames(Y2)[which(Oligo.cor2<(-0.5) & Oligo.fdr2<0.05)],rownames(Y)[which(Oligo.cor<(-0.5) & Oligo.fdr<0.05)])))

Tri2OligLindata=cbind(Oligo.t,t(Oligo.Y[Tri2OligLingenes,]))
Tri2OligLindata=as.data.frame(Tri2OligLindata)
colnames(Tri2OligLindata) = c('pseudotime',Tri2OligLingenes)

Tri2OligLindata$cluster = AllTri2.integrated@meta.data[colnames(AllTri2_SCE),'Maj_Clust_From_Ref']

Tri2OligLindataM = reshape2::melt(data=Tri2OligLindata, id.vars = "pseudotime", measure.vars = Tri2OligLingenes)

pdf('./TestPlots/Tri2_Oligo_TestPseudotime_SmoothLines.pdf',width=12,height=8)
print(ggplot(Tri2OligLindataM, aes(x = pseudotime, y = value, colour = variable)) + geom_smooth())
dev.off()


save(Tri2Oligo.integrated,file='./SeuratObj/Tri2Oligo_ds.integrated.rda',compress=TRUE) ## this was from old version of script...
save(Tri2Oligo.anchors,file='./SeuratObj/Tri2Oligo_ds.anchors.rda',compress=TRUE)
save(Tri2Oligo.monocle,file='./SeuratObj/Tri2Oligo_ds.monocle.rda',compress=TRUE)
save(Tri2Oligo_List,file='./SeuratObj/Tri2Oligo_ds_List.rda',compress=TRUE)
save(Tri2Oligo.markers.NMs,Tri2Oligo.markers.NPs,Tri2Oligo.markers.NRs,Tri2Oligo.markers.TFs,file='./SeuratObj/Tri2Oligo_ds_Markers.rda',compress=TRUE)
save(Oligo.cor,Oligo.fdr,Oligo.pval,Oligo.cor2,Oligo.fdr2,Oligo.pval2,file='./SeuratObj/Tri2Oligo_pseudotime_Correlations.rda',compress=TRUE)


