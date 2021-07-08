library(Seurat)
library(AUCell)
library(biomaRt)
library(Category)
library(org.Hs.eg.db)
library(GOstats)
library(GO.db)

setwd('/local/projects-t3/idea/bherb/Hypothalamus/PubRes')

##########################
#### Load Gene modules that were generated from UMAP relationships in full dataset 

#cl = readRDS('./Analysis/AllSamp_k100Try2.rds') #6677


MEs = readRDS('./Analysis/AllSamp_mergeColors.k100.MEs_CoreCellGene.rds') 

## modules for whole dataset, excluding fibroblast, microglia, and epithelial – found 54 modules across all samples with 10441 genes, 7637 of which are in valid modules and 2804 are in group 0. 


##################
## Load datasets

load('./SeuratObj/AllSamples_integrated.rda') #All.integrated

load('./SeuratObj/Tri2Samples_integrated.rda') #AllTri2.integrated

load('./SeuratObj/Tri1Samples_integrated.rda') #AllTri1.integrated

load('./SeuratObj/Tri2Neuron_ds.integrated.rda')  ## Tri2Neuron_ds.integrated


############################
### GO gene enrichments of tissue specific genes in Radial glia populations 

totgene = data.frame(gene=names(MEs$validColors),stringsAsFactors=FALSE,ME = paste('ME',MEs$validColors,sep=''))

mart = useMart( "ensembl" )
mart = useDataset("hsapiens_gene_ensembl", mart = mart )

gid = getBM(mart = mart, attributes = c("hgnc_symbol","entrezgene_id"),values=totgene) #entrez to refseq

totgene$entrez = gid$entrezgene_id[match(totgene$gene,gid$hgnc_symbol)]

modules = paste('ME',sort(unique(MEs$validColors)),sep='')[-1] ## no ME0
GOcat = c('BP','MF','CC')

### biological process
ME54GO = vector(mode='list',length=length(modules)*length(GOcat))

names(ME54GO) = paste(rep(modules,each=3),rep(GOcat,length(modules)),sep='_')

for(i in modules){
for(k in GOcat){

testGenes = na.omit(totgene$entrez[which(totgene$ME==i)])
universeGenes = unique(na.omit(totgene$entrez))

params <- new('GOHyperGParams',
              geneIds = testGenes,
              universeGeneIds = universeGenes,
              ontology = k,
              pvalueCutoff = 0.05,
              conditional = FALSE,
              testDirection = 'over',
              annotation = "org.Hs.eg.db"
             )
 
 ME54GO[[paste(i,k,sep='_')]] <- hyperGTest(params)
cat(paste(i,k,' done', '\n',sep='_'))

}
}

save(ME54GO,file='./Analysis/ME54_GO.rda',compress=TRUE)




for(k in GOcat){
TMPnames = names(ME54GO)[grep(paste('_',k,sep=''),names(ME54GO))]

for(i in TMPnames) {
if(i==TMPnames[1]){
tmp = summary(ME54GO[[i]])
tmp$FDR = p.adjust(tmp$Pvalue)
tmp=tmp[,c(1,2,8,3:7)]
tmp=tmp[which(tmp$FDR<=0.05),]

write.xlsx2(tmp, file=paste('./Analysis/ME54_',k,'_GO.xlsx',sep=''), sheetName=i,col.names=TRUE, row.names=FALSE, append=FALSE)
} else {
tmp = summary(ME54GO[[i]])
tmp$FDR = p.adjust(tmp$Pvalue)
tmp=tmp[,c(1,2,8,3:7)]
tmp=tmp[which(tmp$FDR<=0.05),]

write.xlsx2(tmp, file=paste('./Analysis/ME54_',k,'_GO.xlsx',sep=''), sheetName=i,col.names=TRUE, row.names=FALSE, append=TRUE)
}
cat(i)
}
cat(k)
}

### quick look at genes in modules by GO category


tmpUni = geneIdUniverse(ME54GO[['ME65_BP']])[['GO:0044782']] ## cilium organization

tmpGenes = totgene[match(tmpUni,totgene$entrez),]



 write.csv(as.data.frame(MEs$validColors),file='./Analysis/ME54_Module_Assignments.csv')






######################
### Cell Type Enrichments

## Major Class All Cells 

## rank sum? t-test?!? 


Eigen = MEs[['eigengenes']]

CellList = All.integrated@meta.data[row.names(Eigen),'MajorClass']

CellTypes = names(table(CellList))[which(table(CellList)>10)]

MajRes = as.data.frame(matrix(1,nrow=ncol(Eigen),ncol=length(CellTypes)),stringsAsFactors=FALSE)

row.names(MajRes)=colnames(Eigen)

colnames(MajRes) = CellTypes

for(j in colnames(Eigen)){
for(i in CellTypes){

MajRes[j,i] = t.test(x=Eigen[which(CellList==i),j],y=Eigen[,j],alternative="greater")$p.value

}
}

save(MajRes,file='./Analysis/KmeansClust_Enrich_MajorCelltypes.rda',compress=TRUE)

### works, but blunt instrument...

CellList = All.integrated@meta.data[row.names(Eigen),'SubClass']

CellTypes = names(table(CellList))[which(table(CellList)>10)]

SubRes = as.data.frame(matrix(1,nrow=ncol(Eigen),ncol=length(CellTypes)),stringsAsFactors=FALSE)

row.names(SubRes)=colnames(Eigen)

colnames(SubRes) = CellTypes

for(j in colnames(Eigen)){
for(i in CellTypes){

SubRes[j,i] = t.test(x=Eigen[which(CellList==i),j],y=Eigen[,j],alternative="greater")$p.value

}
}

save(SubRes,file='./Analysis/KmeansClust_Enrich_SubCelltypes.rda',compress=TRUE)

## Seth's recommendation was to use FindMarkers 

EigenSO = CreateSeuratObject(counts = t(Eigen))

EigenSO@meta.data = cbind(EigenSO@meta.data, All.integrated@meta.data[colnames(EigenSO),c('sample','MajorClass','SubClass','Neuron','Region','Oligo_pseudotime','Dev_pseudotime')])


### plot corresponding dot plot
EigenSO@active.ident=as.factor(EigenSO$MajorClass)
EigenSO.markers.Major = FindAllMarkers(EigenSO,features=row.names(EigenSO),assay='RNA',logfc.threshold=0)

#TopCellTypeAll = rep(NA,length=)


## 

EigenSO@active.ident=as.factor(EigenSO$SubClass)
EigenSO.markers.Sub = FindAllMarkers(EigenSO,features=row.names(EigenSO),assay='RNA',logfc.threshold=0)

EigenSOregion = subset(EigenSO,cells=colnames(EigenSO)[which(!is.na(EigenSO@meta.data$Region))])

EigenSOregion@active.ident=as.factor(EigenSO$Region)
EigenSO.markers.Region = FindAllMarkers(EigenSOregion,features=row.names(EigenSOregion),assay='RNA',logfc.threshold=0)






### look at 3D of all cells 

All3D.integrated = All.integrated


All3D.integrated <- RunUMAP(All3D.integrated, dims = 1:50, n.components = 3L)

plot.data <- FetchData(object = All3D.integrated, vars = c("UMAP_1", "UMAP_2", "UMAP_3"))

plot.data=plot.data[order(plot.data$UMAP_3,decreasing=TRUE),]

plot.data = cbind(plot.data,All.integrated@meta.data[rownames(plot.data),c('sample','MajorClass','SubClass','Neuron','Region','Oligo_pseudotime','Dev_pseudotime','Phase')],Eigen[rownames(plot.data),])

write.csv(plot.data,file='./Analysis/Allintegrated_3D_UMAP.csv')

tmpTF = t(All.integrated@assays$RNA[TFs,rownames(plot.data)])


#tmpTF = tmpTF[,which(colSums(tmpTF)>10)]

write.csv(All.integrated@assays$RNA[rownames(plot.data),file='./Analysis/Allintegrated_3D_UMAP.csv')

write.csv(tmpTF,file='./Analysis/Allintegrated_TFexp_3D_UMAP.csv')

## on laptop
library(plotly)
library(htmlwidgets)
library("rmarkdown")
# Make a column of row name identities (these will be your cell/barcode names)

plot.data=read.csv("~/Documents/UMB_work/Ament/BRAIN_init/Kriegstein500k/Hypo/Paper/Analysis/Allintegrated_3D_UMAP.csv")

tmpTF = read.csv("~/Documents/UMB_work/Ament/BRAIN_init/Kriegstein500k/Hypo/Paper/Analysis/Allintegrated_TFexp_3D_UMAP.csv")

plot.data$label <- paste(rownames(plot.data))

tmpind = sample(1:nrow(plot.data),1000)

p=plot_ly(data = plot.data[tmpind,], x = ~UMAP_1[tmpind], y = ~UMAP_2[tmpind], z = ~UMAP_3[tmpind], color = ~Phase,type = "scatter3d", mode= "markers",  marker = list(size = 2, width=2))

htmlwidgets::saveWidget(p, file="~/Documents/UMB_work/Ament/BRAIN_init/Kriegstein500k/Hypo/Paper/Figure_Drafts/Phase_All_sampled_3D.html",selfcontained=TRUE)

p=plot_ly(data = plot.data, x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, color = ~Oligo_pseudotime,type = "scatter3d", mode= "markers",  marker = list(size = 2, width=2))

htmlwidgets::saveWidget(p, file="~/Documents/UMB_work/Ament/BRAIN_init/Kriegstein500k/Hypo/Paper/Figure_Drafts/All_Oligo_pseudotime_3D.html",selfcontained=TRUE)




plot_ly(data = plot.data, x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, color = ~sample,type = "scatter3d", mode= "markers",  marker = list(size = 2, width=2))

plot_ly(data = plot.data, x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, color = ~Phase,type = "scatter3d", mode= "markers",  marker = list(size = 2, width=2))

plot_ly(data = plot.data, x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, color = ~MajorClass,type = "scatter3d", mode= "markers",  marker = list(size = 2, width=2))

plot_ly(data = plot.data, x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, color = ~Oligo_pseudotime,type = "scatter3d", mode= "markers",  marker = list(size = 2, width=2))

plot_ly(data = plot.data, x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, color = ~ME2,type = "scatter3d", mode= "markers",  marker = list(size = 2, width=2))

plot_ly(data = plot.data, x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, color = ~Neuron,type = "scatter3d", mode= "markers",  marker = list(size = 2, width=2))

plot.data$NR0B1 = tmpTF$NR0B1

plot_ly(data = plot.data, x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, color = ~NR0B1,type = "scatter3d", mode= "markers",  marker = list(size = 2, width=2))

plot.data$OLIG1 = tmpTF$OLIG1

plot_ly(data = plot.data, x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, color = ~OLIG1,type = "scatter3d", mode= "markers",  marker = list(size = 2, width=2))


#### simple look at top expressed TF in mature population to look back on lineage

## take top cor TF's of a module, then plot avg expression of top 20?!? - extend to GRN's? 



for(i in names(Eigen)[20:54]){
tmpTF = intersect(HSTF$Name,names(MEs$validColors)[MEs$validColors==as.numeric(gsub('ME','',i))])

if(length(tmpTF)>0){
pdf(paste('./TestPlots/All_test_',i,'_TFavgExp.pdf',sep=''),width=12,height=8)
for(k in c("sample","MajorClass")){
print(DimPlot(All.integrated, reduction = "umap", group.by = k, label = TRUE, repel = TRUE))
}
print(FeaturePlot(All.integrated, features = i, pt.size = 0.2, ncol = 2,order=TRUE,cols=c('yellow','red')))
All.integrated@meta.data$meanTF = colMeans(All.integrated@assays$RNA[tmpTF,])
print(FeaturePlot(All.integrated, features = 'meanTF', pt.size = 0.3, ncol = 1,order=TRUE))
dev.off()
}
cat(paste(i, ', ',sep=''))
}


load('./Analysis/All.integrated.MajorClass.DEG.rda')

### what do top TF's look like?  - should just do top few 


for(i in unique(All.integrated.MajorClass.DEG$cluster)){
tmpRes = All.integrated.MajorClass.DEG[which(All.integrated.MajorClass.DEG$cluster==i & All.integrated.MajorClass.DEG$avg_log2FC>0.5 & All.integrated.MajorClass.DEG$p_val_adj<=0.05),]
tmpTF = tmpRes$gene

if(length(tmpTF)>0){
pdf(paste('./TestPlots/All_test_',i,'_TFavgExp.pdf',sep=''),width=12,height=8)
for(k in c("sample","MajorClass")){
print(DimPlot(All.integrated, reduction = "umap", group.by = k, label = TRUE, repel = TRUE))
}
All.integrated@meta.data$meanTF = colMeans(All.integrated@assays$RNA[tmpTF,])
print(FeaturePlot(All.integrated, features = 'meanTF', pt.size = 0.3, ncol = 1,order=TRUE))
for(m in tmpTF){
	print(FeaturePlot(All.integrated, features = m, pt.size = 0.4,order=TRUE))
}

dev.off()
}
cat(paste(i, ', ',sep=''))
}


### coloring by top module? 

#ModLookup = data.frame(ME=c("ME2","ME54","ME3","ME62","ME100","ME65","ME80","ME44","ME9","ME47","ME37","ME11","ME12","ME43","ME76","ME92"),Celltype=c("Oligodendrocyte","Oligodendrocyte","Astrocyte","Astrocyte","Astrocyte","Ependymal","Ependymal","Ependymal","Ependymal","Ependymal","Ependymal","ImmatureNeuron","MatureNeuron","MatureNeuron","MatureNeuron","MatureNeuron"),stringsAsFactors=FALSE)


ModLookup = data.frame(ME=c("ME2","ME54","ME3","ME62","ME80","ME9","ME11","ME43","ME92"),Celltype=c("Oligodendrocyte","Oligodendrocyte","Astrocyte","Astrocyte","Ependymal","Ependymal","ImmatureNeuron","MatureNeuron","MatureNeuron"),stringsAsFactors=FALSE)

TopEigen = colnames(Eigen[,ModLookup$ME])[apply(Eigen[,ModLookup$ME],1,which.max)]

TopCT = ModLookup$Celltype[match(TopEigen,ModLookup$ME)]

names(TopCT) = rownames(Eigen)

MaxVal = apply(Eigen[,ModLookup$ME],1,max)

TopCT[which(MaxVal<0.005 )] = NA #

All.integrated@meta.data$ColorByMod = TopCT[colnames(All.integrated)]

All.integrated@meta.data$MaxMod = MaxVal[colnames(All.integrated)]



pdf(paste('./TestPlots/All_test_Color_By_Module_2per.pdf',sep=''),width=12,height=8)
for(k in c("sample","MajorClass","ColorByMod")){
print(DimPlot(All.integrated, reduction = "umap", group.by = k, label = TRUE, repel = TRUE))
}
print(FeaturePlot(All.integrated, features = 'MaxMod', pt.size = 0.4,order=TRUE))
dev.off()


#All.integrated@meta.data$ColorByMod[which(is.na(All.integrated@meta.data$ColorByMod))] = 'None'

#subplots(All.integrated,colName='ColorByMod',filename='./TestPlots/Color_By_Modules_IndvPlots.pdf')
### id cells by TF cor? 


## subplots with proper colors? 

ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}


color_list <- ggplotColours(n=5)

TestLin = sort(unique(All.integrated@meta.data$ColorByMod))

TmpObj1 = All.integrated[,which(!is.na(match(All.integrated@meta.data$MajorClass,c('Astrocyte','Dividing_Progenitor','Neurons','Oligodendrocyte','RadialGlia','Tanycytes'))))]

XLIM = range(Embeddings(All.integrated[['umap']])[,1])
YLIM = range(Embeddings(All.integrated[['umap']])[,2])
## just cells of one type 

pdf('./TestPlots/Color_By_Modules_IndvPlots_Isolated.pdf',width=12,height=8)
TmpObj2 = TmpObj1[,which(!is.na(TmpObj1@meta.data$ColorByMod))]

print(DimPlot(TmpObj2, reduction = "umap", group.by = 'ColorByMod', cols=color_list, label = TRUE, repel = TRUE) + xlim(XLIM) +ylim(YLIM))
for(j in 1:length(TestLin)){
TmpObj3 = TmpObj1[,which(TmpObj1@meta.data$ColorByMod==TestLin[j])]

print(DimPlot(TmpObj3, reduction = "umap", group.by = 'ColorByMod', cols=color_list[j], label = TRUE, repel = TRUE)  + xlim(XLIM) +ylim(YLIM) )

#print(FeaturePlot(TmpObj3, features = 'MaxMod', pt.size = 0.2,order=TRUE)  + xlim(XLIM) +ylim(YLIM) )

}
dev.off()

### marker genes per group 

DefaultAssay(TmpObj2) = 'RNA'

TmpObj2@active.ident=as.factor(TmpObj2$ColorByMod)

TmpObj2.markers.Major = FindAllMarkers(TmpObj2)


table(All.integrated@meta.data$sample[which(All.integrated@meta.data$ColorByMod=='ImmatureNeuron')])/table(All.integrated@meta.data$sample)


### just check, location of NP expression 

pdf(paste('./TestPlots/All_NP_expression.pdf',sep=''),width=12,height=8)
for(k in c("sample","MajorClass","ColorByMod")){
print(DimPlot(All.integrated, reduction = "umap", group.by = k, label = TRUE, repel = TRUE))
}
for(i in NPlist2){
print(FeaturePlot(All.integrated, features = i, pt.size = 0.4,order=TRUE))	
}

dev.off()


## heatmap? 

## also need to eliminate NP's expressed in other cell types - DEG? but also capture progenitors 







##################
### Enrichments with Tri2 Genie 

Tri2Genie7425knn5 = read.csv('./Networks/GENIE_Tri2_7425knn5__100G_16T_seed123.csv')

Tri2Genie7425knn5top = Tri2Genie7425knn5[1:500000,]

MEgenes = names(MEs[['validColors']])[which(MEs[['validColors']]!=0)] #7637

Tri2Genie7425knn5top = Tri2Genie7425knn5top[which(!is.na(match(Tri2Genie7425knn5top$TF,MEgenes))),]

Tri2Genie7425knn5top = Tri2Genie7425knn5top[which(!is.na(match(Tri2Genie7425knn5top$target,MEgenes))),]

Modules = MEs[['validColors']][which(MEs[['validColors']]!=0)]

Modules = Modules[which(!is.na(match(names(Modules),unique(c(Tri2Genie7425knn5top$TF,Tri2Genie7425knn5top$target)))))] #5302

kclusts = sort(unique(Modules)) # still 54

GenieTF = names(table(Tri2Genie7425knn5top$TF))[which(as.numeric(table(Tri2Genie7425knn5top$TF))>20)] ## mean   - 645 TF
 
#fisher.test(matrix(c(A,B,C,D),nrow=2,ncol=2),alternative="greater")

#A = cells expressing and in cluster
#B = cells not expressing but in cluster
#C = cells expressing but not in cluster
#D = cells not expressing and not in cluster

GKperOv = data.frame(matrix(0,nrow=length(GenieTF),ncol=length(kclusts)*6))

### recording - fisher p val, fisher FDR, Per expressing, Avg exp 

rownames(GKperOv) = GenieTF
colnames(GKperOv) = paste(rep(kclusts,each=6),c('FT_pval','FT_FDR','PerClust','TotClust','PerMod','TotMod'),sep='_')
 
 for(i in GenieTF){
 	Gmod = as.character(Tri2Genie7425knn5top$target[which(Tri2Genie7425knn5top$TF==i)])
 for(j in kclusts){
clust=names(Modules)[Modules==j]
AA=length(intersect(Gmod,clust))
BB=length(clust)-AA
CC=length(Gmod)-AA
DD = 5302-AA-BB-CC
GKperOv[i,paste(j,'FT_pval',sep='_')]=fisher.test(matrix(c(AA,BB,CC,DD),nrow=2,ncol=2),alternative="greater")$p.value
GKperOv[i,paste(j,'PerClust',sep='_')]=round((AA/length(clust))*100,2)
GKperOv[i,paste(j,'PerMod',sep='_')]=round((AA/length(Gmod))*100,2)
GKperOv[i,paste(j,'TotClust',sep='_')]=length(clust)
GKperOv[i,paste(j,'TotMod',sep='_')]=length(Gmod)
 	}
 	cat(paste(i,'; ',which(GenieTF==i),' of 1170 done, ',sep=''))
 }

for(j in kclusts){
GKperOv[,paste(j,'FT_FDR',sep='_')] = p.adjust(GKperOv[,paste(j,'FT_pval',sep='_')])
}

ModperOvSpearknn5 = GKperOv

colnames(ModperOvSpearknn5) = paste('ME',colnames(ModperOvSpearknn5),sep='')

save(ModperOvSpearknn5,file='./Analysis/Tri2_GKperOvSpearknn5_GlobalModules.rda',compress=TRUE) 

write.csv(ModperOvSpearknn5,file='./Analysis/Tri2_GKperOvSpearknn5_GlobalModules.csv')

#### plot kmeans average and genie GRN modules 

### All cells (Tri1 + Tri2) - Plot Modules and Tri2 GRNs 


### calculate AUC for GENIE3 results for Trimester 2

TFs = as.character(unique(Tri2Genie7425knn5top$TF))

Targets = unique(Tri2Genie7425knn5top$target)

DefaultAssay(All.integrated) = 'RNA'
All.integrated.ExpRank <- AUCell_buildRankings(All.integrated@assays$RNA[Targets])

Cells = colnames(All.integrated)

AllGenieK5AUC = matrix(0,nrow=length(Cells),ncol=length(TFs))

rownames(AllGenieK5AUC) = Cells
colnames(AllGenieK5AUC) = TFs

for(i in TFs){

tmpTarget = Tri2Genie7425knn5top$target[which(Tri2Genie7425knn5top$TF==i)]

AllGenieK5AUC[,i] <- as.numeric(getAUC(AUCell_calcAUC(as.character(tmpTarget), All.integrated.ExpRank)))
cat(paste(i,', ',sep=''))
}

save(AllGenieK5AUC,file='./Analysis/AllGenieK5AUC.rda',compress=TRUE)	


### plotting 

for(n in colnames(MEs[['eigengenes']])){

All.integrated@meta.data[,n] <- MEs[['eigengenes']][colnames(All.integrated),n]
cat(paste(n,', ',sep=''))
} 


DefaultAssay(All.integrated) = 'RNA'


for(clust in colnames(MEs[['eigengenes']])){

TFs = rownames(ModperOvSpearknn5[order(ModperOvSpearknn5[,paste(clust,'_FT_FDR',sep='')],decreasing=FALSE)[1:10],])

pdf(paste('./TestPlots/All_GRN_Spear_Clust_',clust,'_GlobalModules.pdf',sep=''),width=12,height=8)
#print(DimPlot(All.integrated, label = TRUE) + NoLegend())
for(k in c("sample","MajorClass","SubClass","Neuron","Region")){
print(DimPlot(All.integrated, reduction = "umap", group.by = k, label = TRUE, repel = TRUE))
}
print(FeaturePlot(All.integrated, features = clust, pt.size = 0.2, ncol = 2,order=TRUE,cols=c('yellow','red')))
for(i in TFs){
print(FeaturePlot(All.integrated, features = i, pt.size = 0.3, ncol = 1,order=TRUE))
All.integrated@meta.data$tmp = AllGenieK5AUC[,i]
print(FeaturePlot(All.integrated, features = 'tmp', pt.size = 0.3, ncol = 1,order=TRUE,cols=c('yellow','red')))
}
dev.off()
cat(clust)
}





### calculate AUC for GENIE3 results for Trimester 2

TFs = as.character(unique(Tri2Genie7425knn5top$TF))

Targets = unique(Tri2Genie7425knn5top$target)

DefaultAssay(AllTri2.integrated) = 'RNA'
AllTri2.integrated.ExpRank <- AUCell_buildRankings(AllTri2.integrated@assays$RNA[Targets])

Cells = colnames(AllTri2.integrated)

GenieK5AUC = matrix(0,nrow=length(Cells),ncol=length(TFs))

rownames(GenieK5AUC) = Cells
colnames(GenieK5AUC) = TFs

for(i in TFs){

tmpTarget = Tri2Genie7425knn5top$target[which(Tri2Genie7425knn5top$TF==i)]

GenieK5AUC[,i] <- as.numeric(getAUC(AUCell_calcAUC(as.character(tmpTarget), AllTri2.integrated.ExpRank)))
cat(paste(i,', ',sep=''))
}

save(GenieK5AUC,file='./Analysis/GenieK5AUC.rda',compress=TRUE)







### plotting in Tri2 first

for(n in colnames(MEs[['eigengenes']])){

AllTri2.integrated@meta.data[,n] <- MEs[['eigengenes']][colnames(AllTri2.integrated),n]
cat(paste(n,', ',sep=''))
} 


DefaultAssay(AllTri2.integrated) = 'RNA'


for(clust in colnames(MEs[['eigengenes']])){

TFs = rownames(ModperOvSpearknn5[order(ModperOvSpearknn5[,paste(clust,'_FT_FDR',sep='')],decreasing=FALSE)[1:10],])

pdf(paste('./TestPlots/AllTri2_GRN_Spear_Clust_',clust,'_GlobalModules.pdf',sep=''),width=12,height=8)
#print(DimPlot(AllTri2.integrated, label = TRUE) + NoLegend())
for(k in c("sample","MajorClass","SubClass","Neuron","Region")){
print(DimPlot(AllTri2.integrated, reduction = "umap", group.by = k, label = TRUE, repel = TRUE))
}
print(FeaturePlot(AllTri2.integrated, features = clust, pt.size = 0.2, ncol = 2,order=TRUE,cols=c('yellow','red')))
for(i in TFs){
print(FeaturePlot(AllTri2.integrated, features = i, pt.size = 0.3, ncol = 1,order=TRUE))
AllTri2.integrated@meta.data$tmp = GenieK5AUC[,i]
print(FeaturePlot(AllTri2.integrated, features = 'tmp', pt.size = 0.3, ncol = 1,order=TRUE,cols=c('yellow','red')))
}
dev.off()
cat(clust)
}


EigenSOTri2 = subset(EigenSO,cells=intersect(colnames(EigenSO),colnames(AllTri2.integrated)))


EigenSOTri2.markers.Major = FindAllMarkers(EigenSOTri2,features=row.names(EigenSOTri2),assay='RNA',logfc.threshold=0)

EigenSOTri2@active.ident=as.factor(EigenSOTri2$SubClass)

EigenSOTri2.markers.Sub = FindAllMarkers(EigenSOTri2,features=row.names(EigenSOTri2),assay='RNA',logfc.threshold=0)

## confirms what I see in plots 

### Neurons 


### calculate AUC for GENIE3 results for Trimester 2

TFs = as.character(unique(Tri2Genie7425knn5top$TF))

Targets = unique(Tri2Genie7425knn5top$target)

DefaultAssay(Tri2Neuron_ds.integrated) = 'RNA'
Tri2Neuron_ds.integrated.ExpRank <- AUCell_buildRankings(Tri2Neuron_ds.integrated@assays$RNA[Targets])

Cells = colnames(Tri2Neuron_ds.integrated)

NeuroGenieK5AUC = matrix(0,nrow=length(Cells),ncol=length(TFs))

rownames(NeuroGenieK5AUC) = Cells
colnames(NeuroGenieK5AUC) = TFs

for(i in TFs){

tmpTarget = Tri2Genie7425knn5top$target[which(Tri2Genie7425knn5top$TF==i)]

NeuroGenieK5AUC[,i] <- as.numeric(getAUC(AUCell_calcAUC(as.character(tmpTarget), Tri2Neuron_ds.integrated.ExpRank)))
cat(paste(i,', ',sep=''))
}

save(NeuroGenieK5AUC,file='./Analysis/NeuroGenieK5AUC.rda',compress=TRUE)











CellList = Tri2Neuron_ds.integrated@meta.data[row.names(Eigen),'Neuron']

CellTypes = names(table(CellList))[which(table(CellList)>5)]

NeuroRes = as.data.frame(matrix(1,nrow=ncol(Eigen),ncol=length(CellTypes)),stringsAsFactors=FALSE)

row.names(NeuroRes)=colnames(Eigen)

colnames(NeuroRes) = CellTypes

for(j in colnames(Eigen)){
for(i in CellTypes){

NeuroRes[j,i] = t.test(x=Eigen[which(CellList==i),j],y=Eigen[,j],alternative="greater")$p.value

}
}

save(NeuronRes,file='./Analysis/KmeansClust_Enrich_NeuronSubtypes.rda',compress=TRUE)




for(n in colnames(MEs[['eigengenes']])){

Tri2Neuron_ds.integrated@meta.data[,n] <- MEs[['eigengenes']][colnames(Tri2Neuron_ds.integrated),n]
cat(paste(n,', ',sep=''))
} 

DefaultAssay(Tri2Neuron_ds.integrated) = 'RNA'





pdf(paste('./TestPlots/Tri2Neuron_ds_All_GlobalModules.pdf',sep=''),width=12,height=8)
for(k in c("sample","MajorClass","SubClass","Neuron","Region")){
print(DimPlot(Tri2Neuron_ds.integrated, reduction = "umap", group.by = k, label = TRUE, repel = TRUE))
}
for(clust in colnames(MEs[['eigengenes']])){
print(FeaturePlot(Tri2Neuron_ds.integrated, features = clust, pt.size = 0.2, ncol = 2,order=TRUE,cols=c('yellow','red')))
}
dev.off()










for(clust in colnames(MEs[['eigengenes']])){

TFs = rownames(ModperOvSpearknn5[order(ModperOvSpearknn5[,paste(clust,'_FT_FDR',sep='')],decreasing=FALSE)[1:10],])

pdf(paste('./TestPlots/Tri2Neuron_ds_GRN_Spear_Clust_',clust,'_GlobalModules.pdf',sep=''),width=12,height=8)
print(DimPlot(Tri2Neuron_ds.integrated, label = TRUE) + NoLegend())
for(k in c("sample","MajorClass","SubClass","Neuron","Region")){
print(DimPlot(Tri2Neuron_ds.integrated, reduction = "umap", group.by = k, label = TRUE, repel = TRUE))
}
print(FeaturePlot(Tri2Neuron_ds.integrated, features = clust, pt.size = 0.2, ncol = 2,order=TRUE,cols=c('yellow','red')))
for(i in TFs){
print(FeaturePlot(Tri2Neuron_ds.integrated, features = i, pt.size = 0.3, ncol = 1,order=TRUE))
Tri2Neuron_ds.integrated@meta.data$tmp = GenieK5AUC[colnames(Tri2Neuron_ds.integrated),i]
print(FeaturePlot(Tri2Neuron_ds.integrated, features = 'tmp', pt.size = 0.3, ncol = 1,order=TRUE,cols=c('yellow','red')))
}
dev.off()
cat(clust)
}



EigenSONeuro = subset(EigenSO,cells=intersect(colnames(EigenSO),colnames(Tri2Neuron_ds.integrated)))

EigenSONeuro@active.ident=as.factor(EigenSONeuro$MajorClass)

EigenSONeuro.markers.Major = FindAllMarkers(EigenSONeuro,features=row.names(EigenSONeuro),assay='RNA',logfc.threshold=0)

EigenSONeuro@active.ident=as.factor(EigenSONeuro$SubClass)

EigenSONeuro.markers.Sub = FindAllMarkers(EigenSONeuro,features=row.names(EigenSONeuro),assay='RNA',logfc.threshold=0)

EigenSONeuro@active.ident=as.factor(EigenSONeuro$Region)

EigenSONeuro.markers.Region = FindAllMarkers(EigenSONeuro,features=row.names(EigenSONeuro),assay='RNA',logfc.threshold=0)

EigenSONeuro@active.ident=as.factor(EigenSONeuro$Neuron)

EigenSONeuro.markers.Neuron = FindAllMarkers(EigenSONeuro,features=row.names(EigenSONeuro),assay='RNA',logfc.threshold=0)




