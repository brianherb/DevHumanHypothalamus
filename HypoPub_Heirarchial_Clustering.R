library(Seurat)

setwd('/local/projects-t3/idea/bherb/Hypothalamus/PubRes')

load('./SeuratObj/Tri2Neuron_ds.integrated.rda') ##Tri2Neuron_ds.integrated

HSTF = read.csv('/local/projects-t3/idea/bherb/annotation/Hsapiens/Human_TF.csv')

NPlist = read.csv('/local/projects-t3/idea/bherb/Hypothalamus/Neuropeptide_list.csv',header=FALSE) ## recieved from Hannah on 7/15/20

NPlist = MMtoHS(as.character(NPlist[,1]))

NPlist2 = unique(c(na.omit(as.character(TF_NP$Neuropeptides)),NPlist))


## Neuron is subtypes and Region is Nuclei

### calculate heirarchical clustering and top genes for whole nuclei
DefaultAssay(Tri2Neuron_ds.integrated) = 'integrated' 
Tri2Neuron_ds.integrated@active.ident=as.factor(Tri2Neuron_ds.integrated$Region)
Tri2Neuron_ds.integrated<- FindNeighbors(Tri2Neuron_ds.integrated, dims = 1:10, verbose = FALSE)
Tri2Neuron_ds.integrated= BuildClusterTree(object=Tri2Neuron_ds.integrated,dims=1:10)
nucleitree=Tool(object = Tri2Neuron_ds.integrated, slot = 'BuildClusterTree')

save(nucleitree,file='./Analysis/Neuron_WholeNuclei_Tree_Feb2021.rda')

pdf('./TestPlots/Tri2Neuron_ds_HGnuclei_Feb2021_10PCs_NODELABELS_dendrogram.pdf',width=12,height=8)
plot(nucleitree, edge.width = 2)
ape::nodelabels()
dev.off()


### bifurcation 

DefaultAssay(Tri2Neuron_ds.integrated) = 'RNA' 

Nodes = as.character(sort(unique(nucleitree$edge[,1])))

NodeTFdiff = vector(mode='list',length=nucleitree$Nnode)
names(NodeTFdiff) = Nodes

NodeNPdiff = vector(mode='list',length=nucleitree$Nnode)
names(NodeNPdiff) = Nodes

NodeclustsNuclei = vector(mode='list',length=nucleitree$Nnode)
names(NodeclustsNuclei) = Nodes

for(i in Nodes){

tmpNodes = nucleitree$edge[which(nucleitree$edge[,1]==as.numeric(i)),2]


if(tmpNodes[1]>length(nucleitree$tip.label)){
clusts1 = nucleitree$tip.label[getTips(nucleitree,tmpNodes[1])]
} else {
clusts1 = nucleitree$tip.label[tmpNodes[1]]
}

if(tmpNodes[2]>length(nucleitree$tip.label)){
clusts2 = nucleitree$tip.label[getTips(nucleitree,tmpNodes[2])]
} else {
clusts2 = nucleitree$tip.label[tmpNodes[2]]
}

cells1 = colnames(Tri2Neuron_ds.integrated)[which(!is.na(match(Tri2Neuron_ds.integrated@meta.data$Region,clusts1)))] 
cells2 = colnames(Tri2Neuron_ds.integrated)[which(!is.na(match(Tri2Neuron_ds.integrated@meta.data$Region,clusts2)))] 

Tri2Neuron_ds.integrated@meta.data$tmpGroups = NA
Tri2Neuron_ds.integrated@meta.data[cells1,'tmpGroups'] = 'A'
Tri2Neuron_ds.integrated@meta.data[cells2,'tmpGroups'] = 'B'
Tri2Neuron_ds.integrated@active.ident=as.factor(Tri2Neuron_ds.integrated$tmpGroups)

NodeTFdiff[[i]] = try(FindMarkers(Tri2Neuron_ds.integrated,ident.1='A',ident.2='B',features=intersect(rownames(Tri2Neuron_ds.integrated),HSTF$Name)))

NodeNPdiff[[i]] = try(FindMarkers(Tri2Neuron_ds.integrated,ident.1='A',ident.2='B',features=intersect(rownames(Tri2Neuron_ds.integrated),NPlist2)))

tmp=list(clusts1,clusts2)
names(tmp) = tmpNodes

NodeclustsNuclei[[i]] = tmp

cat('\n\n')
cat(paste(tmpNodes,collapse='_VS_'))
cat('\n\n')

}

NodeTFdiffnuclei = NodeTFdiff
NodeNPdiffnuclei = NodeNPdiff

save(NodeTFdiffnuclei,file='./Analysis/NodeTFdiff_nuclei_10PCs_Feb2021.rda',compress=TRUE)
save(NodeNPdiffnuclei,file='./Analysis/NodeNPdiff_nuclei_10PCs_Feb2021.rda',compress=TRUE)


## plot gene names 
pdf('./TestPlots/Tri2Neuron_ds_HGnuclei_Feb2021_10PCs_TFgenes_dendrogram.pdf',width=12,height=12)
plot(nucleitree, edge.width = 1)
ape::nodelabels()
#ape::edgelabels()
for(i in Nodes){

tmpNodes = nucleitree$edge[which(nucleitree$edge[,1]==as.numeric(i)),2]

edge1 = which(nucleitree$edge[,1]==i & nucleitree$edge[,2]==tmpNodes[1])
edge2 = which(nucleitree$edge[,1]==i & nucleitree$edge[,2]==tmpNodes[2])

tmpDEG = NodeTFdiff[[i]]

sigEdge1 = tmpDEG[which(tmpDEG$avg_log2FC>0.5 & tmpDEG$p_val_adj<0.000000001),]
allEdge1 = na.omit(rownames(sigEdge1[order(sigEdge1$avg_log2FC,decreasing=TRUE),])[1:5])

if(length(allEdge1)>3){
edge1genes = paste(paste(allEdge1[1:3],collapse=','),'\n',paste(allEdge1[4:length(allEdge1)],collapse=','))
} else {
edge1genes = paste(allEdge1,collapse=',')
}

sigEdge2 = tmpDEG[which(tmpDEG$avg_log2FC<(-0.5) & tmpDEG$p_val_adj<0.000000001),]
allEdge2 = na.omit(rownames(sigEdge2[order(sigEdge2$avg_log2FC,decreasing=FALSE),])[1:5])

if(length(allEdge2)>3){
edge2genes = paste(paste(allEdge2[1:3],collapse=','),'\n',paste(allEdge2[4:length(allEdge2)],collapse=','))
} else {
edge2genes = paste(allEdge2,collapse=',')
}

ape::edgelabels(text=edge1genes,edge=edge1,adj=c(1,-1),frame='none',bg='white',cex=0.5)
ape::edgelabels(text=edge2genes,edge=edge2,adj=c(1,-1),frame='none',bg='white',cex=0.5)
}
dev.off()


### plot corresponding dot plot
Tri2Neuron_ds.integrated@active.ident=as.factor(Tri2Neuron_ds.integrated$Region)
Tri2Neuron_ds.markers.nuclei.HSTFs = FindAllMarkers(Tri2Neuron_ds.integrated,features=intersect(rownames(Tri2Neuron_ds.integrated),HSTF$Name))

write.csv(Tri2Neuron_ds.markers.nuclei.HSTFs,file='./Analysis/Tri2Neuron_ds.markers_RNA_nuclei_HSTFs_Feb2021.csv')

## evaluate 
ClustOrder = rev(unique(as.character(Tri2Neuron_ds.integrated@active.ident))[c(1,10,2,3,15,6,11,4,5,8,9,12,7,13,14)])
Tri2Neuron_ds.integrated@active.ident=factor(Tri2Neuron_ds.integrated@active.ident,levels=ClustOrder)
Tri2Neuron_ds.markers.nuclei.HSTFs$perDif=Tri2Neuron_ds.markers.nuclei.HSTFs$pct.1-Tri2Neuron_ds.markers.nuclei.HSTFs$pct.2

#MajTypes = paste('clust_',seq(42,0,-1),sep='')
MajTypes = ClustOrder

for(i in MajTypes){
if(i==MajTypes[1]){
    tmp = Tri2Neuron_ds.markers.nuclei.HSTFs[which(Tri2Neuron_ds.markers.nuclei.HSTFs$cluster==i ),]
    tmp = tmp[which(tmp$p_val_adj<=0.001 & tmp$pct.2<0.2),]
    tmp = tmp[order(tmp$perDif,decreasing=TRUE),]
    #tmp = tmp[which(tmp$avg_logFC>1),]
topPerDiftf <-  tmp$gene[1:3]
} else {
    tmp = Tri2Neuron_ds.markers.nuclei.HSTFs[which(Tri2Neuron_ds.markers.nuclei.HSTFs$cluster==i ),]
    tmp = tmp[which(tmp$p_val_adj<=0.001 & tmp$pct.2<0.2),]
    tmp = tmp[order(tmp$perDif,decreasing=TRUE),]
    #tmp = tmp[which(tmp$avg_logFC>1),]
topPerDiftf <-  c(topPerDiftf,tmp$gene[1:3])
}
}

topPerDiftf = na.omit(topPerDiftf)
topPerDiftf = unique(topPerDiftf)
topPerDiftf = rev(topPerDiftf)

pdf('./TestPlots/Tri2Neuron_ds_Nuclei_DotPlot_TopTF_RNA_Dend_Order_Feb2021.pdf',width=12,height=12)

DotPlot(Tri2Neuron_ds.integrated,features=topPerDiftf,assay='RNA') + RotatedAxis()
dev.off()


### Neuropeptide tree 



## plot gene names 
pdf('./TestPlots/Tri2Neuron_ds_HGnuclei_Feb2021_10PCs_NPgenes_dendrogram.pdf',width=12,height=12)
plot(nucleitree, edge.width = 1)
ape::nodelabels()
#ape::edgelabels()
for(i in Nodes){

tmpNodes = nucleitree$edge[which(nucleitree$edge[,1]==as.numeric(i)),2]

edge1 = which(nucleitree$edge[,1]==i & nucleitree$edge[,2]==tmpNodes[1])
edge2 = which(nucleitree$edge[,1]==i & nucleitree$edge[,2]==tmpNodes[2])

tmpDEG = NodeNPdiff[[i]]

sigEdge1 = tmpDEG[which(tmpDEG$avg_log2FC>0.5 & tmpDEG$p_val_adj<0.01),]
allEdge1 = na.omit(rownames(sigEdge1[order(sigEdge1$avg_log2FC,decreasing=TRUE),])[1:5])

if(length(allEdge1)>3){
edge1genes = paste(paste(allEdge1[1:3],collapse=','),'\n',paste(allEdge1[4:length(allEdge1)],collapse=','))
} else {
edge1genes = paste(allEdge1,collapse=',')
}

sigEdge2 = tmpDEG[which(tmpDEG$avg_log2FC<(-0.5) & tmpDEG$p_val_adj<0.01),]
allEdge2 = na.omit(rownames(sigEdge2[order(sigEdge2$avg_log2FC,decreasing=FALSE),])[1:5])

if(length(allEdge2)>3){
edge2genes = paste(paste(allEdge2[1:3],collapse=','),'\n',paste(allEdge2[4:length(allEdge2)],collapse=','))
} else {
edge2genes = paste(allEdge2,collapse=',')
}

ape::edgelabels(text=edge1genes,edge=edge1,adj=c(1,-1),frame='none',bg='white',cex=0.5)
ape::edgelabels(text=edge2genes,edge=edge2,adj=c(1,-1),frame='none',bg='white',cex=0.5)
}
dev.off()


### plot corresponding dot plot
Tri2Neuron_ds.integrated@active.ident=as.factor(Tri2Neuron_ds.integrated$Region)
Tri2Neuron_ds.markers.nuclei.HSNPs = FindAllMarkers(Tri2Neuron_ds.integrated,features=intersect(rownames(Tri2Neuron_ds.integrated),NPlist2))

write.csv(Tri2Neuron_ds.markers.nuclei.HSNPs,file='./Analysis/Tri2Neuron_ds.markers_RNA_nuclei_HSNPs_Feb2021.csv')

## evaluate 
ClustOrder = rev(unique(as.character(Tri2Neuron_ds.integrated@active.ident))[c(1,10,2,3,15,6,11,4,5,8,9,12,7,13,14)])
Tri2Neuron_ds.integrated@active.ident=factor(Tri2Neuron_ds.integrated@active.ident,levels=ClustOrder)
Tri2Neuron_ds.markers.nuclei.HSNPs$perDif=Tri2Neuron_ds.markers.nuclei.HSNPs$pct.1-Tri2Neuron_ds.markers.nuclei.HSNPs$pct.2

#MajTypes = paste('clust_',seq(42,0,-1),sep='')
MajTypes = ClustOrder

for(i in MajTypes){
if(i==MajTypes[1]){
    tmp = Tri2Neuron_ds.markers.nuclei.HSNPs[which(Tri2Neuron_ds.markers.nuclei.HSNPs$cluster==i ),]
    tmp = tmp[which(tmp$p_val_adj<=0.001 & tmp$pct.2<0.2),]
    tmp = tmp[order(tmp$perDif,decreasing=TRUE),]
    #tmp = tmp[which(tmp$avg_logFC>1),]
topPerDifNP <-  tmp$gene[1:3]
} else {
    tmp = Tri2Neuron_ds.markers.nuclei.HSNPs[which(Tri2Neuron_ds.markers.nuclei.HSNPs$cluster==i ),]
    tmp = tmp[which(tmp$p_val_adj<=0.001 & tmp$pct.2<0.2),]
    tmp = tmp[order(tmp$perDif,decreasing=TRUE),]
    #tmp = tmp[which(tmp$avg_logFC>1),]
topPerDifNP <-  c(topPerDifNP,tmp$gene[1:3])
}
}

topPerDifNP = na.omit(topPerDifNP)
topPerDifNP = unique(topPerDifNP)
topPerDifNP = rev(topPerDifNP)

pdf('./TestPlots/Tri2Neuron_ds_Nuclei_DotPlot_TopNP_RNA_Dend_Order_Feb2021.pdf',width=12,height=12)

DotPlot(Tri2Neuron_ds.integrated,features=topPerDifNP,assay='RNA') + RotatedAxis()
dev.off()

#### 

### Modules - color by eigengene? Color node and tip? 

### aggregate average eigengene per node / gene 

MEs = readRDS('./Analysis/AllSamp_mergeColors.k100.MEs_CoreCellGene.rds') 

Eigen = MEs[['eigengenes']]

Nuclei = nucleitree$tip.label

NucleiEigenAvg = matrix(0,nrow=ncol(Eigen),ncol=length(Nuclei))

colnames(NucleiEigenAvg) = Nuclei
rownames(NucleiEigenAvg) = colnames(Eigen)

for(i in Nuclei){

tmpEigen = Eigen[colnames(Tri2Neuron_ds.integrated)[which(Tri2Neuron_ds.integrated@meta.data$Region==i)],]

NucleiEigenAvg[,i] = colMeans(tmpEigen,na.rm=TRUE)

}


## base colors on min and max values of NucleiEigenAvg


ii <- cut(as.vector(NucleiEigenAvg), breaks = seq(min(NucleiEigenAvg), max(NucleiEigenAvg), len = 100), include.lowest = TRUE)
## Use bin indices, ii, to select color from vector of n-1 equally spaced colors
colors <- colorRampPalette(c("blue", "red"))(99)[ii]

## This call then also produces the plot below
image(seq_along(values), 1, as.matrix(seq_along(values)), col = colors,
      axes = F)



Module = 'ME12'

### modify code chunk below 
## plot gene names 
pdf('./TestPlots/Tri2Neuron_ds_HGnuclei_Feb2021_10PCs_ME12_dendrogram.pdf',width=12,height=12)
plot(nucleitree, edge.width = 1)
ape::nodelabels()
#ape::edgelabels()
for(i in Nodes){

tmpNodes = nucleitree$edge[which(nucleitree$edge[,1]==as.numeric(i)),2]

edge1 = which(nucleitree$edge[,1]==i & nucleitree$edge[,2]==tmpNodes[1])
edge2 = which(nucleitree$edge[,1]==i & nucleitree$edge[,2]==tmpNodes[2])

tmpDEG = NodeNPdiff[[i]]

sigEdge1 = tmpDEG[which(tmpDEG$avg_log2FC>0.5 & tmpDEG$p_val_adj<0.01),]
allEdge1 = na.omit(rownames(sigEdge1[order(sigEdge1$avg_log2FC,decreasing=TRUE),])[1:5])

if(length(allEdge1)>3){
edge1genes = paste(paste(allEdge1[1:3],collapse=','),'\n',paste(allEdge1[4:length(allEdge1)],collapse=','))
} else {
edge1genes = paste(allEdge1,collapse=',')
}

sigEdge2 = tmpDEG[which(tmpDEG$avg_log2FC<(-0.5) & tmpDEG$p_val_adj<0.01),]
allEdge2 = na.omit(rownames(sigEdge2[order(sigEdge2$avg_log2FC,decreasing=FALSE),])[1:5])

if(length(allEdge2)>3){
edge2genes = paste(paste(allEdge2[1:3],collapse=','),'\n',paste(allEdge2[4:length(allEdge2)],collapse=','))
} else {
edge2genes = paste(allEdge2,collapse=',')
}

ape::edgelabels(text=edge1genes,edge=edge1,adj=c(1,-1),frame='none',bg='white',cex=0.5)
ape::edgelabels(text=edge2genes,edge=edge2,adj=c(1,-1),frame='none',bg='white',cex=0.5)
}
dev.off()












if(tmpNodes[1]>length(nucleitree$tip.label)){
clusts1 = nucleitree$tip.label[getTips(nucleitree,tmpNodes[1])]
} else {
clusts1 = nucleitree$tip.label[tmpNodes[1]]
}

if(tmpNodes[2]>length(nucleitree$tip.label)){
clusts2 = nucleitree$tip.label[getTips(nucleitree,tmpNodes[2])]
} else {
clusts2 = nucleitree$tip.label[tmpNodes[2]]
}

cells1 = colnames(Tri2Neuron_ds.integrated)[which(!is.na(match(Tri2Neuron_ds.integrated@meta.data$Region,clusts1)))] 
cells2 = colnames(Tri2Neuron_ds.integrated)[which(!is.na(match(Tri2Neuron_ds.integrated@meta.data$Region,clusts2)))] 







### non neuronal 

load('./SeuratObj/AllSamples_integrated.rda') #All.integrated

### calculate heirarchical clustering and top genes for whole nuclei

NonNeu.integrated = subset(All.integrated,cells=colnames(All.integrated)[which(!is.na(match(All.integrated@meta.data$MajorClass,c('Astrocyte','Ependymal','Oligodendrocyte','RadialGlia','Tanycytes'))))])                                 

NonNeu.integrated = subset(NonNeu.integrated,cells=colnames(NonNeu.integrated)[which(NonNeu.integrated@meta.data$SubClass!='EarlyMicro')])

NonNeu.integrated = subset(NonNeu.integrated,cells=colnames(NonNeu.integrated)[which(NonNeu.integrated@meta.data$SubClass!='MigratingNeuro1')])



DefaultAssay(NonNeu.integrated) = 'integrated' 
NonNeu.integrated@active.ident=as.factor(NonNeu.integrated$SubClass)
NonNeu.integrated<- FindNeighbors(NonNeu.integrated, dims = 1:10, verbose = FALSE)
NonNeu.integrated= BuildClusterTree(object=NonNeu.integrated,dims=1:10)
nucleitree=Tool(object = NonNeu.integrated, slot = 'BuildClusterTree')

save(nucleitree,file='./Analysis/NonNeuron_SubClass_Tree_Feb2021.rda')

pdf('./TestPlots/NonNeuron_SubClass_Feb2021_10PCs_NODELABELS_dendrogram.pdf',width=12,height=8)
plot(nucleitree, edge.width = 2)
ape::nodelabels()
dev.off()


### bifurcation 

DefaultAssay(NonNeu.integrated) = 'RNA' 

Nodes = as.character(sort(unique(nucleitree$edge[,1])))

NodeTFdiff = vector(mode='list',length=nucleitree$Nnode)
names(NodeTFdiff) = Nodes

NodeNPdiff = vector(mode='list',length=nucleitree$Nnode)
names(NodeNPdiff) = Nodes

NodeclustsSubClass = vector(mode='list',length=nucleitree$Nnode)
names(NodeclustsSubClass) = Nodes

for(i in Nodes){

tmpNodes = nucleitree$edge[which(nucleitree$edge[,1]==as.numeric(i)),2]


if(tmpNodes[1]>length(nucleitree$tip.label)){
clusts1 = nucleitree$tip.label[getTips(nucleitree,tmpNodes[1])]
} else {
clusts1 = nucleitree$tip.label[tmpNodes[1]]
}

if(tmpNodes[2]>length(nucleitree$tip.label)){
clusts2 = nucleitree$tip.label[getTips(nucleitree,tmpNodes[2])]
} else {
clusts2 = nucleitree$tip.label[tmpNodes[2]]
}

cells1 = colnames(NonNeu.integrated)[which(!is.na(match(NonNeu.integrated@meta.data$SubClass,clusts1)))] 
cells2 = colnames(NonNeu.integrated)[which(!is.na(match(NonNeu.integrated@meta.data$SubClass,clusts2)))] 

NonNeu.integrated@meta.data$tmpGroups = NA
NonNeu.integrated@meta.data[cells1,'tmpGroups'] = 'A'
NonNeu.integrated@meta.data[cells2,'tmpGroups'] = 'B'
NonNeu.integrated@active.ident=as.factor(NonNeu.integrated$tmpGroups)

NodeTFdiff[[i]] = try(FindMarkers(NonNeu.integrated,ident.1='A',ident.2='B',features=intersect(rownames(NonNeu.integrated),HSTF$Name)))

NodeNPdiff[[i]] = try(FindMarkers(NonNeu.integrated,ident.1='A',ident.2='B',features=intersect(rownames(NonNeu.integrated),NPlist2)))

tmp=list(clusts1,clusts2)
names(tmp) = tmpNodes

NodeclustsNuclei[[i]] = tmp

cat('\n\n')
cat(paste(tmpNodes,collapse='_VS_'))
cat('\n\n')

}

NodeTFdiffnonneu = NodeTFdiff
NodeNPdiffnonneu = NodeNPdiff

save(NodeTFdiffnonneu,file='./Analysis/NodeTFdiff_nonneu_10PCs_Feb2021.rda',compress=TRUE)
save(NodeNPdiffnonneu,file='./Analysis/NodeNPdiff_nonneu_10PCs_Feb2021.rda',compress=TRUE)


## plot gene names 
pdf('./TestPlots/All_NonNeu_Feb2021_10PCs_TFgenes_dendrogram.pdf',width=12,height=12)
plot(nucleitree, edge.width = 1)
ape::nodelabels()
#ape::edgelabels()
for(i in Nodes){

tmpNodes = nucleitree$edge[which(nucleitree$edge[,1]==as.numeric(i)),2]

edge1 = which(nucleitree$edge[,1]==i & nucleitree$edge[,2]==tmpNodes[1])
edge2 = which(nucleitree$edge[,1]==i & nucleitree$edge[,2]==tmpNodes[2])

tmpDEG = NodeTFdiffnonneu[[i]]

sigEdge1 = tmpDEG[which(tmpDEG$avg_log2FC>0.5 & tmpDEG$p_val_adj<0.000000001),]
allEdge1 = na.omit(rownames(sigEdge1[order(sigEdge1$avg_log2FC,decreasing=TRUE),])[1:5])

if(length(allEdge1)>3){
edge1genes = paste(paste(allEdge1[1:3],collapse=','),'\n',paste(allEdge1[4:length(allEdge1)],collapse=','))
} else {
edge1genes = paste(allEdge1,collapse=',')
}

sigEdge2 = tmpDEG[which(tmpDEG$avg_log2FC<(-0.5) & tmpDEG$p_val_adj<0.000000001),]
allEdge2 = na.omit(rownames(sigEdge2[order(sigEdge2$avg_log2FC,decreasing=FALSE),])[1:5])

if(length(allEdge2)>3){
edge2genes = paste(paste(allEdge2[1:3],collapse=','),'\n',paste(allEdge2[4:length(allEdge2)],collapse=','))
} else {
edge2genes = paste(allEdge2,collapse=',')
}

ape::edgelabels(text=edge1genes,edge=edge1,adj=c(1,-1),frame='none',bg='white',cex=0.5)
ape::edgelabels(text=edge2genes,edge=edge2,adj=c(1,-1),frame='none',bg='white',cex=0.5)
}
dev.off()


### plot corresponding dot plot
NonNeu.integrated@active.ident=as.factor(NonNeu.integrated$SubClass)
All.markers.nonneu.HSTFs = FindAllMarkers(NonNeu.integrated,features=intersect(rownames(NonNeu.integrated),HSTF$Name))

write.csv(All.markers.nonneu.HSTFs,file='./Analysis/All.markers_RNA_nonneu_HSTFs_Feb2021.csv')


## evaluate 
ClustOrder = rev(unique(as.character(NonNeu.integrated@active.ident))[c(5,3,6,7,19,14,20,17,9,8,4,1,18,10,13,12,2,16,11,15)])
NonNeu.integrated@active.ident=factor(NonNeu.integrated@active.ident,levels=ClustOrder)
All.markers.nonneu.HSTFs$perDif=All.markers.nonneu.HSTFs$pct.1-All.markers.nonneu.HSTFs$pct.2

#MajTypes = paste('clust_',seq(42,0,-1),sep='')
MajTypes = ClustOrder

for(i in MajTypes){
if(i==MajTypes[1]){
    tmp = All.markers.nonneu.HSTFs[which(All.markers.nonneu.HSTFs$cluster==i ),]
    tmp = tmp[which(tmp$p_val_adj<=0.001 & tmp$pct.2<0.2),]
    tmp = tmp[order(tmp$perDif,decreasing=TRUE),]
    #tmp = tmp[which(tmp$avg_logFC>1),]
topPerDiftf <-  tmp$gene[1:3]
} else {
    tmp = All.markers.nonneu.HSTFs[which(All.markers.nonneu.HSTFs$cluster==i ),]
    tmp = tmp[which(tmp$p_val_adj<=0.001 & tmp$pct.2<0.2),]
    tmp = tmp[order(tmp$perDif,decreasing=TRUE),]
    #tmp = tmp[which(tmp$avg_logFC>1),]
topPerDiftf <-  c(topPerDiftf,tmp$gene[1:3])
}
}

topPerDiftf = na.omit(topPerDiftf)
topPerDiftf = unique(topPerDiftf)
topPerDiftf = rev(topPerDiftf)

pdf('./TestPlots/NonNeu_DotPlot_TopTF_RNA_Dend_Order_Feb2021.pdf',width=12,height=12)

DotPlot(NonNeu.integrated,features=topPerDiftf,assay='RNA') + RotatedAxis()
dev.off()






































### Junque

### did not redo subgroups 

### find node markers 

### try pc's 1:10

DefaultAssay(Tri2Neuron_ds.integrated)='integrated'

Tri2Neuron_ds.integrated@active.ident=as.factor(Tri2Neuron_ds.integrated$Neuron)

Tri2Neuron_ds.integrated<- FindNeighbors(Tri2Neuron_ds.integrated, dims = 1:10, verbose = FALSE)

Tri2Neuron_ds.integrated= BuildClusterTree(object=Tri2Neuron_ds.integrated,dims=1:10)

Neurontree=Tool(object = Tri2Neuron_ds.integrated, slot = 'BuildClusterTree')

pdf('./TestPlots/Tri2Neuron_ds_Neuron_Feb2021_10PCs_dendrogram.pdf',width=12,height=8)
plot(Neurontree, edge.width = 2)
ape::nodelabels()
#ape::tiplabels()
dev.off()

### bifurcation 

DefaultAssay(Tri2Neuron_ds.integrated) = 'RNA'

Nodes = as.character(sort(unique(Neurontree$edge[,1])))

NodeTFdiff = vector(mode='list',length=Neurontree$Nnode)
names(NodeTFdiff) = Nodes

NodeNPdiff = vector(mode='list',length=Neurontree$Nnode)
names(NodeNPdiff) = Nodes

NodeNeurons = vector(mode='list',length=Neurontree$Nnode)
names(NodeNeurons) = Nodes

for(i in Nodes){

tmpNodes = Neurontree$edge[which(Neurontree$edge[,1]==as.numeric(i)),2]


if(tmpNodes[1]>length(Neurontree$tip.label)){
clusts1 = Neurontree$tip.label[getTips(Neurontree,tmpNodes[1])]
} else {
clusts1 = Neurontree$tip.label[tmpNodes[1]]
}

if(tmpNodes[2]>length(Neurontree$tip.label)){
clusts2 = Neurontree$tip.label[getTips(Neurontree,tmpNodes[2])]
} else {
clusts2 = Neurontree$tip.label[tmpNodes[2]]
}

cells1 = colnames(Tri2Neuron_ds.integrated)[which(!is.na(match(Tri2Neuron_ds.integrated@meta.data$Neuron,clusts1)))] 

cells2 = colnames(Tri2Neuron_ds.integrated)[which(!is.na(match(Tri2Neuron_ds.integrated@meta.data$Neuron,clusts2)))] 

Tri2Neuron_ds.integrated@meta.data$tmpGroups = NA

Tri2Neuron_ds.integrated@meta.data[cells1,'tmpGroups'] = 'A'

Tri2Neuron_ds.integrated@meta.data[cells2,'tmpGroups'] = 'B'

Tri2Neuron_ds.integrated@active.ident=as.factor(Tri2Neuron_ds.integrated$tmpGroups)

NodeTFdiff[[i]] = try(FindMarkers(Tri2Neuron_ds.integrated,ident.1='A',ident.2='B',features=intersect(rownames(Tri2Neuron_ds.integrated),HSTF$Name)))

NodeNPdiff[[i]] = try(FindMarkers(Tri2Neuron_ds.integrated,ident.1='A',ident.2='B',features=intersect(rownames(Tri2Neuron_ds.integrated),NPlist2)))

tmp=list(clusts1,clusts2)
names(tmp) = tmpNodes

NodeNeurons[[i]] = tmp

cat('\n\n')
cat(paste(tmpNodes,collapse='_VS_'))
cat('\n\n')

}

NodeTFdiffNeuron = NodeTFdiff
NodeNPdiffNeuron = NodeNPdiff

save(NodeTFdiffNeuron,file='./Analysis/NodeTFdiff_Neuron_10PCs.rda',compress=TRUE)

save(NodeNPdiffNeuron,file='./Analysis/NodeNPdiff_Neuron_10PCs.rda',compress=TRUE)


#### node vs all else

NodeTFvsAll = vector(mode='list',length=Neurontree$Nnode)
names(NodeTFvsAll) = Nodes

NodeNPvsAll = vector(mode='list',length=Neurontree$Nnode)
names(NodeNPvsAll) = Nodes

for(i in Nodes[-1]){

clusts1 = Neurontree$tip.label[getTips(Neurontree,i)]

clusts2 = setdiff(Neurontree$tip.label,clusts1)

cells1 = colnames(Tri2Neuron_ds.integrated)[which(!is.na(match(Tri2Neuron_ds.integrated@meta.data$Neuron,clusts1)))] 

cells2 = colnames(Tri2Neuron_ds.integrated)[which(!is.na(match(Tri2Neuron_ds.integrated@meta.data$Neuron,clusts2)))] 

Tri2Neuron_ds.integrated@meta.data$tmpGroups = NA

Tri2Neuron_ds.integrated@meta.data[cells1,'tmpGroups'] = 'A'

Tri2Neuron_ds.integrated@meta.data[cells2,'tmpGroups'] = 'B'

Tri2Neuron_ds.integrated@active.ident=as.factor(Tri2Neuron_ds.integrated$tmpGroups)

NodeTFvsAll[[i]] = try(FindMarkers(Tri2Neuron_ds.integrated,ident.1='A',ident.2='B',features=intersect(rownames(Tri2Neuron_ds.integrated),HSTF$Name)))

NodeNPvsAll[[i]] = try(FindMarkers(Tri2Neuron_ds.integrated,ident.1='A',ident.2='B',features=intersect(rownames(Tri2Neuron_ds.integrated),NPlist2)))

cat('\n\n')
cat(i)
cat('\n\n')

}

NodeTFvsAllNeuron = NodeTFvsAll
NodeNPvsAllNeuron = NodeNPvsAll

save(NodeTFvsAllNeuron,file='./Analysis/NodeTFvsAll_Neuron_10PCs.rda',compress=TRUE)

save(NodeNPvsAllNeuron,file='./Analysis/NodeNPvsAll_Neuron_10PCs.rda',compress=TRUE)







