library(Seurat)

setwd('/local/projects-t3/idea/bherb/Hypothalamus/PubRes')

source('./Code/HypoPub_LoadFunctions.R')

HSTF = read.csv('/local/projects-t3/idea/bherb/annotation/Hsapiens/Human_TF.csv')

NPlist = read.csv('/local/projects-t3/idea/bherb/Hypothalamus/Neuropeptide_list.csv',header=FALSE) ## recieved from Hannah on 7/15/20

NPlist = MMtoHS(as.character(NPlist[,1]))

NPlist2 = unique(c(na.omit(as.character(TF_NP$Neuropeptides)),NPlist))

TF_NP = read.csv('/local/projects-t3/idea/bherb/Hypothalamus/TF_NP_list_from_Moffitt.csv') #direct from Moffitt suppmental materials 

for(i in 1:ncol(TF_NP)){
    TF_NP[,i] = MMtoHS(TF_NP[,i])
} ## human convention 


load('./SeuratObj/Tri2Samples_PostSCT.rda') #Tri2NormDat
load('./SeuratObj/GW18CtxSamples_PostSCT.rda') #GW18CtxNormDat
load('./SeuratObj/GW18GESamples_PostSCT.rda') #GW18GENormDat


#### more cells? The sample with the least number of cells at GW18 is Hypo, with 5201 - most others have ~8000. Try double - 4000 each, 

GW18AllGE4KGE4KCombDat = c(Tri2NormDat[['GW18']],GW18CtxNormDat,GW18GENormDat)
names(GW18AllGE4KCombDat)[1] = 'GW18_Hypo'

GW18AllGE4KSamples = names(GW18AllGE4KCombDat)

## new random sampling 

for(i in GW18AllGE4KSamples){
GW18AllGE4KCombDat[[i]] = subset(GW18AllGE4KCombDat[[i]],cells = colnames(GW18AllGE4KCombDat[[i]])[sample(c(1:ncol(GW18AllGE4KCombDat[[i]])),4000)])
}

for (i in GW18AllGE4KSamples) {
    GW18AllGE4KCombDat[[i]] <- NormalizeData(GW18AllGE4KCombDat[[i]], verbose = FALSE)
    GW18AllGE4KCombDat[[i]] <- FindVariableFeatures(GW18AllGE4KCombDat[[i]], selection.method = "vst", nfeatures = 2000, verbose = FALSE)
}

options(future.globals.maxSize = 4000 * 1024^2)

GW18AllGE4K.anchors <- FindIntegrationAnchors(object.list = GW18AllGE4KCombDat, dims = 1:30)
GW18AllGE4K.integrated <- IntegrateData(anchorset = GW18AllGE4K.anchors, dims = 1:30)
DefaultAssay(GW18AllGE4K.integrated) <- "integrated"
GW18AllGE4K.integrated <- ScaleData(GW18AllGE4K.integrated, verbose = FALSE)
GW18AllGE4K.integrated <- RunPCA(GW18AllGE4K.integrated, npcs = 30, verbose = FALSE)
GW18AllGE4K.integrated <- RunUMAP(GW18AllGE4K.integrated, reduction = "pca", dims = 1:30)

pdf('./TestPlots/GW18AllGE4K_Merged.pdf',width=12,height=8)
print(DimPlot(GW18AllGE4K.integrated, label = TRUE) + NoLegend())
for(i in c("sample","Cell_Class","Cell_State","Cell_Type","Maj_Clust_From_Ref","neuron_name")){
print(DimPlot(GW18AllGE4K.integrated, reduction = "umap", group.by = i, label = TRUE, repel = TRUE))
}
print(FeaturePlot(GW18AllGE4K.integrated, features = c("SLC17A6","GAD1"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(GW18AllGE4K.integrated, features = c("RAX","CHGB"), pt.size = 0.2, ncol = 2))
dev.off()

pdf('./TestPlots/GW18AllGE4K_Merged_SOX10.pdf',width=12,height=8)
print(DimPlot(GW18AllGE4K.integrated, label = TRUE) + NoLegend())
for(i in c("sample","Cell_Class","Cell_State","Cell_Type","Maj_Clust_From_Ref","neuron_name")){
print(DimPlot(GW18AllGE4K.integrated, reduction = "umap", group.by = i, label = TRUE, repel = TRUE))
}
print(FeaturePlot(GW18AllGE4K.integrated, features = c("SOX10"), pt.size = 0.4))
dev.off()

subplots(GW18AllGE4K.integrated,'sample',file='./TestPlots/GW18AllGE4K_Samples_IndvPlots.pdf')


save(GW18AllGE4K.integrated,file='./SeuratObj/GW18AllGE4K.integrated.rda',compress=TRUE)

save(GW18AllGE4KCombDat,file='./SeuratObj/GW18AllGE4KCombDat.rda',compress=TRUE)

### clustering 


GW18AllGE4K.integrated  <- CellCycleScoring(GW18AllGE4K.integrated, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)


### cluster cells

GW18AllGE4K.integrated <- FindNeighbors(GW18AllGE4K.integrated, dims = 1:30, verbose = FALSE)
GW18AllGE4K.integrated <- FindClusters(GW18AllGE4K.integrated, verbose = FALSE, resolution=0.5)

pdf('./TestPlots/GW18AllGE4K_Clusters_05.pdf',width=12,height=8)
print(DimPlot(GW18AllGE4K.integrated, label = TRUE) + NoLegend())
for(i in c("seurat_clusters","sample","Cell_Class","Cell_State","Cell_Type","Maj_Clust_From_Ref","neuron_name")){
print(DimPlot(GW18AllGE4K.integrated, reduction = "umap", group.by = i, label = TRUE, repel = TRUE))
}
dev.off()


for(i in seq(from=0.5, to=2, by=0.1)){
GW18AllGE4K.integrated <- FindClusters(GW18AllGE4K.integrated, verbose = FALSE, resolution=i)
}

pdf('./TestPlots/GW18AllGE4K_Check_clusters.pdf',width=12,height=8)
print(DimPlot(GW18AllGE4K.integrated, label = TRUE) + NoLegend())
for(i in c("seurat_clusters","sample","Cell_Class","Cell_State","Cell_Type","Maj_Clust_From_Ref","neuron_name")){
print(DimPlot(GW18AllGE4K.integrated, reduction = "umap", group.by = i, label = TRUE, repel = TRUE))
}
for (k in seq(from=0.5,to=2,by=0.1)){
print(DimPlot(GW18AllGE4K.integrated, reduction = "umap", group.by = paste("integrated_snn_res.",k,sep=""), label = TRUE, repel = TRUE))
}
dev.off()




## res of 1.5 looks good - 


## res = 1.5 captures major cell types 
subplots(GW18AllGE4K.integrated,'integrated_snn_res.1.5',file='./TestPlots/GW18AllGE4K_Res15_Indiv_Clusters.pdf')

subplots(GW18AllGE4K.integrated,'sample',file='./TestPlots/GW18AllGE4K_Indiv_Samples.pdf')

### explore DEGs between various cell types 


#TF markers - res 1.5



GW18AllGE4K.integrated@meta.data$seurat_clusters = paste('clust_',GW18AllGE4K.integrated@meta.data$integrated_snn_res.1.5,sep='')

GW18AllGE4K.integrated@active.ident=as.factor(GW18AllGE4K.integrated$seurat_clusters)

#GW18AllGE4K.integrated@active.ident=as.factor(paste('Cluster_',as.character(GW18AllGE4K.integrated$integrated_snn_res.1.5),'_Yup',sep=''))


DefaultAssay(GW18AllGE4K.integrated) = 'RNA'
Markers.TFs = FindAllMarkers(GW18AllGE4K.integrated,features=intersect(toupper(unique(hstf$Name)),rownames(GW18AllGE4K.integrated@assays$RNA)))

DefaultAssay(GW18AllGE4K.integrated) = 'integrated'
VarGenes = VariableFeatures(GW18AllGE4K.integrated)
DefaultAssay(GW18AllGE4K.integrated) = 'RNA'
Markers.VARs = FindAllMarkers(GW18AllGE4K.integrated,features=VarGenes)




Markers.NPs = FindAllMarkers(GW18AllGE4K.integrated,features=intersect(toupper(unique(TF_NP$Neuropeptides)),rownames(GW18AllGE4K.integrated@assays$RNA)))

Markers.NRs = FindAllMarkers(GW18AllGE4K.integrated,features=intersect(toupper(unique(TF_NP$Neuromodulator_receptors)),rownames(GW18AllGE4K.integrated@assays$RNA)))

Markers.NMs = FindAllMarkers(GW18AllGE4K.integrated,features=intersect(toupper(unique(TF_NP$Neuromodulator_production_and_transport)),rownames(GW18AllGE4K.integrated@assays$RNA)))


Markers.Neuro = FindAllMarkers(GW18AllGE4K.integrated,features=intersect(toupper(unique(c(TF_NP[,1],TF_NP[,3],TF_NP[,4]))),rownames(GW18AllGE4K.integrated@assays$RNA)))



write.csv(Markers.VARs,file='./Analysis/GW18AllGE4K_ClusterMarkers_VarGenes.csv')


##MAdeg = read.table('./Analysis/MiniAtlas_Sup6_DEGs.txt',header=TRUE,sep=',')
## presented as all the same cluster name?!? 

DefaultAssay(GW18AllGE4K.integrated) = 'RNA'
Markers.Layer = FindAllMarkers(GW18AllGE4K.integrated,features=intersect(c('SLC30A3','OTOF','RORB','RSPO1','FEZF2','SULF1','CAR3','FAM84N','SLA2','FOXP2','NXPH4'),rownames(GW18AllGE4K.integrated@assays$RNA)))



### need to redo the 2D umap 

GW18AllGE4K_2Dumap.integrated = GW18AllGE4K.integrated

#FindVariableFeatures(GW18AllGE4K_2Dumap.integrated)

#GW18AllGE4K_2Dumap.integrated <- RunPCA(GW18AllGE4K_2Dumap.integrated, npcs = 30, verbose = FALSE)
GW18AllGE4K_2Dumap.integrated <- RunUMAP(GW18AllGE4K_2Dumap.integrated, reduction = "pca", dims = 1:30)

pdf('./TestPlots/GW18AllGE4K_2Dumap_Merged.pdf',width=12,height=8)
print(DimPlot(GW18AllGE4K_2Dumap.integrated, label = TRUE) + NoLegend())
for(i in c("sample","seurat_clusters","Cell_Class","Cell_State","Cell_Type","Maj_Clust_From_Ref","neuron_name")){
print(DimPlot(GW18AllGE4K_2Dumap.integrated, reduction = "umap", group.by = i, label = TRUE, repel = TRUE))
}
print(FeaturePlot(GW18AllGE4K_2Dumap.integrated, features = c("SLC17A6","GAD1"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(GW18AllGE4K_2Dumap.integrated, features = c("RAX","CHGB"), pt.size = 0.2, ncol = 2))
dev.off()

pdf('./TestPlots/GW18AllGE4K_2Dumap_Merged_LayerMarkers.pdf',width=12,height=8)
print(DimPlot(GW18AllGE4K_2Dumap.integrated, label = TRUE) + NoLegend())
for(i in c("sample","seurat_clusters","Cell_Class","Cell_State","Cell_Type","Maj_Clust_From_Ref","neuron_name")){
print(DimPlot(GW18AllGE4K_2Dumap.integrated, reduction = "umap", group.by = i, label = TRUE, repel = TRUE))
}
for(i in intersect(c('CUX2','SLC30A3','OTOF','RORB','RSPO1','FEZF2','SULF1','CAR3','FAM84N','SLA2','FOXP2','NXPH4'),rownames(GW18AllGE4K.integrated@assays$RNA))){
	print(FeaturePlot(GW18AllGE4K_2Dumap.integrated, features = i, pt.size = 0.4))
}

dev.off()

### trying human IT cell markers by layer from cross species paper 

pdf('./TestPlots/GW18AllGE4K_2Dumap_Merged_LayerMarkersCrossSpeciesPaper.pdf',width=12,height=8)
print(DimPlot(GW18AllGE4K_2Dumap.integrated, label = TRUE) + NoLegend())
for(i in c("sample","seurat_clusters","Cell_Class","Cell_State","Cell_Type","Maj_Clust_From_Ref","neuron_name")){
print(DimPlot(GW18AllGE4K_2Dumap.integrated, reduction = "umap", group.by = i, label = TRUE, repel = TRUE))
}
for(i in intersect(c('CUX2','FAM19A1','LOC105376457','OR51E2','LOC105374971','RORB','NPSR1_AS1','EGFEM1P','HS3ST4','SEMA5A','ATP10A','RGS12','LINC00343','LOC1019282784','THEMIS','MDFIC','SEMA3D'),rownames(GW18AllGE4K_2Dumap.integrated@assays$RNA))){
	print(FeaturePlot(GW18AllGE4K_2Dumap.integrated, features = i, pt.size = 0.4))
}

dev.off()


pdf('./TestPlots/GW18AllGE4K_2Dumap_Merged_LayerMarkersWiki.pdf',width=12,height=8)
print(DimPlot(GW18AllGE4K_2Dumap.integrated, label = TRUE) + NoLegend())
for(i in c("sample","seurat_clusters","Cell_Class","Cell_State","Cell_Type","Maj_Clust_From_Ref","neuron_name")){
print(DimPlot(GW18AllGE4K_2Dumap.integrated, reduction = "umap", group.by = i, label = TRUE, repel = TRUE))
}
for(i in intersect(c('CUX2','SATB2','TBR1','OTX1','CTIP2'),rownames(GW18AllGE4K_2Dumap.integrated@assays$RNA))){
	print(FeaturePlot(GW18AllGE4K_2Dumap.integrated, features = i, pt.size = 0.4))
}

dev.off()




pdf('./TestPlots/GW18AllGE4K_2Dumap_Merged_Inhib.pdf',width=12,height=8)
print(DimPlot(GW18AllGE4K_2Dumap.integrated, label = TRUE) + NoLegend())
for(i in c("sample","seurat_clusters","Cell_Class","Cell_State","Cell_Type","Maj_Clust_From_Ref","neuron_name")){
print(DimPlot(GW18AllGE4K_2Dumap.integrated, reduction = "umap", group.by = i, label = TRUE, repel = TRUE))
}
for(i in intersect(c('GAD1','PROX1','LAMP5','NDNF','SNCG','VIP','LHX6','SOX6','SIP1','SATB1','SST','CHODL','PVALB','VIPR2'),rownames(GW18AllGE4K.integrated@assays$RNA))){
	print(FeaturePlot(GW18AllGE4K_2Dumap.integrated, features = i, pt.size = 0.4))
}

dev.off()


pdf('./TestPlots/GW18AllGE4K_2Dumap_BrainRegions_Replot.pdf',width=12,height=8)
print(DimPlot(GW18AllGE4K_2Dumap.integrated, reduction = "umap", group.by = 'BrainRegion', label = TRUE, repel = TRUE))
dev.off()


## port major class cell types to object 

GW18AllGE4K_2Dumap.integrated@meta.data$MajorClass = NA

GW18AllGE4K_2Dumap.integrated@meta.data$MajorClass = All.integrated@meta.data[colnames(GW18AllGE4K_2Dumap.integrated),'MajorClass']


GW18AllGE4K_2Dumap.integrated@meta.data$SubClass = NA

GW18AllGE4K_2Dumap.integrated@meta.data$SubClass = All.integrated@meta.data[colnames(GW18AllGE4K_2Dumap.integrated),'SubClass']



pdf('./TestPlots/GW18AllGE4K_2Dumap_CellTypes_Replot.pdf',width=12,height=8)
print(DimPlot(GW18AllGE4K_2Dumap.integrated, label = TRUE) + NoLegend())
for(i in c("sample","seurat_clusters","Cell_Class","Cell_State","Cell_Type","MajorClass","SubClass","Phase")){
print(DimPlot(GW18AllGE4K_2Dumap.integrated, reduction = "umap", group.by = i, label = TRUE, repel = TRUE))
}
dev.off()






save(GW18AllGE4K_2Dumap.integrated,file='./SeuratObj/GW18AllGE4K_2Dumap.integrated.rda',compress=TRUE)



genes=c('NEUROD2','NEUROD6')

pdf(file='./TestPlots/GW18AllGE4K_2Dumap_ExcGenes.pdf',width=12,height=8)
for(i in genes){
print(FeaturePlot(GW18AllGE4K_2Dumap.integrated, features = i, pt.size = 0.2,order=TRUE))
}
dev.off()


### follow up on genes driving GPC3 expression


genes=c('GPC3',as.character(genieGW18$TF[which(genieGW18$target=='GPC3')][1:10]))

DefaultAssay(GW18AllGE4K_2Dumap.integrated)='RNA'

pdf(file='./TestPlots/GW18AllGE4K_2Dumap_GPC3_Top10TF.pdf',width=12,height=8)
for(i in c("sample","seurat_clusters","Cell_Type","Maj_Clust_From_Ref")){
print(DimPlot(GW18AllGE4K_2Dumap.integrated, reduction = "umap", group.by = i, label = TRUE, repel = TRUE))
}
for(i in genes){
print(FeaturePlot(GW18AllGE4K_2Dumap.integrated, features = i, pt.size = 0.2,order=TRUE))
}
dev.off()



#### sample differences 

## inhib - clusters 15, 17, 21 


MajTypes=unique(GW18AllGE4K.integrated$seurat_clusters)

MajTypes=sort(MajTypes)

SpeciesMarkers = vector(mode='list',length=length(MajTypes))
names(SpeciesMarkers) = MajTypes

for(i in MajTypes){

GW18AllGE_Tmp = subset(GW18AllGE4K.integrated,idents=i)

GW18AllGE_Tmp@active.ident=as.factor(GW18AllGE_Tmp$sample)

SpeciesMarkers[[i]] = FindAllMarkers(GW18AllGE_Tmp,features=VarGenes)
}

save(SpeciesMarkers,file='./Analysis/Hypo_Ctx_GE_Species_DEG.rda',compress=TRUE) ### this is tissue (sample level) comparisons within each cluster 

## heatmaps within clusters across samples 



for(j in MajTypes){
tmpDEG = SpeciesMarkers[[j]]
MajTypes2 = unique(tmpDEG$cluster)
for(i in MajTypes2){
if(i==MajTypes2[1]){
    tmp = tmpDEG[which(tmpDEG$cluster==i ),]
    tmp = tmp[which(tmp$p_val_adj<=0.01 & tmp$avg_logFC>0 & tmp$pct.2<0.1),]
    tmp$perDif = tmp$pct.1-tmp$pct.2
    tmp = tmp[order(tmp$perDif,decreasing=TRUE),]
    #tmp = tmp[which(tmp$avg_logFC>1),]
topPerDif <-  as.character(na.omit(tmp$gene[1:10]))
} else {
    tmp = tmpDEG[which(tmpDEG$cluster==i ),]
    tmp = tmp[which(tmp$p_val_adj<=0.01),]
    tmp$perDif = tmp$pct.1-tmp$pct.2
    tmp = tmp[order(tmp$perDif,decreasing=TRUE),]
    #tmp = tmp[which(tmp$avg_logFC>1),]
topPerDif <-  c(topPerDif,as.character(na.omit(tmp$gene[1])))
}
}

DefaultAssay(GW18AllGE4K.integrated) = 'integrated'

GW18AllGE_Tmp = subset(GW18AllGE4K.integrated,idents=j)

pdf(paste('./TestPlots/GW18AllGE4K_',j,'_topVAR_heatmap.pdf',sep=''),width=12,height=8)

print(DoHeatmap(GW18AllGE_Tmp,features = topPerDif,group.by = 'sample',label=FALSE,raster=FALSE) + NoLegend())

dev.off()

cat(paste('\n\n',j,' Done','\n\n',sep=''))

}

#topPerDif = na.omit(topPerDif)


#### group by region - Hypo, GE, Ctx 

GW18AllGE4K.integrated@meta.data$BrainRegion = NA

GW18AllGE4K.integrated@meta.data$BrainRegion[which(GW18AllGE4K.integrated@meta.data$sample == 'GW18')] = 'Hypo'

GW18AllGE4K.integrated@meta.data$BrainRegion[grep('GW18_',GW18AllGE4K.integrated@meta.data$sample)] =  'Cortex'

GW18AllGE4K.integrated@meta.data$BrainRegion[grep('GE',GW18AllGE4K.integrated@meta.data$sample)] =  'GE'



GW18AllGE4K_2Dumap.integrated@meta.data$BrainRegion = NA

GW18AllGE4K_2Dumap.integrated@meta.data$BrainRegion[which(GW18AllGE4K_2Dumap.integrated@meta.data$sample == 'GW18')] = 'Hypo'

GW18AllGE4K_2Dumap.integrated@meta.data$BrainRegion[grep('GW18_',GW18AllGE4K_2Dumap.integrated@meta.data$sample)] =  'Cortex'

GW18AllGE4K_2Dumap.integrated@meta.data$BrainRegion[grep('GE',GW18AllGE4K_2Dumap.integrated@meta.data$sample)] =  'GE'


### individual plots for publication

pdf('./TestPlots/GW18AllGE4K_2Dumap_Merged_Regions.pdf',width=12,height=8)
print(DimPlot(GW18AllGE4K_2Dumap.integrated, reduction = "umap", group.by = 'BrainRegion', label = TRUE, repel = TRUE))
dev.off()

pdf('./TestPlots/GW18AllGE4K_2Dumap_Merged_Cortex_CellType.pdf',width=12,height=8)
print(DimPlot(GW18AllGE4K_2Dumap.integrated, reduction = "umap", group.by = 'Cell_Type', label = TRUE, repel = TRUE))
dev.off()




MajTypes=unique(GW18AllGE4K.integrated$seurat_clusters)

MajTypes=sort(MajTypes)

RegionMarkers = vector(mode='list',length=length(MajTypes))
names(RegionMarkers) = MajTypes

for(i in MajTypes){

GW18AllGE_Tmp = subset(GW18AllGE4K.integrated,idents=i)

GW18AllGE_Tmp@active.ident=as.factor(GW18AllGE_Tmp$BrainRegion)

RegionMarkers[[i]] = FindAllMarkers(GW18AllGE_Tmp,features=VarGenes)
}

save(RegionMarkers,file='./Analysis/Hypo_Ctx_GE_BrainRegion_DEG.rda',compress=TRUE) #major region comp within cluster 


clusts = paste('clust_',c(0:29),sep='')
for(i in clusts) {
if(i==clusts[1]){
tmp = RegionMarkers[[i]]
tmp$FDR = p.adjust(tmp$p_val)
tmp=tmp[,c(7,6,2:4,1,5,8)]
tmp=tmp[which(tmp$FDR<=0.05),]

write.xlsx2(tmp, file='./Analysis/Hypo_Ctx_GE_BrainRegion_DEG.xlsx', sheetName=i,col.names=TRUE, row.names=FALSE, append=FALSE)
} else {
tmp = RegionMarkers[[i]]
tmp$FDR = p.adjust(tmp$p_val)
tmp=tmp[,c(7,6,2:4,1,5,8)]
tmp=tmp[which(tmp$FDR<=0.05),]

write.xlsx2(tmp, file='./Analysis/Hypo_Ctx_GE_BrainRegion_DEG.xlsx', sheetName=i,col.names=TRUE, row.names=FALSE, append=TRUE)
}
cat(i)
gc()
}














for(j in MajTypes){
tmpDEG = RegionMarkers[[j]]
MajTypes2 = unique(tmpDEG$cluster)
for(i in rev(MajTypes2)){
if(i==MajTypes2[3]){
    tmp = tmpDEG[which(tmpDEG$cluster==i ),]
    tmp = tmp[which(tmp$p_val_adj<=0.01 & tmp$avg_logFC>0 & tmp$pct.2<0.1),]
    tmp$perDif = tmp$pct.1-tmp$pct.2
    tmp = tmp[order(tmp$perDif,decreasing=TRUE),]
    #tmp = tmp[which(tmp$avg_logFC>1),]
topPerDif <-  as.character(na.omit(tmp$gene[1:10]))
} else {
    tmp = tmpDEG[which(tmpDEG$cluster==i ),]
    tmp = tmp[which(tmp$p_val_adj<=0.01 & tmp$avg_logFC>0 & tmp$pct.2<0.1),]
    tmp$perDif = tmp$pct.1-tmp$pct.2
    tmp = tmp[order(tmp$perDif,decreasing=TRUE),]
    #tmp = tmp[which(tmp$avg_logFC>1),]
topPerDif <-  c(topPerDif,as.character(na.omit(tmp$gene[1:10])))
}
}

DefaultAssay(GW18AllGE4K.integrated) = 'integrated'

GW18AllGE_Tmp = subset(GW18AllGE4K.integrated,idents=j)

pdf(paste('./TestPlots/GW18AllGE4K_',j,'_topVAR_BrainRegion_heatmap.pdf',sep=''),width=12,height=8)

print(DoHeatmap(GW18AllGE_Tmp,features = topPerDif,group.by = 'sample',label=FALSE,raster=FALSE) + NoLegend())

dev.off()

cat(paste('\n\n',j,' Done','\n\n',sep=''))

}


head(RegionMarkers[['clust_16']][which(RegionMarkers[['clust_16']]$avg_logFC>0 & RegionMarkers[['clust_16']]$pct.2<0.1 & RegionMarkers[['clust_16']]$cluster=='Cortex'),],50)

CtxGene = RegionMarkers[['clust_16']]$gene[which(RegionMarkers[['clust_16']]$avg_logFC>0 & RegionMarkers[['clust_16']]$pct.2<0.1 & RegionMarkers[['clust_16']]$cluster=='Cortex')]


rownames(head(RegionMarkers[['clust_16']][which(RegionMarkers[['clust_16']]$avg_logFC>0 & RegionMarkers[['clust_16']]$pct.2<0.1 & RegionMarkers[['clust_16']]$cluster=='GE'),],50))

GEGene = RegionMarkers[['clust_16']]$gene[which(RegionMarkers[['clust_16']]$avg_logFC>0 & RegionMarkers[['clust_16']]$pct.2<0.1 & RegionMarkers[['clust_16']]$cluster=='GE')]

HypoGene = RegionMarkers[['clust_16']]$gene[which(RegionMarkers[['clust_16']]$avg_logFC>0 & RegionMarkers[['clust_16']]$pct.2<0.1 & RegionMarkers[['clust_16']]$cluster=='Hypo')]



GW18AllGE_Tmp = subset(GW18AllGE4K.integrated,idents='clust_16')


CtxGene2 = intersect(CtxGene[1:10],rownames(GW18AllGE_Tmp))

GEGene2 = intersect(GEGene[1:10],rownames(GW18AllGE_Tmp))

HypoGene2 = intersect(HypoGene[1:10],rownames(GW18AllGE_Tmp))

totGenes = c(CtxGene2,GEGene2,HypoGene2)

pdf(paste('./TestPlots/GW18AllGE4K_clust_16_topDif_BrainRegion_heatmap.pdf',sep=''),width=12,height=8)

print(DoHeatmap(GW18AllGE_Tmp,features = totGenes,group.by = 'BrainRegion',label=FALSE,raster=FALSE) + NoLegend())

dev.off()







### individual cluster comparisons 

GW18AllGE4K.integrated.test = GW18AllGE4K.integrated

DefaultAssay(GW18AllGE4K.integrated.test) = 'RNA'

#GW18AllGE4K.integrated.test <- FindVariableFeatures(GW18AllGE4K.integrated.test, selection.method = "vst", nfeatures = 10000)

#VarGenes = VariableFeatures(GW18AllGE4K.integrated.test)

GW18AllGE4K.integrated.test <- ScaleData(object = GW18AllGE4K.integrated.test)


#DefaultAssay(GW18AllGE4K.integrated.test) = 'integrated'

GW18AllGE4K.integrated@active.ident=as.factor(GW18AllGE4K.integrated$BrainRegion)

GW18AllGE_Tmp = subset(GW18AllGE4K.integrated,idents=unique(GW18AllGE4K.integrated@active.ident)[2])

GW18AllGE_Tmp@active.ident=as.factor(GW18AllGE_Tmp$seurat_clusters)

TMPdeg = FindMarkers(GW18AllGE_Tmp,ident.1='clust_14',ident.2='clust_19')


TMPdegSig = TMPdeg[which(TMPdeg$p_val_adj<0.001),]
TMPdegSig = TMPdegSig[order(TMPdegSig$avg_logFC,decreasing=TRUE),]


FindMarkers(ident.1,ident.2)



### Trying slingshot for linage analysis - seems straightforward and can accept precomputed PCAs



data("slingshotExample")
rd <- slingshotExample$rd
cl <- slingshotExample$cl

dim(rd)


GW18AllGE4K.integrated@active.ident=as.factor(GW18AllGE4K.integrated$seurat_clusters)

GW18AllGE_Tmp = subset(GW18AllGE4K.integrated,idents=unique(GW18AllGE4K.integrated@active.ident)[c(14,1,24,8,28,26)])

GW18AllGE_PCA = GW18AllGE_Tmp@reductions[["pca"]]@cell.embeddings[,1:30]

GW18AllGE_Clust = GW18AllGE_Tmp$seurat_clusters

GW18AllGE_SCE = Seurat::as.SingleCellExperiment(GW18AllGE_Tmp) ## just 2000 var genes - try for now? - transfered PCA


GW18AllGE_SCE <- slingshot(GW18AllGE_SCE, clusterLabels = 'seurat_clusters', reducedDim = 'PCA')



colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(GW18AllGE_SCE$slingPseudotime_1, breaks=100)]

pdf(file='./TestPlots/GW18AllGE4K_EarlyNeuroLineage_PC1_2.pdf',width=12,height=8)
plot(reducedDims(GW18AllGE_SCE)$PCA, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(GW18AllGE_SCE), lwd=2, col='black')
plot(reducedDims(GW18AllGE_SCE)$PCA, col = brewer.pal(9,'Set1')[GW18AllGE_SCE$seurat_clusters], pch=16, asp = 1)
lines(SlingshotDataSet(GW18AllGE_SCE), lwd=2, type = 'lineages', col = 'black')
dev.off()


lin1 <- getLineages(GW18AllGE_PCA, GW18AllGE_Clust, start.clus = 'clust_14')


### lineages: 2 
#Lineage1: clust_14  clust_19  clust_7  clust_6  clust_10  
#Lineage2: clust_14  clust_19  clust_7  clust_6  clust_0 

## separate datasets 


GW18AllGE_Ctx = subset(GW18AllGE4K.integrated,idents=unique(GW18AllGE4K.integrated@active.ident)[c(14,24,8,28,26)])

GW18AllGE_PCA_Ctx = GW18AllGE_Ctx@reductions[["pca"]]@cell.embeddings[,1:30]

GW18AllGE_Clust_Ctx = GW18AllGE_Ctx$seurat_clusters

linCtx <- getLineages(GW18AllGE_PCA_Ctx, GW18AllGE_Clust_Ctx, start.clus = 'clust_14')

GW18AllGE_Ctx_SCE = Seurat::as.SingleCellExperiment(GW18AllGE_Ctx) ## just 2000 var genes - try for now? - transfered PCA

GW18AllGE_Ctx_SCE <- slingshot(GW18AllGE_Ctx_SCE, clusterLabels = 'seurat_clusters', reducedDim = 'PCA')

t <- GW18AllGE_Ctx_SCE$slingPseudotime_1

# for time, only look at the 100 most variable genes
#Y <- log1p(assays(GW18AllGE_Ctx_SCE)$norm)
#var100 <- names(sort(apply(Y,1,var),decreasing = TRUE))[1:100]
#Y <- Y[var100,]
Y = GW18AllGE_Ctx_SCE@assays@data[[1]]


# fit a GAM with a loess term for pseudotime
Ctx.pval <- apply(Y,1,function(z){
    d <- data.frame(z=z, t=t)
    suppressWarnings({
      tmp <- suppressWarnings(gam(z ~ lo(t), data=d))
    })
    p <- summary(tmp)[3][[1]][2,3]
    p
})

##FDR? 

Ctx.fdr = p.adjust(gam.pval,method='fdr')

Ctx.cor = NA

for(i in 1:nrow(Y)){
Ctx.cor[i] = cor(t,Y[i,])
if(i%%100==0) cat(paste(i,', ',sep=''))
}





## also calculate cor 

## check pseudotime

GW18AllGE4K_2Dumap.integrated@meta.data$Ctx_pseudotime = 0

GW18AllGE4K_2Dumap.integrated@meta.data[colnames(GW18AllGE_Ctx_SCE),'Ctx_pseudotime'] = t


pdf('./TestPlots/GW18AllGE4K_2D_CTX_TestPseudotime.pdf')
print(FeaturePlot(GW18AllGE4K_2Dumap.integrated, features = 'Ctx_pseudotime', pt.size = 0.4))
for(i in c('TUBA1A','STMN2','MEF2C','SYT4','TUBB','SLA','PTN','CENPF','UBE2C','TOP2A','HMGB2','PTTG1')){
print(FeaturePlot(GW18AllGE4K_2Dumap.integrated, features = i, pt.size = 0.4))
}
dev.off()

## does look good



#### Full cortex IT lineage 


GW18AllGE_CtxIT = subset(GW18AllGE4K.integrated,idents=unique(GW18AllGE4K.integrated@active.ident)[c(14,24,8,28,26,19,16,18)])


GW18AllGE_CtxIT = subset(GW18AllGE_CtxIT,cells = colnames(GW18AllGE_CtxIT)[which(GW18AllGE_CtxIT@meta.data$BrainRegion=='Cortex')])


GW18AllGE_PCA_CtxIT = GW18AllGE_CtxIT@reductions[["pca"]]@cell.embeddings[,1:30]

GW18AllGE_Clust_CtxIT = GW18AllGE_CtxIT$seurat_clusters

linCtxIT <- getLineages(GW18AllGE_PCA_CtxIT, GW18AllGE_Clust_CtxIT, start.clus = 'clust_14')

GW18AllGE_CtxIT_SCE = Seurat::as.SingleCellExperiment(GW18AllGE_CtxIT) ## just 2000 var genes - try for now? - transfered PCA

GW18AllGE_CtxIT_SCE <- slingshot(GW18AllGE_CtxIT_SCE, clusterLabels = 'seurat_clusters', reducedDim = 'PCA',start.clus = 'clust_14', end.clus = 'clust_1')

CtxIT.t <- GW18AllGE_CtxIT_SCE$slingPseudotime_1

CtxIT.Y = GW18AllGE_CtxIT_SCE@assays@data[[1]]

GW18AllGE4K_2Dumap.integrated@meta.data$CtxIT_pseudotime = 0

GW18AllGE4K_2Dumap.integrated@meta.data[colnames(GW18AllGE_CtxIT_SCE),'CtxIT_pseudotime'] = CtxIT.t


pdf('./TestPlots/GW18AllGE4K_2D_CTXIT_TestPseudotime.pdf')
print(FeaturePlot(GW18AllGE4K_2Dumap.integrated, features = 'CtxIT_pseudotime', pt.size = 0.4))
dev.off()

# fit a GAM with a loess term for pseudotime
CtxIT.pval <- apply(CtxIT.Y,1,function(z){
    d <- data.frame(z=z, t=CtxIT.t)
    suppressWarnings({
      tmp <- suppressWarnings(gam(z ~ lo(CtxIT.t), data=d))
    })
    p <- summary(tmp)[3][[1]][2,3]
    p
})

##FDR
CtxIT.fdr = p.adjust(CtxIT.pval,method='fdr')
CtxIT.cor = NA

for(i in 1:nrow(CtxIT.Y)){
CtxIT.cor[i] = cor(CtxIT.t,CtxIT.Y[i,])
if(i%%100==0) cat(paste(i,', ',sep=''))
}

Ctx.ITres=data.frame(gene=names(CtxIT.fdr),pval=CtxIT.pval,fdr=CtxIT.fdr,cor=CtxIT.cor)

rownames(Ctx.ITres)[which(Ctx.ITres$fdr<0.05 & Ctx.ITres$cor>0.5)]
"TUBA1A" "STMN2"  "MEF2C"  "SYT4"   "TUBB"   "SLA" 


pdf('./TestPlots/GW18AllGE4K_2D_CTX_IT_TestPseudotime.pdf')
print(FeaturePlot(GW18AllGE4K_2Dumap.integrated, features = 'CtxIT_pseudotime', pt.size = 0.4))
for(i in c(rownames(Ctx.ITres)[which(Ctx.ITres$fdr<0.05 & Ctx.ITres$cor>0.5)],rownames(Ctx.ITres)[which(Ctx.ITres$fdr<0.05 & Ctx.ITres$cor<(-0.5))])){
print(FeaturePlot(GW18AllGE4K_2Dumap.integrated, features = i, pt.size = 0.4,order=TRUE))
}
dev.off()

Ctx.ITgenes = c(rownames(Ctx.ITres)[which(Ctx.ITres$fdr<0.05 & Ctx.ITres$cor>0.5)],rownames(Ctx.ITres)[which(Ctx.ITres$fdr<0.05 & Ctx.ITres$cor<(-0.5))])

Ctx.ITdata=cbind(CtxIT.t,t(CtxIT.Y[Ctx.ITgenes,]))
Ctx.ITdata=as.data.frame(Ctx.ITdata)
colnames(Ctx.ITdata) = c('pseudotime',Ctx.ITgenes)

Ctx.ITdata$cluster = GW18AllGE4K_2Dumap.integrated@meta.data[colnames(GW18AllGE_CtxIT_SCE),'seurat_clusters']

Ctx.ITdataM = reshape2::melt(data=Ctx.ITdata, id.vars = "pseudotime", measure.vars = Ctx.ITgenes)

pdf('./TestPlots/GW18AllGE4K_2D_CTX_IT_TestPseudotime_SmoothLines.pdf',width=12,height=8)
print(ggplot(Ctx.ITdataM, aes(x = pseudotime, y = value, colour = variable)) + geom_smooth())
dev.off()

Ctx.ITclusts = unique(Ctx.ITdata$cluster)[c(4,7,2,8,5,6,3,1)]

Ctx.ITmins = tapply(Ctx.ITdata$pseudotime,INDEX=Ctx.ITdata$cluster,FUN=min)[Ctx.ITclusts]
Ctx.ITmaxes = tapply(Ctx.ITdata$pseudotime,INDEX=Ctx.ITdata$cluster,FUN=max)[Ctx.ITclusts]

Ctx.ITd=data.frame(x1=Ctx.ITmins, x2=Ctx.ITmaxes, y1=seq(from=(-0.2),to=(-0.2875),by=-0.0125), y2=seq(from=(-0.2125),to=(-0.3),by=-0.0125), t=Ctx.ITclusts, r=c(1:length(Ctx.ITclusts)))

if(FALSE){ ## error
pdf('./TestPlots/GW18AllGE4K_2D_CTX_IT_TestPseudotime_SmoothLines_clusts.pdf',width=12,height=8)
print(ggplot(data=Ctx.ITdataM, aes(x = pseudotime, y = value, colour = variable)) + geom_smooth() +  geom_rect(data=Ctx.ITd, inherit.aes = FALSE, aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2, fill=t), color="black", alpha=0.5) + geom_text(data=Ctx.ITd, aes(x=x1+(x2-x1)/2, y=y1+(y2-y1)/2, label=t), size=4))
dev.off()
}

pdf('./TestPlots/GW18AllGE4K_2D_CTX_IT_TestPseudotime_SmoothLines_clusts.pdf',width=16,height=8)
print(ggplot(data=Ctx.ITdataM, aes(x = pseudotime, y = value, colour = variable)) + geom_smooth() +  annotate("text", x = rowMeans(Ctx.ITd[,1:2]), y = -0.2, label = Ctx.ITd$t))
dev.off() ## gene expression over pseudotime 





### GE and Hypo  - repeat if needed 

GW18AllGE_GEH_C3 = subset(GW18AllGE4K.integrated,idents=unique(GW18AllGE4K.integrated@active.ident)[c(14,1,3)])

##sample 

GW18AllGE_GEH_C3 = subset(GW18AllGE_GEH_C3,cells = colnames(GW18AllGE_GEH_C3)[which(GW18AllGE_GEH_C3@meta.data$BrainRegion!='Cortex')])


GW18AllGE_PCA_GEH_C3 = GW18AllGE_GEH_C3@reductions[["pca"]]@cell.embeddings[,1:30]

GW18AllGE_Clust_GEH_C3 = GW18AllGE_GEH_C3$seurat_clusters

linGEH_C3 <- getLineages(GW18AllGE_PCA_GEH_C3, GW18AllGE_Clust_GEH_C3, start.clus = 'clust_14')

GW18AllGE_GEH_C3_SCE = Seurat::as.SingleCellExperiment(GW18AllGE_GEH_C3) ## just 2000 var genes - try for now? - transfered PCA

GW18AllGE_GEH_C3_SCE <- slingshot(GW18AllGE_GEH_C3_SCE, clusterLabels = 'seurat_clusters', reducedDim = 'PCA',start.clus = 'clust_14', end.clus = 'clust_3')


GEH_C3.t <- GW18AllGE_GEH_C3_SCE$slingPseudotime_1

GW18AllGE4K_2Dumap.integrated@meta.data$GEH_C3_pseudotime = 0

GW18AllGE4K_2Dumap.integrated@meta.data[colnames(GW18AllGE_GEH_C3_SCE),'GEH_C3_pseudotime'] = GEH_C3.t

pdf('./TestPlots/GW18AllGE4K_2D_GEH_C3_TestPseudotime.pdf')
print(FeaturePlot(GW18AllGE4K_2Dumap.integrated, features = 'GEH_C3_pseudotime', pt.size = 0.4))
dev.off()



# for time, only look at the 100 most variable genes
#Y <- log1p(assays(GW18AllGE_GEH_C3_SCE)$norm)
#var100 <- names(sort(apply(Y,1,var),decreasing = TRUE))[1:100]
#Y <- Y[var100,]
GEH_C3.Y = GW18AllGE_GEH_C3_SCE@assays@data[[1]]

plan(multiprocess)

# fit a GAM with a loess term for pseudotime
GEH_C3.pval <- future_apply(GEH_C3.Y,1,function(z){
    d <- data.frame(z=z, t=GEH_C3.t)
    suppressWarnings({
      tmp <- suppressWarnings(gam(z ~ lo(GEH_C3.t), data=d))
    p <- summary(tmp)[3][[1]][2,3]
    p
})
    })
## multithread did not work 


##FDR? 

GEH_C3.fdr = p.adjust(GEH_C3.pval,method='fdr')

GEH_C3.cor = NA

for(i in 1:nrow(GEH_C3.Y)){
GEH_C3.cor[i] = cor(GEH_C3.t,GEH_C3.Y[i,])
if(i%%100==0) cat(paste(i,', ',sep=''))
}

names(GEH_C3.cor) = names(GEH_C3.fdr)


GEH_C3res=data.frame(gene=names(GEH_C3.fdr),pval=GEH_C3.pval,fdr=GEH_C3.fdr,cor=GEH_C3.cor)

rownames(GEH_C3res)[which(GEH_C3res$fdr<0.05 & GEH_C3res$cor>0.5)]
"TUBA1A" "STMN2"  "MEF2C"  "SYT4"   "TUBB"   "SLA" 


pdf('./TestPlots/GW18AllGE4K_2D_GEH_C3_TestPseudotime.pdf')
print(FeaturePlot(GW18AllGE4K_2Dumap.integrated, features = 'GEH_C3_pseudotime', pt.size = 0.4))
for(i in c(rownames(GEH_C3res)[which(GEH_C3res$fdr<0.05 & GEH_C3res$cor>0.5)],rownames(GEH_C3res)[which(GEH_C3res$fdr<0.05 & GEH_C3res$cor<(-0.5))])){
print(FeaturePlot(GW18AllGE4K_2Dumap.integrated, features = i, pt.size = 0.4,order=TRUE))
}
dev.off()

GEH_C3.genes = c(rownames(GEH_C3res)[which(GEH_C3res$fdr<0.05 & GEH_C3res$cor>0.5)],rownames(GEH_C3res)[which(GEH_C3res$fdr<0.05 & GEH_C3res$cor<(-0.5))])

GEH_C3.data=cbind(GEH_C3.t,t(GEH_C3.Y[GEH_C3.genes,]))
GEH_C3.data=as.data.frame(GEH_C3.data)
colnames(GEH_C3.data) = c('pseudotime',GEH_C3.genes)

GEH_C3.data$cluster = GW18AllGE4K_2Dumap.integrated@meta.data[colnames(GW18AllGE_GEH_C3_SCE),'seurat_clusters']

GEH_C3.dataM = reshape2::melt(data=GEH_C3.data, id.vars = "pseudotime", measure.vars = GEH_C3.genes)

pdf('./TestPlots/GW18AllGE4K_2D_GEH_C3_TestPseudotime_SmoothLines.pdf',width=12,height=8)
print(ggplot(GEH_C3.dataM, aes(x = pseudotime, y = value, colour = variable)) + geom_smooth())
dev.off()



GEH_C3.clusts = unique(GEH_C3.data$cluster)[c(3,1,2)]

GEH_C3.mins = tapply(GEH_C3.data$pseudotime,INDEX=GEH_C3.data$cluster,FUN=min)[GEH_C3.clusts]
GEH_C3.maxes = tapply(GEH_C3.data$pseudotime,INDEX=GEH_C3.data$cluster,FUN=max)[GEH_C3.clusts]

GEH_C3.d=data.frame(x1=GEH_C3.mins, x2=GEH_C3.maxes, y1=seq(from=(-0.2),to=(-0.225),by=-0.0125), y2=seq(from=(-0.2125),to=(-0.2375),by=-0.0125), t=GEH_C3.clusts, r=c(1:length(GEH_C3.clusts)))

if(FALSE){ ## error
pdf('./TestPlots/GW18AllGE4K_2D_GEH_C3_TestPseudotime_SmoothLines_clusts.pdf',width=12,height=8)
print(ggplot(data=dataM, aes(x = pseudotime, y = value, colour = variable)) + geom_smooth() +  geom_rect(data=d, inherit.aes = FALSE, aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2, fill=t), color="black", alpha=0.5) + geom_text(data=d, aes(x=x1+(x2-x1)/2, y=y1+(y2-y1)/2, label=t), size=4))
dev.off()
}

pdf('./TestPlots/GW18AllGE4K_2D_GEH_C3_TestPseudotime_SmoothLines_clusts.pdf',width=12,height=8)
print(ggplot(data=GEH_C3.dataM, aes(x = pseudotime, y = value, colour = variable)) + geom_smooth() +  annotate("text", x = rowMeans(GEH_C3.d[,1:2]), y = -0.2, label = GEH_C3.d$t))
dev.off()



## C2 

GW18AllGE_GEH_C2 = subset(GW18AllGE4K.integrated,idents=unique(GW18AllGE4K.integrated@active.ident)[c(14,1,5)])

GW18AllGE_GEH_C2 = subset(GW18AllGE_GEH_C2,cells = colnames(GW18AllGE_GEH_C2)[which(GW18AllGE_GEH_C2@meta.data$BrainRegion!='Cortex')])

GW18AllGE_GEH_C2_SCE = Seurat::as.SingleCellExperiment(GW18AllGE_GEH_C2) ## just 2000 var genes - try for now? - transfered PCA

GW18AllGE_GEH_C2_SCE <- slingshot(GW18AllGE_GEH_C2_SCE, clusterLabels = 'seurat_clusters', reducedDim = 'PCA',start.clus = 'clust_14', end.clus = 'clust_2')

GEH_C2.t <- GW18AllGE_GEH_C2_SCE$slingPseudotime_1

GW18AllGE4K_2Dumap.integrated@meta.data$GEH_C2_pseudotime = 0

GW18AllGE4K_2Dumap.integrated@meta.data[colnames(GW18AllGE_GEH_C2_SCE),'GEH_C2_pseudotime'] = GEH_C2.t

pdf('./TestPlots/GW18AllGE4K_2D_GEH_C2_TestPseudotime.pdf')
print(FeaturePlot(GW18AllGE4K_2Dumap.integrated, features = 'GEH_C2_pseudotime', pt.size = 0.4))
dev.off()

GEH_C2.Y = GW18AllGE_GEH_C2_SCE@assays@data[[1]]

# fit a GAM with a loess term for pseudotime
GEH_C2.pval <- apply(GEH_C2.Y,1,function(z){
    d <- data.frame(z=z, t=GEH_C2.t)
    suppressWarnings({
      tmp <- suppressWarnings(gam(z ~ lo(GEH_C2.t), data=d))
    p <- summary(tmp)[3][[1]][2,3]
    p
})
    })

GEH_C2.fdr = p.adjust(GEH_C2.pval,method='fdr')

GEH_C2.cor = NA

for(i in 1:nrow(GEH_C2.Y)){
GEH_C2.cor[i] = cor(GEH_C2.t,GEH_C2.Y[i,])
if(i%%100==0) cat(paste(i,', ',sep=''))
}

names(GEH_C2.cor) = names(GEH_C2.fdr)

GEH_C2res=data.frame(gene=names(GEH_C2.fdr),pval=GEH_C2.pval,fdr=GEH_C2.fdr,cor=GEH_C2.cor)

rownames(GEH_C2res)[which(GEH_C2res$fdr<0.05 & GEH_C2res$cor>0.5)]
"TUBA1A" "STMN2"  "MEF2C"  "SYT4"   "TUBB"   "SLA" 


pdf('./TestPlots/GW18AllGE4K_2D_GEH_C2_TestPseudotime.pdf')
print(FeaturePlot(GW18AllGE4K_2Dumap.integrated, features = 'GEH_C2_pseudotime', pt.size = 0.4))
for(i in c(rownames(GEH_C2res)[which(GEH_C2res$fdr<0.05 & GEH_C2res$cor>0.5)],rownames(GEH_C2res)[which(GEH_C2res$fdr<0.05 & GEH_C2res$cor<(-0.5))])){
print(FeaturePlot(GW18AllGE4K_2Dumap.integrated, features = i, pt.size = 0.4,order=TRUE))
}
dev.off()

GEH_C2.genes = c(rownames(GEH_C2res)[which(GEH_C2res$fdr<0.05 & GEH_C2res$cor>0.5)],rownames(GEH_C2res)[which(GEH_C2res$fdr<0.05 & GEH_C2res$cor<(-0.5))])

GEH_C2.data=cbind(GEH_C2.t,t(GEH_C2.Y[GEH_C2.genes,]))
GEH_C2.data=as.data.frame(GEH_C2.data)
colnames(GEH_C2.data) = c('pseudotime',GEH_C2.genes)

GEH_C2.data$cluster = GW18AllGE4K_2Dumap.integrated@meta.data[colnames(GW18AllGE_GEH_C2_SCE),'seurat_clusters']

GEH_C2.dataM = reshape2::melt(data=GEH_C2.data, id.vars = "pseudotime", measure.vars = GEH_C2.genes)

GEH_C2.clusts = unique(GEH_C2.data$cluster)[c(3,1,2)]

GEH_C2.mins = tapply(GEH_C2.data$pseudotime,INDEX=GEH_C2.data$cluster,FUN=min)[GEH_C2.clusts]
GEH_C2.maxes = tapply(GEH_C2.data$pseudotime,INDEX=GEH_C2.data$cluster,FUN=max)[GEH_C2.clusts]

GEH_C2.d=data.frame(x1=GEH_C2.mins, x2=GEH_C2.maxes, y1=seq(from=(-0.2),to=(-0.225),by=-0.0125), y2=seq(from=(-0.2125),to=(-0.2375),by=-0.0125), t=GEH_C2.clusts, r=c(1:length(GEH_C2.clusts)))

pdf('./TestPlots/GW18AllGE4K_2D_GEH_C2_TestPseudotime_SmoothLines_clusts.pdf',width=12,height=8)
print(ggplot(data=GEH_C2.dataM, aes(x = pseudotime, y = value, colour = variable)) + geom_smooth() +  annotate("text", x = rowMeans(GEH_C2.d[,1:2]), y = -0.2, label = GEH_C2.d$t))
dev.off()

### try plotting two genes that typify lineage 


pdf('./TestPlots/GW18AllGE4K_2D_GEH_C2_TestPseudotime_PTTG1_STMN2.pdf')

print(FeaturePlot(object = GW18AllGE4K_2Dumap.integrated,features = c("PTTG1", "STMN2"),cols  = c("grey", "blue", "red", "pink"),reduction = "umap", blend = TRUE))
	dev.off()

pdf('./TestPlots/GW18AllGE4K_2D_GEH_C2_TestPseudotime_PTTG1_STMN2.pdf')

print(FeaturePlot(object = GW18AllGE4K_2Dumap.integrated,features = "PTTG1",reduction = "umap"))
print(FeaturePlot(object = GW18AllGE4K_2Dumap.integrated,features =  "STMN2",reduction = "umap"))
	dev.off()




## C14 - C0 - C9 - C2 

GW18AllGE_GEH_C9C2 = subset(GW18AllGE4K.integrated,idents=unique(GW18AllGE4K.integrated@active.ident)[c(14,1,11,5)])

GW18AllGE_GEH_C9C2 = subset(GW18AllGE_GEH_C9C2,cells = colnames(GW18AllGE_GEH_C9C2)[which(GW18AllGE_GEH_C9C2@meta.data$BrainRegion!='Cortex')])

GW18AllGE_GEH_C9C2_SCE = Seurat::as.SingleCellExperiment(GW18AllGE_GEH_C9C2) ## just 2000 var genes - try for now? - transfered PCA

GW18AllGE_GEH_C9C2_SCE <- slingshot(GW18AllGE_GEH_C9C2_SCE, clusterLabels = 'seurat_clusters', reducedDim = 'PCA',start.clus = 'clust_14', end.clus = 'clust_2')

GEH_C9C2.t <- GW18AllGE_GEH_C9C2_SCE$slingPseudotime_1

GW18AllGE4K_2Dumap.integrated@meta.data$GEH_C9C2_pseudotime = 0

GW18AllGE4K_2Dumap.integrated@meta.data[colnames(GW18AllGE_GEH_C9C2_SCE),'GEH_C9C2_pseudotime'] = GEH_C9C2.t

pdf('./TestPlots/GW18AllGE4K_2D_GEH_C9C2_TestPseudotime.pdf')
print(FeaturePlot(GW18AllGE4K_2Dumap.integrated, features = 'GEH_C9C2_pseudotime', pt.size = 0.4))
dev.off()

GEH_C9C2.Y = GW18AllGE_GEH_C9C2_SCE@assays@data[[1]]

# fit a GAM with a loess term for pseudotime
GEH_C9C2.pval <- apply(GEH_C9C2.Y,1,function(z){
    d <- data.frame(z=z, t=GEH_C9C2.t)
    suppressWarnings({
      tmp <- suppressWarnings(gam(z ~ lo(GEH_C9C2.t), data=d))
    p <- summary(tmp)[3][[1]][2,3]
    p
})
    })

GEH_C9C2.fdr = p.adjust(GEH_C9C2.pval,method='fdr')

GEH_C9C2.cor = NA

for(i in 1:nrow(GEH_C9C2.Y)){
GEH_C9C2.cor[i] = cor(GEH_C9C2.t,GEH_C9C2.Y[i,])
if(i%%100==0) cat(paste(i,', ',sep=''))
}

names(GEH_C9C2.cor) = names(GEH_C9C2.fdr)

GEH_C9C2res=data.frame(gene=names(GEH_C9C2.fdr),pval=GEH_C9C2.pval,fdr=GEH_C9C2.fdr,cor=GEH_C9C2.cor)

rownames(GEH_C9C2res)[which(GEH_C9C2res$fdr<0.05 & GEH_C9C2res$cor>0.5)]
"TUBA1A" "STMN2"  "MEF2C"  "SYT4"   "TUBB"   "SLA" 


pdf('./TestPlots/GW18AllGE4K_2D_GEH_C9C2_TestPseudotime.pdf')
print(FeaturePlot(GW18AllGE4K_2Dumap.integrated, features = 'GEH_C9C2_pseudotime', pt.size = 0.4))
for(i in c(rownames(GEH_C9C2res)[which(GEH_C9C2res$fdr<0.05 & GEH_C9C2res$cor>0.5)],rownames(GEH_C9C2res)[which(GEH_C9C2res$fdr<0.05 & GEH_C9C2res$cor<(-0.5))])){
print(FeaturePlot(GW18AllGE4K_2Dumap.integrated, features = i, pt.size = 0.4,order=TRUE))
}
dev.off()

GEH_C9C2.genes = c(rownames(GEH_C9C2res)[which(GEH_C9C2res$fdr<0.05 & GEH_C9C2res$cor>0.5)],rownames(GEH_C9C2res)[which(GEH_C9C2res$fdr<0.05 & GEH_C9C2res$cor<(-0.5))])

GEH_C9C2.data=cbind(GEH_C9C2.t,t(GEH_C9C2.Y[GEH_C9C2.genes,]))
GEH_C9C2.data=as.data.frame(GEH_C9C2.data)
colnames(GEH_C9C2.data) = c('pseudotime',GEH_C9C2.genes)

GEH_C9C2.data$cluster = GW18AllGE4K_2Dumap.integrated@meta.data[colnames(GW18AllGE_GEH_C9C2_SCE),'seurat_clusters']

GEH_C9C2.dataM = reshape2::melt(data=GEH_C9C2.data, id.vars = "pseudotime", measure.vars = GEH_C9C2.genes)

GEH_C9C2.clusts = unique(GEH_C9C2.data$cluster)[c(4,1,3,2)]

GEH_C9C2.mins = tapply(GEH_C9C2.data$pseudotime,INDEX=GEH_C9C2.data$cluster,FUN=min)[GEH_C9C2.clusts]
GEH_C9C2.maxes = tapply(GEH_C9C2.data$pseudotime,INDEX=GEH_C9C2.data$cluster,FUN=max)[GEH_C9C2.clusts]

GEH_C9C2.d=data.frame(x1=GEH_C9C2.mins, x2=GEH_C9C2.maxes, y1=seq(from=(-0.2),to=(-0.2375),by=-0.0125), y2=seq(from=(-0.2125),to=(-0.25),by=-0.0125), t=GEH_C9C2.clusts, r=c(1:length(GEH_C9C2.clusts)))

pdf('./TestPlots/GW18AllGE4K_2D_GEH_C9C2_TestPseudotime_SmoothLines_clusts.pdf',width=12,height=8)
print(ggplot(data=GEH_C9C2.dataM, aes(x = pseudotime, y = value, colour = variable)) + geom_smooth() +  annotate("text", x = rowMeans(GEH_C9C2.d[,1:2]), y = -0.2, label = GEH_C9C2.d$t))
dev.off()


### hypo and GE separate? repeat C9C2

## C2C9 GE 


GW18AllGE_GE_C9C2 = subset(GW18AllGE4K.integrated,idents=unique(GW18AllGE4K.integrated@active.ident)[c(14,1,11,5)])

GW18AllGE_GE_C9C2 = subset(GW18AllGE_GE_C9C2,cells = colnames(GW18AllGE_GE_C9C2)[which(GW18AllGE_GE_C9C2@meta.data$BrainRegion=='GE')])

GW18AllGE_GE_C9C2_SCE = Seurat::as.SingleCellExperiment(GW18AllGE_GE_C9C2) ## just 2000 var genes - try for now? - transfered PCA

GW18AllGE_GE_C9C2_SCE <- slingshot(GW18AllGE_GE_C9C2_SCE, clusterLabels = 'seurat_clusters', reducedDim = 'PCA',start.clus = 'clust_14', end.clus = 'clust_2')

GE_C9C2.t <- GW18AllGE_GE_C9C2_SCE$slingPseudotime_1

GW18AllGE4K_2Dumap.integrated@meta.data$GE_C9C2_pseudotime = 0

GW18AllGE4K_2Dumap.integrated@meta.data[colnames(GW18AllGE_GE_C9C2_SCE),'GE_C9C2_pseudotime'] = GE_C9C2.t

pdf('./TestPlots/GW18AllGE4K_2D_GE_C9C2_TestPseudotime.pdf')
print(FeaturePlot(GW18AllGE4K_2Dumap.integrated, features = 'GE_C9C2_pseudotime', pt.size = 0.4))
dev.off()

GE_C9C2.Y = GW18AllGE_GE_C9C2_SCE@assays@data[[1]]

# fit a GAM with a loess term for pseudotime
GE_C9C2.pval <- apply(GE_C9C2.Y,1,function(z){
    d <- data.frame(z=z, t=GE_C9C2.t)
    suppressWarnings({
      tmp <- suppressWarnings(gam(z ~ lo(GE_C9C2.t), data=d))
    p <- summary(tmp)[3][[1]][2,3]
    p
})
    })

GE_C9C2.fdr = p.adjust(GE_C9C2.pval,method='fdr')

GE_C9C2.cor = NA

for(i in 1:nrow(GE_C9C2.Y)){
GE_C9C2.cor[i] = cor(GE_C9C2.t,GE_C9C2.Y[i,])
if(i%%100==0) cat(paste(i,', ',sep=''))
}

names(GE_C9C2.cor) = names(GE_C9C2.fdr)

GE_C9C2res=data.frame(gene=names(GE_C9C2.fdr),pval=GE_C9C2.pval,fdr=GE_C9C2.fdr,cor=GE_C9C2.cor)

rownames(GE_C9C2res)[which(GE_C9C2res$fdr<0.05 & GE_C9C2res$cor>0.5)]
"TUBA1A" "STMN2"  "MEF2C"  "SYT4"   "TUBB"   "SLA" 


pdf('./TestPlots/GW18AllGE4K_2D_GE_C9C2_TestPseudotime.pdf')
print(FeaturePlot(GW18AllGE4K_2Dumap.integrated, features = 'GE_C9C2_pseudotime', pt.size = 0.4))
for(i in c(rownames(GE_C9C2res)[which(GE_C9C2res$fdr<0.05 & GE_C9C2res$cor>0.5)],rownames(GE_C9C2res)[which(GE_C9C2res$fdr<0.05 & GE_C9C2res$cor<(-0.5))])){
print(FeaturePlot(GW18AllGE4K_2Dumap.integrated, features = i, pt.size = 0.4,order=TRUE))
}
dev.off()

GE_C9C2.genes = c(rownames(GE_C9C2res)[which(GE_C9C2res$fdr<0.05 & GE_C9C2res$cor>0.6)],rownames(GE_C9C2res)[which(GE_C9C2res$fdr<0.06 & GE_C9C2res$cor<(-0.5))])

GE_C9C2.data=cbind(GE_C9C2.t,t(GE_C9C2.Y[GE_C9C2.genes,]))
GE_C9C2.data=as.data.frame(GE_C9C2.data)
colnames(GE_C9C2.data) = c('pseudotime',GE_C9C2.genes)

GE_C9C2.data$cluster = GW18AllGE4K_2Dumap.integrated@meta.data[colnames(GW18AllGE_GE_C9C2_SCE),'seurat_clusters']

GE_C9C2.dataM = reshape2::melt(data=GE_C9C2.data, id.vars = "pseudotime", measure.vars = GE_C9C2.genes)

GE_C9C2.clusts = unique(GE_C9C2.data$cluster)[c(4,1,3,2)]

GE_C9C2.mins = tapply(GE_C9C2.data$pseudotime,INDEX=GE_C9C2.data$cluster,FUN=min)[GE_C9C2.clusts]
GE_C9C2.maxes = tapply(GE_C9C2.data$pseudotime,INDEX=GE_C9C2.data$cluster,FUN=max)[GE_C9C2.clusts]

GE_C9C2.d=data.frame(x1=GE_C9C2.mins, x2=GE_C9C2.maxes, y1=seq(from=(-0.2),to=(-0.2375),by=-0.0125), y2=seq(from=(-0.2125),to=(-0.25),by=-0.0125), t=GE_C9C2.clusts, r=c(1:length(GE_C9C2.clusts)))

pdf('./TestPlots/GW18AllGE4K_2D_GE_C9C2_TestPseudotime_SmoothLines_clusts.pdf',width=12,height=8)
print(ggplot(data=GE_C9C2.dataM, aes(x = pseudotime, y = value, colour = variable)) + geom_smooth() +  annotate("text", x = rowMeans(GE_C9C2.d[,1:2]), y = -0.2, label = GE_C9C2.d$t))
dev.off()


## C2C9 Hypo 


GW18AllGE_Hypo_C9C2 = subset(GW18AllGE4K.integrated,idents=unique(GW18AllGE4K.integrated@active.ident)[c(14,1,11,5)])

GW18AllGE_Hypo_C9C2 = subset(GW18AllGE_Hypo_C9C2,cells = colnames(GW18AllGE_Hypo_C9C2)[which(GW18AllGE_Hypo_C9C2@meta.data$BrainRegion=='Hypo')])

GW18AllGE_Hypo_C9C2_SCE = Seurat::as.SingleCellExperiment(GW18AllGE_Hypo_C9C2) ## just 2000 var genes - try for now? - transfered PCA

GW18AllGE_Hypo_C9C2_SCE <- slingshot(GW18AllGE_Hypo_C9C2_SCE, clusterLabels = 'seurat_clusters', reducedDim = 'PCA',start.clus = 'clust_14', end.clus = 'clust_2')

Hypo_C9C2.t <- GW18AllGE_Hypo_C9C2_SCE$slingPseudotime_1

GW18AllGE4K_2Dumap.integrated@meta.data$Hypo_C9C2_pseudotime = 0

GW18AllGE4K_2Dumap.integrated@meta.data[colnames(GW18AllGE_Hypo_C9C2_SCE),'Hypo_C9C2_pseudotime'] = Hypo_C9C2.t

pdf('./TestPlots/GW18AllGE4K_2D_Hypo_C9C2_TestPseudotime.pdf')
print(FeaturePlot(GW18AllGE4K_2Dumap.integrated, features = 'Hypo_C9C2_pseudotime', pt.size = 0.4,order=TRUE))
dev.off()

Hypo_C9C2.Y = GW18AllGE_Hypo_C9C2_SCE@assays@data[[1]]

# fit a GAM with a loess term for pseudotime
Hypo_C9C2.pval <- apply(Hypo_C9C2.Y,1,function(z){
    d <- data.frame(z=z, t=Hypo_C9C2.t)
    suppressWarnings({
      tmp <- suppressWarnings(gam(z ~ lo(Hypo_C9C2.t), data=d))
    p <- summary(tmp)[3][[1]][2,3]
    p
})
    })

Hypo_C9C2.fdr = p.adjust(Hypo_C9C2.pval,method='fdr')

Hypo_C9C2.cor = NA

for(i in 1:nrow(Hypo_C9C2.Y)){
Hypo_C9C2.cor[i] = cor(Hypo_C9C2.t,Hypo_C9C2.Y[i,])
if(i%%100==0) cat(paste(i,', ',sep=''))
}

names(Hypo_C9C2.cor) = names(Hypo_C9C2.fdr)

Hypo_C9C2res=data.frame(gene=names(Hypo_C9C2.fdr),pval=Hypo_C9C2.pval,fdr=Hypo_C9C2.fdr,cor=Hypo_C9C2.cor)

rownames(Hypo_C9C2res)[which(Hypo_C9C2res$fdr<0.05 & Hypo_C9C2res$cor>0.5)]
"TUBA1A" "STMN2"  "MEF2C"  "SYT4"   "TUBB"   "SLA" 


pdf('./TestPlots/GW18AllGE4K_2D_Hypo_C9C2_TestPseudotime.pdf')
print(FeaturePlot(GW18AllGE4K_2Dumap.integrated, features = 'Hypo_C9C2_pseudotime', pt.size = 0.4))
for(i in c(rownames(Hypo_C9C2res)[which(Hypo_C9C2res$fdr<0.05 & Hypo_C9C2res$cor>0.5)],rownames(Hypo_C9C2res)[which(Hypo_C9C2res$fdr<0.05 & Hypo_C9C2res$cor<(-0.5))])){
print(FeaturePlot(GW18AllGE4K_2Dumap.integrated, features = i, pt.size = 0.4,order=TRUE))
}
dev.off()

Hypo_C9C2.genes = c(rownames(Hypo_C9C2res)[which(Hypo_C9C2res$fdr<0.05 & Hypo_C9C2res$cor>0.5)],rownames(Hypo_C9C2res)[which(Hypo_C9C2res$fdr<0.05 & Hypo_C9C2res$cor<(-0.5))])

Hypo_C9C2.data=cbind(Hypo_C9C2.t,t(Hypo_C9C2.Y[Hypo_C9C2.genes,]))
Hypo_C9C2.data=as.data.frame(Hypo_C9C2.data)
colnames(Hypo_C9C2.data) = c('pseudotime',Hypo_C9C2.genes)

Hypo_C9C2.data$cluster = GW18AllGE4K_2Dumap.integrated@meta.data[colnames(GW18AllGE_Hypo_C9C2_SCE),'seurat_clusters']

Hypo_C9C2.dataM = reshape2::melt(data=Hypo_C9C2.data, id.vars = "pseudotime", measure.vars = Hypo_C9C2.genes)

Hypo_C9C2.clusts = unique(Hypo_C9C2.data$cluster)[c(4,1,3,2)]

Hypo_C9C2.mins = tapply(Hypo_C9C2.data$pseudotime,INDEX=Hypo_C9C2.data$cluster,FUN=min)[Hypo_C9C2.clusts]
Hypo_C9C2.maxes = tapply(Hypo_C9C2.data$pseudotime,INDEX=Hypo_C9C2.data$cluster,FUN=max)[Hypo_C9C2.clusts]

Hypo_C9C2.d=data.frame(x1=Hypo_C9C2.mins, x2=Hypo_C9C2.maxes, y1=seq(from=(-0.2),to=(-0.2375),by=-0.0125), y2=seq(from=(-0.2125),to=(-0.25),by=-0.0125), t=Hypo_C9C2.clusts, r=c(1:length(Hypo_C9C2.clusts)))

pdf('./TestPlots/GW18AllGE4K_2D_Hypo_C9C2_TestPseudotime_SmoothLines_clusts.pdf',width=12,height=8)
print(ggplot(data=Hypo_C9C2.dataM, aes(x = pseudotime, y = value, colour = variable)) + geom_smooth() +  annotate("text", x = rowMeans(Hypo_C9C2.d[,1:2]), y = -0.2, label = Hypo_C9C2.d$t))
dev.off()



### top gene dif? 













save(Hypo_C9C2res,GE_C9C2res,GEH_C9C2res,GEH_C2res,GEH_C3res,Ctx.ITres,file='./Analysis/GW18AllGE4K_Pseudotime_Var_Gene_Cor.rda',compress=TRUE)




### heatmaps of pseudotime genes 

## bin by pseudotime? 



ComGenes = intersect(Ctx.ITgenes,intersect(GE_C9C2.genes,Hypo_C9C2.genes))

JointGenes = c(Ctx.ITgenes,GE_C9C2.genes,Hypo_C9C2.genes)


### need to pull out data for genes 

### translate pseudotime into quartiles? - add as index 

CtxQuant = as.numeric(cut(Ctx.ITdata$pseudotime,quantile(Ctx.ITdata$pseudotime, probs=seq(from=0, to=1,by=0.05)), include.lowest=TRUE))

GEQuant = as.numeric(cut(GE_C9C2.data$pseudotime,quantile(GE_C9C2.data$pseudotime, probs=seq(from=0, to=1,by=0.05)), include.lowest=TRUE))

HypoQuant = as.numeric(cut(Hypo_C9C2.data$pseudotime,quantile(Hypo_C9C2.data$pseudotime, probs=seq(from=0, to=1,by=0.05)), include.lowest=TRUE))

for(i in 1:length(JointGenes)){
if(i==1){
ComDatCtx = tapply(X=as.numeric(t(CtxIT.Y[JointGenes[i],])),INDEX=CtxQuant,FUN=mean)
ComDatGE = tapply(X=as.numeric(t(GE_C9C2.Y[JointGenes[i],])),INDEX=GEQuant,FUN=mean)
ComDatHypo = tapply(X=as.numeric(t(Hypo_C9C2.Y[JointGenes[i],])),INDEX=HypoQuant,FUN=mean)
} else {
ComDatCtx = rbind(ComDatCtx,tapply(X=as.numeric(t(CtxIT.Y[JointGenes[i],])),INDEX=CtxQuant,FUN=mean))
ComDatGE = rbind(ComDatGE,tapply(X=as.numeric(t(GE_C9C2.Y[JointGenes[i],])),INDEX=GEQuant,FUN=mean))
ComDatHypo = rbind(ComDatHypo,tapply(X=as.numeric(t(Hypo_C9C2.Y[JointGenes[i],])),INDEX=HypoQuant,FUN=mean))
}
}


ComDat = rbind(ComDatCtx,ComDatGE,ComDatHypo)

pdf('./TestPlots/GW18AllGE4K_ExcitLineage_Heatmap.pdf',width=12,height=8)

heatmap(ComDat,Rowv=NA,Colv=NA)

dev.off()



### pairwise comparisons 


comps = data.frame(AA = c(2,2,3,2,3,9,3,3,3,2,2,2,2,2,2,1,1,4,9,10,7,6,0), BB=c(3,9,9,10,10,10,15,17,21,8,13,26,27,4,11,2,9,9,11,19,10,10,10),stringsAsFactors=FALSE)


DefaultAssay(GW18AllGE4K_2Dumap.integrated)='RNA'

GW18AllGE4K_2Dumap.integrated@active.ident=as.factor(GW18AllGE4K_2Dumap.integrated$seurat_clusters)


clustComps = vector(mode='list',length=nrow(comps))


for(k in 1:nrow(comps)){

id1 = paste('clust_',comps[k,1],sep='')
id2 = paste('clust_',comps[k,2],sep='')


tmpDEG = FindMarkers(GW18AllGE4K_2Dumap.integrated,ident.1=id1,ident.2=id2)

tmpDEG$gene= rownames(tmpDEG)

clustComps[[k]] = tmpDEG

names(clustComps)[k] = paste(id1,'_VS_',id2,sep='')
}

save(clustComps,file='/Analysis/GW18AllGE4K_2Dumap_Cluster_DEGs.rda',compress=TRUE)

for(k in 1:nrow(comps)){
id1 = paste('clust_',comps[k,1],sep='')
id2 = paste('clust_',comps[k,2],sep='')
tmpDEG=clustComps[[k]]
tmp=tmpDEG
tmp = tmp[which(tmp$p_val_adj<=0.01 & tmp$avg_logFC>0 & tmp$pct.2<0.1),]
tmp$perDif = tmp$pct.1-tmp$pct.2
tmp = tmp[order(tmp$perDif,decreasing=TRUE),]
topPerDif1 <-  as.character(na.omit(tmp$gene[1:5]))
tmp=tmpDEG
tmp = tmp[which(tmp$p_val_adj<=0.01 & tmp$avg_logFC<0 & tmp$pct.1<0.1),]
tmp$perDif = tmp$pct.2-tmp$pct.1
tmp = tmp[order(tmp$perDif,decreasing=TRUE),]
topPerDif2 <-  as.character(na.omit(tmp$gene[1:5]))
topPerDif = c(topPerDif1,topPerDif2)

if(length(topPerDif)>0){
pdf(paste('./TestPlots/GW18AllGE4K_2Dumap_',id1,'_VS_',id2,'.pdf',sep=''),width=12,height=8)
for(j in c("sample","seurat_clusters")){
print(DimPlot(GW18AllGE4K_2Dumap.integrated, reduction = "umap", group.by = j, label = TRUE, repel = TRUE))
}
for(g in topPerDif){
print(FeaturePlot(GW18AllGE4K_2Dumap.integrated,features=g, pt.size = 0.4))
}
dev.off()
}
cat(paste(id1,'_VS_',id2,'\n',sep=''))
}




DefaultAssay(GW18AllGE4K.integrated) = 'integrated'

GW18AllGE_Tmp = subset(GW18AllGE4K.integrated,idents=j)

pdf(paste('./TestPlots/GW18AllGE4K_',j,'_topVAR_heatmap.pdf',sep=''),width=12,height=8)

print(DoHeatmap(GW18AllGE_Tmp,features = topPerDif,group.by = 'sample',label=FALSE,raster=FALSE) + NoLegend())

dev.off()

cat(paste('\n\n',j,' Done','\n\n',sep=''))

}

### focus on what makes cluster 2 special / why it is there 


GW18AllGE4K_2Dumap.integrated@meta.data$Myst = 'Out'

GW18AllGE4K_2Dumap.integrated@meta.data$Myst[which(!is.na(match(GW18AllGE4K_2Dumap.integrated@meta.data$seurat_clusters,c('clust_2','clust_26','clust_13','clust_8','clust_27'))))] = 'In'


pdf('./TestPlots/GW18AllGE4K_2Dumap_Check_Mystery_region.pdf',width=12,height=8)
for(j in c("sample","seurat_clusters","Myst")){
print(DimPlot(GW18AllGE4K_2Dumap.integrated, reduction = "umap", group.by = j, label = TRUE, repel = TRUE))
}
dev.off()


DefaultAssay(GW18AllGE4K_2Dumap.integrated)='RNA'

GW18AllGE4K_2Dumap.integrated@active.ident=as.factor(GW18AllGE4K_2Dumap.integrated$Myst)

MystDEG = FindMarkers(GW18AllGE4K_2Dumap.integrated,ident.1='In',ident.2='Out')

MystDEGsig = MystDEG[which(MystDEG$p_val_adj<0.001),]

MystDEGsig = MystDEGsig[order(MystDEGsig$avg_logFC,decreasing=TRUE),]

rownames(MystDEGsig[which(MystDEGsig$pct.2<0.1),])[1:25]

pdf('./TestPlots/GW18AllGE4K_2Dumap_High_Myst_genes.pdf',width=12,height=8)
for(j in c("sample","seurat_clusters")){
print(DimPlot(GW18AllGE4K_2Dumap.integrated, reduction = "umap", group.by = j, label = TRUE, repel = TRUE))
}
for(g in rownames(MystDEGsig[which(MystDEGsig$pct.2<0.1),])[1:25]){
print(FeaturePlot(GW18AllGE4K_2Dumap.integrated,features=g, pt.size = 0.4,order=TRUE))
}
dev.off()


MystDEGsig = MystDEGsig[order(MystDEGsig$avg_logFC,decreasing=FALSE),]


pdf('./TestPlots/GW18AllGE4K_2Dumap_Low_Myst_genes.pdf',width=12,height=8)
for(j in c("sample","seurat_clusters")){
print(DimPlot(GW18AllGE4K_2Dumap.integrated, reduction = "umap", group.by = j, label = TRUE, repel = TRUE))
}
for(g in rownames(MystDEGsig[which(MystDEGsig$pct.1<0.1),])[1:25]){
print(FeaturePlot(GW18AllGE4K_2Dumap.integrated,features=g, pt.size = 0.4))
}
dev.off()


GW18AllGE4K_2Dumap.integrated@meta.data$MajGroups = 'Out'

GW18AllGE4K_2Dumap.integrated@meta.data$MajGroups[which(!is.na(match(GW18AllGE4K_2Dumap.integrated@meta.data$seurat_clusters,c('clust_2','clust_26','clust_13','clust_8','clust_27'))))] = 'ExcET'

GW18AllGE4K_2Dumap.integrated@meta.data$MajGroups[which(!is.na(match(GW18AllGE4K_2Dumap.integrated@meta.data$seurat_clusters,c('clust_19','clust_7','clust_6','clust_0','clust_11','clust_4','clust_1','clust_10','clust_9'))))] = 'ExcIT'

GW18AllGE4K_2Dumap.integrated@meta.data$MajGroups[which(!is.na(match(GW18AllGE4K_2Dumap.integrated@meta.data$seurat_clusters,c('clust_3','clust_15','clust_17','clust_21'))))] = 'Inhib'



pdf('./TestPlots/GW18AllGE4K_2Dumap_Check_MajGroups_region.pdf',width=12,height=8)
for(j in c("sample","seurat_clusters","MajGroups")){
print(DimPlot(GW18AllGE4K_2Dumap.integrated, reduction = "umap", group.by = j, label = TRUE, repel = TRUE))
}
dev.off()


DefaultAssay(GW18AllGE4K_2Dumap.integrated)='RNA'

GW18AllGE4K_2Dumap.integrated@active.ident=as.factor(GW18AllGE4K_2Dumap.integrated$MajGroups)

MajGroupsDEG = FindMarkers(GW18AllGE4K_2Dumap.integrated,ident.1='ExcET',ident.2='Inhib')

MajGroupsDEGsig = MajGroupsDEG[which(MajGroupsDEG$p_val_adj<0.001),]

MajGroupsDEGsig = MajGroupsDEGsig[order(MajGroupsDEGsig$avg_logFC,decreasing=TRUE),]

rownames(MajGroupsDEGsig[which(MajGroupsDEGsig$pct.2<0.1),])[1:25]

pdf('./TestPlots/GW18AllGE4K_2Dumap_High_MajGroups_genes.pdf',width=12,height=8)
for(j in c("sample","seurat_clusters")){
print(DimPlot(GW18AllGE4K_2Dumap.integrated, reduction = "umap", group.by = j, label = TRUE, repel = TRUE))
}
for(g in rownames(MajGroupsDEGsig[which(MajGroupsDEGsig$pct.2<0.1),])[1:25]){
print(FeaturePlot(GW18AllGE4K_2Dumap.integrated,features=g, pt.size = 0.4,order=TRUE))
}
dev.off()


MajGroupsDEGsig = MajGroupsDEGsig[order(MajGroupsDEGsig$avg_logFC,decreasing=FALSE),]


pdf('./TestPlots/GW18AllGE4K_2Dumap_Low_MajGroups_genes.pdf',width=12,height=8)
for(j in c("sample","seurat_clusters")){
print(DimPlot(GW18AllGE4K_2Dumap.integrated, reduction = "umap", group.by = j, label = TRUE, repel = TRUE))
}
for(g in rownames(MajGroupsDEGsig[which(MajGroupsDEGsig$pct.1<0.1),])[1:25]){
print(FeaturePlot(GW18AllGE4K_2Dumap.integrated,features=g, pt.size = 0.4))
}
dev.off()


### checking in on a few genes

samples = unique(GW18AllGE4K_2Dumap.integrated@meta.data$sample)
gene='NEUROD6'
for(i in samples){
	cat(i)
	cat('\n')
	cat(table(GW18AllGE4K_2Dumap.integrated@assays$RNA[gene,][which(GW18AllGE4K_2Dumap.integrated@meta.data$sample==i)]))
		cat('\n')
}

gene='NKX2-1'
for(i in samples){
	cat(i)
	cat('\n')
	cat(table(GW18AllGE4K_2Dumap.integrated@assays$RNA[gene,][which(GW18AllGE4K_2Dumap.integrated@meta.data$sample==i & GW18AllGE4K_2Dumap.integrated@meta.data$seurat_clusters=='clust_16')]))
		cat('\n')
}





GW18AllGE4K_2Dumap.integrated@meta.data$Tissue = 'Hypo'
GW18AllGE4K_2Dumap.integrated@meta.data$Tissue[grep('GW18_',GW18AllGE4K_2Dumap.integrated@meta.data$sample)]='Cortex'
GW18AllGE4K_2Dumap.integrated@meta.data$Tissue[grep('GE',GW18AllGE4K_2Dumap.integrated@meta.data$sample)]='GE'


pdf('./TestPlots/GW18AllGE4K_2Dumap_GE_Hypo_Ctx_genes.pdf',width=12,height=8)
for(j in c("sample","seurat_clusters","Maj_Clust_From_Ref")){
print(DimPlot(GW18AllGE4K_2Dumap.integrated, reduction = "umap", group.by = j, label = TRUE, repel = TRUE))
}
for(g in c("PPP1R17","PPP2R2B","NRN1","NEUROD2","NEUROD6","SSTR2","RND3","SLA")){
print(FeaturePlot(GW18AllGE4K_2Dumap.integrated,features=g, pt.size = 0.4,split.by='Tissue',combine=TRUE,order=TRUE))
}
dev.off()

pdf('./TestPlots/GW18AllGE4K_2Dumap_GE_Maturation_Genes.pdf',width=12,height=8)
for(j in c("sample","seurat_clusters","Maj_Clust_From_Ref")){
print(DimPlot(GW18AllGE4K_2Dumap.integrated, reduction = "umap", group.by = j, label = TRUE, repel = TRUE))
}
for(g in c("NES","CCND2","ASCL1","DCX")){
print(FeaturePlot(GW18AllGE4K_2Dumap.integrated,features=g, pt.size = 0.4,split.by='Tissue',combine=TRUE,order=TRUE))
}
dev.off() ### these are very broad  - from Mayer 2018



GW18AllGE4K_2Dumap.integrated@active.ident=as.factor(GW18AllGE4K_2Dumap.integrated$seurat_clusters)

pdf('./TestPlots/GW18AllGE4K_2Dumap_GE_Hypo_Ctx_genes_VLN.pdf',width=12,height=8)

for(i in c("PPP1R17","PPP2R2B","NRN1","NEUROD2","NEUROD6","SSTR2","RND3","SLA")){
print(VlnPlot(GW18AllGE4K_2Dumap.integrated, features=i,split.by='Tissue'))
}
dev.off()

pdf('./TestPlots/GW18AllGE4K_2Dumap_GE_Hypo_Ctx_genes_DOT.pdf',width=12,height=8)

for(i in c("PPP1R17","PPP2R2B","NRN1","NEUROD2","NEUROD6","SSTR2","RND3","SLA")){
print(DotPlot(GW18AllGE4K_2Dumap.integrated, features=i,cols = c("red","green","blue"), group.by ='seurat_clusters' ,split.by="Tissue"))
}
dev.off()




## heatmap 

pdf('./TestPlots/GW18AllGE4K_2Dumap_GE_Hypo_Ctx_genes_heatmap.pdf',width=12,height=8)
DoHeatmap(GW18AllGE4K_2Dumap.integrated, features = "PPP1R17",group.by = "seurat_clusters",label=FALSE,raster=FALSE) + NoLegend()
dev.off()



### solidify layering in cortex result? 

pdf('./TestPlots/GW18AllGE4K_2Dumap_Projection_neuron_genes.pdf',width=12,height=8)
for(j in c("sample","seurat_clusters")){
print(DimPlot(GW18AllGE4K_2Dumap.integrated, reduction = "umap", group.by = j, label = TRUE, repel = TRUE))
}
for(g in c('OTOF','RSPO1','FEZF2','RXFP1','OSR1','BCL11B','TLE4','SATB2')){
print(FeaturePlot(GW18AllGE4K_2Dumap.integrated,features=g, pt.size = 0.4,order=TRUE))
}
dev.off()


### tissue DEGs per cluster 

GW18AllGE4K_2Dumap.integrated@active.ident=as.factor(GW18AllGE4K_2Dumap.integrated$Tissue)

DefaultAssay(GW18AllGE4K_2Dumap.integrated) = 'SCT'

clusts = unique(GW18AllGE4K_2Dumap.integrated@meta.data$seurat_clusters)

TissueClusterDEG = vector(mode='list',length=length(clusts))
names(TissueClusterDEG) = clusts

for(i in clusts){
 tmp = FindAllMarkers(subset(GW18AllGE4K_2Dumap.integrated,cells=colnames(GW18AllGE4K_2Dumap.integrated)[which(GW18AllGE4K_2Dumap.integrated@meta.data$seurat_clusters==i)]))

 tmp$seurat_cluster=i
TissueClusterDEG[[i]] = tmp
cat(i)
}

save(TissueClusterDEG,file='./Analysis/TissueClusterDEG_HypoGEcrtx.rda',compress=TRUE) # maj region dif within cluster - uses all genes by SCT 


### major tissue markers of each cluster (have this?)



### cluster  DEGs per tissue

GW18AllGE4K_2Dumap.integrated@active.ident=as.factor(GW18AllGE4K_2Dumap.integrated$seurat_clusters)

DefaultAssay(GW18AllGE4K_2Dumap.integrated) = 'SCT'

tissues = unique(GW18AllGE4K_2Dumap.integrated@meta.data$Tissue)

ClusterDEGbyTissue = vector(mode='list',length=length(tissues))
names(ClusterDEGbyTissue) = tissues

for(i in tissues){
 tmp = FindAllMarkers(subset(GW18AllGE4K_2Dumap.integrated,cells=colnames(GW18AllGE4K_2Dumap.integrated)[which(GW18AllGE4K_2Dumap.integrated@meta.data$Tissue==i)]))

 tmp$Tissue=i
ClusterDEGbyTissue[[i]] = tmp
cat(i)
}

save(ClusterDEGbyTissue,file='./Analysis/ClusterDEGbyTissue_HypoGEcrtx.rda',compress=TRUE) # maj region dif within cluster - uses all genes by SCT 


### evaluate

tissue = c("Hypo","Cortex","GE")[3]
cluster = 'clust_14'

test = ClusterDEGbyTissue[[tissue]][which(ClusterDEGbyTissue[[tissue]]$cluster==cluster & ClusterDEGbyTissue[[tissue]]$avg_logFC>0),]

head(test,50)



ClusterDEGbyTissue[[tissue]][which(ClusterDEGbyTissue[[tissue]]$gene=='DLX1' & ClusterDEGbyTissue[[tissue]]$avg_logFC>0),]

## hypo specific 

clust 14 - MFNG (all express this - but is c14 spec) , EGFR - precurors and 14, TFDP2, SPARCL1 is mostly RG, but hypo specific  - CFAP126 similar (but cilia, and RG exp), GPC3 specific, but what does it do? 

EOMES is in clust 19...how much? contaimination? - like 5 cells... 

pdf('./TestPlots/GW18AllGE4K_2Dumap_Hypo_clust14_Genes.pdf',width=12,height=8)
for(j in c("sample","seurat_clusters","Maj_Clust_From_Ref")){
print(DimPlot(GW18AllGE4K_2Dumap.integrated, reduction = "umap", group.by = j, label = TRUE, repel = TRUE))
}
for(g in c('SPARCL1','KCNQ1OT1','PTN','CFAP126','ABCA5','PCP4','GPC3','MFNG','EGFR','TFDP2','SRI','EGR1','EEF1D','BTG2','FOS','ASCL1','CHD7','NOB1','NES','CKS2','GTF2F2','NFIX')){
print(FeaturePlot(GW18AllGE4K_2Dumap.integrated,features=g, pt.size = 0.4,split.by='Tissue',combine=TRUE,order=TRUE))
}
dev.off() 


### genie 3 result merge 
genes=c('GPC3',as.character(genieGW18$TF[which(genieGW18$target=='GPC3')][1:10]))


test[match(genes,test$gene),] #  - FOS, HES1, EGR1 all also elevated in Hypo. - checked in TissueClusterDEG

## focus on FOS, HES1 and EGR1 

TissueClusterDEG[['clust_14']][which(
TissueClusterDEG[['clust_14']]$gene=='NFIB'),]






## cortex

clust 14 - EOMES - precursor and 14, 

#PPP1R17 (that's the odd one - RG in GE/Hypo), NEUROD4, PAX6 (overlaps with GE), BTG(exp in all), NEUROG1 - Ctx specific, and roughly c14 specific, PENK is interesting - c14 enriched in ctx, turned on very late in GE (ET) and not at all in Hypo, NHLH1 (TF) is both c14 and ctx specific - not much lit, but proven to be essential to neuro dev(kruger 2002) - CA12 is similar, but some expression in RG, carbonic anhydrase?!?, ) 

pdf('./TestPlots/GW18AllGE4K_2Dumap_Cortex_clust14_Genes.pdf',width=12,height=8)
for(j in c("sample","seurat_clusters","Maj_Clust_From_Ref")){
print(DimPlot(GW18AllGE4K_2Dumap.integrated, reduction = "umap", group.by = j, label = TRUE, repel = TRUE))
}
for(g in test$gene[1:25]){
print(FeaturePlot(GW18AllGE4K_2Dumap.integrated,features=g, pt.size = 0.4,split.by='Tissue',combine=TRUE,order=TRUE))
}
dev.off() 



## GE

DLX1, DLX2, DLX5, ZEB2, CADM1, COL1A2, HAT1, MOB3B

DLX - turned on earl in GE, and not at all in cortex

ZEB2 - on all the time in GE, late in cortex and early in hypo(RG) - TF linked to neuro dev

COL1A2 - very early in GE, but only in GE - collagen





Tissue comp - 

tissue = c("Hypo","Cortex","GE")[2]
cluster = 'clust_16'

test2 = TissueClusterDEG[[cluster]][which(TissueClusterDEG[[cluster]]$cluster==tissue & TissueClusterDEG[[cluster]]$avg_logFC>0),]

head(test2,50)



### top interesting genes for the paper 

pdf('./TestPlots/GW18AllGE4K_2Dumap_GE_clust14_Genes.pdf',width=12,height=8)
for(j in c("sample","seurat_clusters","Maj_Clust_From_Ref")){
print(DimPlot(GW18AllGE4K_2Dumap.integrated, reduction = "umap", group.by = j, label = TRUE, repel = TRUE))
}
for(g in c("DLX1","DLX2","DLX5","GLCCI1","ZEB2","CADM1","COL1A2","HAT1","MOB3B","UBE2C","PCNA","RAB3IP","EPHA4","MYO1B","CENPU","DLX6","FGD3","ARX","HSPA2","PDZRN3","TAC3")){
print(FeaturePlot(GW18AllGE4K_2Dumap.integrated,features=g, pt.size = 0.4,split.by='Tissue',combine=TRUE,order=TRUE))
}
dev.off() 










gene='NEUROD6'

tmplist = vector(mode='list',length=length(clusts))
names(tmplist)=clusts
for(i in clusts){
tmplist[[i]] = TissueClusterDEG[[i]][gene,]
}
tmpres=do.call(rbind,tmplist)
tmpres = tmpres[which(!is.na(tmpres[,1])),]

## super simple - grab pos FC, and low pct.2 - then table genes with 2 or more tissues - this is to find 'flips'

totres = do.call(rbind,TissueClusterDEG)

## sig / pos

totres = totres[which(!is.na(match(totres$seurat_cluster,paste('clust_',c(9, 14, 15, 16, 17, 21, 22, 23),sep='')))),] #20754

totres = totres[which(totres$avg_logFC>0 & totres$p_val_adj<=0.05),]
totres = totres[which(totres$pct.2<=0.1),]
geneTiss = tapply(totres$cluster,INDEX=totres$gene,FUN=table)

geneTiss2 = do.call(rbind,geneTiss)

geneTiss2 = as.data.frame(geneTiss2)

geneTiss2[which(geneTiss2$Cortex>0 & geneTiss2$Hypo>0),]

geneTiss2[which(geneTiss2$Cortex>0 & geneTiss2$GE>0),]

geneTiss2[which(geneTiss2$Hypo>0 & geneTiss2$GE>0),]


rownames(geneTiss2)[which(geneTiss2$Cortex>0 & geneTiss2$Hypo>0)]








pdf('./TestPlots/GW18AllGE4K_2Dumap_GE_Hypo_Ctx_QuickTestgenes.pdf',width=12,height=8)
for(j in c("sample","seurat_clusters","Maj_Clust_From_Ref")){
print(DimPlot(GW18AllGE4K_2Dumap.integrated, reduction = "umap", group.by = j, label = TRUE, repel = TRUE))
}
for(g in rownames(geneTiss2)[which(geneTiss2$Cortex>0 & geneTiss2$Hypo>0)]){
print(FeaturePlot(GW18AllGE4K_2Dumap.integrated,features=g, pt.size = 0.4,split.by='Tissue',combine=TRUE,order=TRUE))
}
dev.off()


### differences in RG - what happens 

## focus on clust 5 , 16 22 

if(FALSE){
## using topGO - looking at tutorial https://hms-dbmi.github.io/qmb-2016/exercises/10-gene-set-enrichment-and-pathways/

allGenes = TissueClusterDEG[['clust_16']]$p_val_adj
names(allGenes) = TissueClusterDEG[['clust_16']]$gene

selFun = TissueClusterDEG[['clust_16']]$cluster=='Cortex' & TissueClusterDEG[['clust_16']]$avg_logFC>0.5 & TissueClusterDEG[['clust_16']]$p_val_adj<=0.05
names(selFun) = TissueClusterDEG[['clust_16']]$gene

GOdata <- new('topGOdata',ontology='BP',allGenes=allGenes, geneSel = selFun, annot=annFUN.org,mapping = 'org.Hs.eg.db',ID = 'symbol')
### topGO seems like a real pain in the butt, trying GOstat
}

### looking at Dave Tang's blog: https://davetang.org/muse/2010/11/10/gene-ontology-enrichment-analysis/

## of course, compute nodes can't access internets 

ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
 
my_chr <- c(1:22, 'M', 'X', 'Y')
my_ensembl_gene <- getBM(attributes='ensembl_gene_id',filters = 'chromosome_name', values = my_chr, mart = ensembl)




### GO gene enrichments of tissue specific genes in Radial glia populations 

totDEG = do.call(rbind,TissueClusterDEG)

totgene = unique(totDEG$gene)

mart = useMart( "ensembl" )
mart = useDataset("hsapiens_gene_ensembl", mart = mart )

gid = getBM(mart = mart, attributes = c("hgnc_symbol","entrezgene_id"),values=totgene) #entrez to refseq

totDEG$entrez = gid$entrezgene_id[match(totDEG$gene,gid$hgnc_symbol)]

#entrez_object <- org.Hs.egGO

tissues = unique(totDEG$cluster)
clusters = c('clust_5','clust_16','clust_22')

### biological process
nonNeuGOBP = vector(mode='list',length=length(tissues)*length(clusters))

names(nonNeuGOBP) = paste(rep(tissues,each=3),rep(clusters,3),sep='_')

nonNeuGOBPgenes = vector(mode='list',length=length(tissues)*length(clusters))

names(nonNeuGOBPgenes) = paste(rep(tissues,each=3),rep(clusters,3),sep='_')

for(i in tissues){
for(j in clusters){

testGenes = na.omit(totDEG$entrez[which(totDEG$cluster==i & totDEG$avg_logFC>0.5 & totDEG$p_val_adj<=0.05 & totDEG$seurat_cluster==j)])
universeGenes = unique(na.omit(totDEG$entrez))

params <- new('GOHyperGParams',
              geneIds = testGenes,
              universeGeneIds = universeGenes,
              ontology = 'BP',
              pvalueCutoff = 0.05,
              conditional = FALSE,
              testDirection = 'over',
              annotation = "org.Hs.eg.db"
             )
 
 nonNeuGOBP[[paste(i,j,sep='_')]] <- hyperGTest(params)

### add in catch for genes per GO id 

goTerms = summary(nonNeuGOBP[[paste(i,j,sep='_')]])[,1]

goGenes = data.frame(GOID=goTerms,entrezGenes=NA,ncbiGenes=NA)

rownames(goGenes) = goTerms

for(k in goTerms){
tmpUni = geneIdUniverse(nonNeuGOBP[[paste(i,j,sep='_')]])[[k]]
tmpGenes = intersect(testGenes,tmpUni)

goGenes[k,'entrezGenes'] = paste(tmpGenes,collapse=',')

tmpGeneId = gid$hgnc_symbol[match(tmpGenes,gid$entrezgene_id)]
goGenes[k,'ncbiGenes'] = paste(tmpGeneId,collapse=',')

}

tmpRes = cbind(summary(nonNeuGOBP[[paste(i,j,sep='_')]]),goGenes$ncbiGenes)
colnames(tmpRes)[8] = 'Genes'

nonNeuGOBPgenes[[paste(i,j,sep='_')]] = tmpRes
cat(paste(i,j,' done', '\n',sep='_'))

}
}






### most of the GO makes sense - central nervous system growth for cortex, some cillia for hypothalamus 


### molecular function
nonNeuGOMF = vector(mode='list',length=length(tissues)*length(clusters))

names(nonNeuGOMF) = paste(rep(tissues,each=3),rep(clusters,3),sep='_')

nonNeuGOMFgenes = vector(mode='list',length=length(tissues)*length(clusters))

names(nonNeuGOMFgenes) = paste(rep(tissues,each=3),rep(clusters,3),sep='_')

for(i in tissues){
for(j in clusters){

testGenes = na.omit(totDEG$entrez[which(totDEG$cluster==i & totDEG$avg_logFC>0.5 & totDEG$p_val_adj<=0.05 & totDEG$seurat_cluster==j)])
universeGenes = unique(na.omit(totDEG$entrez))

params <- new('GOHyperGParams',
              geneIds = testGenes,
              universeGeneIds = universeGenes,
              ontology = 'MF',
              pvalueCutoff = 0.05,
              conditional = FALSE,
              testDirection = 'over',
              annotation = "org.Hs.eg.db"
             )
 
 nonNeuGOMF[[paste(i,j,sep='_')]] <- hyperGTest(params)
#cat(paste(i,j,' done', '\n',sep='_'))

### add in catch for genes per GO id 

goTerms = summary(nonNeuGOMF[[paste(i,j,sep='_')]])[,1]

goGenes = data.frame(GOID=goTerms,entrezGenes=NA,ncbiGenes=NA)

rownames(goGenes) = goTerms

for(k in goTerms){
tmpUni = geneIdUniverse(nonNeuGOMF[[paste(i,j,sep='_')]])[[k]]
tmpGenes = intersect(testGenes,tmpUni)

goGenes[k,'entrezGenes'] = paste(tmpGenes,collapse=',')

tmpGeneId = gid$hgnc_symbol[match(tmpGenes,gid$entrezgene_id)]
goGenes[k,'ncbiGenes'] = paste(tmpGeneId,collapse=',')

}

tmpRes = cbind(summary(nonNeuGOMF[[paste(i,j,sep='_')]]),goGenes$ncbiGenes)
colnames(tmpRes)[8] = 'Genes'

nonNeuGOMFgenes[[paste(i,j,sep='_')]] = tmpRes
cat(paste(i,j,' done', '\n',sep='_'))
}
}


## CC cell compartment
nonNeuGOCC = vector(mode='list',length=length(tissues)*length(clusters))

names(nonNeuGOCC) = paste(rep(tissues,each=3),rep(clusters,3),sep='_')

nonNeuGOCCgenes = vector(mode='list',length=length(tissues)*length(clusters))

names(nonNeuGOCCgenes) = paste(rep(tissues,each=3),rep(clusters,3),sep='_')

for(i in tissues){
for(j in clusters){

testGenes = na.omit(totDEG$entrez[which(totDEG$cluster==i & totDEG$avg_logFC>0.5 & totDEG$p_val_adj<=0.05 & totDEG$seurat_cluster==j)])
universeGenes = unique(na.omit(totDEG$entrez))

params <- new('GOHyperGParams',
              geneIds = testGenes,
              universeGeneIds = universeGenes,
              ontology = 'CC',
              pvalueCutoff = 0.05,
              conditional = FALSE,
              testDirection = 'over',
              annotation = "org.Hs.eg.db"
             )
 
 nonNeuGOCC[[paste(i,j,sep='_')]] <- hyperGTest(params)
#cat(paste(i,j,' done', '\n',sep='_'))
goTerms = summary(nonNeuGOCC[[paste(i,j,sep='_')]])[,1]

goGenes = data.frame(GOID=goTerms,entrezGenes=NA,ncbiGenes=NA)

rownames(goGenes) = goTerms

for(k in goTerms){
tmpUni = geneIdUniverse(nonNeuGOCC[[paste(i,j,sep='_')]])[[k]]
tmpGenes = intersect(testGenes,tmpUni)

goGenes[k,'entrezGenes'] = paste(tmpGenes,collapse=',')

tmpGeneId = gid$hgnc_symbol[match(tmpGenes,gid$entrezgene_id)]
goGenes[k,'ncbiGenes'] = paste(tmpGeneId,collapse=',')

}

tmpRes = cbind(summary(nonNeuGOCC[[paste(i,j,sep='_')]]),goGenes$ncbiGenes)
colnames(tmpRes)[8] = 'Genes'

nonNeuGOCCgenes[[paste(i,j,sep='_')]] = tmpRes
cat(paste(i,j,' done', '\n',sep='_'))
}
}


save(nonNeuGOBP,nonNeuGOMF,nonNeuGOCC,nonNeuGOBPgenes,nonNeuGOMFgenes,nonNeuGOCCgenes,file='/local/projects-t3/idea/bherb/Hypothalamus/PubRes/Analysis/RG_CtxHypoGE_GOterms.rda',compress=TRUE)

### example genes from GO - for paper 


head(summary(nonNeuGOBP[["Cortex_clust_5"]]))


pdf('./TestPlots/GW18AllGE4K_2Dumap_Cortex_clust_5_GO_teledev.pdf',width=12,height=8)
for(j in c("sample","seurat_clusters","Maj_Clust_From_Ref")){
print(DimPlot(GW18AllGE4K_2Dumap.integrated, reduction = "umap", group.by = j, label = TRUE, repel = TRUE))
}
for(g in strsplit(as.character(nonNeuGOBPgenes[["Cortex_clust_5"]][2,'Genes']),',')[[1]]){
print(FeaturePlot(GW18AllGE4K_2Dumap.integrated,features=g, pt.size = 0.4,split.by='Tissue',combine=TRUE,order=TRUE))
}
dev.off()

### these genes are not as specific to RG as I would like... use anova to separate effect of sample and cell type? 




### GO enrichments among clusters along neuro lineages 
tissues = unique(totDEG$cluster)
clusters = c('clust_14','clust_10','clust_2','clust_3','clust_9')
GOcat = c('BP','MF','CC')

### biological process
NeuGO = vector(mode='list',length=length(tissues)*length(clusters)*length(GOcat))

names(NeuGO) = paste(rep(tissues,each=15),rep(clusters,each=3),rep(GOcat,15),sep='_')

for(i in tissues){
for(j in clusters){
for(k in GOcat){

testGenes = na.omit(totDEG$entrez[which(totDEG$cluster==i & totDEG$avg_logFC>0.5 & totDEG$p_val_adj<=0.05 & totDEG$seurat_cluster==j)])
universeGenes = unique(na.omit(totDEG$entrez))

params <- new('GOHyperGParams',
              geneIds = testGenes,
              universeGeneIds = universeGenes,
              ontology = k,
              pvalueCutoff = 0.05,
              conditional = FALSE,
              testDirection = 'over',
              annotation = "org.Hs.eg.db"
             )
 
 NeuGO[[paste(i,j,k,sep='_')]] <- hyperGTest(params)
cat(paste(i,j,k,' done', '\n',sep='_'))

}
}
}


### most of the GO makes sense - but more dev genes activated in GE vs hypothalamus?



save(NeuGO,file='/local/projects-t3/idea/bherb/Hypothalamus/PubRes/Analysis/Neu_CtxHypoGE_GOterms.rda',compress=TRUE)



## pull out genes only over expressed in one cluster?

sigDEG = totDEG[which(totDEG$avg_logFC>0 & totDEG$p_val_adj<=0.05),]


tableDEG = table(sigDEG$gene)

intersect(sigDEG$gene[which(sigDEG$cluster=='Cortex' & sigDEG$seurat_cluster=='clust_5')],names(tableDEG)[which(tableDEG==1)])

intersect(sigDEG$gene[which(sigDEG$cluster=='Hypo' & sigDEG$seurat_cluster=='clust_5')],names(tableDEG)[which(tableDEG==1)])

intersect(sigDEG$gene[which(sigDEG$cluster=='GE' & sigDEG$seurat_cluster=='clust_5')],names(tableDEG)[which(tableDEG==1)])

intersect(sigDEG$gene[which(sigDEG$cluster=='Cortex' & sigDEG$seurat_cluster=='clust_16')],names(tableDEG)[which(tableDEG==1)])

intersect(sigDEG$gene[which(sigDEG$cluster=='Hypo' & sigDEG$seurat_cluster=='clust_16')],names(tableDEG)[which(tableDEG==1)])

intersect(sigDEG$gene[which(sigDEG$cluster=='GE' & sigDEG$seurat_cluster=='clust_16')],names(tableDEG)[which(tableDEG==1)])

### so some mixture....increase regional specificity of gene expression 


### Good enough? Direct GE to Hypo comp? 

test = TissueClusterDEG[['clust_10']]

head(test[which(test$cluster=='Hypo'),],40)


pdf('./TestPlots/GW18AllGE4K_2Dumap_Hypo_clust10_QuickTestgenes.pdf',width=12,height=8)
for(j in c("sample","seurat_clusters","Maj_Clust_From_Ref")){
print(DimPlot(GW18AllGE4K_2Dumap.integrated, reduction = "umap", group.by = j, label = TRUE, repel = TRUE))
}
for(g in c('BEX5','BEX2','PEG10','SCG2','SCG5','RIT2','PCP4','ANK3','CLU')){
print(FeaturePlot(GW18AllGE4K_2Dumap.integrated,features=g, pt.size = 0.4,split.by='Tissue',combine=TRUE,order=TRUE))
}
dev.off()




pdf('./TestPlots/GW18AllGE4K_2Dumap_GE_clust10_QuickTestgenes.pdf',width=12,height=8)
for(j in c("sample","seurat_clusters","Maj_Clust_From_Ref")){
print(DimPlot(GW18AllGE4K_2Dumap.integrated, reduction = "umap", group.by = j, label = TRUE, repel = TRUE))
}
for(g in c('DLX6-AS1','PFN2','DLX5','DLX1')){
print(FeaturePlot(GW18AllGE4K_2Dumap.integrated,features=g, pt.size = 0.4,split.by='Tissue',combine=TRUE,order=TRUE))
}
dev.off()



### overlap with TF's? 

test = TissueClusterDEG[['clust_14']]

test = test[na.omit(match(HSTF$Name,test$gene)),]

test[which(test$cluster=='GE' & test$avg_logFC>0 & test$p_val_adj<=0.05),]

test[which(test$cluster=='Hypo'  & test$p_val_adj<=0.05& test$avg_logFC>0),]

test[which(test$cluster=='Cortex'  & test$p_val_adj<=0.05& test$avg_logFC>0),]

x

## common markers in RG across samples 

GW18AllGE4K_2Dumap.integrated@meta.data$RadialGlia = 'N'

GW18AllGE4K_2Dumap.integrated@meta.data$RadialGlia[which(GW18AllGE4K_2Dumap.integrated@meta.data$seurat_clusters=='clust_16')] = 'Y'

GW18AllGE4K_2Dumap.integrated@meta.data$RadialGlia[which(GW18AllGE4K_2Dumap.integrated@meta.data$seurat_clusters=='clust_5')] = 'Y'


GW18AllGE4K_2Dumap.integrated@active.ident=as.factor(GW18AllGE4K_2Dumap.integrated$RadialGlia)

DefaultAssay(GW18AllGE4K_2Dumap.integrated) = 'SCT'

RGmarkers = FindMarkers(GW18AllGE4K_2Dumap.integrated,ident.1='Y',ident.2='N')

RGmarkersTF = FindMarkers(GW18AllGE4K_2Dumap.integrated,ident.1='Y',ident.2='N',features=intersect(rownames(GW18AllGE4K_2Dumap.integrated),HSTF$Name))


pdf('./TestPlots/GW18AllGE4K_2Dumap_RG_shared_markers.pdf',width=12,height=8)
for(j in c("sample","seurat_clusters","Maj_Clust_From_Ref")){
print(DimPlot(GW18AllGE4K_2Dumap.integrated, reduction = "umap", group.by = j, label = TRUE, repel = TRUE))
}
for(g in rownames(RGmarkersTF)[1:10]){
print(FeaturePlot(GW18AllGE4K_2Dumap.integrated,features=g, pt.size = 0.4,split.by='Tissue',combine=TRUE,order=TRUE))
}
dev.off()




## spearman....

GW18AllGE4K_2Dumap.integrated@assays$SCT[,which(GW18AllGE4K_2Dumap.integrated@meta.data$sample=='GW18')]


PPPexp = as.numeric(GW18AllGE4K_2Dumap.integrated@assays$SCT["PPP1R17",which(GW18AllGE4K_2Dumap.integrated@meta.data$sample=='GW18')])

HypoMat = t(as.matrix(GW18AllGE4K_2Dumap.integrated@assays$SCT[,which(GW18AllGE4K_2Dumap.integrated@meta.data$sample=='GW18')]))

ColSum = colSums(HypoMat)

ZeroInd = 


testCor = cor(PPPexp,HypoMat,method='pearson')



#### trying tradeSeq

## for fitGAM function - need singlecellexperiment, pseudotime, and cellWeights is lineage assignments  - can add 


### comparisons - 1. inib lineage hypo vs GE (including mature Ctx cells) - maybe can reassign clust 15,17,21 as one group  2. ctx excite vs GE/Hypo lineage (early? 14, 19,7,6,0/10 )


### inhib lineage 

### using guide 

RNGversion("3.5.0")
palette(brewer.pal(8, "Dark2"))
data(countMatrix, package = "tradeSeq")
counts <- as.matrix(countMatrix) ## row=gene, col=cell - raw counts 
rm(countMatrix)
data(crv, package = "tradeSeq") ## output from slingshot
data(celltype, package = "tradeSeq") ## vector of celltype assignments - 

### run slingshot

GW18AllGE_Hypo_Inhib = subset(GW18AllGE4K.integrated,idents=unique(GW18AllGE4K.integrated@active.ident)[c(14,1,3,10,29,15)])


GW18AllGE_Hypo_Inhib@meta.data$Lin_Group=NA

GW18AllGE_Hypo_Inhib@meta.data$Lin_Group[which(GW18AllGE_Hypo_Inhib@meta.data$BrainRegion=='Hypo' & GW18AllGE_Hypo_Inhib@meta.data$seurat_clusters=='clust_14')] = "Hypo_1"

GW18AllGE_Hypo_Inhib@meta.data$Lin_Group[which(GW18AllGE_Hypo_Inhib@meta.data$BrainRegion=='GE' & GW18AllGE_Hypo_Inhib@meta.data$seurat_clusters=='clust_14')] = "GE_1"

GW18AllGE_Hypo_Inhib@meta.data$Lin_Group[which(GW18AllGE_Hypo_Inhib@meta.data$BrainRegion=='Hypo' & GW18AllGE_Hypo_Inhib@meta.data$seurat_clusters=='clust_10')] = "Hypo_2"

GW18AllGE_Hypo_Inhib@meta.data$Lin_Group[which(GW18AllGE_Hypo_Inhib@meta.data$BrainRegion=='GE' & GW18AllGE_Hypo_Inhib@meta.data$seurat_clusters=='clust_10')] = "GE_2"

GW18AllGE_Hypo_Inhib@meta.data$Lin_Group[which(GW18AllGE_Hypo_Inhib@meta.data$BrainRegion=='Hypo' & GW18AllGE_Hypo_Inhib@meta.data$seurat_clusters=='clust_3')] = "Hypo_3"

GW18AllGE_Hypo_Inhib@meta.data$Lin_Group[which(GW18AllGE_Hypo_Inhib@meta.data$BrainRegion=='GE' & GW18AllGE_Hypo_Inhib@meta.data$seurat_clusters=='clust_3')] = "GE_3"

GW18AllGE_Hypo_Inhib@meta.data$Lin_Group[which(GW18AllGE_Hypo_Inhib@meta.data$BrainRegion=='Hypo' & GW18AllGE_Hypo_Inhib@meta.data$seurat_clusters=='clust_15')] = "Hypo_4"

GW18AllGE_Hypo_Inhib@meta.data$Lin_Group[which(GW18AllGE_Hypo_Inhib@meta.data$BrainRegion=='GE' & GW18AllGE_Hypo_Inhib@meta.data$seurat_clusters=='clust_15')] = "GE_4"

GW18AllGE_Hypo_Inhib@meta.data$Lin_Group[which(GW18AllGE_Hypo_Inhib@meta.data$BrainRegion=='Cortex' & GW18AllGE_Hypo_Inhib@meta.data$seurat_clusters=='clust_15')] = "GE_4"

GW18AllGE_Hypo_Inhib@meta.data$Lin_Group[which(GW18AllGE_Hypo_Inhib@meta.data$BrainRegion=='Hypo' & GW18AllGE_Hypo_Inhib@meta.data$seurat_clusters=='clust_17')] = "Hypo_4"

GW18AllGE_Hypo_Inhib@meta.data$Lin_Group[which(GW18AllGE_Hypo_Inhib@meta.data$BrainRegion=='GE' & GW18AllGE_Hypo_Inhib@meta.data$seurat_clusters=='clust_17')] = "GE_4"

GW18AllGE_Hypo_Inhib@meta.data$Lin_Group[which(GW18AllGE_Hypo_Inhib@meta.data$BrainRegion=='Cortex' & GW18AllGE_Hypo_Inhib@meta.data$seurat_clusters=='clust_17')] = "GE_4"

GW18AllGE_Hypo_Inhib@meta.data$Lin_Group[which(GW18AllGE_Hypo_Inhib@meta.data$BrainRegion=='Hypo' & GW18AllGE_Hypo_Inhib@meta.data$seurat_clusters=='clust_21')] = "Hypo_4"

GW18AllGE_Hypo_Inhib@meta.data$Lin_Group[which(GW18AllGE_Hypo_Inhib@meta.data$BrainRegion=='GE' & GW18AllGE_Hypo_Inhib@meta.data$seurat_clusters=='clust_21')] = "GE_4"

GW18AllGE_Hypo_Inhib@meta.data$Lin_Group[which(GW18AllGE_Hypo_Inhib@meta.data$BrainRegion=='Cortex' & GW18AllGE_Hypo_Inhib@meta.data$seurat_clusters=='clust_21')] = "GE_4"


GW18AllGE_Hypo_Inhib = subset(GW18AllGE_Hypo_Inhib,cells = colnames(GW18AllGE_Hypo_Inhib)[-which(is.na(GW18AllGE_Hypo_Inhib@meta.data$Lin_Group))])


GW18AllGE_Hypo_Inhib@meta.data$Lin_Group2 = GW18AllGE_Hypo_Inhib@meta.data$Lin_Group

 GW18AllGE_Hypo_Inhib@meta.data$Lin_Group2 = gsub('GE','Stage',GW18AllGE_Hypo_Inhib@meta.data$Lin_Group2)

 GW18AllGE_Hypo_Inhib@meta.data$Lin_Group2 = gsub('Hypo','Stage',GW18AllGE_Hypo_Inhib@meta.data$Lin_Group2)

GW18AllGE_Hypo_Inhib_SCE = Seurat::as.SingleCellExperiment(GW18AllGE_Hypo_Inhib) ## just 2000 var genes - try for now? - transfered PCA

GW18AllGE_Hypo_Inhib_SCE <- slingshot(GW18AllGE_Hypo_Inhib_SCE, clusterLabels = 'Lin_Group2', reducedDim = 'PCA',start.clus = 'Stage_1', end.clus = 'Stage_4')



GW18AllGE4K_2Dumap.integrated@meta.data$GE_Hypo_Inhib_pseudotime = 0

GW18AllGE4K_2Dumap.integrated@meta.data[colnames(GW18AllGE_Hypo_Inhib_SCE),'GE_Hypo_Inhib_pseudotime'] = GW18AllGE_Hypo_Inhib_SCE$slingPseudotime_1

GW18AllGE4K_2Dumap.integrated@meta.data$GE_Hypo_Inhib_Groups1= NA

GW18AllGE4K_2Dumap.integrated@meta.data[colnames(GW18AllGE_Hypo_Inhib_SCE),'GE_Hypo_Inhib_Groups1'] = GW18AllGE_Hypo_Inhib_SCE$Lin_Group

GW18AllGE4K_2Dumap.integrated@meta.data$GE_Hypo_Inhib_Groups2= NA

GW18AllGE4K_2Dumap.integrated@meta.data[colnames(GW18AllGE_Hypo_Inhib_SCE),'GE_Hypo_Inhib_Groups2'] = GW18AllGE_Hypo_Inhib_SCE$Lin_Group2

pdf('./TestPlots/GW18AllGE4K_2D_GE_Hypo_Inhib_TestPseudotime.pdf')
print(DimPlot(GW18AllGE4K_2Dumap.integrated, reduction = "umap", group.by = 'GE_Hypo_Inhib_Groups1', label = TRUE, repel = TRUE))
print(DimPlot(GW18AllGE4K_2Dumap.integrated, reduction = "umap", group.by = 'GE_Hypo_Inhib_Groups2', label = TRUE, repel = TRUE))
print(FeaturePlot(GW18AllGE4K_2Dumap.integrated, features = 'GE_Hypo_Inhib_pseudotime', pt.size = 0.4,order=TRUE))

dev.off()



GW18AllGE_Hypo_Inhib_PCA = GW18AllGE_Hypo_Inhib@reductions[["pca"]]@cell.embeddings[,1:30]

#GW18AllGE_Hypo_Inhib_Lin = GW18AllGE_Hypo_Inhib$Lin_Group2
GW18AllGE_Hypo_Inhib_Lin = GW18AllGE_Hypo_Inhib$Lin_Group


linGeHypo <- getLineages(GW18AllGE_Hypo_Inhib_PCA, GW18AllGE_Hypo_Inhib_Lin)

### forcing lineages


linGeHypo@lineages[[1]] = c('GE_1','GE_2','GE_3','GE_4')
linGeHypo@lineages[[2]] = c('Hypo_1','Hypo_2','Hypo_3','Hypo_4')

### change adjacency by hand
linGeHypo@adjacency[1,] = c(0,0,0,1,0,0,0,0)
linGeHypo@adjacency[2,] = c(1,0,1,0,0,0,0,0)
linGeHypo@adjacency[3,] = c(0,1,0,0,0,0,0,0)
linGeHypo@adjacency[4,] = c(1,0,0,0,0,0,0,0)
linGeHypo@adjacency[5,] = c(0,0,0,0,0,1,0,0)
linGeHypo@adjacency[6,] = c(0,0,0,0,1,0,0,1)
linGeHypo@adjacency[7,] = c(0,0,0,0,0,0,0,1)
linGeHypo@adjacency[8,] = c(0,0,0,0,0,1,1,0)

## slingParams might need to be changed...

linGeHypo = getCurves(linGeHypo)

pdf('./TestPlots/GW18AllGE_Hypo_Inhib_test_slingshot_lin.pdf',width=12,height=8)
plot(GW18AllGE_Hypo_Inhib_PCA[,1:2], col = brewer.pal(9,"Set1")[as.numeric(as.factor(GW18AllGE_Hypo_Inhib_Lin))], asp = 1, pch = 16)
lines(linGeHypo, lwd = 3, col = 'black')
dev.off()

VarGenes = VariableFeatures(GW18AllGE_Hypo_Inhib)

GeHypoCounts = as.matrix(GW18AllGE_Hypo_Inhib@assays$RNA[VarGenes,])

#cw.geh = matrix(0,nrow=ncol(GeHypoCounts),ncol=2)
#cw.geh[grep('Hypo_',GW18AllGE_Hypo_Inhib$Lin_Group),1]=1
#cw.geh[grep('GE_',GW18AllGE_Hypo_Inhib$Lin_Group),2]=1
cw.geh = rep(0,ncol(GeHypoCounts))
cw.geh[grep('Hypo_',GW18AllGE_Hypo_Inhib$Lin_Group)]=1
cw.geh[grep('GE_',GW18AllGE_Hypo_Inhib$Lin_Group)]=2



set.seed(5)
icMat <- evaluateK(counts = GeHypoCounts, sds = linGeHypo, k = 3:10, nGenes = 200, verbose = T)

sce <- fitGAM(counts = GeHypoCounts, pseudotime = GW18AllGE_Hypo_Inhib_SCE$slingPseudotime_1, cellWeights = cw.geh,nknots = 6, verbose = FALSE)


startRes <- startVsEndTest(sce)

oStart <- order(startRes$waldStat, decreasing = TRUE)
sigGeneStart <- names(sce)[oStart[3]]

pdf('./TestPlots/GW18AllGE4K_2D_GE_Hypo_Inhib_TestTradeSeqProgenitor.pdf')
plotSmoothers(sce[[sigGeneStart]])
dev.off()









Hypo_C9C2.t <- GW18AllGE_Hypo_Inhib_SCE$slingPseudotime_1


set.seed(5)
icMat <- evaluateK(counts = counts, sds = crv, k = 3:10, 
                   nGenes = 200, verbose = T)






### do metaneighbor on 4K clustering 








#### full projections 

GW18AllGE4K_2Dumap.integrated.data = as.matrix(GW18AllGE4K_2Dumap.integrated@assays$SCT[1:20311,1:36000])

proRes = vector(mode='list',length=length(proRef))

metaLength=ncol(cbind(GW18AllGE4K_2Dumap.integrated@meta.data))

#for(i in 1:length(proMat)){
for(i in c(1,4,6,7,8)){

tmpMat = proMat[[i]]
na.ind = which(is.na(tmpMat[,1]))
if(length(na.ind)>0){
tmpMat = tmpMat[-na.ind,] #22498 genes 
}
uniqGenes = names(table(tmpMat[,1]))[which(table(tmpMat[,1])==1)]
tmpMat = tmpMat[match(uniqGenes,tmpMat[,1]),]
rownames(tmpMat) = tmpMat[,1]
tmpMat = as.matrix(tmpMat[,-1])

proRes[[i]]=projectR(data=GW18AllGE4K_2Dumap.integrated.data, loadings = tmpMat, full = FALSE) ##
GW18AllGE4K_2Dumap.integrated@meta.data = cbind(GW18AllGE4K_2Dumap.integrated@meta.data,t(proRes[[i]]))

clrs00=c("black","darkred","red","orange","yellow")# "originColor" - 1st to last correspond to lo to hi values

pdf(paste('./TestPlots/GW18AllGE4K_2Dumap_Cortex_projection_',proNames[i],'.pdf',sep=''),width=12,height=8)
for(j in c("sample","seurat_clusters","Cell_Type","Maj_Clust_From_Ref")){
print(DimPlot(GW18AllGE4K_2Dumap.integrated, reduction = "umap", group.by = j, label = TRUE, repel = TRUE))
}
print(FeaturePlot(GW18AllGE4K_2Dumap.integrated,features='nCount_RNA', pt.size = 0.4))
for(k in rownames(proRes[[i]])){
clr=color.scale(proRes[[i]][k,],extremes=clrs00)
print(FeaturePlot(GW18AllGE4K_2Dumap.integrated,features=k, pt.size = 0.4,order=TRUE))
}
dev.off()
GW18AllGE4K_2Dumap.integrated@meta.data = GW18AllGE4K_2Dumap.integrated@meta.data[,1:39]
cat(paste(i,'/8: ',proNames[i],', ',sep=''))
}

GW18AllGE4K_2DumapproRes = proRes

### Decon projection neuron data

tmpMat = HypoLMDpat

tmpRes=projectR(data=GW18AllGE4K_2Dumap.integrated.data, loadings = tmpMat, full = FALSE) ##
GW18AllGE4K_2Dumap.integrated@meta.data = cbind(GW18AllGE4K_2Dumap.integrated@meta.data,t(tmpRes))

clrs00=c("black","darkred","red","orange","yellow")# "originColor" - 1st to last correspond to lo to hi values

pdf(paste('./TestPlots/GW18AllGE4K_2Dumap_','DeCon_Projection','.pdf',sep=''),width=12,height=8)
for(j in c("sample","seurat_clusters","Cell_Type","Maj_Clust_From_Ref")){
print(DimPlot(GW18AllGE4K_2Dumap.integrated, reduction = "umap", group.by = j, label = TRUE, repel = TRUE))
}
print(FeaturePlot(GW18AllGE4K_2Dumap.integrated,features='nCount_RNA', pt.size = 0.4))
for(k in rownames(tmpRes)){
clr=color.scale(tmpRes[k,],extremes=clrs00)
print(FeaturePlot(GW18AllGE4K_2Dumap.integrated,features=k, pt.size = 0.4,order=TRUE))
}
dev.off()





GW18AllGE4K_2Dumap.integrated@meta.data = GW18AllGE4K_2Dumap.integrated@meta.data[,1:39]
cat(paste(i,'/8: ',proNames[i],', ',sep=''))
}




#save(GW18AllGE4K_2Dumap.integrated,file='./SeuratObj/GW18AllGE4K_2Dumap.integrated.rda',compress=TRUE)

save(GW18AllGE4K_2DumapproRes,file='./Analysis/GW18AllGE4K_2DumapprojectionResults.rda',compress=T)

















pdf('./TestPlots/GW18AllGE4K_2D_GEH_TestPseudotime.pdf')
print(FeaturePlot(GW18AllGE4K_2Dumap.integrated, features = 'GEH_pseudotime', pt.size = 0.4,order=TRUE))
for(i in c('CENPF','HMGB2','PTTG1','MKI67','HMGN2','BIRC5','PAX6','TAGLN2','TUBA1A','MALAT1','STMN2','MEF2C','GPR22','SLA')){
print(FeaturePlot(GW18AllGE4K_2Dumap.integrated, features = i, pt.size = 0.4,order=TRUE))
}
dev.off()












assays(GW18AllGE_SCE)$norm 

reducedDims(GW18AllGE_SCE)











## if scaling is needed 
Tri2Neuron_ds.integrated.test = Tri2Neuron_ds.integrated

 DefaultAssay(Tri2Neuron_ds.integrated.test) = 'RNA'

Tri2Neuron_ds.integrated.test <- ScaleData(object = Tri2Neuron_ds.integrated.test, features = rownames(Tri2Neuron_ds.integrated.test))


pdf('./TestPlots/Tri2Neuron_ds_topDifTF_heatmap.pdf',width=12,height=8)
DoHeatmap(Tri2Neuron_ds.integrated.test, features = topDiftf,group.by="cluster_name",label=FALSE,raster=FALSE) + NoLegend()
dev.off()

pdf('./TestPlots/Tri2Neuron_ds_topPerDifTF_heatmap.pdf',width=12,height=8)
DoHeatmap(Tri2Neuron_ds.integrated.test, features = topPerDiftf,group.by="cluster_name",label=FALSE,raster=FALSE) + NoLegend()
dev.off()













### also group by ge, hypo and ctx





### top genes per cluster 

MajTypes=unique(GW18AllGE4K.integrated$seurat_clusters)

MajTypes=sort(MajTypes)

ClustMarkers = vector(mode='list',length=length(MajTypes))
names(ClustMarkers) = MajTypes


for(i in MajTypes){
    tmp = Markers.VARs[which(Markers.VARs$cluster==i & Markers.VARs$p_val_adj<=0.01),]
    tmp = tmp[order(tmp$avg_logFC,decreasing=TRUE),]
    tmp = tmp[which(tmp$avg_logFC>0.5),]
ClustMarkers[[i]] <-  as.character(na.omit(tmp$gene[1:10]))
}



for(i in MajTypes){
    tmp = Markers.TFs[which(Markers.TFs$cluster==i & Markers.TFs$p_val_adj<=0.01),]
    tmp = tmp[order(tmp$avg_logFC,decreasing=TRUE),]
    tmp = tmp[which(tmp$avg_logFC>0.5),]
ClustMarkers[[i]] <- unique(c(ClustMarkers[[i]],as.character(na.omit(tmp$gene[1:10]))))
}





### also check to see if these markers are expressed in each sample 









### proportion of sample in each cluster 

GW18AllGE4K.integrated@meta.data$samp_by_clust = 

data2= table(GW18AllGE4K.integrated@meta.data$sample,GW18AllGE4K.integrated@meta.data$seurat_clusters)

data3 = data.table::melt(data2)

# Stacked

pdf(file='./TestPlots/GW18AllGE4K_Cluster_content.pdf',width=20,height=8)
ggplot(data3, aes(fill=Var1, y=value, x=Var2)) + 
    geom_bar(position="fill", stat="identity")
dev.off()


### try monocle on this dataset 



DefaultAssay(GW18AllGE4K.integrated) = 'integrated' ## integrated has zeros

geneMeta = data.frame(gene_short_name=rownames(as.matrix(GetAssayData(GW18AllGE4K.integrated))))       
rownames(geneMeta) = rownames(as.matrix(GetAssayData(GW18AllGE4K.integrated)))

DefaultAssay(GW18AllGE4K.integrated) = 'RNA' ## integrated has zeros

dat = as.matrix(GetAssayData(GW18AllGE4K.integrated))
dat = dat[as.character(geneMeta[,1]),]
rownames(dat) = rownames(geneMeta)

GW18AllGE4KSeu <- new_cell_data_set(dat,cell_metadata = GW18AllGE4K.integrated@meta.data,gene_metadata = geneMeta)

GW18AllGE4KCDS <- preprocess_cds(GW18AllGE4KSeu, num_dim = 30)

DefaultAssay(GW18AllGE4K.integrated) = 'integrated' ## integrated has zeros

GW18AllGE4KCDS@preprocess_aux$gene_loadings = GW18AllGE4K.integrated@reductions[["pca"]]@feature.loadings[,1:30] ## PC by gene 

reducedDims(GW18AllGE4KCDS)[['PCA']] = GW18AllGE4K.integrated@reductions[["pca"]]@cell.embeddings[,1:30] ## from code - PC by cell 

GW18AllGE4KCDS <- reduce_dimension(GW18AllGE4KCDS,reduction_method='UMAP',preprocess_method='PCA')

GW18AllGE4KCDS <- cluster_cells(GW18AllGE4KCDS,reduction_method='UMAP',k=5)

Project = 'GW18AllGE4K'


pdf(paste('./TestPlots/',Project,'_Monocle.pdf',sep=''),width=12,height=8)

plot_cells(GW18AllGE4KCDS,reduction_method='UMAP',color_cells_by='partition',show_trajectory_graph = FALSE,cell_size=0.7,group_label_size=5)

plot_cells(GW18AllGE4KCDS,reduction_method='UMAP',color_cells_by='cluster',show_trajectory_graph = FALSE,cell_size=0.7,group_label_size=5)

plot_cells(GW18AllGE4KCDS,reduction_method='UMAP',color_cells_by='sample',show_trajectory_graph = FALSE,cell_size=0.7,group_label_size=5)

plot_cells(GW18AllGE4KCDS,reduction_method='UMAP',color_cells_by='seurat_clusters',show_trajectory_graph = FALSE,cell_size=0.7,group_label_size=5)

plot_cells(GW18AllGE4KCDS,reduction_method='UMAP',color_cells_by='Maj_Clust_From_Ref',show_trajectory_graph = FALSE,cell_size=0.7,group_label_size=5)


plot_cells(GW18AllGE4KCDS,reduction_method='UMAP',color_cells_by='Cell_Type',show_trajectory_graph = FALSE,cell_size=0.7,group_label_size=5)

dev.off()


GW18AllGE4KCDS <- learn_graph(GW18AllGE4KCDS,use_partition=FALSE) ## set use_partition=FALSE to connect branches 

pdf(paste('./TestPlots/',Project,'_Monocle_Trajectories.pdf',sep=""),width=12,height=8)

plot_cells(GW18AllGE4KCDS,color_cells_by = "sample",label_groups_by_cluster=FALSE,label_leaves=FALSE,label_branch_points=FALSE,cell_size=0.7)

plot_cells(GW18AllGE4KCDS, color_cells_by = "sample",label_cell_groups=FALSE,label_leaves=TRUE,label_branch_points=TRUE,graph_label_size=1.5,cell_size=0.7)

plot_cells(GW18AllGE4KCDS,color_cells_by = 'Maj_Clust_From_Ref',label_groups_by_cluster=FALSE,label_leaves=FALSE,label_branch_points=FALSE,cell_size=0.7)

plot_cells(GW18AllGE4KCDS, color_cells_by = 'Maj_Clust_From_Ref',label_cell_groups=FALSE,label_leaves=TRUE,label_branch_points=TRUE,graph_label_size=1.5,cell_size=0.7)

plot_cells(GW18AllGE4KCDS,color_cells_by = 'Cell_Type',label_groups_by_cluster=FALSE,label_leaves=FALSE,label_branch_points=FALSE,cell_size=0.7)

plot_cells(GW18AllGE4KCDS, color_cells_by = 'Cell_Type',label_cell_groups=FALSE,label_leaves=TRUE,label_branch_points=TRUE,graph_label_size=1.5,cell_size=0.7)

plot_cells(GW18AllGE4KCDS,color_cells_by = "seurat_clusters",label_groups_by_cluster=FALSE,label_leaves=FALSE,label_branch_points=FALSE,cell_size=0.7)

plot_cells(GW18AllGE4KCDS, color_cells_by = "seurat_clusters",label_cell_groups=FALSE,label_leaves=TRUE,label_branch_points=TRUE,graph_label_size=1.5,cell_size=0.7)

dev.off()

save(GW18AllGE4KCDS,file='./SeuratObj/GW18AllGE4KCDS_monocle.rda',compress=TRUE)




## 3D umaps? 

GW18AllGE4K.integrated <- RunUMAP(GW18AllGE4K.integrated, reduction = "pca", dims = 1:30,n.components = 3L)

umap_1 <- GW18AllGE4K.integrated[["umap"]]@cell.embeddings[,1]
umap_2 <- GW18AllGE4K.integrated[["umap"]]@cell.embeddings[,2]
umap_3 <- GW18AllGE4K.integrated[["umap"]]@cell.embeddings[,3]

# Visualize what headings are called so that you can extract them to form a dataframe
head(Embeddings(object = GW18AllGE4K.integrated, reduction = "umap"))

# Prepare a dataframe for cell plotting
plot.data <- FetchData(object = GW18AllGE4K.integrated, vars = c("UMAP_1", "UMAP_2", "UMAP_3", "seurat_clusters"))


#plot.data=plot.data[order(as.numeric(gsub('clust_','',plot.data$seurat_clusters)),decreasing=FALSE),]

plot.data=plot.data[order(plot.data$UMAP_3,decreasing=TRUE),]

pdf(file='./TestPlots/GW18AllGE4K_3D_Umap.pdf')
scatter3D(x = plot.data$UMAP_1, y = plot.data$UMAP_2, z =plot.data$UMAP_3,col.var = as.numeric(gsub('clust_','',plot.data$seurat_clusters))+1, pch='.', bty = "b2",add = FALSE,col = base::sample(gg.col(1000),30),phi = 0)
dev.off()



pdf(file='./TestPlots/GW18AllGE4K_3D_Umap.pdf')
scatter3D(x = plot.data$UMAP_1, y = plot.data$UMAP_2, z =plot.data$UMAP_3, pch='.', bty = "b2",add = FALSE,col = base::sample(gg.col(1000),30)[as.numeric(gsub('clust_','',plot.data$seurat_clusters))+1],phi = 0)
dev.off()


pdf(file='./TestPlots/GW18AllGE4K_3D_Umap_P20T90.pdf')
scatter3D(x = plot.data$UMAP_1, y = plot.data$UMAP_2, z =plot.data$UMAP_3, pch='.', bty = "b2",add = FALSE,col = base::sample(gg.col(32),30)[as.numeric(gsub('clust_','',plot.data$seurat_clusters))+1],theta = 90, phi = 20)
dev.off()





pdf(file='./TestPlots/3D_Test.pdf')

scatter3D(x, y, z, bty = "g", pch = 18, 
          col.var = as.integer(iris$Species), 
          col = c("#1B9E77", "#D95F02", "#7570B3"),
          pch = 18, ticktype = "detailed",
          colkey = list(at = c(2, 3, 4), side = 1, 
          addlines = TRUE, length = 0.5, width = 0.5,
          labels = c("setosa", "versicolor", "virginica")) )

dev.off()


plot.test=plot.data[which(plot.data$seurat_clusters=='clust_5' |plot.data$seurat_clusters=='clust_6' |plot.data$seurat_clusters=='clust_8' ),]


pdf(file='./TestPlots/3D_TestPlot.pdf')

scatter3D(x = plot.test$UMAP_1, y = plot.test$UMAP_2, z =plot.test$UMAP_3, pch='.', bty = "g", col = plot.test$color)

dev.off()







### try scatterplot3D
colors <- c("#999999", "#E69F00", "#56B4E9")
colors <- colors[as.numeric(iris$Species)]


pdf(file='./TestPlots/3D_Test.pdf')

scatterplot3d(iris[,1:3], pch = 16, color = colors,
              grid=TRUE, box=FALSE)

dev.off()


pdf(file='./TestPlots/3D_Test.pdf')


s3d <- scatterplot3d(iris[, 1:3], pch = "", grid=FALSE, box=FALSE)
# 3. Add grids
addgrids3d(iris[, 1:3], grid = c("xy", "xz", "yz"))
# 4. Add points
s3d$points3d(iris[, 1:3], pch = 16,color=colors)

dev.off()



colors <- c("#999999", "#E69F00", "#56B4E9")
colors <- colors[as.numeric(as.factor(plot.test$seurat_clusters))]

pdf(file='./TestPlots/3D_TestPlot.pdf')

scatterplot3d(x = plot.test$UMAP_1, y = plot.test$UMAP_2, z =plot.test$UMAP_3, pch='.', color = colors, grid=TRUE, box=FALSE )

dev.off()


library(scales)

GW18AllGE4K.integrated@active.ident=as.factor(GW18AllGE4K.integrated$seurat_clusters)

identities <- levels(GW18AllGE4K.integrated$integrated_snn_res.1.5)

# Create vector of default ggplot2 colors
my_color_palette <- hue_pal()(length(identities))

colors <- my_color_palette
colors <- colors[as.numeric(GW18AllGE4K.integrated$integrated_snn_res.1.5)]


plot.data$color = colors


plot.data = cbind(plot.data,GW18AllGE4K.integrated@meta.data[,c('sample','Maj_Clust_From_Ref','neuron_name','Cell_Type','Phase')])


write.csv(plot.data,file='./Analysis/Hypo_Ctx_GE_3D_UMAP.csv')

pdf(file='./TestPlots/3D_TestPlot3.pdf')

scatterplot3d(x = plot.data$UMAP_1, y = plot.data$UMAP_2, z =plot.data$UMAP_3, pch='.', color = colors, grid=TRUE, box=FALSE ,angle=90)

dev.off()





library(plotly)
library(htmlwidgets)
library("rmarkdown")
# Make a column of row name identities (these will be your cell/barcode names)

plot.data=read.csv("~/Documents/UMB_work/Ament/BRAIN_init/Kriegstein500k/Hypo/Paper/Analysis/Hypo_Ctx_GE_3D_UMAP.csv")

plot.data$label <- paste(rownames(plot.data))

#pdf(file='./TestPlots/GW18AllGE4K_3D_Umap.pdf')


### did on laptop



p=plot_ly(data = plot.data, x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, color = ~Phase,type = "scatter3d", mode= "markers",  marker = list(size = 2, width=2))

htmlwidgets::saveWidget(p, "~/Documents/UMB_work/Ament/BRAIN_init/Kriegstein500k/Hypo/Paper/Figure_Drafts/Phase_CxGeH_3D.html")


plot_ly(data = plot.data, x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, color = ~sample,type = "scatter3d", mode= "markers",  marker = list(size = 2, width=2))

plot_ly(data = plot.data, x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, color = ~Maj_Clust_From_Ref,type = "scatter3d", mode= "markers",  marker = list(size = 2, width=2))


plot_ly(data = plot.data, x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, color = ~seurat_clusters,type = "scatter3d", mode= "markers",  marker = list(size = 2, width=2))


plot_ly(data = plot.data, x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, color = ~Cell_Type,type = "scatter3d", mode= "markers",  marker = list(size = 2, width=2))

#orca(fig, './TestPlots/GW18AllGE4K_3D_Umap.pdf')


plot.test = plot.data[which(plot.data$sample == 'GW18' | plot.data$sample == 'GW18_CGE' | plot.data$sample == 'GW18_LGE' | plot.data$sample == 'GW18_MGE'),]

plot_ly(data = plot.test, x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, color = ~sample,type = "scatter3d", mode= "markers",  marker = list(size = 2, width=2))





GW18AllGE4K.integrated[["UMAP"]]@cell.embeddings


### GE and Hypo dev markers 

DefaultAssay(GW18AllGE4K.integrated) = 'RNA'

genes = c('NKX2-1','SHH','SOX2','TBX2','TBX3','BMP2','BMP4','BMP7','FGF8','FGF10','WNT8B','TCF7L1','RAX','OTP','SIX6','LHX2','POU3F2','SIM2','DLX1','DBX1','BSX','ASCL1','NEUROG3','NEUROD1','NHLH2','FEZF1','NR4A1','NHLH2','GSX1','IKZF1','HMX2','HMX3','NR5A2','FOXB1','SP8','SP9')

pdf(file='./TestPlots/GW18AllGE4K_NeuroDevGenes.pdf',width=12,height=8)
print(DimPlot(GW18AllGE4K.integrated, reduction = "umap", group.by = 'seurat_clusters', label = TRUE, repel = TRUE))
print(DimPlot(GW18AllGE4K.integrated, reduction = "umap", group.by = 'seurat_clusters', label = TRUE, repel = TRUE,dims=c(1,3)))
print(DimPlot(GW18AllGE4K.integrated, reduction = "umap", group.by = 'seurat_clusters', label = TRUE, repel = TRUE,dims=c(2,3)))
for(i in genes){
print(FeaturePlot(GW18AllGE4K.integrated, features = i, pt.size = 0.2,order=TRUE))
print(FeaturePlot(GW18AllGE4K.integrated, features = i, pt.size = 0.2,order=TRUE,dims=c(1,3)))
print(FeaturePlot(GW18AllGE4K.integrated, features = i, pt.size = 0.2,order=TRUE,dims=c(2,3)))
}
dev.off()



genes = c('SLC17A7','SLC30A3','OTOF','RORB','RSPO1','FEZF2','SULF1','CAR3','FAM84B','SLA2','FOXP2','NXPH4','GAD1','PROX1','LAMP5','NDNF','SNCG','VIP','LHX6','SST','CHODL','PVALB','VIPR2')

genes = intersect(rownames(GW18AllGE4K.integrated),genes)


pdf(file='./TestPlots/GW18AllGE4K_EIGenes.pdf',width=12,height=8)
for(i in genes){
print(FeaturePlot(GW18AllGE4K.integrated, features = i, pt.size = 0.2,order=TRUE))
}
dev.off()




## save sparse matrix instead?!?

write.csv(as.matrix(GW18AllGE4K.integrated@assays$RNA[1:33694,1:36000]),file='./Analysis/GW18AllGE4K_counts.csv')







#### just Hypo and GE's 




GW18HypoGECombDat = c(Tri2NormDat[['GW18']],GW18GENormDat)
names(GW18HypoGECombDat)[1] = 'GW18_Hypo'

GW18HypoGESamples = names(GW18HypoGECombDat)

### just redoing the subsampling for time

for(i in GW18HypoGESamples){
GW18HypoGECombDat[[i]] = subset(GW18HypoGECombDat[[i]],cells = colnames(GW18HypoGECombDat[[i]])[sample(c(1:ncol(GW18HypoGECombDat[[i]])),2000)])
}

for (i in GW18HypoGESamples) {
    GW18HypoGECombDat[[i]] <- NormalizeData(GW18HypoGECombDat[[i]], verbose = FALSE)
    GW18HypoGECombDat[[i]] <- FindVariableFeatures(GW18HypoGECombDat[[i]], selection.method = "vst", nfeatures = 2000, verbose = FALSE)
}

options(future.globals.maxSize = 4000 * 1024^2)

GW18HypoGE.anchors <- FindIntegrationAnchors(object.list = GW18HypoGECombDat, dims = 1:30)
GW18HypoGE.integrated <- IntegrateData(anchorset = GW18HypoGE.anchors, dims = 1:30)
DefaultAssay(GW18HypoGE.integrated) <- "integrated"
GW18HypoGE.integrated <- ScaleData(GW18HypoGE.integrated, verbose = FALSE)
GW18HypoGE.integrated <- RunPCA(GW18HypoGE.integrated, npcs = 30, verbose = FALSE)
GW18HypoGE.integrated <- RunUMAP(GW18HypoGE.integrated, reduction = "pca", dims = 1:30)

pdf('./TestPlots/GW18HypoGE_Merged.pdf',width=12,height=8)
print(DimPlot(GW18HypoGE.integrated, label = TRUE) + NoLegend())
for(i in c("sample","Maj_Clust_From_Ref","neuron_name")){
print(DimPlot(GW18HypoGE.integrated, reduction = "umap", group.by = i, label = TRUE, repel = TRUE))
}
print(FeaturePlot(GW18HypoGE.integrated, features = c("SLC17A6","GAD1"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(GW18HypoGE.integrated, features = c("RAX","CHGB"), pt.size = 0.2, ncol = 2))
dev.off()

pdf('./TestPlots/GW18HypoGE_Merged_SOX10.pdf',width=12,height=8)
print(DimPlot(GW18HypoGE.integrated, label = TRUE) + NoLegend())
for(i in c("sample","Maj_Clust_From_Ref","neuron_name")){
print(DimPlot(GW18HypoGE.integrated, reduction = "umap", group.by = i, label = TRUE, repel = TRUE))
}
print(FeaturePlot(GW18HypoGE.integrated, features = c("SOX10"), pt.size = 0.4))
dev.off()

subplots(GW18HypoGE.integrated,'sample',file='./TestPlots/GW18HypoGE_Samples_IndvPlots.pdf')

subplots(GW18HypoGE.integrated,'Maj_Clust_From_Ref',file='./TestPlots/GW18HypoGE_Maj_Clust_From_Ref_IndvPlots.pdf')

### GE projection 



GW18HypoGE.integrated.data = as.matrix(GW18HypoGE.integrated@assays$SCT[1:19211,1:8000])

proRes = vector(mode='list',length=length(proRef))

metaLength=ncol(cbind(GW18HypoGE.integrated@meta.data))

#for(i in 1:length(proMat)){
for(i in c(1,4,6,7,8)){

tmpMat = proMat[[i]]
na.ind = which(is.na(tmpMat[,1]))
if(length(na.ind)>0){
tmpMat = tmpMat[-na.ind,] #22498 genes 
}
uniqGenes = names(table(tmpMat[,1]))[which(table(tmpMat[,1])==1)]
tmpMat = tmpMat[match(uniqGenes,tmpMat[,1]),]
rownames(tmpMat) = tmpMat[,1]
tmpMat = as.matrix(tmpMat[,-1])

proRes[[i]]=projectR(data=GW18HypoGE.integrated.data, loadings = tmpMat, full = FALSE) ##
GW18HypoGE.integrated@meta.data = cbind(GW18HypoGE.integrated@meta.data,t(proRes[[i]]))

clrs00=c("black","darkred","red","orange","yellow")# "originColor" - 1st to last correspond to lo to hi values

pdf(paste('./TestPlots/GW18HypoGE_Cortex_projection_',proNames[i],'.pdf',sep=''),width=12,height=8)
for(j in c("sample","neuron_name","Maj_Clust_From_Ref")){
print(DimPlot(GW18HypoGE.integrated, reduction = "umap", group.by = j, label = TRUE, repel = TRUE))
}
print(FeaturePlot(GW18HypoGE.integrated,features='nCount_RNA', pt.size = 0.4))
for(k in rownames(proRes[[i]])){
clr=color.scale(proRes[[i]][k,],extremes=clrs00)
print(FeaturePlot(GW18HypoGE.integrated,features=k, pt.size = 0.4))
}
dev.off()
GW18HypoGE.integrated@meta.data = GW18HypoGE.integrated@meta.data[,1:11]
cat(paste(i,'/8: ',proNames[i],', ',sep=''))
}

GW18HypoGEproRes = proRes

save(GW18HypoGE.integrated,file='./SeuratObj/GW18HypoGE.integrated.rda',compress=TRUE)

save(GW18HypoGEproRes,file='./Analysis/GW18HypoGEprojectionResults.rda',compress=T)


### color GE+Hypo plot with clusters from GE+Hypo+Ctx merge 

GW18HypoGE.integrated@meta.data$CtxGEH_SC = GW18AllGE4K_2Dumap.integrated@meta.data[colnames(GW18HypoGE.integrated),'seurat_clusters']

GW18HypoGE.integrated@meta.data$CtxGEH_SC_Lin = NA

for(i in paste('clust_',c(14,10,9,3,2),sep='')){
GW18HypoGE.integrated@meta.data$CtxGEH_SC_Lin[which(GW18HypoGE.integrated@meta.data$CtxGEH_SC==i)] = i
}

