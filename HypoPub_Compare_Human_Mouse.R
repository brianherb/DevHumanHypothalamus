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



load('./SeuratObj/BlackSamples_PostSCT.rda')


## Trimester 2 Arc 









### CS22 



##CS22 with region info

load('./SeuratObj/Tri1Samples_integrated.rda')
load('./SeuratObj/Tri1Samples_PostSCT.rda') 

CombDat = list(BlackNormDat[['E11']],BlackNormDat[['E12']],BlackNormDat[['E13']],Tri1NormDat[['CS22_1']],Tri1NormDat[['CS22_2']])

names(CombDat) = c('Black_E11','Black_E12','Black_E13','Human_CS22_1','Human_CS22_2')

## add human region info - (hannah just did tri2 so far ) 


tmpregion = data.frame(cell=colnames(AllTri1.integrated),region=AllTri1.integrated@meta.data$region)
tmpregion = tmpregion[-which(is.na(tmpregion$region)),]

tmpregion1 = tmpregion[grep('CS22_1',tmpregion$cell),]
tmpregion2 = tmpregion[grep('CS22_2',tmpregion$cell),]


CombDat[['Human_CS22_1']]@meta.data$region = NA
CombDat[['Human_CS22_1']]@meta.data[as.character(tmpregion1$cell),'region'] = as.character(tmpregion1$region)

CombDat[['Human_CS22_2']]@meta.data$region = NA
CombDat[['Human_CS22_2']]@meta.data[as.character(tmpregion2$cell),'region'] = as.character(tmpregion2$region)

## since CS22 has less cells, no subsetting - and use same cells in Tri2 merge 
for(i in names(CombDat)[1:3]){
CombDat[[i]] = subset(CombDat[[i]],cells=colnames(E11_E13_GW18.integrated)[which(E11_E13_GW18.integrated$sample==i)])
}

DefaultAssay(CombDat[['Human_CS22_1']])='RNA'
DefaultAssay(CombDat[['Human_CS22_2']])='RNA'

for(i in names(CombDat)){
CombDat[[i]]@meta.data$sample = i
}

for (i in 1:length(CombDat)) {
    CombDat[[i]] <- NormalizeData(CombDat[[i]], verbose = FALSE)
    CombDat[[i]] <- FindVariableFeatures(CombDat[[i]], selection.method = "vst", nfeatures = 2000, verbose = FALSE)
}

options(future.globals.maxSize = 4000 * 1024^2)

E11_E13_CS22.anchors <- FindIntegrationAnchors(object.list = CombDat, dims = 1:30)
E11_E13_CS22.integrated <- IntegrateData(anchorset = E11_E13_CS22.anchors, dims = 1:30)
DefaultAssay(E11_E13_CS22.integrated) <- "integrated"
E11_E13_CS22.integrated <- ScaleData(E11_E13_CS22.integrated, verbose = FALSE)
E11_E13_CS22.integrated <- RunPCA(E11_E13_CS22.integrated, npcs = 30, verbose = FALSE)
E11_E13_CS22.integrated <- RunUMAP(E11_E13_CS22.integrated, reduction = "pca", dims = 1:30)
E11_E13_CS22.integrated <- RunTSNE(E11_E13_CS22.integrated, reduction = "pca", dims = 1:30)




### blackshaw E11 - E13 metadata 
BlackMetaEarly = read.csv('/local/projects-t3/idea/bherb/Hypothalamus/Blackshaw/data/GSE132355_E11-E13_umap_data.csv',row.names=1)

E11_E13_CS22.integrated@meta.data$BlackMetaEarly = NA
E11_E13_CS22.integrated@meta.data$BlackMetaEarly = BlackMetaEarly[colnames(E11_E13_CS22.integrated),'Cluster']

pdf('./TestPlots/E11_E13_CS22_Major_Sub_clusters_tSNE.pdf',width=12,height=8)
print(DimPlot(E11_E13_CS22.integrated, label = TRUE, reduction = "tsne") + NoLegend())
for(i in c("sample","Cluster","BlackMetaEarly","region")){
print(DimPlot(E11_E13_CS22.integrated, reduction = "tsne", group.by = i, label = TRUE, repel = TRUE))
}
print(FeaturePlot(E11_E13_CS22.integrated, features = c("SLC17A6","GAD1"), pt.size = 0.2, ncol = 2, reduction = "tsne"))
print(FeaturePlot(E11_E13_CS22.integrated, features = c("RAX","CHGB"), pt.size = 0.2, ncol = 2, reduction = "tsne"))
dev.off()


pdf('./TestPlots/E11_E13_CS22_Major_Sub_clusters_uMAP.pdf',width=12,height=8)
print(DimPlot(E11_E13_CS22.integrated, label = TRUE, reduction = "umap") + NoLegend())
for(i in c("sample","Cluster","BlackMetaEarly","region")){
print(DimPlot(E11_E13_CS22.integrated, reduction = "umap", group.by = i, label = TRUE, repel = TRUE))
}
print(FeaturePlot(E11_E13_CS22.integrated, features = c("SLC17A6","GAD1"), pt.size = 0.2, ncol = 2, reduction = "umap"))
print(FeaturePlot(E11_E13_CS22.integrated, features = c("RAX","CHGB"), pt.size = 0.2, ncol = 2, reduction = "umap"))
dev.off()

### try metaneighbor to compare mouse to human regions? 


E11_E13_CS22.integrated@meta.data$species = 'mouse'
E11_E13_CS22.integrated@meta.data$species[grep('CS22',E11_E13_CS22.integrated@meta.data$sample)] = 'human'

E11_E13_CS22.integrated@meta.data$MergedRegion = E11_E13_CS22.integrated@meta.data$Cluster
E11_E13_CS22.integrated@meta.data$MergedRegion[grep('CS22',E11_E13_CS22.integrated@meta.data$sample)] = E11_E13_CS22.integrated@meta.data$region[grep('CS22',E11_E13_CS22.integrated@meta.data$sample)]


### common regions 

#RegionLookup = dataframe(Mouse=c('ARC','DMH','MMN'),MouseCom=c('ARC_MM','DMH_MM','MMN_MM'),Human=c('ARC','DMH'),HumanCom=c('ARC_HS','DMH_HS'))



E11_E13_CS22.integrated@meta.data$MouseRegion = NA

E11_E13_CS22.integrated@meta.data$MouseRegion[which(E11_E13_CS22.integrated@meta.data$BlackMetaEarly=='ARC')] = 'ARC'
E11_E13_CS22.integrated@meta.data$MouseRegion[which(E11_E13_CS22.integrated@meta.data$BlackMetaEarly=='DMH')] = 'DMH'
E11_E13_CS22.integrated@meta.data$MouseRegion[which(E11_E13_CS22.integrated@meta.data$BlackMetaEarly=='MMH')] = 'MMH'
E11_E13_CS22.integrated@meta.data$MouseRegion[which(E11_E13_CS22.integrated@meta.data$BlackMetaEarly=='PMN')] = 'PMN'
E11_E13_CS22.integrated@meta.data$MouseRegion[which(E11_E13_CS22.integrated@meta.data$BlackMetaEarly=='POA')] = 'POA'
E11_E13_CS22.integrated@meta.data$MouseRegion[which(E11_E13_CS22.integrated@meta.data$BlackMetaEarly=='Prethalamus')] = 'Prethalamus'
E11_E13_CS22.integrated@meta.data$MouseRegion[which(E11_E13_CS22.integrated@meta.data$BlackMetaEarly=='PVH & SON')] = 'PVH & SON'
E11_E13_CS22.integrated@meta.data$MouseRegion[which(E11_E13_CS22.integrated@meta.data$BlackMetaEarly=='SCN')] = 'SCN'
E11_E13_CS22.integrated@meta.data$MouseRegion[which(E11_E13_CS22.integrated@meta.data$BlackMetaEarly=='SMN')] = 'SMN'
E11_E13_CS22.integrated@meta.data$MouseRegion[which(E11_E13_CS22.integrated@meta.data$BlackMetaEarly=='VMH')] = 'VMH'


E11_E13_CS22.integrated@meta.data$HumanRegion = NA


E11_E13_CS22.integrated@meta.data$HumanRegion[which(E11_E13_CS22.integrated@meta.data$region=='ARC')] = 'ARC'
E11_E13_CS22.integrated@meta.data$HumanRegion[which(E11_E13_CS22.integrated@meta.data$region=='DMH')] = 'DMH'
E11_E13_CS22.integrated@meta.data$HumanRegion[which(E11_E13_CS22.integrated@meta.data$region=='MMH')] = 'MMH'
E11_E13_CS22.integrated@meta.data$HumanRegion[which(E11_E13_CS22.integrated@meta.data$region=='PMN')] = 'PMN'
E11_E13_CS22.integrated@meta.data$HumanRegion[which(E11_E13_CS22.integrated@meta.data$region=='POA')] = 'POA'
E11_E13_CS22.integrated@meta.data$HumanRegion[which(E11_E13_CS22.integrated@meta.data$region=='PreThal')] = 'Prethalamus'
E11_E13_CS22.integrated@meta.data$HumanRegion[which(E11_E13_CS22.integrated@meta.data$region=='PVH_SON')] = 'PVH & SON'
E11_E13_CS22.integrated@meta.data$HumanRegion[which(E11_E13_CS22.integrated@meta.data$region=='SCN')] = 'SCN'
E11_E13_CS22.integrated@meta.data$HumanRegion[which(E11_E13_CS22.integrated@meta.data$region=='SMN')] = 'SMN'
E11_E13_CS22.integrated@meta.data$HumanRegion[which(E11_E13_CS22.integrated@meta.data$region=='VMH')] = 'VMH'



pdf('./TestPlots/E11_E13_CS22_Major_Nuclei_uMAP.pdf',width=12,height=8)
print(DimPlot(E11_E13_CS22.integrated, label = TRUE, reduction = "umap") + NoLegend())
for(i in c("sample","Cluster","BlackMetaEarly","region","MouseRegion","HumanRegion")){
print(DimPlot(E11_E13_CS22.integrated, reduction = "umap", group.by = i, label = TRUE, repel = TRUE))
}
dev.off()


E11_E13_CS22.integrated$MergedEarlyRegion = NA

E11_E13_CS22.integrated$MergedEarlyRegion[which(!is.na(E11_E13_CS22.integrated@meta.data$HumanRegion))] = paste('HS_',E11_E13_CS22.integrated@meta.data$HumanRegion[which(!is.na(E11_E13_CS22.integrated@meta.data$HumanRegion))],sep='')

E11_E13_CS22.integrated$MergedEarlyRegion[which(!is.na(E11_E13_CS22.integrated@meta.data$MouseRegion))] = paste('MM_',E11_E13_CS22.integrated@meta.data$MouseRegion[which(!is.na(E11_E13_CS22.integrated@meta.data$MouseRegion))],sep='')


E11_E13_CS22.integrated@active.ident=as.factor(E11_E13_CS22.integrated$MergedEarlyRegion)

cellTypes = unique(na.omit(E11_E13_CS22.integrated@meta.data$MouseRegion))

BlackCS22DEGrna = vector(mode='list',length=length(cellTypes))
names(BlackCS22DEGrna) = cellTypes

DefaultAssay(E11_E13_CS22.integrated) = 'RNA'


nonzeroMouseGenes = names(which(rowMeans(E11_E13_CS22.integrated@assays$RNA[,which(E11_E13_CS22.integrated@meta.data$species=='mouse')])>0))

nonzeroHumanGenes = names(which(rowMeans(E11_E13_CS22.integrated@assays$RNA[,which(E11_E13_CS22.integrated@meta.data$species=='human')])>0))

okGenes = intersect(intersect(nonzeroMouseGenes,nonzeroHumanGenes),rownames(E11_E13_CS22.integrated@assays$RNA))

for(i in cellTypes){

BlackCS22DEGrna[[i]] = FindMarkers(E11_E13_CS22.integrated, ident.1 = paste("HS_",i,sep=''), ident.2 = paste("MM_",i,sep=''),features=okGenes)
cat(paste('\n',i,' Done','\n',sep=''))

}

for(i in cellTypes){

BlackCS22DEGrna[[i]]$gene = rownames(BlackCS22DEGrna[[i]])
BlackCS22DEGrna[[i]]$region = i

}



head(BlackCS22DEGrna[['ARC']][intersect(rownames(BlackCS22DEGrna[['ARC']]),TF_NP$Transcription_factors),],60) ## looks good 

head(BlackCS22DEGrna[[sort(cellTypes)[9]]][intersect(rownames(BlackCS22DEGrna[[sort(cellTypes)[9]]]),NPlist2),],30)

## 

Arc  M - POMC ,    H - SCG3, SCG2, NUCB2

DMH. M - SST  ,     H - NUCB2, SCG3, CHGB

PMN  M - ,         H - GHRH, SCG3, NUCB2, CHGB, SST, PNOC

POA M- UBL5, SCG5, PNOC,  H- SCG3, NUCB2

Prethal  M- UBL5       h- SCG2 , SCG3, TAC1, SST, NUCB2

PVH SON   M  - CARTPT, SCG5, CBLN1, UBL5    H- SCG3, TRH SCG2 NUCB2, CHCB

SCN. M -  PNOC DBI      H - SCG3, NXPH1, SCG2, NTS,

SMN. M -  ADCYAP1, NXPH1, NXPH4, UBL5, CBLN2, CHGB, SCG5, DBI     H- SCG3, SCG2, NUCB2, TAC1 

VMH M- SCG5 CBLN1. H - SCG3, NUCB2 


NPdif = c()


Totdeg = do.call(rbind,BlackCS22DEGrna)

TotdegNP = Totdeg[which(!is.na(match(Totdeg$gene,intersect(Totdeg$gene,NPlist2)))),]


TotdegNPsig = TotdegNP[which(TotdegNP$p_val_adj<=0.05),]


TotdegNPsigMM = TotdegNPsig[which(TotdegNPsig$avg_log2FC<0),]
TotdegNPsigMM = TotdegNPsigMM[order(TotdegNPsigMM$region),]

TotdegNPsigHS = TotdegNPsig[which(TotdegNPsig$avg_log2FC>0),]
TotdegNPsigHS = TotdegNPsigHS[order(TotdegNPsigHS$region),]

NPdif = unique(c(rev(unique(TotdegNPsigHS$gene)),rev(unique(TotdegNPsigMM$gene))))

NPdif = NA

for(i in cellTypes){

NPdif=c(NPdif,intersect(rownames(BlackCS22DEGrna[[i]]),NPlist2)[1:5])
}

NPdif = unique(NPdif[-1])


## dot plot not working? 


pdf(file=paste("./TestPlots/E11_E13_CS22_NP_Dif_genes_VlnPlot.pdf",sep=''),width=12,height=8)
print(VlnPlot(E11_E13_CS22.integrated, features = rev(NPdif)[1])
)
dev.off()

E11_E13_CS22.nuclei = subset(E11_E13_CS22.integrated,cells=colnames(E11_E13_CS22.integrated)[which(!is.na(E11_E13_CS22.integrated@meta.data$MergedEarlyRegion))])


regions = paste(rep(c('HS_','MM_'),9), rep(sort(unique(na.omit(E11_E13_CS22.nuclei@meta.data$MouseRegion))),each=2),sep='')

E11_E13_CS22.nuclei@active.ident=factor(E11_E13_CS22.nuclei$MergedEarlyRegion,levels=rev(regions)) 

pdf(file=paste("./TestPlots/E11_E13_CS22_NP_Dif_genes_DotPlot.pdf",sep=''),width=12,height=8)
print(DotPlot(E11_E13_CS22.nuclei, features = rev(NPdif))
 + RotatedAxis() )
dev.off()

## try color by species 

E11_E13_CS22.nuclei@meta.data$JustRegion = flexsplit(E11_E13_CS22.nuclei@meta.data$MergedEarlyRegion,'_')[,2]


E11_E13_CS22.nuclei@active.ident=factor(E11_E13_CS22.nuclei$JustRegion,levels=rev(sort(unique(na.omit(E11_E13_CS22.nuclei$JustRegion))))) 

pdf(file=paste("./TestPlots/E11_E13_CS22_NP_Dif_genes_SpeciesColor_DotPlot.pdf",sep=''),width=12,height=8)
print(DotPlot(E11_E13_CS22.nuclei, group.by='JustRegion',split.by='species',features = NPdif,cols=c('red','blue'))
 + RotatedAxis() )
dev.off()


#head(BlackCS22DEGsct[['ARC']][intersect(rownames(BlackCS22DEGsct[['ARC']]),TF_NP$Transcription_factors),],60) ## looks good 



Prethalamus  PreThal

PVH & SON   PVH_SON


PVH_SON_DEG = FindMarkers(E11_E13_CS22.integrated, ident.1 = paste("Black_",'PVH & SON',sep=''), ident.2 = paste("CS22_",'PVH_SON',sep=''),features=okGenes)

head(PVH_SON_DEG[intersect(rownames(PVH_SON_DEG),TF_NP$Transcription_factors),],60) 

PreThal_DEG = FindMarkers(E11_E13_CS22.integrated, ident.1 = paste("Black_",'Prethalamus',sep=''), ident.2 = paste("CS22_",'PreThal',sep=''),features=okGenes)

head(PreThal_DEG[intersect(rownames(PreThal_DEG),TF_NP$Transcription_factors),],20) 


BlackCS22DEGrna[['PVH_SON']] = PVH_SON_DEG
BlackCS22DEGrna[['PreThal']] = PreThal_DEG


for(i in 1:length(BlackCS22DEGrna)){
BlackCS22DEGrna[[i]]$gene = rownames(BlackCS22DEGrna[[i]])
BlackCS22DEGrna[[i]]$region = names(BlackCS22DEGrna)[i]
}

save(BlackCS22DEGrna,file='./Analysis/BlackCS22DEGrna.rda',compress=TRUE)

### save into excel


TMPnames = names(BlackCS22DEGrna)
for(i in TMPnames) {
if(i==TMPnames[1]){
tmp = BlackCS22DEGrna[[i]]

write.xlsx2(tmp, file='./Analysis/E11_E13_CS22_DEGs_by_region.xlsx', sheetName=i,col.names=TRUE, row.names=FALSE, append=FALSE)
} else {
tmp = BlackCS22DEGrna[[i]]

write.xlsx2(tmp, file='./Analysis/E11_E13_CS22_DEGs_by_region.xlsx', sheetName=i,col.names=TRUE, row.names=FALSE, append=TRUE)
}
cat(i)
}


save(E11_E13_CS22.integrated,file='./SeuratObj/E11_E13_CS22_integrated.rda',compress=TRUE)




 DefaultAssay(E11_E13_CS22.integrated) = 'RNA'

pdf('./TestPlots/E11_E13_CS22_Check_ARC_DMH_genes.pdf',width=12,height=8)

for(i in c("sample","MergedRegion","MouseRegion","HumanRegion")){
print(DimPlot(E11_E13_CS22.integrated, reduction = "umap", group.by = i, label = TRUE, repel = TRUE))
}

for(i in c("POMC","CHCHD10","S100A10","SIX6","AGRP","NPY","HDC","LHX8","GPC3","ISL1","TBX3","TBX2","PCP4","JUNB")){
print(FeaturePlot(E11_E13_CS22.integrated, features = i, pt.size = 0.4, reduction = "umap"))
}
dev.off()



#### Trimester 2 Arc comp


load('/local/projects-t3/idea/bherb/Hypothalamus/Blackshaw/Hypo_E10_P45_updated_Neurons.Robj') #object name is Hypo - but just seems to be subseted neurons 

E10P45Neuron = UpdateSeuratObject(Hypo)

pdf('./TestPlots/E10P45black_Check_Neuron_clusters.pdf',width=12,height=8)
print(DimPlot(E10P45Neuron, label = TRUE) + NoLegend())
for(i in c("orig.ident","Cluster")){
print(DimPlot(E10P45Neuron, reduction = "umap", group.by = i, label = TRUE, repel = TRUE))
}
dev.off()


### integrate with Trimester 2 neurons 

HumanSamples = unique(Tri2Neuron_ds.integrated@meta.data$sample)

MouseSamples = unique(E10P45Neuron@meta.data$orig.ident)

CombDat = vector(mode='list',length=length(HumanSamples)+length(MouseSamples))

names(CombDat) = c(HumanSamples,MouseSamples)


for(i in HumanSamples){

tmp = subset(Tri2Neuron_ds.integrated,cells=colnames(Tri2Neuron_ds.integrated)[which(Tri2Neuron_ds.integrated@meta.data$sample==i)])

DefaultAssay(tmp)='RNA'
tmp@meta.data$species='human'
tmp@meta.data$NeuronCluster = tmp@meta.data$hgclust
tmp <- NormalizeData(tmp, verbose = FALSE)
tmp <- FindVariableFeatures(tmp, selection.method = "vst", nfeatures = 1000, verbose = FALSE)
CombDat[[i]] = tmp
cat(paste(i,', ',sep=''))
}


### Gene names! - consider downsampling 
for(i in MouseSamples){

tmp = subset(E10P45Neuron,cells=colnames(E10P45Neuron)[which(E10P45Neuron@meta.data$orig.ident==i)])

tmp2 = subset(BlackAllDat[[i]],cells=colnames(tmp))

DefaultAssay(tmp2)='RNA'
tmp2@meta.data$species='mouse'
tmp2@meta.data$NeuronCluster = tmp@meta.data$Cluster

tmp2 <- NormalizeData(tmp2, verbose = FALSE)
tmp2 <- FindVariableFeatures(tmp2, selection.method = "vst", nfeatures = 1000, verbose = FALSE)
CombDat[[i]] = tmp2
cat(paste(i,', ',sep=''))
}



for(i in c(6:length(CombDat))){
if(length(colnames(CombDat[[i]]))>500) {

CombDat[[i]] = subset(CombDat[[i]],cells = colnames(CombDat[[i]])[sample(c(1:ncol(CombDat[[i]])),500)])
cat(i)
} else {
  next
}
}

save(CombDat,file='./SeuratObj/Blackshaw_Tri2_Neuron_Obj.rda',compress=TRUE)


options(future.globals.maxSize = 4000 * 1024^2)

HM_Neurons.anchors <- FindIntegrationAnchors(object.list = CombDat[c(1:5,12,14,15,16)], dims = 1:20,k.filter=50) # 13 has 10cells
save(HM_Neurons.anchors,file='./Analysis/HM_Neurons_anchors_Tri2_E16_E18_P4_P8.rda',compress=TRUE)
HM_Neurons.integrated <- IntegrateData(anchorset = HM_Neurons.anchors, dims = 1:20)
save(HM_Neurons.integrated,file='./SeuratObj/Blackshaw_Tri2_Neuron_Integrated_Obj.rda',compress=TRUE)

DefaultAssay(HM_Neurons.integrated) <- "integrated"
HM_Neurons.integrated <- ScaleData(HM_Neurons.integrated, verbose = FALSE)
HM_Neurons.integrated <- RunPCA(HM_Neurons.integrated,npcs = 30, verbose = FALSE)
#HM_Neurons.integrated <- RunUMAP(HM_Neurons.integrated, reduction = "pca", dims = 1:30)

HM_Neurons.integrated <- RunTSNE(HM_Neurons.integrated, reduction = "pca", dims = 1:30)




HM_Neurons.integrated@meta.data$hgclust = Tri2Neuron_ds.integrated@meta.data[colnames(HM_Neurons.integrated),'hgclust']



pdf('./TestPlots/Tri2_Black_Neuron_clusters.pdf',width=12,height=8)
print(DimPlot(HM_Neurons.integrated, label = TRUE) + NoLegend())
for(i in c("species","hgclust","NeuronCluster")){
print(DimPlot(HM_Neurons.integrated, reduction = "tsne", group.by = i, label = TRUE, repel = TRUE))
}
dev.off()



save(HM_Neurons.integrated,file='./SeuratObj/Blackshaw_Tri2_Neuron_Integrated_Obj.rda',compress=TRUE)


### not the greatest merging, but if I do again in the future, don't include postnatal samples


## Direct comparison of clusters - metaneighbor? 


#E11_E13_GW18_sce = as.SingleCellExperiment(E11_E13_GW18.integrated,assay='SCT')

#mnresTri2 =  MetaNeighborUS(var_genes = VariableFeatures(E11_E13_GW18.integrated),dat = E11_E13_GW18_sce,study_id = E11_E13_GW18_sce$species,cell_type = E11_E13_GW18_sce$MergedRegion)

### repeat with new cluster Feb25 

Tri2ArcInd = which(Tri2Neuron_ds.integrated@meta.data$nuclei=='ARC')
Tri2ArcMat = Tri2Neuron_ds.integrated@assays$RNA[,Tri2ArcInd]

Tri2ArcPheno = data.frame(Sample_ID=colnames(Tri2Neuron_ds.integrated)[Tri2ArcInd],Celltype = Tri2Neuron_ds.integrated@meta.data[Tri2ArcInd,'hgclust'],Study_ID = "human_Feb25")


BlackArcInd = which(as.numeric(as.character(E10P45Neuron@meta.data$Cluster))>8 & as.numeric(as.character(E10P45Neuron@meta.data$Cluster))<18)

BlackArcMat = E10P45Neuron@assays$RNA[,BlackArcInd]

tmpGene = rownames(BlackArcMat)
tmpGene2 = MMtoHS(tmpGene)
naInd = which(is.na(tmpGene2))
dupInd = which(dupMatch!=c(1:length(tmpGene2)))
rownames(BlackArcMat) = tmpGene2
BlackArcMat = BlackArcMat[-dupInd,]
BlackArcMat = BlackArcMat[-1,] ## first gene is NA

BlackArcPheno = data.frame(Sample_ID=colnames(E10P45Neuron)[BlackArcInd],Celltype = E10P45Neuron@meta.data[BlackArcInd,'Cluster'],Study_ID = "mouse")



compMat2(xmat=Tri2ArcMat,ymat=BlackArcMat,xpheno=Tri2ArcPheno,ypheno=BlackArcPheno,cexRow=0.3,cexCol=0.3)


### Follow - up on metaneighbor result 

## Arc_7 only seems similar to Arc_4? 

### DEG list is muddled, checking with dot plot 


Tri2Arc = subset(Tri2Neuron_ds.integrated,cells = colnames(Tri2Neuron_ds.integrated)[Tri2ArcInd])


Tri2Arc@active.ident=as.factor(Tri2Arc$hgclust)

Arc_NP = sort(unique(c("SST","NPY","CORT","PNOC","AGRP","NPY","SST","SCG2","SCG5","TAC1","POMC","SCG2","PROK2","VIP","TAC1","POMC","UBL5","SCG5","SST","TAC1","DBI","TAC3","DBI","CHGB","TAC1","CHGA","SCG2","UBL5","SCG5","DBI","TAC1","SCG2","NUCB2","SST")))

pdf(file=paste("./TestPlots/Tri2_Arc_NP_genes.pdf",sep=''),width=12,height=12)
print(VlnPlot(Tri2Arc, features = Arc_NP, ncol = 4)
)
dev.off()


pdf(file=paste("./TestPlots/Tri2_Arc_NP_genes_DotPlot.pdf",sep=''),width=12,height=8)
print(DotPlot(Tri2Arc, features = rev(Arc_NP))
)
dev.off()

## make similar plots for mouse

BlackArc = subset(E10P45Neuron,cells = colnames(E10P45Neuron)[BlackArcInd])

BlackArc@active.ident=as.factor(BlackArc$Cluster)

DefaultAssay(BlackArc)='RNA'



pdf(file=paste("./TestPlots/Black_Arc_NP_genes.pdf",sep=''),width=12,height=12)
print(VlnPlot(BlackArc, features = HStoMM(Arc_NP), ncol = 4)
)
dev.off()


pdf(file=paste("./TestPlots/Black_Arc_NP_genes_DotPlot.pdf",sep=''),width=12,height=8)
print(DotPlot(BlackArc, features = HStoMM(rev(Arc_NP))))
dev.off()





