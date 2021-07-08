## changes since RHEL8 install - created a conda environment for R-4.0.3: conda activate r_4.0.3

library(Seurat)

setwd('/local/projects-t3/idea/bherb/Hypothalamus/PubRes')

load('./SeuratObj/AllSamples_integrated.rda') #All.integrated

source('./Code/HypoPub_LoadFunctions.R')

## cell type assignments from Hannah

ctBroad = read.csv('./Analysis/AllCells_RadialGliaBroad_19APR21.csv',row.names=1)

## collapse Astro, Oligo, Micro, VLMC

ctBroad$Pop[grep('Astro_',ctBroad$Pop,fixed=TRUE)] = 'Astrocyte'

ctBroad$Pop[grep('Olig',ctBroad$Pop,fixed=TRUE)] = 'Oligodendrocyte'

ctBroad$Pop[grep('Micro_',ctBroad$Pop,fixed=TRUE)] = 'Microglia'

ctBroad$Pop[grep('VLMC_',ctBroad$Pop,fixed=TRUE)] = 'VLMC'

ctBroad$Pop[grep('DividingProgen',ctBroad$Pop,fixed=TRUE)] = 'Dividing_Progenitor'


ctSub = read.csv('./Analysis/AllCells_RadialGliaClustered_19APR21.csv',row.names=1)

ctNeuron = read.csv('./Analysis/FetalAssignments_24FEB21_AprUpdate_HG.csv',row.names=1)

## Region - currently $region in meta.data has CS22 assignments

ctTot = data.frame(MajorClass = ctBroad$Pop,SubClass=ctSub$Pop)
row.names(ctTot) = row.names(ctBroad)

ctTot$Neuron = NA

ctTot[row.names(ctNeuron),'Neuron'] = ctNeuron$Pop

ctTot$Region = flexsplit(ctTot$Neuron,'_')[,1]

ctTot$Region[grep('Preoptic',ctTot$Region,fixed=TRUE)] = 'POA'

cs22ind = which(!is.na(All.integrated@meta.data$region))

ctTot[colnames(All.integrated)[cs22ind],'Region'] = as.character(All.integrated@meta.data$region[cs22ind])

ctTot$Region[grep('PVH_SON',ctTot$Region,fixed=TRUE)] = 'PVH'


write.csv(ctTot,file='./Analysis/CellTypeAssignments.csv')

All.integrated@meta.data$MajorClass = as.character(ctTot[colnames(All.integrated),'MajorClass'])

All.integrated@meta.data$SubClass = as.character(ctTot[colnames(All.integrated),'SubClass'])

All.integrated@meta.data$Neuron = as.character(ctTot[colnames(All.integrated),'Neuron'])

All.integrated@meta.data$Region = as.character(ctTot[colnames(All.integrated),'Region'])

## Clean set eliminates doublets and unknown cells, Skipping this and using Hannah's assignments. 
#All.integrated = subset(All.integrated, cells = colnames(All.integrated)[which(All.integrated@meta.data$Clean_Set=='Yes')])

## test plots 

pdf('./TestPlots/All_Major_Sub_clusters_Hannah.pdf',width=12,height=8)
print(DimPlot(All.integrated, label = TRUE) + NoLegend())
for(i in c("sample","MajorClass","SubClass","Neuron","Region","Phase","region","Maj_Clust_From_Ref","Consensus_Cell_Type","neuron_name")){
print(DimPlot(All.integrated, reduction = "umap", group.by = i, label = TRUE, repel = TRUE))
}
dev.off()

save(All.integrated,file='./SeuratObj/AllSamples_integrated.rda',compress=TRUE)


## apply cell type assignments to Tri2, TRi1 and neurons 

load('./SeuratObj/Tri2Samples_integrated.rda') #AllTri2.integrated


AllTri2.integrated@meta.data$MajorClass = as.character(ctTot[colnames(AllTri2.integrated),'MajorClass'])

AllTri2.integrated@meta.data$SubClass = as.character(ctTot[colnames(AllTri2.integrated),'SubClass'])

AllTri2.integrated@meta.data$Neuron = as.character(ctTot[colnames(AllTri2.integrated),'Neuron'])

AllTri2.integrated@meta.data$Region = as.character(ctTot[colnames(AllTri2.integrated),'Region'])

pdf('./TestPlots/AllTri2_Major_Sub_clusters_Hannah.pdf',width=12,height=8)
#print(DimPlot(AllTri2.integrated, label = TRUE) + NoLegend())
for(i in c("sample","MajorClass","SubClass","Neuron","Region","Maj_Clust_From_Ref","neuron_name")){
print(DimPlot(AllTri2.integrated, reduction = "umap", group.by = i, label = TRUE, repel = TRUE))
}
dev.off()

save(AllTri2.integrated,file='./SeuratObj/Tri2Samples_integrated.rda',compress=TRUE)


##  Tri1 

load('./SeuratObj/Tri1Samples_integrated.rda') #AllTri1.integrated

AllTri1.integrated@meta.data$MajorClass = as.character(ctTot[colnames(AllTri1.integrated),'MajorClass'])

AllTri1.integrated@meta.data$SubClass = as.character(ctTot[colnames(AllTri1.integrated),'SubClass'])

AllTri1.integrated@meta.data$Region = as.character(ctTot[colnames(AllTri1.integrated),'Region'])

pdf('./TestPlots/AllTri1_Major_Sub_clusters_Hannah.pdf',width=12,height=8)
#print(DimPlot(AllTri1.integrated, label = TRUE) + NoLegend())
for(i in c("sample","MajorClass","SubClass","Region","Consensus_Cell_Type","Phase")){
print(DimPlot(AllTri1.integrated, reduction = "umap", group.by = i, label = TRUE, repel = TRUE))
}
dev.off()

save(AllTri1.integrated,file='./SeuratObj/Tri1Samples_integrated.rda',compress=TRUE)

## Tri2 example Neurons 

load('./SeuratObj/Tri2Neuron_ds.integrated.rda')  ## Tri2Neuron_ds.integrated


Tri2Neuron_ds.integrated@meta.data$MajorClass = as.character(ctTot[colnames(Tri2Neuron_ds.integrated),'MajorClass'])

Tri2Neuron_ds.integrated@meta.data$SubClass = as.character(ctTot[colnames(Tri2Neuron_ds.integrated),'SubClass'])

Tri2Neuron_ds.integrated@meta.data$Neuron = as.character(ctTot[colnames(Tri2Neuron_ds.integrated),'Neuron'])

Tri2Neuron_ds.integrated@meta.data$Region = as.character(ctTot[colnames(Tri2Neuron_ds.integrated),'Region'])

pdf('./TestPlots/Tri2Neuron_ds_Major_Sub_clusters_Hannah.pdf',width=12,height=8)
#print(DimPlot(Tri2Neuron_ds.integrated, label = TRUE) + NoLegend())
for(i in c("sample","MajorClass","SubClass","Neuron","Region","Neuron_Class","cluster_name_TF")){
print(DimPlot(Tri2Neuron_ds.integrated, reduction = "umap", group.by = i, label = TRUE, repel = TRUE))
}
dev.off()

save(Tri2Neuron_ds.integrated,file='./SeuratObj/Tri2Neuron_ds.integrated.rda',compress=TRUE)



## CS22 

load('/local/projects-t3/idea/bherb/Hypothalamus/Seurat_CS22_combined_with_regions.rda')

## test plot, then use hannah major class 

## to match cell names, dump _1  and _2 in CS22 object and replace with CS22_1 and CS22_2 at begining

## region in CS22.integrated is old 


newCellNames = paste(as.character(CS22.integrated@meta.data$sample),as.character(flexsplit(colnames(CS22.integrated),'_')[,1]),sep='_')
 
CS22.integrated@meta.data$MajorClass = All.integrated@meta.data[newCellNames,'MajorClass']

CS22.integrated@meta.data$SubClass = All.integrated@meta.data[newCellNames,'SubClass']

CS22.integrated@meta.data$region = NA

CS22.integrated@meta.data$region = All.integrated@meta.data[newCellNames,'region']

CS22.integrated@meta.data$MergedClasses = CS22.integrated@meta.data$MajorClass

EIind = which(CS22.integrated@meta.data$predicted.id=='Excitatory' | CS22.integrated@meta.data$predicted.id=='Inhibitory')

EIind = intersect(EIind,which(CS22.integrated@meta.data$MajorClass=='Neurons'))

CS22.integrated@meta.data$MergedClasses[EIind] = CS22.integrated@meta.data$predicted.id[EIind] 

## get rid of Ambiguius and NA, and other 'predicted.id' terms - Astrocytes, 

pdf('./TestPlots/CS22_Major_Sub_clusters_Hannah.pdf',width=12,height=8)
#print(DimPlot(AllTri1.integrated, label = TRUE) + NoLegend())
for(i in c("sample","MajorClass","SubClass","region","predicted.id","MergedClasses")){
print(DimPlot(CS22.integrated, reduction = "umap", group.by = i, label = TRUE, repel = TRUE, cells=colnames(CS22.integrated)[which(!is.na(CS22.integrated@meta.data$MergedClasses))]))
}
dev.off()


### regions 


## double check VMH - Blackshaw - FEZF1, DEGs- RPS2, NEGR1, RPLP1, HNRNPA1

pdf('./TestPlots/Check_CS22_VMH.pdf')
print(DimPlot(CS22.integrated, reduction = "umap", group.by = 'region', label = TRUE, repel = TRUE))
for(i in c('FEZF1','RPS2','NEGR1','RPLP1','HNRNPA1')){
print(FeaturePlot(CS22.integrated, features = i, pt.size = 0.4,order=TRUE))
}
dev.off()

## redo top region markers 




