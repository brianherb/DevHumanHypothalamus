library(Seurat)
library(xlsx)

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

load('./SeuratObj/AllSamples_integrated.rda') #All.integrated

All.integrated@active.ident=factor(All.integrated$MajorClass)

DefaultAssay(All.integrated) = 'RNA'
All.integrated.MajorClass.DEG = FindAllMarkers(All.integrated,features=rownames(All.integrated))

save(All.integrated.MajorClass.DEG,file='./Analysis/All.integrated.MajorClass.DEG.rda',compress=T)

totDEG = All.integrated.MajorClass.DEG
tissues = unique(totDEG$cluster)

for(i in tissues) {
if(i==tissues[1]){
tmp = totDEG[which(totDEG$cluster==i),]
tmp=tmp[,c(7,2:4,1,5)]
tmp=tmp[which(tmp$p_val_adj<=0.05),]

write.xlsx2(tmp, file='./Analysis/Tri12_CellTypes_DEG.xlsx', sheetName=i,col.names=TRUE, row.names=FALSE, append=FALSE)
} else {
tmp = totDEG[which(totDEG$cluster==i),]
tmp=tmp[,c(7,2:4,1,5)]
tmp=tmp[which(tmp$p_val_adj<=0.05),]

write.xlsx2(tmp, file='./Analysis/Tri12_CellTypes_DEG.xlsx', sheetName=i,col.names=TRUE, row.names=FALSE, append=TRUE)
}
cat(i)
gc()
}


## check top genes

sigDEG = totDEG[which(totDEG$p_val_adj<=0.05),]

sigTF = sigDEG[which(!is.na(match(HSTF$Name,sigDEG$gene))),]



### GO analysis 
library(biomaRt)
library(GOstats)
library(xlsx)

totgene = data.frame(gene=totDEG$gene,celltype = totDEG$cluster,stringsAsFactors=FALSE)

mart = useMart( "ensembl" )
mart = useDataset("hsapiens_gene_ensembl", mart = mart )

gid = getBM(mart = mart, attributes = c("hgnc_symbol","entrezgene_id"),values=totgene) #entrez to refseq

totgene$entrez = gid$entrezgene_id[match(totgene$gene,gid$hgnc_symbol)]

celltypes = unique(totgene$celltype)
GOcat = c('BP','MF','CC')

### biological process
degGO = vector(mode='list',length=length(celltypes)*length(GOcat))

names(degGO) = paste(rep(celltypes,each=3),rep(GOcat,length(celltypes)),sep='_')

for(i in celltypes){
for(k in GOcat){

testGenes = na.omit(totgene$entrez[which(totgene$celltype==i)])
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
 
degGO[[paste(i,k,sep='_')]] <- hyperGTest(params)
cat(paste(i,k,' done', '\n',sep='_'))

}
}

save(degGO,file='./Analysis/Tri12_CellTypes_GO.rda',compress=TRUE)



for(k in GOcat){
TMPnames = names(degGO)[grep(paste('_',k,sep=''),names(degGO))]

for(i in TMPnames) {
if(i==TMPnames[1]){
tmp = summary(degGO[[i]])
tmp$FDR = p.adjust(tmp$Pvalue)
tmp=tmp[,c(1,2,8,3:7)]
tmp=tmp[which(tmp$FDR<=0.05),]

write.xlsx2(tmp, file=paste('./Analysis/Tri12_CellTypes_',k,'_GO.xlsx',sep=''), sheetName=i,col.names=TRUE, row.names=FALSE, append=FALSE)
} else {
tmp = summary(degGO[[i]])
tmp$FDR = p.adjust(tmp$Pvalue)
tmp=tmp[,c(1,2,8,3:7)]
tmp=tmp[which(tmp$FDR<=0.05),]

write.xlsx2(tmp, file=paste('./Analysis/Tri12_CellTypes_',k,'_GO.xlsx',sep=''), sheetName=i,col.names=TRUE, row.names=FALSE, append=TRUE)
}
cat(i)
}
cat(k)
}





library(future)
plan(strategy = "multicore", workers = 24)



All.integrated@active.ident=factor(All.integrated$SubClass)

DefaultAssay(All.integrated) = 'RNA'
All.integrated.SubClass.DEG = FindAllMarkers(All.integrated,features=rownames(All.integrated))


## quick look at oligos 






save(All.integrated.SubClass.DEG,file='./Analysis/All.integrated.SubClass.DEG.rda',compress=T)


All.integrated@active.ident=factor(All.integrated$MajorClass)

DefaultAssay(All.integrated) = 'RNA'

Markers.TFs = FindAllMarkers(All.integrated,features=intersect(HSTF$Name,rownames(All.integrated@assays$RNA)))

Markers.NPs = FindAllMarkers(All.integrated,features=intersect(NPlist2,rownames(All.integrated@assays$RNA)))

Markers.NRs = FindAllMarkers(All.integrated,features=intersect(unique(TF_NP$Neuromodulator_receptors),rownames(All.integrated@assays$RNA)))

Markers.NMs = FindAllMarkers(All.integrated,features=intersect(unique(TF_NP$Neuromodulator_production_and_transport),rownames(All.integrated@assays$RNA)))


MajTypes = sort(unique(All.integrated$MajorClass))

Markers.TFs$perDif = Markers.TFs$pct.1 - Markers.TFs$pct.2

for(i in MajTypes){
if(i==MajTypes[1]){
    tmp = Markers.TFs[which(Markers.TFs$cluster==i ),]
    tmp = tmp[which(tmp$p_val_adj<=0.001 & tmp$pct.2<0.3),]
    tmp = tmp[order(tmp$perDif,decreasing=TRUE),]
    #tmp = tmp[which(tmp$avg_logFC>1),]
topPerDiftf <-  tmp$gene[1:5]
} else {
    tmp = Markers.TFs[which(Markers.TFs$cluster==i ),]
    tmp = tmp[which(tmp$p_val_adj<=0.001 & tmp$pct.2<0.3),]
    tmp = tmp[order(tmp$perDif,decreasing=TRUE),]
    #tmp = tmp[which(tmp$avg_logFC>1),]
topPerDiftf <-  c(topPerDiftf,tmp$gene[1:5])
}
}

topPerDiftf = na.omit(topPerDiftf)

tmpTableTF = table(topPerDiftf)

topPerDiftf2 = unique(topPerDiftf)

topPerDiftf3 = intersect(topPerDiftf2,names(tmpTableTF)[tmpTableTF<3])

All.integrated@active.ident=factor(All.integrated$MajorClass,levels=rev(MajTypes))

pdf('./TestPlots/All_Tri12_MajorClass_DotPlot_TopTF.pdf',width=12,height=8)

DotPlot(All.integrated,features=topPerDiftf3) + RotatedAxis()
dev.off()


### plotting cell type proportion per sample 
data2= table(All.integrated$MajorClass,All.integrated@meta.data$sample)

data3 = data.table::melt(data2)

# Stacked

pdf(file='./TestPlots/All_Tri12_MajorClass_CellType_content.pdf',width=20,height=8)
ggplot(data3, aes(fill=Var1, y=value, x=Var2)) + 
    geom_bar(position="fill", stat="identity")
dev.off()


### check when NP's are turned on 

### only look in DivPro, RG and Neurons? 

### initial Dot plot?  - which NP's are of interest? 

NeuroUniv = subset(All.integrated,cells=colnames(All.integrated)[which(!is.na(match(All.integrated@meta.data$MajorClass,c("Neurons","RadialGlia","Dividing_Progenitor"))))])


NPs = c(intersect(NPlist2,rownames(NeuroUniv@assays$RNA)),'HDC','TSHB','B2M','OTP')

NPclusterPer = data.frame(matrix(0,nrow=length(NPs),ncol=length(unique(NeuroUniv@meta.data$sample))))

rownames(NPclusterPer) = as.character(NPs)
colnames(NPclusterPer) = as.character(unique(NeuroUniv@meta.data$sample))

for(i in colnames(NPclusterPer)){
	for(j in NPs){
	tmpdat = NeuroUniv@assays$RNA[j,]
	tmpdat[which(tmpdat>0)] = 1
tmp = table(tmpdat[which(NeuroUniv@meta.data$sample==i)])
tmpPer = round((tmp/length(which(NeuroUniv@meta.data$sample==i)))*100,2)
if(length(tmpPer)==2){
NPclusterPer[j,i] = tmpPer[2]
}
}
cat(paste(i,', ',sep=""))
}


## mentioned in the paper: 

PubNP = c('POMC','NPY','GHRH','TRH','AVP','GAL','AGRP','KISS1','OXT','CRH','TSHB','SST','HDC','OTP','B2M')


## Heatmap? 

NeuroUniv.scale = NeuroUniv

NeuroUniv.scale <- ScaleData(object = NeuroUniv.scale, features = rownames(NeuroUniv.scale))

pdf('./TestPlots/NP_Expression_acrossDevTime_heatmap.pdf',width=12,height=8)

DoHeatmap(NeuroUniv,features = PubNP,group.by = 'sample',label=FALSE,raster=FALSE) + NoLegend()

dev.off()




## simple Tri1 vs Tri2 - based on NPclusterPer[PubNP,]

Tri1 - 'POMC','GHRH','TRH','AVP','GAL','SST', HDC, OTP,

Tri2 - 'AGRP','KISS1','OXT','CRH','NPY'

Never = 'TSHB',




binaryAll = 

NonZeroCounts = 

NPlistExp = rownames(All.integrated)[which()]


DefaultAssay(Tri2Neuron_ds.integrated)='RNA'

Tri2Neuron_ds.integrated@active.ident=as.factor(Tri2Neuron_ds.integrated$hgclust)

pdf('./TestPlots/Tri2Neuron_ds_HGclust_DotPlot_AllNP.pdf',width=20,height=8)

DotPlot(Tri2Neuron_ds.integrated,features=rev(sort(NPlist2))) + RotatedAxis()
dev.off()






### marker genes from references


load("./SeuratObj/Mus_Ref_12000cells_Indv_Objects.rda")


## Mof

MofSub = RefNormDat[['Mof']]

MofSub=subset(MofSub,cells=colnames(MofSub)[which(MofSub@meta.data$Neuronal_cluster!='')])

DefaultAssay(MofSub) = 'RNA'

MofSub=FindVariableFeatures(MofSub, selection.method = "vst", 
        nfeatures = 5000, verbose = FALSE)

MofSub$Neuronal_cluster = as.character(MofSub$Neuronal_cluster)

MofSub@active.ident=factor(MofSub$Neuronal_cluster)

MofSubDEG = FindAllMarkers(MofSub,features=VariableFeatures(MofSub))

## Chen

ChenSub = RefNormDat[['Chen']]

ChenSub=subset(ChenSub,cells=c(colnames(ChenSub)[grep('GABA',ChenSub@meta.data$Major_Cell_Type)],colnames(ChenSub)[grep('Glu',ChenSub@meta.data$Major_Cell_Type)]))

DefaultAssay(ChenSub) = 'RNA'

ChenSub=FindVariableFeatures(ChenSub, selection.method = "vst", 
        nfeatures = 5000, verbose = FALSE)

ChenSub$Major_Cell_Type = as.character(ChenSub$Major_Cell_Type)

ChenSub@active.ident=factor(ChenSub$Major_Cell_Type)

ChenSubDEG = FindAllMarkers(ChenSub,features=VariableFeatures(ChenSub))



## Mick.male

Mick.maleSub = RefNormDat[['Mick.male']]

Mick.maleSub=subset(Mick.maleSub,cells=colnames(Mick.maleSub)[grep('Neurons',Mick.maleSub@meta.data$Major_Cell_Type)])

DefaultAssay(Mick.maleSub) = 'RNA'

Mick.maleSub=FindVariableFeatures(Mick.maleSub, selection.method = "vst", 
        nfeatures = 5000, verbose = FALSE)

Mick.maleSub$Sub_Cell_Type = as.character(Mick.maleSub$Sub_Cell_Type)

Mick.maleSub@active.ident=factor(Mick.maleSub$Sub_Cell_Type)

Mick.maleSubDEG = FindAllMarkers(Mick.maleSub,features=VariableFeatures(Mick.maleSub))



## Mick.female

Mick.femaleSub = RefNormDat[['Mick.female']]

Mick.femaleSub=subset(Mick.femaleSub,cells=colnames(Mick.femaleSub)[grep('Neurons',Mick.femaleSub@meta.data$Major_Cell_Type)])

DefaultAssay(Mick.femaleSub) = 'RNA'

Mick.femaleSub=FindVariableFeatures(Mick.femaleSub, selection.method = "vst", 
        nfeatures = 5000, verbose = FALSE)

Mick.femaleSub$Sub_Cell_Type = as.character(Mick.femaleSub$Sub_Cell_Type)

Mick.femaleSub@active.ident=factor(Mick.femaleSub$Sub_Cell_Type)

Mick.femaleSubDEG = FindAllMarkers(Mick.femaleSub,features=VariableFeatures(Mick.femaleSub))



## Kim_FC1

Kim_FC1Sub = RefNormDat[['Kim_FC1']]

Kim_FC1Sub=subset(Kim_FC1Sub,cells=c(colnames(Kim_FC1Sub)[grep('GABA',Kim_FC1Sub@meta.data$Major_Cell_Type)],colnames(Kim_FC1Sub)[grep('Gluta',Kim_FC1Sub@meta.data$Major_Cell_Type)]))

DefaultAssay(Kim_FC1Sub) = 'RNA'

Kim_FC1Sub=FindVariableFeatures(Kim_FC1Sub, selection.method = "vst", 
        nfeatures = 5000, verbose = FALSE)

Kim_FC1Sub$Major_Cell_Type = as.character(Kim_FC1Sub$Major_Cell_Type)

Kim_FC1Sub@active.ident=factor(Kim_FC1Sub$Major_Cell_Type)

Kim_FC1SubDEG = FindAllMarkers(Kim_FC1Sub,features=VariableFeatures(Kim_FC1Sub))


## Kim_MC1

Kim_MC1Sub = RefNormDat[['Kim_MC1']]

Kim_MC1Sub=subset(Kim_MC1Sub,cells=c(colnames(Kim_MC1Sub)[grep('GABA',Kim_MC1Sub@meta.data$Major_Cell_Type)],colnames(Kim_MC1Sub)[grep('Gluta',Kim_MC1Sub@meta.data$Major_Cell_Type)]))

DefaultAssay(Kim_MC1Sub) = 'RNA'

Kim_MC1Sub=FindVariableFeatures(Kim_MC1Sub, selection.method = "vst", 
        nfeatures = 5000, verbose = FALSE)

Kim_MC1Sub$Major_Cell_Type = as.character(Kim_MC1Sub$Major_Cell_Type)

Kim_MC1Sub@active.ident=factor(Kim_MC1Sub$Major_Cell_Type)

Kim_MC1SubDEG = FindAllMarkers(Kim_MC1Sub,features=VariableFeatures(Kim_MC1Sub))



## Kim_MC2

Kim_MC2Sub = RefNormDat[['Kim_MC2']]

Kim_MC2Sub=subset(Kim_MC2Sub,cells=c(colnames(Kim_MC2Sub)[grep('GABA',Kim_MC2Sub@meta.data$Major_Cell_Type)],colnames(Kim_MC2Sub)[grep('Gluta',Kim_MC2Sub@meta.data$Major_Cell_Type)]))

DefaultAssay(Kim_MC2Sub) = 'RNA'

Kim_MC2Sub=FindVariableFeatures(Kim_MC2Sub, selection.method = "vst", 
        nfeatures = 5000, verbose = FALSE)

Kim_MC2Sub$Major_Cell_Type = as.character(Kim_MC2Sub$Major_Cell_Type)

Kim_MC2Sub@active.ident=factor(Kim_MC2Sub$Major_Cell_Type)

Kim_MC2SubDEG = FindAllMarkers(Kim_MC2Sub,features=VariableFeatures(Kim_MC2Sub))



## Kim_MC3

Kim_MC3Sub = RefNormDat[['Kim_MC3']]

Kim_MC3Sub=subset(Kim_MC3Sub,cells=c(colnames(Kim_MC3Sub)[grep('GABA',Kim_MC3Sub@meta.data$Major_Cell_Type)],colnames(Kim_MC3Sub)[grep('Gluta',Kim_MC3Sub@meta.data$Major_Cell_Type)]))

DefaultAssay(Kim_MC3Sub) = 'RNA'

Kim_MC3Sub=FindVariableFeatures(Kim_MC3Sub, selection.method = "vst", 
        nfeatures = 5000, verbose = FALSE)

Kim_MC3Sub$Major_Cell_Type = as.character(Kim_MC3Sub$Major_Cell_Type)

Kim_MC3Sub@active.ident=factor(Kim_MC3Sub$Major_Cell_Type)

Kim_MC3SubDEG = FindAllMarkers(Kim_MC3Sub,features=VariableFeatures(Kim_MC3Sub))


## Wen

WenSub = RefNormDat[['Wen']]

WenSub=subset(WenSub,cells=colnames(WenSub)[grep('Neuron',WenSub@meta.data$Major_Cell_Type)])

DefaultAssay(WenSub) = 'RNA'

WenSub=FindVariableFeatures(WenSub, selection.method = "vst", 
        nfeatures = 5000, verbose = FALSE)

WenSub$sampleClust = as.character(WenSub$sampleClust)

WenSub@active.ident=factor(WenSub$sampleClust)

WenSubDEG = FindAllMarkers(WenSub,features=VariableFeatures(WenSub))



## Cam

CamSub = RefNormDat[['Cam']]

CamSub=subset(CamSub,cells=colnames(CamSub)[grep('Neuron',CamSub@meta.data$Major_Cell_Type)])

DefaultAssay(CamSub) = 'RNA'

CamSub=FindVariableFeatures(CamSub, selection.method = "vst", 
        nfeatures = 5000, verbose = FALSE)

CamSub$Sub_Cell_Type = as.character(CamSub$Sub_Cell_Type)

CamSub@active.ident=factor(CamSub$Sub_Cell_Type)

CamSubDEG = FindAllMarkers(CamSub,features=VariableFeatures(CamSub))



## blackshaw

load('/local/projects-t3/idea/bherb/Hypothalamus/Blackshaw/Hypo_E10_P45_updated_Neurons.Robj') #object name is Hypo - but just seems to be subseted neurons 

E10P45Neuron = UpdateSeuratObject(Hypo)

DefaultAssay(E10P45Neuron) = 'RNA'

E10P45Neuron=FindVariableFeatures(E10P45Neuron, selection.method = "vst", 
        nfeatures = 5000, verbose = FALSE)

E10P45Neuron$Cluster = as.character(E10P45Neuron$Cluster)

E10P45Neuron@active.ident=factor(E10P45Neuron$Cluster)

E10P45NeuronDEG = FindAllMarkers(E10P45Neuron,features=VariableFeatures(E10P45Neuron))




save(MofSubDEG,ChenSubDEG,Mick.maleSubDEG,Mick.femaleSubDEG,Kim_FC1SubDEG,Kim_MC1SubDEG,Kim_MC2SubDEG,Kim_MC3SubDEG,WenSubDEG,CamSubDEG,file='./Analysis/Reference_DEGs_1.rda',compress=TRUE)


save(E10P45NeuronDEG,file='./Analysis/Reference_DEGs_2.rda',compress=TRUE)




## writing out tables 

library(xlsx)
options(java.parameters = "-Xmx100g")

tmp=MofSubDEG[which(MofSubDEG$p_val_adj<=0.05),]

write.xlsx2(tmp, file='./Analysis/Mouse_Reference_Neurons_DEG_Mof.xlsx', sheetName="Moffitt_2018",col.names=TRUE, row.names=FALSE, append=FALSE)

tmp=ChenSubDEG[which(ChenSubDEG$p_val_adj<=0.05),]

write.xlsx2(tmp, file='./Analysis/Mouse_Reference_Neurons_DEG_Chen.xlsx', sheetName="Chen_2017",col.names=TRUE, row.names=FALSE, append=FALSE)

tmp=Mick.maleSubDEG[which(Mick.maleSubDEG$p_val_adj<=0.05),]

write.xlsx2(tmp, file='./Analysis/Mouse_Reference_Neurons_DEG_MickMale.xlsx', sheetName="Mickelsen_Male_2019",col.names=TRUE, row.names=FALSE, append=FALSE)

tmp=Mick.femaleSubDEG[which(Mick.femaleSubDEG$p_val_adj<=0.05),]

write.xlsx2(tmp, file='./Analysis/Mouse_Reference_Neurons_DEG_MickFemale.xlsx', sheetName="Mickelsen_Female_2019",col.names=TRUE, row.names=FALSE, append=FALSE)

tmp=Kim_FC1SubDEG[which(Kim_FC1SubDEG$p_val_adj<=0.05),]

write.xlsx2(tmp, file='./Analysis/Mouse_Reference_Neurons_DEG_Kim_FC1.xlsx', sheetName="Kim_FC1_2019",col.names=TRUE, row.names=FALSE, append=FALSE)

tmp=Kim_MC1SubDEG[which(Kim_MC1SubDEG$p_val_adj<=0.05),]

write.xlsx2(tmp, file='./Analysis/Mouse_Reference_Neurons_DEG_Kim_MC1.xlsx', sheetName="Kim_MC1_2019",col.names=TRUE, row.names=FALSE, append=FALSE)

tmp=Kim_MC2SubDEG[which(Kim_MC2SubDEG$p_val_adj<=0.05),]

write.xlsx2(tmp, file='./Analysis/Mouse_Reference_Neurons_DEG_Kim_MC2.xlsx', sheetName="Kim_MC2_2019",col.names=TRUE, row.names=FALSE, append=FALSE)

tmp=Kim_MC3SubDEG[which(Kim_MC3SubDEG$p_val_adj<=0.05),]

write.xlsx2(tmp, file='./Analysis/Mouse_Reference_Neurons_DEG_Kim_MC3.xlsx', sheetName="Kim_MC3_2019",col.names=TRUE, row.names=FALSE, append=FALSE)


tmp=WenSubDEG[which(WenSubDEG$p_val_adj<=0.05),]

write.xlsx2(tmp, file='./Analysis/Mouse_Reference_Neurons_DEG_Wen.xlsx', sheetName="Wen_2020",col.names=TRUE, row.names=FALSE, append=FALSE)

tmp=CamSubDEG[which(CamSubDEG$p_val_adj<=0.05),]

write.xlsx2(tmp, file='./Analysis/Mouse_Reference_Neurons_DEG_Campbell.xlsx', sheetName="Campbell_2017",col.names=TRUE, row.names=FALSE, append=FALSE)


tmp=E10P45NeuronDEG[which(E10P45NeuronDEG$p_val_adj<=0.05),]

write.xlsx2(tmp, file='./Analysis/Mouse_Reference_Neurons_DEG_Blackshaw.xlsx', sheetName="Kim_Blackshaw_2020",col.names=TRUE, row.names=FALSE, append=FALSE)


