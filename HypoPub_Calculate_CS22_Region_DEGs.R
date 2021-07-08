
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




##CS22 with region info

load('./SeuratObj/Tri1Samples_integrated.rda') #AllTri1.integrated

load('/local/projects-t3/idea/bherb/Hypothalamus/Seurat_CS22_combined_with_regions.rda') #CS22.integrated


## check genes 

 DefaultAssay(CS22.integrated) = 'RNA'

pdf('./TestPlots/CS22_Check_ARC_DMH_genes.pdf',width=12,height=8)

for(i in c("sample","MajorClass","region")){
print(DimPlot(CS22.integrated, reduction = "umap", group.by = i, label = TRUE, repel = TRUE))
}

for(i in c("VIP","VIPR2","POMC","CHCHD10","S100A10","SOX14","NR5A1","SIX6","AGRP","NPY","HDC","LHX8","GPC3","ISL1","TBX3","TBX2","SST","CITED1","OTP","PCP4","JUNB","ONECUT2","ONECUT3","CARTPT")){
print(FeaturePlot(CS22.integrated, features = i, pt.size = 0.4, reduction = "umap"))
}
dev.off()

pdf('./TestPlots/CS22_GAD2_SLC17A6_genes.pdf',width=12,height=8)

for(i in c("sample","MajorClass","region")){
print(DimPlot(CS22.integrated, reduction = "umap", group.by = i, label = TRUE, repel = TRUE))
}

for(i in c("GAD2","SLC17A6")){
print(FeaturePlot(CS22.integrated, features = i, pt.size = 0.4, reduction = "umap",order=TRUE))
}
dev.off()