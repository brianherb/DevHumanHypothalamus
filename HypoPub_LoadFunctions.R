
source('/home/sament/knn-smoothing/knn_smooth.R')

## custom functions 


### old metaneighbor function 
source("/home/bherb/Rlib/MetaNeighbor-master/2017-08-28-runMN-US.R")
load("/home/bherb/Rlib/MetaNeighbor-master/MetaNeighbor_US_data.Rdata")

#######################
### functions

compMat2 <- function(xmat,ymat,xpheno,ypheno,hits=TRUE,fullResults=TRUE,...){

CombGenes <- intersect(rownames(xmat),rownames(ymat))

CombMatrix <- cbind(xmat[CombGenes,],ymat[CombGenes,])

CombPheno <- rbind(xpheno[,c("Sample_ID","Celltype","Study_ID")],ypheno[,c("Sample_ID","Celltype","Study_ID")])

Combcelltypes <- unique(CombPheno$Celltype)

Comb.var.genes=get_variable_genes(data=CombMatrix, pheno=CombPheno)

Comb.celltype.NV=run_MetaNeighbor_US(vargenes=Comb.var.genes, data=CombMatrix, pheno=CombPheno, celltypes=Combcelltypes)

cols=rev(colorRampPalette(brewer.pal(11,"RdYlBu"))(100))
breaks=seq(0,1,length=101)

    pdf(paste(xpheno$Study_ID[1],ypheno$Study_ID[1],"heatmap.pdf",sep="_"))
    #par(mar=c(10, 2, 2, 6))
heatmap.2(Comb.celltype.NV,trace="none",density.info="none",col=cols,breaks=breaks,margins=c(12,12),srtCol=45)
dev.off()

if(hits==TRUE){
    top_hits=get_top_hits(Comb.celltype.NV,pheno=CombPheno,threshold=0.9,filename=paste(xpheno$Study_ID[1],ypheno$Study_ID[1],"top_hits.txt",sep="_"))
    }

if(fullResults==TRUE){
    all_hits=get_top_hits(Comb.celltype.NV,pheno=CombPheno,threshold=0,filename=paste(xpheno$Study_ID[1],ypheno$Study_ID[1],"all_results.txt",sep="_"))
    }
}




flexsplit <- function(dat,char){
test=strsplit(as.character(dat),char,fixed=TRUE)
n.obs <- sapply(test, length)
seq.max <- seq_len(max(n.obs))
mat <- t(sapply(test, "[", i = seq.max))
mat
}

subplots <- function(sobj,colName,filename){
classes = unique(sobj@meta.data[colName])[,1]
pdf(filename)
for(i in classes){
sobj@meta.data[i]="Other"
sobj@meta.data[which(sobj@meta.data[colName]==i),i]=i
print(DimPlot(sobj, reduction = "umap", group.by = i,order=c(i,'other')) + ggplot2::labs(title=i))
}
dev.off()
}

#HStoMM <- function(x){
#    substr(x,2,nchar(x)) <- tolower(substr(x,2,nchar(x)))
#    x[1]
#}

MM2HSref = read.csv('/local/projects-t3/idea/bherb/annotation/MM2HS_EG100.csv')
colnames(MM2HSref) = c('MM_ID','MM_Gene','HS_ID','HS_Gene')
HS2MMref = read.csv('/local/projects-t3/idea/bherb/annotation/HS2MM_EG100.csv')
colnames(HS2MMref) = c('HS_ID','HS_Gene','MM_ID','MM_Gene')


HStoMM <- function(x){
    tmpind = match(as.character(x),HS2MMref$HS_Gene)
return(as.character(HS2MMref$MM_Gene[tmpind]))
}

MMtoHS <- function(x){
    tmpind = match(as.character(x),MM2HSref$MM_Gene)
return(as.character(MM2HSref$HS_Gene[tmpind]))
}


addMofMeta = function(x,neuronName='CellType'){
MofInhib = read.csv('/local/projects-t3/idea/bherb/Hypothalamus/Moffitt/Inhib_hierarchy.csv',fill=T,header=T)
MofExcit = read.csv('/local/projects-t3/idea/bherb/Hypothalamus/Moffitt/Excit_hierarchy.csv',fill=T,header=T)
for(k in c('level1','level2','level3','level4','level5')){
MofInhib[,k] = paste('Inhib_',MofInhib[,k],sep="")
MofExcit[,k] = paste('Excit_',MofExcit[,k],sep="")
}
TotNeuro = rbind(MofInhib,MofExcit)
tmpmeta = x@meta.data
newCol = c('Neuron_level_1','Neuron_level_2','Neuron_level_3','Neuron_level_4','Neuron_level_5','Neuron_Region','Neuron_Driver')
refCol = colnames(TotNeuro)[2:8]
for(i in newCol){
    tmpmeta[,i] <- NA
}
tmpind = base::match(tmpmeta[,neuronName],TotNeuro$NeuronType)
for(j in 1:7){
tmpmeta[,newCol[j]] = TotNeuro[tmpind,refCol[j]]
}
x@meta.data = tmpmeta
return(x)
}

getDescendants<-function(tree,node,curr=NULL){
  if(is.null(curr)) curr<-vector()
  daughters<-tree$edge[which(tree$edge[,1]==node),2]
  curr<-c(curr,daughters)
  w<-which(daughters>=length(tree$tip))
  if(length(w)>0) for(i in 1:length(w))
    curr<-getDescendants(tree,daughters[w[i]],curr)
  return(curr)
}

getTips = function(tree=NA,node=NA){
tmp = getDescendants(tree,node)
tmp2 = tmp[which(tmp<=length(tree$tip.label))]
return(tmp2)
}

## 

TF_NP = read.csv('/local/projects-t3/idea/bherb/Hypothalamus/TF_NP_list_from_Moffitt.csv') #direct from Moffitt suppmental materials 

for(i in 1:ncol(TF_NP)){
    TF_NP[,i] = MMtoHS(TF_NP[,i])
} ## human convention 

HSTF = read.csv('/local/projects-t3/idea/bherb/annotation/Hsapiens/Human_TF.csv')

NPlist = read.csv('/local/projects-t3/idea/bherb/Hypothalamus/Neuropeptide_list.csv',header=FALSE) ## recieved from Hannah on 7/15/20

NPlist = MMtoHS(as.character(NPlist[,1]))

NPlist2 = unique(c(na.omit(as.character(TF_NP$Neuropeptides)),NPlist))

## cell cycle genes loaded in with seurat 
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
