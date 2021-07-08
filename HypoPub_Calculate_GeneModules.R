## changes since RHEL8 install - created a conda environment for R-4.0.3: conda activate r_4.0.3

options(java.parameters = "-Xmx80g")

library( monocle3 )
library( cicero )
library( Seurat )
library( edgeR )
library(xlsx)

setwd('/local/projects-t3/idea/bherb/Hypothalamus/PubRes')

Tri1Genie7003knn5 = read.csv('./Networks/GENIE_Tri1_7003knn5__100G_24T_seed123.csv')

GenieResTri1 = Tri1Genie7003knn5[1:500000,]


Tri2Genie7425knn5 = read.csv('./Networks/GENIE_Tri2_7425knn5__100G_16T_seed123.csv')

GenieResTri2 = Tri2Genie7425knn5[1:500000,]

totVar=unique(c(GenieResTri1$TF,GenieResTri1$target,GenieResTri2$TF,GenieResTri2$target)) #10441

rm(Tri1Genie7003knn5,GenieResTri1,Tri2Genie7425knn5,GenieResTri2)


## capture grns in excel 

totTF = unique(GenieResTri2$TF)

for(i in totTF) {
if(i==totTF[1]){
tmp = GenieResTri2[which(GenieResTri2$TF==i),]
tmp=tmp[,c(2:4)]

write.xlsx2(tmp, file='./Analysis/Tri2_GENIE3_GRN.xlsx', sheetName=i,col.names=TRUE, row.names=FALSE, append=FALSE)
} else {
tmp = GenieResTri2[which(GenieResTri2$TF==i),]
tmp=tmp[,c(2:4)]

write.xlsx2(tmp, file='./Analysis/Tri2_GENIE3_GRN.xlsx', sheetName=i,col.names=TRUE, row.names=FALSE, append=TRUE)
}
cat(i)
gc()
}





gc()

### All samples 
#setwd('/local/projects-t3/idea/bherb/Hypothalamus/PubRes')
load('./SeuratObj/AllSamples_integrated.rda') #All.integrated
obj = All.integrated

rm(All.integrated)
gc()

tmpCells = colnames(obj)[intersect(which(obj@meta.data$Clean_Set=='Yes'),which(is.na(match(obj@meta.data$Consensus_Cell_Type,c('Fibroblast','Microglia','Doublet?','UKN','Epithelial','Blood')))))]


obj = subset(obj,cells=tmpCells,features=totVar)

DefaultAssay(obj) = 'RNA'
cellinfo = obj@meta.data[,c(4,13)]
id='AllSamp'
k=100

metaplot=c('sample','Consensus_Cell_Type') ## metadata to plot from the obj



studies = as.character( unique( cellinfo$sample ))
mtx_full = GetAssayData(obj)
pct.study = matrix(NA,nrow = nrow(mtx_full) , ncol = length(studies) )
for( i in 1:length(studies) ) {
  idx = which( cellinfo$sample == studies[i] )
  sub = mtx_full[ , idx ]
  pct.i = rowSums( sub > 0 ) / ncol(sub)
  pct.study[ , i ] = pct.i
}

#min.pct.study = apply( pct.study , 1 , min ) ## could be a number of Tri2 that are not expressed in Tri1? Skip this step? 

#g = rownames(obj)[ min.pct.study > 0 ]
#obj = subset( obj , features = g )
geneinfo = data.frame( genes = rownames(obj) )
rownames(geneinfo) = rownames(obj)

## fix cell info 
cellinfo = cellinfo[colnames(obj),]

#input_cds <-  suppressWarnings(new_cell_data_set(
 # GetAssayData(obj,slot='counts'),
 # cell_metadata = cellinfo ,
 # gene_metadata = geneinfo  ))

#input_cds <- monocle3::detect_genes(input_cds)
#Ensure there are no peaks (genes?) included with zero reads - according to function, the default setting is just looking for genes that have zero expression across all genes, which should have been previously filtered out.      
#input_cds <- input_cds[Matrix::rowSums(exprs(input_cds)) != 0,]
### make cicero cds (aggregated, k=50)                                          
umap = Embeddings( obj , reduction = 'umap' )

#input_cds=newCellDataSet(cellData = GetAssayData(obj,slot='counts'),phenoData=cellinfo,featureData=geneinfo)


#input_cds = as.CellDataSet(input_cds)
input_cds = as.CellDataSet(obj) ## needs more memory 

cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = umap )
saveRDS( cicero_cds , file = paste('./Analysis/',id,'_cicero_cds_CoreCellGene.rds',sep='' )) ## 10441 genes X 4418 

#counts = exprs(cicero_cds)

counts = cicero_cds@assayData

saveRDS( counts , file = paste('./Analysis/',id,'_cicero_agg_counts_CoreCellGene.rds',sep='' ))

test=ExpressionSet(cicero_cds@assayData)

storageMode(test)='list'

counts = as.matrix(test@assayData$exprs)

#features=unique(c(GenieRes$TF,GenieRes$target))
features = totVar


#counts = counts[ intersect(features,rownames(counts)), ]

logcpm = cpm( counts , log = T )
datExpr=sweep(logcpm,1,apply(logcpm,1,mean),"-")
indx.sd=(apply(logcpm,1,sd))==0 # these will produce NAs
datExpr=sweep(datExpr,1,apply(datExpr,1,sd),"/")
datExpr[indx.sd,]=0
if(sum(is.na(datExpr))!=0){print("NAs in exprsMTX.Z Zscores!")}
#kmeans

cl=kmeans(x=datExpr,centers=k,iter.max=10000,nstart=10,alg="Lloyd")
saveRDS(cl , file = paste('./Analysis/',id,'_k',k,'CoreCellGene.rds',sep='') )

clorig=cl

minModSize = 10
kAEtoStay = 0.3
cutHeight = 0.25

#clusters = readRDS('k25.5perc.rds')
clusters = cl

#datExpr = readRDS('cicero_agg_counts.rds')
datExpr = counts
datExpr = as.matrix( datExpr )
datExpr = datExpr[ names(clusters$cluster) , ]

kAE = cor( t(clusters$centers) , t(datExpr) )

colors = clusters$cluster

for( i in 1:length(colors) ) {
  r.i = kAE[ colors[i] , i ]
  if( r.i < kAEtoStay ) colors[i] = 0
}

size = table( colors )
too.small = as.numeric(names(which(size<minModSize)))
colors[ colors %in% too.small ] = 0


centers = sapply( sort(unique(colors)) , function(i)
  colMeans(datExpr[ colors == i , ]) )

colnames(centers) = paste( 'AE' , sort(unique(colors)) , sep = '' )

r = cor( centers )

d = as.dist( 1 - r )

cl = cutree( hclust( d , method = 'average' ) , h = cutHeight )

mergeColors = rep(NA,length(colors))
for( i in 1:max(cl) ) {
  idx = as.numeric( gsub( 'AE','', names(cl)[ cl == i ] ))
  mergeColors[ colors %in% idx ] = i
}
mergeColors = mergeColors - 1
names(mergeColors) = names(colors)

saveRDS( mergeColors , file = paste('./Analysis/',id,'_mergeColors.k',k,'.cutHeight0.25_CoreCellGene.rds',sep='') )

AE = sapply( sort(unique(mergeColors)) , function(i)
  colMeans(datExpr[ mergeColors == i , ]) )
colnames(AE) = paste( 'AE' , 0:max(mergeColors) , sep = '' )

kAE = cor( t(datExpr) , AE )
sort( kAE[ mergeColors == 1 , 'AE1' ] , decreasing = T )[1:50]

saveRDS( kAE , file = paste('./Analysis/',id,'_mergeColors.k',k,'.cutHeight0.25.kAE_CoreCellGene.rds',sep='') )


## Restart in R 3.6. Some WGCNA dependencies not available in 4.0.2 - not needed now, since using conda 4.0.3

#.libPaths('/local/projects-t3/idea/bherb/software/R_lib/R_3_6')
#library(Seurat,lib.loc='/usr/local/packages/r-3.6.0/lib64/R/library') ## 3.0.1
#library(WGCNA)

library(WGCNA)


## reload definitions at top 

#colors = readRDS(paste('./Analysis/',id,'_mergeColors.k',k,'.cutHeight0.25.rds',sep=''))

mtx = as.matrix(GetAssayData(obj))
mtx = mtx[ names(colors) , ] ##genes intersecting with genie3 (essentially is genie3 genes )
datExpr = t(mtx) ## ok that this will be written over? 

MEs = moduleEigengenes( datExpr , colors , softPower = 1 , excludeGrey = T )

saveRDS( MEs , file = paste('./Analysis/',id,'_mergeColors.k',k,'.MEs_CoreCellGene.rds',sep='') )


for(n in colnames(MEs[['eigengenes']])){

obj@meta.data[,n] <- MEs[['eigengenes']][colnames(obj),n]
cat(paste(n,', ',sep=''))
} 


pdf(paste('./TestPlots/',id,'_Kmean',ncol(MEs[['eigengenes']]),'groupsEigengenes_CoreCellGene.pdf',sep=''),width=12,height=8)
#print(DimPlot(obj, label = TRUE) + NoLegend())
for(x in metaplot){
print(DimPlot(obj, reduction = "umap", group.by = x, label = TRUE, repel = TRUE))
}
for(i in colnames(MEs[['eigengenes']])){
print(FeaturePlot(obj, features = i, pt.size = 0.2, ncol = 2,order=TRUE,cols=c('yellow','red')))
}
dev.off()

## enrichment 
GenieTF = names(table(GenieRes$TF))[which(as.numeric(table(GenieRes$TF))>20)] ## mean 427  - 1122 TF
 
#fisher.test(matrix(c(A,B,C,D),nrow=2,ncol=2),alternative="greater")

#A = cells expressing and in cluster
#B = cells not expressing but in cluster
#C = cells expressing but not in cluster
#D = cells not expressing and not in cluster

kclusts = sort(unique(colors))

GKperOv = data.frame(matrix(0,nrow=length(GenieTF),ncol=length(kclusts)*6))
