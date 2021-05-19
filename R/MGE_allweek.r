library(Seurat)
library(destiny)
library(slingshot)
library(dplyr)
library(data.table)

##(1)#################read data and clustering##################
project<-subset(project,ident=c("MGE1","MGE2"))
project<- FindVariableFeatures(object =project,selection.method = "vst",nfeatures = 1000)
project <- NormalizeData(object =project, verbose = TRUE,scale.factor = 10000)
project <- ScaleData(object = project, verbose = TRUE,
                     features = project@assays$RNA@var.features)
project <- RunPCA(object = project, npcs =30, 
                  verbose = TRUE,
                  features=VariableFeatures(project))
project$batch<-project$Week
project$batch<-gsub('GW09','V3',project$batch)
project$batch<-gsub('GW10','V2',project$batch)
project$batch<-gsub('GW11','V3',project$batch)
project$batch<-gsub('GW12','V3',project$batch)
project <- RunHarmony(project,group.by.vars = c('Week','batch'),
                      dims.use = 1:30)
project <- FindNeighbors(project, reduction = "harmony", dims = 1:30) 
project <- RunTSNE(object = project, reduction = "harmony",dims=1:30)
project<- FindClusters(resolution = 0.3,project)

diffs <- FindAllMarkers(object = project, logfc.threshold = 0.25, 
                        only.pos = TRUE, min.pct = 0.1,assay="RNA")
diffs2<-subset(diffs[grep("^RP[L|S]",diffs$gene, ignore.case = FALSE,invert=TRUE),],
               subset=avg_logFC>0.75 & (pct.1 > pct.2) & pct.1> 0.5)
##(2)####################slingshot analysis######################
ct <- as.ExpressionSet(as.data.frame(Embeddings(project,'harmony')))
ct$celltype<-project@active.ident
dm <- DiffusionMap(ct,sigma='global',k=30,density_norm = T,n_pcs = NA)
sling<-slingshot(data = dm@eigenvectors[,1:2],clusterLabels = MGE@active.ident,
                 start.clus=5,end.clus=c(3,6))

##(3)###################plot######################
##figS7a
cols=c(brewer.pal(12,"Paired"),brewer.pal(8,"Dark2"),brewer.pal(9,"Set1"),brewer.pal(8,"Set2"),brewer.pal(12,"Set3"),brewer.pal(8,"Pastel2"),brewer.pal(9,"Pastel1"),brewer.pal(8,"Accent"))
col=unique(cols)[c(1,5,33,13,16,34,9,36)]
DimPlot(project,reduction = 'tsne',label = F,cols = col)
ggsave('~/results/MGE_celltype_tsne_allweek.pdf',width = 6.6,height = 4)

##figS7c
cell_cluster<-as.data.frame(project@active.ident)
colnames(cell_cluster)<-'Cluster'
cluster<-unique(cell_cluster$Cluster)
for(i in 1:length(cluster)){
  index=which(cell_cluster$Cluster==cluster[i])
  cell_cluster[index,2]<-col[i]
}
pdf('~/results/sling_shot_MGE_allweek.pdf')
plot(dm@eigenvectors[,c(1,2)], pch=20,cex=0.5,col=cell_cluster$V2)
lines(sling@curves$curve1, lwd=2, col='black')
lines(sling@curves$curve3, lwd=2, col='black')
legend('bottomleft',legend=unique(cell_cluster$Cluster),col=unique(cell_cluster$V2),pch=20,xjust=1,yjust=1,cex = 0.7)
dev.off()

##figS7b
levels(project)<-c('5','4','0','1','3','6','7','2')
heatmap_genes<-c('UBE2C','CENPE','CCNB1','ANGPT2','CRABP1','RPH3A','PDE4DIP',
                 'ZEB2','MAF','ERBB4','POU3F2','PDE1A','CNTNAP2','SHH','GSTP1',
                 'FGF19','ZIC1','ZIC2','ZIC4','SLC5A7','RBMS1','GBX2','HS6ST3',
                 'GABRB2','FGF10','LHX8','NKX2-1','BEST3')
MGE_average<-AverageExpression(project,return.seurat = T)
cm<-GetAssayData(MGE_average,slot = 'data')[heatmap_genes,]
cm<-as.data.frame(cm)
pheatmap::pheatmap(heatmap_genes,scale = 'row',cluster_rows = F,cluster_cols = F,border_color = '#3f3f44',
                   color =colorRampPalette(c('#3da4ab','#fcf9ea','#fe8a71'))(100),
                   filename = '~/results//MGE_heatmap_allweek.pdf')
