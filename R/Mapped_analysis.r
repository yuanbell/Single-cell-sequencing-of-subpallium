library(Seurat)
library(ggplot2)
library(ggalluvial)
library(dplyr)
library(data.table)
library(pbapply)
##(1)##########################read data and preprocess#############################
###set analysis time point
GW='GW09'
#GW='GW11'
##read GE dataset
project<-readRDS('~/results/project.rds')
##pick up GW09 or GW11 cells and clusters
project<-subset(project,cells = row.names(subset(project@meta.data,
                                                 project@meta.data$Week==GW)))
project<-subset(project,idents=c('MGE1','MGE2','CGE1','LGE1','LGE2','LGE3','P1','P2','P3','P4','P5','P6'))
##read M1 dataset
M1_project<-readRDS('~/data/M1_project.rds')
##set cell labels
large_celltype<-c('PVALB','PVALB','PVALB','SST','SST','SST','SST',
                  'VIP','VIP',"ID2",'ID2','NDNF','NDNF')
names(large_celltype)<-levels(M1_project)
M1_project<-RenameIdents(M1_project,large_celltype)

###
options(future.globals.maxSize = 10000 * 1024^2)
###merge cell matrix
common_genes<-intersect(row.names(project@assays$RNA@counts),
                        row.names(M1_project@assays$RNA@counts))
cm<-cbind(project@assays$RNA@counts[common_genes,],
          M1_project@assays$RNA@counts[common_genes,])

##(2)#######################integrated analysis################################
###Make Seurat project
project.integrated<-CreateSeuratObject(counts =cm)
remove(cm)
project.integrated$Stage<-rep(c(GW,'Adults M1'),
                              times=c(ncol(project),ncol(M1_project)))
project.integrated$orig.celltype<-c(as.character(project$celltype),
                                    as.character(M1_project@active.ident))
###SCTransform to normalize
project.integrated <- SplitObject(project.integrated, split.by = "Stage")
for (i in 1:length(project.integrated)) {
  project.integrated[[i]] <- SCTransform(project.integrated[[i]])
}
###integrated dataset
project.features <- SelectIntegrationFeatures(object.list = project.integrated, nfeatures = 3000)
project.integrated <- PrepSCTIntegration(object.list = project.integrated, anchor.features = project.features)
pancreas.anchors <- FindIntegrationAnchors(object.list = project.integrated, normalization.method = "SCT", 
                                           anchor.features = project.features,k.anchor = 5,
                                           k.filter = 200)
project.integrated <- IntegrateData(anchorset = pancreas.anchors, normalization.method = "SCT")
###RunPCA and Umap
project.integrated <- RunPCA(project.integrated, verbose = FALSE)
project.integrated <- RunUMAP(project.integrated, dims = 1:20)

DimPlot(project.integrated,group.by = 'orig.celltype',reduction = 'umap',split.by = 'Stage')
rm(pancreas.anchors)
gc()

##(3)##########################Mapped cells#################################
##this function reference from https://github.com/mayer-lab/Mayer-et-al-2018_IntegratedAnalysis/tree/master/R
MapCells=function(object, unknown.id="GW09", known.id="Adults_M1", num.k=10, thresh.require=0.9,map.col="celltype",new.col="Mapped") {
  input.dims=Embeddings(object,reduction= "pca")[,1:20]
  input.dist=as.matrix(dist(input.dims))
  cells.end <- colnames(object[,object@meta.data$Stage == known.id])
  cells.map=colnames(object[,object@meta.data$Stage== unknown.id])
  map.id=pbsapply(cells.map,function(x)
    names(which(sort(table(object@meta.data[names(sort(input.dist[x,cells.end])[1:num.k]),map.col]),
                     decreasing = T)>=num.k*thresh.require))[1])
  nn.same=pbsapply(cells.map, function(x) length(intersect(names(sort(input.dist[x,])[1:num.k]),cells.map)))
  print(table(map.id))
  map.id[is.na(map.id)]="Unassigned"
  map.id[names(which(nn.same>=num.k*1))]="Unassigned"
  print(table(map.id))
  object@meta.data[new.col]=c(as.character(map.id),as.character(object@meta.data[cells.end,]$orig.celltype))
  return(object)
}

project.integrated<-MapCells(project.integrated,unknown.id = GW,known.id = 'Adults M1',
                             map.col = 'orig.celltype',new.col = 'Mapped')

project.integrated@meta.data$MapStage <- apply(project.integrated@meta.data, MARGIN = 1, FUN = function(x) {
  if (x["Stage"] == GW) {
    paste(x["Stage"], x["Mapped"])
  } else {
    x["Stage"]
  }
})

project.integrated@meta.data$MapStage <- factor(project.integrated@meta.data$MapStage, 
                                      levels = c(paste(GW,c("SST","PVALB", "VIP", 
                                                 "ID2", "NDNF")), paste(GW,"Unassigned"), 
                                                 "Adults M1"))
##find_conserved_genes
Idents(project.integrated)<-project.integrated$Mapped
find_conserved_genes<-function(ident1,ident2,filename){
  tmp<-FindConservedMarkers(project.integrated,ident.1 = ident1,ident.2 = ident2,
                            slot = 'data',assay = 'SCT',grouping.var = 'Stage')
  tmp<-subset(tmp,tmp$max_pval<0.01)
  write.csv(tmp,filename)
  return(tmp)
}

SSTvsPV.conserved<-find_conserved_genes(ident1 = 'SST',ident2 = 'PVALB',
                                        filename = paste0('~/M1_GE_map/',GW,'_SSTvsPV.csv'))
VIPvsID2.conserved<-find_conserved_genes(ident1 = 'VIP',ident2 = 'ID2',
                                         filename = paste0('~/M1_GE_map/',GW,'_VIPvsID2.csv'))
VIPvsNDNF.conserved<-find_conserved_genes(ident1 = 'VIP',ident2 = 'NDNF',
                                          filename =paste0('~/M1_GE_map/',GW,'_VIPvsNDNF.csv'))
ID2vsNDNF.conserved<-find_conserved_genes(ident1 = 'ID2',ident2 = 'NDNF',
                                          filename =paste0('~/M1_GE_map/',GW,'_ID2vsNDNF.csv'))

##map subcelltype
project.integrated$orig.subcelltype<-c(as.character(project$celltype),as.character(M1_project$celltype))

MapCells2=function(object, unknown.id="GW09", known.id="Adults_M1", num.k=10, thresh.require=0.9,map.col="celltype",new.col="Mapped") {
  input.dims=Embeddings(object,reduction= "pca")[,1:20]
  input.dist=as.matrix(dist(input.dims))
  cells.end <- colnames(object[,object@meta.data$Stage == known.id])
  cells.map=colnames(object[,object@meta.data$Stage== unknown.id])
  map.id=pbsapply(cells.map,function(x)
    names(which(sort(table(object@meta.data[names(sort(input.dist[x,cells.end])[1:num.k]),map.col]),
                     decreasing = T)>=num.k*thresh.require))[1])
  nn.same=pbsapply(cells.map, function(x) length(intersect(names(sort(input.dist[x,])[1:num.k]),cells.map)))
  print(table(map.id))
  map.id[is.na(map.id)]="Unassigned"
  map.id[names(which(nn.same>=num.k*1))]="Unassigned"
  print(table(map.id))
  object@meta.data[new.col]=c(as.character(map.id),as.character(object@meta.data[cells.end,]$orig.subcelltype))
  return(object)
}
project.integrated<-MapCells2(project.integrated,unknown.id = GW,known.id = 'Adults M1',
                             map.col = 'orig.subcelltype',new.col = 'Mapped2')
project.integrated@meta.data$MapStage2 <- apply(project.integrated@meta.data, MARGIN = 1, FUN = function(x) {
  if (x["Stage"] == GW) {
    paste(x["Stage"], x["Mapped2"])
  } else {
    x["Stage"]
  }
})

##(4)########################plot##########################
##fig6d or g
unsign_cells<-row.names(subset(project.integrated@meta.data,Mapped=='Unassigned'))
sign_cell<-setdiff(colnames(project),unsign_cells)
pdf(paste0('~/results/',GW,'_integrated1.pdf'),width = 8,height = 6) 
DimPlot(project.integrated,pt.size = 1,group.by ='MapStage',
        cells =sign_cell,reduction = 'umap',
        cols = c('#FFD800','#46B5D1','#8CBA51','#A593E0','#D68189'))
dev.off()


##fig6b or f
project.integrated$M1_celltype<-c(as.character(project@meta.data$Week),
                                  as.character(project.integrated@meta.data$orig.celltype[(ncol(project)+1):ncol(project.integrated)]))
pdf(paste0('~/results/',GW,'_integrated2.pdf'),width = 7.5,height = 6) 
DimPlot(project.integrated, reduction = "umap",pt.size = 0.5,group.by ='M1_celltype',
        order = rev(c("SST", "PVALB", "VIP", "ID2","NDNF",GW)),
        cols = c('#FFD800','#46B5D1','#8CBA51','#A593E0','#D68189','#FFE9DE'))
dev.off()

##import region information
meta.data<-project.integrated@meta.data
meta.data<-subset(meta.data,meta.data$Stage==GW)
MGE_celltype<-read.csv('~/mapping/MGE_celltype.csv',row.names = 1)
LGE_celltype<-read.csv('~/mapping/LGE_celltype.csv',row.names = 1)
MGE_celltype$x<-as.character(MGE_celltype$x)
LGE_celltype$x<-as.character(LGE_celltype$x)
mMGE_cells<-row.names(MGE_celltype)[which(MGE_celltype$x=='mMGE')]
lMGE_cells<-row.names(MGE_celltype)[which(MGE_celltype$x=='lMGE')]
MGE_cells<-row.names(MGE_celltype)[which(MGE_celltype$x=='MGE')]
dLGE_cells<-row.names(LGE_celltype)[which(LGE_celltype$x=='dLGE')]
LGE_cells<-row.names(LGE_celltype)[which(LGE_celltype$x=='LGE')]
meta.data$orig.celltype<-as.character(meta.data$orig.celltype)
meta.data$orig.celltype[row.names(meta.data) %in% mMGE_cells]<-'mMGE'
meta.data$orig.celltype[row.names(meta.data) %in% lMGE_cells]<-'lMGE'
meta.data$orig.celltype[row.names(meta.data) %in% MGE_cells]<-'MGE'
meta.data$orig.celltype[row.names(meta.data) %in% dLGE_cells]<-'dLGE'
meta.data$orig.celltype[row.names(meta.data) %in% LGE_cells]<-'LGE'
meta.data$Mapped<-as.character(meta.data$Mapped)

##calculate cellcluster percent
celltype_table<-table(meta.data$orig.celltype,meta.data$Mapped)[,-5]
for (i in 1:5) {
  sum_count<-sum(celltype_table[,i])
  for (j in 1:12) {
    celltype_table[j,i]=(celltype_table[j,i]/sum_count)*100
  }
}
celltype_table<-as.data.frame(celltype_table)

celltype_table$Var1<-factor(celltype_table$Var1,levels=c('P1','P2','P3','P4','P5','P6','lMGE','MGE',
                               'mMGE','CGE1','LGE','dLGE'))

##fig6i or j
celltype_table<-table(meta.data$orig.celltype,meta.data$Mapped)[,-5]
for (i in 1:ncol(celltype_table)) {
  percent <- round(celltype_table[,i]/sum(celltype_table[,i])*100, 1)
  label <- paste(rownames(celltype_table), "(", percent, "% )")
  pdf(paste0('~/M1_GE_map/',colnames(celltype_table)[i],'_',GW,'_pie.pdf'))
  pie(celltype_table[,i],col = c('#ffbd69','#f6def6','#726a95','#900c3f','#fc9d9d','#ffcac2',
                                 '#f4f7c5','#aacdbe'),labels = label,
      border="white",main = colnames(celltype_table)[i])
  dev.off()
}

##figS15
Avg_project<-AverageExpression(project.integrated,add.ident = 'Stage',return.seurat = T,
                               assays = 'SCT')
data<-as.matrix(Avg_project@assays$SCT@data)
data<-data[,1:6]
data<-t(apply(data,1,scale))
colnames(data)<-colnames(Avg_project@assays$SCT@data)[1:6]
data<-data[,paste0(c('Unassigned_','SST_','PVALB_',
              'VIP_','ID2_','NDNF_'),GW)]
SSTvsPV<-data[row.names(SSTvsPV.conserved[order(SSTvsPV.conserved[,2]),]),c(2,3)]
VIPvsID2<-data[row.names(VIPvsID2.conserved[order(VIPvsID2.conserved[,2]),]),c(4,5)]
VIPvsNDNF<-data[row.names(VIPvsNDNF.conserved[order(VIPvsNDNF.conserved[,2]),]),c(4,6)]
ID2vsNDNF<-data[row.names(ID2vsNDNF.conserved[order(ID2vsNDNF.conserved[,2]),]),c(5,6)]

pheatmap::pheatmap(SSTvsPV,cluster_rows = F,cluster_cols = F,color =colorRampPalette(c('#3ca4ab','#fcf8e8','#fe8a71'))(100),
                   border_color = '#ededed',filename =paste0('~/M1_GE_map/','SSTvsPV_',GW,'_heatmap.pdf'),width = 3,height = 7.3)
pheatmap::pheatmap(VIPvsID2,cluster_rows = F,cluster_cols = F,color =colorRampPalette(c('#3ca4ab','#fcf8e8','#fe8a71'))(100),
                   border_color = '#ededed',filename =paste0('~/M1_GE_map/','VIPvsID2_',GW,'_heatmap.pdf'),width = 3,height = 7.3)
pheatmap::pheatmap(VIPvsNDNF,cluster_rows = F,cluster_cols = F,color =colorRampPalette(c('#3ca4ab','#fcf8e8','#fe8a71'))(100),
                   border_color = '#ededed',filename = paste0('~/M1_GE_map/','VIPvsNDNF_',GW,'_heatmap.pdf'),width = 3,height = 7.3)
pheatmap::pheatmap(ID2vsNDNF,cluster_rows = F,cluster_cols = F,color =colorRampPalette(c('#3ca4ab','#fcf8e8','#fe8a71'))(100),
                   border_color = '#ededed',filename = paste0('~/M1_GE_map/','ID2vsNDNF_',GW,'_heatmap.pdf'),width = 3,height = 7.3)


##figS14
project.integrated$Stage<-factor(project.integrated$Stage,levels = c(GW,'Adults M1'))
features<-c('SST','ERBB4','VIP','NRG1','RELN')
for (i in 1:length(features)) {
  FeaturePlot(project.integrated,features = features[i],
              split.by = 'Stage',order = T)
  ggsave(paste0('~/M1_GE_map/',GW,'_',features[i],'_feature_plot.pdf'),width = 9.03,height = 4.47)
}



##figS13b or c
unsign_cells<-row.names(subset(project.integrated@meta.data,Mapped2=='Unassigned'))
sign_cell<-setdiff(colnames(project),unsign_cells)

project.integrated@meta.data$MapStage2 <- factor(project.integrated@meta.data$MapStage2, 
                                                levels = c(paste(GW,levels(M1_project)), paste(GW,'Unassigned'), 
                                                           "Adults M1"))
pdf(paste0('~/results/',GW,'_integrated_subcelltype1.pdf'),width =8,height = 5) 
DimPlot(project.integrated,pt.size = 2,group.by ='MapStage2',
        cells =sign_cell,reduction = 'umap',
        cols = c('#A5DFF9','#1EC0FF','#6AAFE6','#F6EA8C','#FFB140','#CE6D39','#EC7357',
                 '#8FBC94','#C5E99B','#A593E0','#9055A2','#D499B9','#D68189'))
dev.off()

##fig7e or h
meta.data<-subset(meta.data,meta.data$Stage==GW)
data2<-table(meta.data$orig.celltype,meta.data$Mapped)[,-5]
frac<-vector()
for (i in 1:ncol(data2)) {
  tmp_frac<-round(data2[,i]/sum(data2[,i])*100,1)
  frac<-c(frac,tmp_frac)
}
data2<-as.data.frame(data2)
data2$Frac<-frac
colnames(data2)<-c('celltype','Mapped_celltype','Count','Fraction')
data2$celltype<-factor(data2$celltype,levels = c('P1','P2','P3','P4','P5','P6','MGE','mMGE','lMGE','CGE1','LGE',
                            'dLGE'))
data2$Mapped_celltype<-factor(data2$Mapped_celltype,levels =c('SST','PVALB','VIP','ID2','NDNF'))
data2<-subset(data2,data2$Count>0)
ggplot(data2,
       aes(y = Count, axis1 =celltype, axis2= Mapped_celltype)) +
  geom_alluvium(aes(fill = celltype), width = 1/3,alpha=0.7) +
  geom_stratum(width = 1/3, aes(fill =celltype), color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c(GW, 'M1'), expand = c(.05, .05))+
  scale_fill_manual(values = c('#ffbd69','#f6def6','#726a95','#900c3f','#fc9d9d','#ffcac2',
                               '#f4f7c5','#f4f7c5','#f4f7c5','#f4f7c5','#aacdbe'))+
  theme_bw()
ggsave(paste0('~/resulst/',GW,'_Alluvial_Plots.pdf'))  

##figS13d
table<-table(project.integrated$Mapped2[1:ncol(project)])
frac_donutchart<-function(celltype,filename,colors){
  table<-table[grep(celltype,names(table))]
  table<-as.data.frame(table)
  table$Frac<-round(table$Freq/sum(table$Freq)*100,2)
  table$Week<-GW
  ggplot(table,aes(y=Frac,x=Week,fill=Var1))+
    geom_bar(stat="identity",colour = "white", position = "stack")+
    scale_fill_manual(values = colors)+
    xlab("Samples")+ylab("Relative percent")+
    theme(legend.position="right")+
    theme(panel.background=element_blank(),axis.text.x = element_text(size = 11, angle=0 , hjust = 1,vjust=1, colour = "black"),axis.text.y=element_text(size=11,colour="black"))
  ggsave(filename,width = 6,height = 4)
}
frac_donutchart(celltype = 'SST',filename = paste0('~/results/',GW,'SST_subcelltype_frac_bar.pdf'),
                colors = c('#F6EA8C','#FFB140','#CE6D39','#EC7357'))
frac_donutchart(celltype = 'PVALB',filename =paste0('~/results/',GW,'PVALB_subcelltype_frac_bar.pdf') ,
                colors = c('#A5DFF9','#1EC0FF','#6AAFE6'))
frac_donutchart(celltype = 'VIP',filename = paste0('~/results/',GW,'VIP_subcelltype_frac_bar.pdf'),
                colors = c('#8FBC94','#C5E99B'))                 
frac_donutchart(celltype = 'ID2',filename =paste0('~/results/',GW,'ID2_subcelltype_frac_bar.pdf'),
                colors = c('#A593E0','#9055A2'))                    
frac_donutchart(celltype = 'NDNF',filename = paste0('~/results/',GW,'NDNF_subcelltype_frac_bar.pdf'),
                colors = c('#D499B9','#D68189')) 
