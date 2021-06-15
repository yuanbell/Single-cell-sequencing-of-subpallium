##
library(Seurat)
library(ggplot2)
library(data.table)
library(RColorBrewer)
library(dplyr)
library(cowplot)
library(harmony)
options(warn=-1)

# (1) ######################read data and preprocess############################

##read 10X digital expression data 
sample_name<-c('./9W','./10W','./11W','./12W')
project.data<-list()
for (i in 1:length(sample_name)) {
  setwd('~/data/')
  project.data[[i]]<-Read10X(data.dir=sample_name[i])
}

###add cell barcode index (GW09 -1,GW10 -2,GW11 -3,GW12 -4)
for (i in 1:length(project.data)) {
  colnames(project.data[[i]])<-gsub(pattern = 1,replacement = i,colnames(project.data[[i]]))
}


##filter doublet cells
doublets_filename<-list.files('~/doublets/',pattern = '.txt')
doublets_files<-list()
for (i in 1:length(doublets_filename)) {
  doublets_files[[i]]<-read.table(paste0('~/doublets/',doublets_filename[i]))
}

for (i in 1:length(project.data)) {
  project.data[[i]]<-project.data[[i]][,which(doublets_files[[i]]$V1==0)]
}

###merge cell matrix of 4 samples
cm<-cbind(project.data[[1]],project.data[[2]],project.data[[3]],project.data[[4]])

###Create metadata
md<-rep(c('GW09','GW10','GW11','GW12'),times=c(length(colnames(project.data[[1]])),
                                                     length(colnames(project.data[[2]])),
                                                     length(colnames(project.data[[3]])),
                                                     length(colnames(project.data[[4]]))))
md<-as.data.frame(md)
row.names(md)<-c(colnames(project.data[[1]]),
                       colnames(project.data[[2]]),
                       colnames(project.data[[3]]),
                       colnames(project.data[[4]]))
colnames(md)<-'Week'



# (2) ######################## Seurat Workflow to cluster cells ############################


##Create Seurat object
project <- CreateSeuratObject(counts = cm, project = "project",meta.data = md, min.cells =60, min.features = 0)

##cell stat
project[["percent.mt"]] <- PercentageFeatureSet(project, pattern = "^MT-") 
project[["percent.ribo"]] <- PercentageFeatureSet(project, pattern = "^RP[L|S]") 
red.genes <- c("HBA1","HBA2","HBB",'HBD','HBE1','HBG1','HBG2','HBM','HBQ1','HBZ') 
project[["percent.redcell"]] <- PercentageFeatureSet(project,
                                                     features=red.genes[!is.na(match(red.genes,rownames(project)))]) 

##set filter threshold
nFeature_lowlimit=1000
nFeature_highlimit=10000
nCount_RNA_highlimit=30000
mito_cutoff=10
redcell_cutoff=10
ribo_cutoff=40
pcs_number=30
resolution=0.6
top_num=10
##filter cells
project <- subset(project, subset = nFeature_RNA > nFeature_lowlimit  &  
                    nCount_RNA <  nCount_RNA_highlimit& 
                    percent.mt < mito_cutoff & 
                    percent.redcell < redcell_cutoff & 
                    percent.ribo < ribo_cutoff)
##remove sex genes
sex.genes<-c('DDX3Y','EIF2S3Y','UTY','KDM5D','XIST','TSIX')
project<-subset(project,features = setdiff(row.names(project.data),c(red.genes,mt.genes,sex.genes)))

##find highly variable features
project <- FindVariableFeatures(object = project,selection.method = "vst")

##normalize
project <- NormalizeData(object = project, verbose = TRUE,scale.factor = 10000)

##cell cycle score
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
project <- CellCycleScoring(project, s.features = s.genes, 
                            g2m.features = g2m.genes, set.ident = TRUE)
project$CC.Difference <- project$S.Score - project$G2M.Score               

##scale and regress cell cycle effect
project <- ScaleData(object = project, verbose = TRUE,
                     vars.to.regress = 'CC.Difference',
                     features = project@assays$RNA@var.features)

##PCA dimension reduce
project <- RunPCA(object = project, npcs = pcs_number, 
                  verbose = TRUE,
                  features=VariableFeatures(project))

##harmony to remove batch effect
project<-RunHarmony(project,"Week", plot_convergence = TRUE,dims.use = 1:pcs_number)

##cluster
project <- FindNeighbors(project, reduction = "harmony", dims = 1:pcs_number) 
project <- FindClusters(resolution = resolution, project)

##T-SNE dimension reduce
project <- RunTSNE(object = project, reduction = "harmony",dims=1:pcs_number)

##Calculate cluster markers
diffs <- FindAllMarkers(object = project, logfc.threshold = 0.25, only.pos = TRUE, test.use="bimod", min.pct = 0.1,assay="RNA")
diffs2<-subset(diffs[grep("^RP[L|S]",diffs$gene, ignore.case = FALSE,invert=TRUE),],subset=avg_logFC>0.25 & pct.1 > 0.25 & (pct.1 > pct.2))
diffs2<-subset(diffs2,subset=!duplicated(gene))

##Annotate celltype and Rename clusters
CellTypes <- c("MGE1","P1","N1","LGE1","CGE1","P2","LGE2","N2","MGE2","P6",
               "N3","P3","IPCs","CGE2","P4","P5","N4","OG","EC1","MGE3","MG","EC2","CGE3")
names(CellTypes) <- levels(project)
project <- RenameIdents(project, CellTypes)
levels(project)<-c('P1','P2','P3','P4','P5','P6','IPCs2','pre-IN',
                   'MGE1','MGE2','CGE1','LGE1','LGE2','LGE3','OPC',
                   'MG','EC1','EC2','N1','N2','N3')



# (3) ######################## Plot ############################
##Fig1b
cols=c("#efefb6","#A1D99B","#41AB5D","#238B45","#006D2C","#4292C6","#2171B5","#d8c962","#ff8ba7",
       "#ffc6c7","#ffbd69","#ad62aa","#c295d8","#faafff","#6ba8a9","#844685","#ef962d","#9c5518",
       "#DEEBF7","#9ECAE1","#6BAED6")
pdf(file="/results/tsne_cluster_allsample.pdf",width=8,height=6)
p <- DimPlot(object = project,reduction = "tsne",
             pt.size=0.1, label.size = 4, label = FALSE,cols=cols,
             order=rev(levels(project)))
print(p)
dev.off()

##figS1a
pdf(file="tsne_cluster_eachsample.pdf",width=14,height=10)
p<-DimPlot(project, reduction = "tsne", pt.size = 0.1,
        split.by = 'Week',ncol = 2,cols = cols)
print(p)
dev.off()

##figS1c
pdf(file="tsne_cluster_GW.pdf",width=14,height=10)
p<-DimPlot(project, reduction = "tsne", pt.size = 0.1,
           group.by = 'Week',ncol = 2)
print(p)
dev.off()

##fig1c
gene2plot<-c('TOP2A','HES5','GAD2','PDGFRA','SPP1','IGFBP7')
pdf(file = 'feature_plot_majorcluster.pdf')
FeaturePlot(project,features = gene2plot,pt.size = 0.1,order = T)
dev.off()

##figS2
gene2plot2<-c('EOMES','NEUROD2','HES1','NUSAP1','TTYH1','TYMS','GAD2','ASCL1','DLX2','LHX6','LHX8',
              'MAF','MAFB','ERBB4','NKX2-1','SST','NR2F1','NR2F2','EBF1','ISL1','SIX3','ZNF503')
pdf(file = 'feature_plot_cluster.pdf')
FeaturePlot(project,features = gene2plot2,pt.size = 0.1,order = T)
dev.off()

##figS1d
project<-subset(project,idents=c("P1","P2","P3","P4","P5",'pre-IN','MGE1',"MGE2",'CGE1',"LGE1","LGE2",'LGE3',
                                 'OPC','MG','EC1','EC2'))
genelist<-as.character(read.csv('~/heatmap_features.csv',header = F)$V1)
project_avg<-Seurat::AverageExpression(project,return.seurat = T)
cm<-GetAssayData(project_avg,slot = 'data')[row.names(project) %in% genelist,]
cm<-cm[genelist,]
cm<-t(apply(cm,1,scale))
colnames(cm)<-colnames(project_avg)
cm[cm>1.5]=1.5
cm[cm<(-1.5)]=-1.5
pheatmap::pheatmap(cm,cluster_cols = F,cluster_rows = F,color=colorRampPalette(c("#442e63","#088682","#f4f262"))(100),
                   fontsize_row =7,border_color = NA,filename = '~/pheatmap.pdf',width =5.5,height = 10)

##figs1b
input<-data.frame(Barcode=colnames(project),Sample=project$Week,Cluster=project@active.ident)
data<-melt(table(input$Week,input$Cluster))
colnames(data)<-c("Sample","Cluster","Count")
data$Cluster<-ordered(data$Cluster)
data<-ddply(data,"Sample",transform, Percent = Count / sum(Count))
pdf(file="cluster_dis.pdf",width=5,height=6)
ggplot(data, aes(x=Sample, y=Percent,fill=Cluster)) + geom_bar(stat="identity",colour = "white", position = "stack")+xlab("Samples")+ylab("Relative percent")+theme(legend.position="right")+theme(panel.background=element_blank(),axis.text.x = element_text(size = 11, angle=0 , hjust = 1,vjust=1, colour = "black"),axis.text.y=element_text(size=11,colour="black"))+scale_fill_manual(values=col)
dev.off()


##save data
saveRDS(project,'~/results/project_allsample.rds')
