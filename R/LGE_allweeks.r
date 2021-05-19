library(Seurat)
library(harmony)
library(data.table)
library(dplyr)
library(monocle)
##(1) ############################read data and clustering####################################
project<-readRDS('~/results/project.rds')
LGE<-subset(project,idents=c('LGE1','LGE2','LGE3'))
LGE<-NormalizeData(object = LGE, verbose = TRUE,scale.factor = 10000)
LGE <- FindVariableFeatures(object = LGE,selection.method = "vst")
LGE <- ScaleData(object = LGE , verbose = TRUE)
LGE <- RunPCA(LGE,npcs = 30)
LGE$batch<-LGE$Week
LGE$batch<-gsub('GW09','V3',LGE$batch)
LGE$batch<-gsub('GW10','V2',LGE$batch)
LGE$batch<-gsub('GW11','V3',LGE$batch)
LGE$batch<-gsub('GW12','V3',LGE$batch)
LGE <- RunHarmony(LGE,group.by.vars = c('Week','batch'),dims.use = 1:20)
LGE <- FindNeighbors(LGE,dims = 1:20,reduction = 'harmony')
LGE <- FindClusters(LGE,resolution = 0.6)
LGE <- RunTSNE(LGE,dims = 1:20,reduction = 'harmony')

##calculate cluster markers
diffs<- FindAllMarkers(object = LGE, logfc.threshold = 0.25, only.pos = TRUE, test.use="bimod", min.pct = 0.1,assay="RNA")
diffs2<-subset(diffs[grep("^RP[L|S]",diffs$gene, ignore.case = FALSE,invert=TRUE),],subset=avg_logFC>0.5 & pct.1 > 0.6& pct.2 < 0.4 & (pct.1 > pct.2))
diffs2<-subset(diffs2,subset=!duplicated(gene))

##annoate cluster
newidents<-c('ISL1,EBF1 2','ISL1,EBF1 3','PAX6,ETV1','SIX3,SOX2 2',
             'SIX3,SOX2 3','SIX3,SOX2 4','SIX3,SOX2 1','SIX3,SOX2 5',
             'ISL1,EBF1 4','ISL1,EBF1 1','TOP2A,ASPM')
names(newidents)<-levels(LGE)
LGE<-RenameIdents(LGE,newidents)
levels(LGE)<-c('TOP2A,ASPM',"PAX6,ETV1",'SIX3,SOX2 1','SIX3,SOX2 2','SIX3,SOX2 3','SIX3,SOX2 4','SIX3,SOX2 5','ISL1,EBF1 1',
               'ISL1,EBF1 2','ISL1,EBF1 3','ISL1,EBF1 4')

##(2)########################monocle analysis############################
##create monocle2 project and analysis
cm<-GetAssayData(LGE,slot = 'counts')
LGE$celltype<-LGE@active.ident
md<-LGE@meta.data
gene_annotation <- data.frame(gene_id=rownames(cm), gene_short_name=rownames(cm))
rownames(gene_annotation) <- as.vector(gene_annotation[, 1])
fd <- new("AnnotatedDataFrame", data = gene_annotation)
pd <- new("AnnotatedDataFrame", data = md)
monocl_data <- newCellDataSet(cm, phenoData = pd, 
                              featureData = fd, expressionFamily=negbinomial())
monocl_data <- estimateSizeFactors(monocl_data)
monocl_data <- estimateDispersions(monocl_data)
##using the clusters markers genes to ordering cell and dimension reduce
diffs<- FindAllMarkers(object = LGE,logfc.threshold = 0.25, only.pos = TRUE, test.use="bimod", min.pct = 0.1,assay="RNA")
diffs2<-subset(diffs[grep("^RP[L|S]",diffs$gene, ignore.case = FALSE,invert=TRUE),],subset=avg_logFC>0.75 & (pct.1 > pct.2))
diffs2<-subset(diffs2,subset=!duplicated(gene))
monocl_data <- setOrderingFilter(monocl_data,diffs2$gene)
monocl_data <- reduceDimension(monocl_data, max_components =2
                               , reduction_method = 'DDRTree', 
                               norm_method = c("log"),residualModelFormulaStr="~Week")
monocl_data <- orderCells(monocl_data)

## (3)######################plot########################
##figS10a
DimPlot(LGE,reduction = 'tsne',pt.size = 0.5,label = F)+
  scale_color_manual(values = c('#6a197d','#fa7d09','#00a8cc','#87dfd6',
                                '#086972','#abf0e9','#709fb0',
                                '#9f5f80','#e36387','#ffcac2','#fc9d9d'))
ggsave("~/results/LGE_tsne_allweeks.pdf")

##figS10c
diffs2 %>% group_by(cluster) %>% top_n(-3, p_val_adj) -> top3
top3 %>% group_by(cluster) %>% top_n(3, avg_logFC) -> top3
top3<-top3[order(top3$cluster),]
levels(LGE)<-rev(c('TOP2A,ASPM',"PAX6,ETV1",'SIX3,SOX2 1','SIX3,SOX2 2','SIX3,SOX2 3','SIX3,SOX2 4','SIX3,SOX2 5','ISL1,EBF1 1',
                   'ISL1,EBF1 2','ISL1,EBF1 3','ISL1,EBF1 4'))
Seurat::DotPlot(LGE,features = top3$gene)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  scale_color_gradient2(high = '#fe8a71',mid = '#fcf9ea',low = '#3da4ab',midpoint = 0)
ggsave('~/results/LGE_dotplot_allweeks.pdf')  

##figS10b
monocl_data$celltype<-LGE@active.ident
plot_cell_trajectory(monocl_data,color_by = 'celltype',cell_size =1,show_branch_points = T)+
  scale_color_manual(values = c('#6a197d','#fa7d09','#00a8cc','#87dfd6',
                                    '#086972','#abf0e9','#709fb0',
                                    '#9f5f80','#e36387','#ffcac2','#fc9d9d'))
ggsave('~/results/LGE_monocle_allweeks.pdf')
