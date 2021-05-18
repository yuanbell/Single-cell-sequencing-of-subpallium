library(Seurat)
library(monocle)

##(1) ############################read data and clustering####################################
project<-readRDS('~/results/project.rds')
LGE_project<-subset(project,idents=c('LGE1','LGE2','LGE3'))
LGE_project<-subset(LGE_project,cells = row.names(subset(LGE_project@meta.data,
                                                         LGE_projectE@meta.data$Week=='GW09'|
                                                           LGE_project@meta.data$Week=='GW11')))
LGE_project<-NormalizeData(object = LGE_project, verbose = TRUE,scale.factor = 10000)
LGE_project <- FindVariableFeatures(object = LGE_project,selection.method = "vst")
LGE_project <- ScaleData(object = LGE_project , verbose = TRUE)
LGE_project <- RunPCA(LGE_project,npcs = 20)
LGE_project <- RunHarmony(LGE_project,group.by.vars = 'Week',dims.use = 1:15)
LGE_project <- FindNeighbors(LGE_project,dims = 1:15,reduction = 'harmony')
LGE_project <- FindClusters(LGE_project,resolution = 0.68)
LGE_project <- RunTSNE(LGE_project,dims = 1:15,reduction = 'harmony')

##calculate marker genes
diffs<- FindAllMarkers(object = LGE_project, logfc.threshold = 0.25, only.pos = TRUE, test.use="bimod", min.pct = 0.1,assay="RNA")
diffs2<-subset(diffs[grep("^RP[L|S]",diffs$gene, ignore.case = FALSE,invert=TRUE),],subset=avg_logFC>0.5 & pct.1 > 0.6& pct.2 < 0.4 & (pct.1 > pct.2))
diffs2<-subset(diffs2,subset=!duplicated(gene))

##annoated cluster
new_idents<-c('SIX3+,SOX2+ 1','ISL1+,EBF1+ 3','ISL1+,EBF1+ 2','SIX3+,SOX2+ 2','ISL1+,EBF1+ 1','SIX3+,SOX2+ 3',
              'ISL1+,EBF1+ 4','SIX3+,SOX2+ 4','SIX3+,SOX2+ 5','PAX6+,ETV1+','ASPM+,TOP2A+')
names(new_idents)<-levels(LGE)
LGE<-RenameIdents(LGE,new_idents)
levels(LGE)<-rev(c('ASPM+,TOP2A+',"PAX6+,ETV1+",'SIX3+,SOX2+ 1','SIX3+,SOX2+ 2','SIX3+,SOX2+ 3','SIX3+,SOX2+ 4','SIX3+,SOX2+ 5','ISL1+,EBF1+ 1',
                   'ISL1+,EBF1+ 2','ISL1+,EBF1+ 3','ISL1+,EBF1+ 4'))



## (2)#########################monocle2 analysis#############################

## create monocle2 project and analysis
cm<-GetAssayData(LGE,slot = 'counts')
md<-LGE@meta.data
gene_annotation <- data.frame(gene_id=rownames(cm), gene_short_name=rownames(cm))
rownames(gene_annotation) <- as.vector(gene_annotation[, 1])
fd <- new("AnnotatedDataFrame", data = gene_annotation)
pd <- new("AnnotatedDataFrame", data = md)
monocl_data <- newCellDataSet(cm, phenoData = pd, 
                              featureData = fd, expressionFamily=negbinomial())
monocl_data$celltype<-LGE@active.ident
monocl_data <- estimateSizeFactors(monocl_data)
monocl_data <- estimateDispersions(monocl_data)
##using the clusters markers genes to ordering cell and dimension reduce
monocl_data <- setOrderingFilter(monocl_data,diffs2$gene)
monocl_data <- reduceDimension(monocl_data, max_components =2, reduction_method = 'DDRTree', 
                               norm_method = c("log"),residualModelFormulaStr=(~Week)
)
monocl_data <- orderCells(monocl_data)

##BEAM analyssi of branch_point 2
BEAM_res<-BEAM(monocl_data, branch_point = 2, cores = 4)
BEAM_res<-subset(BEAM_res,BEAM_res$qval<0.001)
BEAN_res<-BEAM_res[order(BEAM_res$qval),]
BEAN_res$gene_short_name<-as.character(BEAN_res$gene_short_name)

##performe Lin1 and Lin2 differential analysis
LGE_new<-LGE
Idents(LGE_new)<-monocl_data$State
State_idents<-c('Lin1','Lin1','Lin2','Pre-Branch','Lin1')
names(State_idents)<-levels(LGE_new)
LGE_new<-RenameIdents(LGE_new,State_idents)
Lin1vsLIn2<-FindMarkers(LGE_new,ident.1 = 'Lin1',ident.2 = 'Lin2',assay = 'RNA',test.use = 'wilcox',logfc.threshold = 0,min.pct = 0)
Lin1vsLIn2$type<-ifelse(Lin1vsLIn2$p_val_adj<0.05&abs(Lin1vsLIn2$avg_logFC)>0.25,
                        ifelse(Lin1vsLIn2$avg_logFC>0,'Lin1','Lin2'),'NO')

##GSEA analysis
lin1_geneid<-bitr(row.names(subset(Lin1vsLIn2,Lin1vsLIn2$type=='Lin1')),fromType = 'SYMBOL',
                  toType = 'ENTREZID',OrgDb = org.Hs.eg.db)
Lin1vsLIn2<-Lin1vsLIn2[order(Lin1vsLIn2$avg_logFC),]
rank_genes<-Lin1vsLIn2$avg_logFC
names(rank_genes)<-row.names(Lin1vsLIn2)
rank_genes <- sort(rank_genes, decreasing = TRUE)
ego3 <- gseGO(geneList     = rank_genes,
              OrgDb        = org.Hs.eg.db,
              ont          = "BP",
              nPerm        = 1000,
              minGSSize    = 100,
              maxGSSize    = 500,
              pvalueCutoff = 0.1,
              pAdjustMethod = 'fdr',
              verbose      = FALSE,
              keyType="SYMBOL")


## (3)####################plot########################
##fig5a
DimPlot(LGE,reduction = 'tsne',pt.size = 0.5)+
  scale_color_manual(values = c('#6a197d','#fa7d09','#00a8cc','#87dfd6','#9f5f80','#086972','#abf0e9','#709fb0','#ffcac2','#fc9d9d','#e36387'))
ggsave("~/results/LGE_tsne_GW9-GW11.pdf")

##figS9b
diffs2 %>% group_by(cluster) %>% top_n(-3, p_val_adj) -> top3
top3 %>% group_by(cluster) %>% top_n(3, avg_logFC) -> top3
top3<-top3[order(top3$cluster),]
Seurat::DotPlot(LGE,features = dot_plot$X1)+theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave('~/results/dotplot_LGE_GW9-GW11.pdf')

##figS9a
features<-c('ASPM','TOP2A','ETV1','PAX6','SIX3','SOX2','EBF1','ISL1')
FeaturePlot(LGE,features = features,reduction = 'tsne',order = T)
ggsave('~/results/featureplot_LGE_GW9-GW11.pdf')

##fig5c
plot_cell_trajectory(monocl_data,color_by = 'celltype',cell_size =1,show_branch_points = T)+
  scale_color_manual(values = c('#6a197d','#fa7d09','#00a8cc','#87dfd6','#9f5f80',
                                           '#086972','#abf0e9','#709fb0','#ffcac2','#fc9d9d','#e36387'))
ggsave('~/results/cell_trajectory_LGE_GW9-GW11.pdf')

##figS9c
plot_complex_cell_trajectory(monocl_data,root_states = 4,color_by = 'celltype',cell_size = 0.5)+
  scale_color_manual(values = c('#6a197d','#fa7d09','#00a8cc','#87dfd6','#9f5f80',
                                '#086972','#abf0e9','#709fb0','#ffcac2','#fc9d9d','#e36387'))
ggsave('~/results/cell_trajectory_complex_LGE_GW9-GW11.pdf')


##fig5d
sig_genes<-c('PAX6','TAC1','EBF1','ISL1','IGFBP5')
for (i in 1:length(sig_genes)) {
  plot_genes_branched_pseudotime(LGE_monocl_data[sig_genes[i],],branch_point = 2,color_by = 'celltype',
                                method = 'loess',branch_labels=c('lin1','lin2'),cell_size = 0.5)+
    scale_color_manual(values = c('#ffcac2','#fc9d9d','#e36387','#6a197d','#ffd31d','#e71414','#fa7d09',
                                  '#00a8cc','#87dfd6','#9f5f80','#086972','#abf0e9','#709fb0'))+
    scale_fill_manual(values = c("#ffd31d","#e71414"))
  ggsave(paste0('~/results/',sig_genes[i],'_genes_branched_pseudotime.pdf'),height =10,width = 4)
}

##fig5h
Lin1vsLIn2$avg_logFC[Lin1vsLIn2$avg_logFC>1.5]=1.5
Lin1vsLIn2$avg_logFC[Lin1vsLIn2$avg_logFC<(-1.5)]=-1.5
Lin1vsLIn2$p_val_adj[Lin1vsLIn2$p_val_adj<0]=0
Lin1vsLIn2$p_val_adj[-log10(Lin1vsLIn2$p_val_adj)>200]=exp(x = -200)
ggplot(data = Lin1vsLin2,aes(y=-log10(p_val_adj),x=avg_logFC,color=type))+geom_point(size=1)+
  theme_classic()+scale_color_manual(values = c("#ffd31d","#e71414",'gray'))+
  xlim(c(-1.6,1.6))+
  geom_hline(aes(yintercept=-log10(0.05)),colour="gray", linetype="dashed",size=1.5)+
  geom_vline(aes(xintercept=-0.25),colour="gray", linetype="dashed",size=1.5)+
  geom_vline(aes(xintercept=0.25),colour="gray", linetype="dashed",size=1.5)
ggsave('~/results/LGE_volcano_GW9-GW11.pdf')  

##fig5i
gseaplot2(ego3, geneSetID = c(5,50), pvalue_table = TRUE,
          color = c("#fb7813","#17706e"),subplots = 1:2)

gseaplot2(ego3, geneSetID = c(7,56), pvalue_table = TRUE,
          color = c("#fb7813","#17706e"),subplots = 1:2)
