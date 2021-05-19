library(Seurat)
library(harmony)
library(monocle)
library(GENIE3)
library(clusterProfiler)
library(org.Hs.eg.db)
##(1)#####################read data and clustering######################
CGE<-subset(project,idents='CGE1')
CGE<-subset(CGE,cells = row.names(subset(CGE@meta.data,
                                         CGE@meta.data$Week=='GW09'|MGE@meta.data$Week=='GW11')))
CGE<-NormalizeData(object = CGE, verbose = TRUE,scale.factor = 10000)
CGE <- FindVariableFeatures(object = CGE,selection.method = "vst")
CGE <- ScaleData(object = CGE , verbose = TRUE)
CGE <- RunPCA(CGE,npcs = 20)
CGE <- RunHarmony(CGE,group.by.vars = 'Week',max.iter.harmony = 20)
CGE <- FindNeighbors(CGE,dims = 1:20,reduction = 'harmony')
CGE <- FindClusters(CGE,resolution = 0.3)
CGE <- RunTSNE(CGE,dims = 1:20,reduction = 'harmony')

##calculate marker genes
diffs<- FindAllMarkers(object = CGE, logfc.threshold = 0.25, only.pos = TRUE, test.use="bimod", min.pct = 0.1,assay="RNA")
diffs2<-subset(diffs[grep("^RP[L|S]",diffs$gene, ignore.case = FALSE,invert=TRUE),],subset=avg_logFC>0.25 & pct.1 > 0.4& pct.2 < 0.4 & (pct.1 > pct.2))
diffs2<-subset(diffs2,subset=!duplicated(gene))

##annoated cluster
cluster_name<-c('CALb2+,SST-','NFIA+','ANK3+,SST-','ANK3+,SST+','CGE_pro')
names(cluster_name)<-levels(CGE)
CGE<-RenameIdents(CGE,cluster_name)

##(2)########################monocle analysis############################
CGE$celltype<-CGE@active.ident
cm<-GetAssayData(CGE,slot = 'counts')
md<-CGE@meta.data
gene_annotation <- data.frame(gene_id=rownames(cm), gene_short_name=rownames(cm))
rownames(gene_annotation) <- as.vector(gene_annotation[, 1])
fd <- new("AnnotatedDataFrame", data = gene_annotation)
pd <- new("AnnotatedDataFrame", data = md)
monocl_data <- newCellDataSet(cm, phenoData = pd, 
                              featureData = fd, expressionFamily=negbinomial())
monocl_data <- estimateSizeFactors(monocl_data)
monocl_data <- estimateDispersions(monocl_data)
monocl_data <- setOrderingFilter(monocl_data, diffs2$gene)
monocl_data <- reduceDimension(monocl_data, max_components = 2, reduction_method = 'DDRTree', 
                               norm_method = c("log"),residualModelFormulaStr=(~Week)
)
monocl_data <- orderCells(monocl_data,reverse = T)

##calculate pseudotime associate TF
diff_test_res2 <- differentialGeneTest(monocl_data,cores = 4,
                                       fullModelFormulaStr = "~sm.ns(Pseudotime,df=3)")
diffs_genes_pseudotimes<-subset(diff_test_res2,diff_test_res2$qval<1e-4)
diffs_genes_pseudotimes<-diffs_genes_pseudotimes[order(diffs_genes_pseudotimes$qval),]
Homo_TF <- read.table('~/data/Homo_sapiens_TF.txt',header = T,sep = '')

##(3)#########################use GENEIE3 to analysis TF coexpression network#################################
regulators<-c('TOX3','MEIS2','ARX','DLX5','ST18','HES6','ZNF608','CSRNP3','MYT1L','NR2F1')
weightMat <- GENIE3(as.matrix(cm), regulators=regulators,nCores = 4)
linkList <- getLinkList(weightMat)
linkList <- getLinkList(weightMat, reportMax=1000)
##export table and visable by Cytoscope
write.table(linkList,'CGE_linlist.txt',row.names = F,col.names = T,sep = ',')

##GO analysis

GO_TOX3<-enrichGO(gene=subset(linkList,linkList$regulatoryGene=='TOX3')[,2],OrgDb =org.Hs.eg.db,keyType = 'SYMBOL',
                  ont = 'BP',pvalueCutoff = 0.05,qvalueCutoff = 0.05)
GO_MEIS2<-enrichGO(gene=subset(linkList,linkList$regulatoryGene=='MEIS2')[,2],OrgDb =org.Hs.eg.db,keyType = 'SYMBOL',
                   ont = 'BP',pvalueCutoff = 0.05,qvalueCutoff = 0.05)
GO_ARX<-enrichGO(gene=subset(linkList,linkList$regulatoryGene=='ARX')[,2],OrgDb =org.Hs.eg.db,keyType = 'SYMBOL',
                 ont = 'BP',pvalueCutoff = 0.05,qvalueCutoff = 0.05)
GO_CSRNP3<-enrichGO(gene=subset(linkList,linkList$regulatoryGene=='CSRNP3')[,2],OrgDb =org.Hs.eg.db,keyType = 'SYMBOL',
                    ont = 'BP',pvalueCutoff = 0.05,qvalueCutoff = 0.05)
GO_DLX5<-enrichGO(gene=subset(linkList,linkList$regulatoryGene=='DLX5')[,2],OrgDb =org.Hs.eg.db,keyType = 'SYMBOL',
                  ont = 'BP',pvalueCutoff = 0.05,qvalueCutoff = 0.05)
GO_HES6<-enrichGO(gene=subset(linkList,linkList$regulatoryGene=='HES6')[,2],OrgDb =org.Hs.eg.db,keyType = 'SYMBOL',
                  ont = 'BP',pvalueCutoff = 0.05,qvalueCutoff = 0.05)
GO_MYT1L<-enrichGO(gene=subset(linkList,linkList$regulatoryGene=='MYT1L')[,2],OrgDb =org.Hs.eg.db,keyType = 'SYMBOL',
                   ont = 'BP',pvalueCutoff = 0.05,qvalueCutoff = 0.05)
GO_NR2F1<-enrichGO(gene = subset(linkList,linkList$regulatoryGene=='NR2F1')[,2],OrgDb =org.Hs.eg.db,keyType = 'SYMBOL',
                   ont = 'BP',pvalueCutoff = 0.05,qvalueCutoff = 0.05)
GO_ST18<-enrichGO(gene = subset(linkList,linkList$regulatoryGene=='ST18')[,2],OrgDb =org.Hs.eg.db,keyType = 'SYMBOL',
                  ont = 'BP',pvalueCutoff = 0.05,qvalueCutoff = 0.05)
GO_ZNF608<-enrichGO(gene = subset(linkList,linkList$regulatoryGene=='ZNF608')[,2],OrgDb =org.Hs.eg.db,keyType = 'SYMBOL',
                    ont = 'BP',pvalueCutoff = 0.05,qvalueCutoff = 0.05)

###GO analysis
GO_TOX3<-enrichGO(gene=subset(linkList,linkList$regulatoryGene=='TOX3')[,2],OrgDb =org.Hs.eg.db,keyType = 'SYMBOL',
                  ont = 'BP',pvalueCutoff = 0.05,qvalueCutoff = 0.05)
GO_MEIS2<-enrichGO(gene=subset(linkList,linkList$regulatoryGene=='MEIS2')[,2],OrgDb =org.Hs.eg.db,keyType = 'SYMBOL',
                   ont = 'BP',pvalueCutoff = 0.05,qvalueCutoff = 0.05)
GO_ARX<-enrichGO(gene=subset(linkList,linkList$regulatoryGene=='ARX')[,2],OrgDb =org.Hs.eg.db,keyType = 'SYMBOL',
                 ont = 'BP',pvalueCutoff = 0.05,qvalueCutoff = 0.05)
GO_CSRNP3<-enrichGO(gene=subset(linkList,linkList$regulatoryGene=='CSRNP3')[,2],OrgDb =org.Hs.eg.db,keyType = 'SYMBOL',
                    ont = 'BP',pvalueCutoff = 0.05,qvalueCutoff = 0.05)
GO_DLX5<-enrichGO(gene=subset(linkList,linkList$regulatoryGene=='DLX5')[,2],OrgDb =org.Hs.eg.db,keyType = 'SYMBOL',
                  ont = 'BP',pvalueCutoff = 0.05,qvalueCutoff = 0.05)
GO_HES6<-enrichGO(gene=subset(linkList,linkList$regulatoryGene=='HES6')[,2],OrgDb =org.Hs.eg.db,keyType = 'SYMBOL',
                  ont = 'BP',pvalueCutoff = 0.05,qvalueCutoff = 0.05)
GO_MYT1L<-enrichGO(gene=subset(linkList,linkList$regulatoryGene=='MYT1L')[,2],OrgDb =org.Hs.eg.db,keyType = 'SYMBOL',
                   ont = 'BP',pvalueCutoff = 0.05,qvalueCutoff = 0.05)
GO_NR2F1<-enrichGO(gene = subset(linkList,linkList$regulatoryGene=='NR2F1')[,2],OrgDb =org.Hs.eg.db,keyType = 'SYMBOL',
                   ont = 'BP',pvalueCutoff = 0.05,qvalueCutoff = 0.05)
GO_ST18<-enrichGO(gene = subset(linkList,linkList$regulatoryGene=='ST18')[,2],OrgDb =org.Hs.eg.db,keyType = 'SYMBOL',
                  ont = 'BP',pvalueCutoff = 0.05,qvalueCutoff = 0.05)
GO_ZNF608<-enrichGO(gene = subset(linkList,linkList$regulatoryGene=='ZNF608')[,2],OrgDb =org.Hs.eg.db,keyType = 'SYMBOL',
                    ont = 'BP',pvalueCutoff = 0.05,qvalueCutoff = 0.05)

##(4)#############################plot###############################
##fig7a
DimPlot(CGE,reduction = 'tsne',label = F,cols = col)
ggsave('~/results/CGE_celltype_tsne_GW9-GW11.pdf',width = 6.6,height = 4)

##fig7b
vinplot_genes<-read.csv('~/data/CGE_vinplot_genes.csv',header = F)
vinplot_cm<-GetAssayData(CGE,slot = 'data')[vinplot_genes$V1,]
vinplot_data<-as.data.frame(t(vinplot_cm))
vinplot_data$celltype<-CGE@active.ident
vinplot_data<-melt(vinplot_data,id.vars = 'celltype')
p1<-ggplot(vinplot_data,aes(x=celltype,y=value,fill=celltype))+
  geom_violin(scale='width')+facet_wrap(~variable,scales = 'free_y',ncol=1,strip.position = 'right',as.table=F,shrink = F)+
  theme_minimal_vgrid()+
  theme(axis.text.y=element_blank(),strip.text.y = element_text(angle = 0,size=10),axis.ticks.y=element_blank())+
  scale_fill_manual(values = col)+
  guides(fill=FALSE)
p1
ggsave('~/results/vinlin_CGE_markers_GW9-GW11.pdf',height = 11,width = 5)

##fig7c
monocl_data$celltype<-CGE@active.ident
plot_cell_trajectory(monocl_data,color_by = 'celltype')+scale_color_manual(values = col)
pdf('cell_complex_trajectory_CGE_GW9-GW11.pdf',width = 6,height = 6)
plot_complex_cell_trajectory(monocl_data,color_by = 'celltype',root_states = 1,cell_size =1)
dev.off()

##fig7d-f
regulators<-c('TOX3','MEIS2','ARX','DLX5','ST18','HES6','ZNF608','CSRNP3','MYT1L','NR2F1')
setwd('~/results/')
for (i in 1:length(regulators)) {
  plot_genes_in_pseudotime(monocl_data[plot_genes_list[i],],color_by ='Pseudotime')+
    scale_colour_continuous(low='#527318',high='#ffe75e')
  ggsave(paste0(plot_genes_list[i],'_CGE_pesudotime.pdf'),height = 3, width = 6)
}

##figS12b
GO_TF<-read.csv('~/data/GO_heatmap.csv',header = T,row.names = 1)
GO_TF_heatmap=-log(GO_TF)
pheatmap::pheatmap(GO_TF_heatmap,cluster_rows = T,cluster_cols = T,color =colorRampPalette(c('#fcf8e8','#d63447'))(100),
                   border_color = '#ededed')
