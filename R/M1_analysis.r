library(Seurat)
library(data.table)
library(dplyr)

##(1)#######################read data and clustering analysis ##########################
##read M1 martix and metadata
cm<-read.csv('~/data/M1_expr.csv')
md<-read.csv('~/metadata_M1.csv')

##select GABAergic cells
use_cells<-row.names(subset(md,md$class_label=='GABAergic'))
cm<-cm[,use_cells]

##SCT Transform
M1_project<-CreateSeuratObject(cm,meta.data = subset(md,md$class_label=='GABAergic'))
M1_project <- SplitObject(M1_project, split.by = "external_donor_name_oredr")
for (i in 1:length(M1_project)) {
  project.integrated[[i]] <- SCTransform(project.integrated[[i]])
}
project.features <- SelectIntegrationFeatures(object.list = M1_project, nfeatures = 3000)
M1_project <- PrepSCTIntegration(object.list = M1_project, anchor.features = project.features)
pancreas.anchors <- FindIntegrationAnchors(object.list =M1_project, normalization.method = "SCT", 
                                           anchor.features = project.features,k.anchor = 5)
M1_project<- IntegrateData(anchorset = pancreas.anchors, normalization.method = "SCT")

##RunPCA and Umap
M1_project <- RunPCA(M1_project, verbose = FALSE)
M1_project <- RunUMAP(M1_project, dims = 1:30)
M1_project <- FindNeighbors(M1_project,dims = 1:30)
M1_project <- FindClusters(M1_project,resolution = 0.15)

##calculate marker genes
diffs<- FindAllMarkers(object = LGE_project, logfc.threshold = 0.25, only.pos = TRUE, test.use="bimod", min.pct = 0.1,assay="RNA")
diffs2<-subset(diffs[grep("^RP[L|S]",diffs$gene, ignore.case = FALSE,invert=TRUE),],subset=avg_logFC>0.5 & pct.1 > 0.6& pct.2 < 0.4 & (pct.1 > pct.2))
diffs2<-subset(diffs2,subset=!duplicated(gene))

##
new_idents<-c('PVALB-1 Basket','SST-1 Martinotti','VIP-1 Bipolar','SST-2 Martinotti','VIP-2 Multipolar (CCK)',
              'ID2-1','ID2-2','NDNF-1','NDNF-2','PVALB-2 Chandelier','PVALB-3 Chandelier','SST-4 NPY')
names(new_idents)<-levels(M1_project)
M1_project<-RenameIdents(M1_project,new_idents)
levels(M1_project)<-c("PVALB-1 Basket","PVALB-2 Chandelier","PVALB-3 Chandelier","SST-1 Martinotti","SST-2 Martinotti",
                      "SST-3 Martinotti","SST-4 NPY","VIP-1 Bipolar","VIP-2 Multipolar (CCK)","ID2-1","ID2-2","NDNF-1","NDNF-2")                

##(2)###############################plot################################
##fig7a
M1_project$celltype<-M1_project@active.ident
large_celltype<-c('PVALB','PVALB','PVALB','SST','SST','SST','SST',
                  'VIP','VIP',"ID2",'ID2','NDNF','NDNF')
names(large_celltype)<-levels(M1_project)
M1_project<-RenameIdents(M1_project,large_celltype)
DimPlot(M1_project,cols = c('#FFD800','#46B5D1','#8CBA51','#A593E0','#D68189'))
ggsave('~/results/M1_tsne_lagrgelabel.pdf')
##fig7b
FeaturePlot(M1_project,features = c('SST','PVALB','VIP','ID2','NDNF'))
ggsave('~/results/featureplot_M1.pdf')

##figS13a
DimPlot(M1_project,
        cols = c('#A5DFF9','#1EC0FF','#6AAFE6','#F6EA8C','#FFB140','#CE6D39','#EC7357',
                          '#8FBC94','#C5E99B','#A593E0','#9055A2','#D499B9','#D68189'))
ggsave('~/results/M1_tsne_sublabel.pdf')
