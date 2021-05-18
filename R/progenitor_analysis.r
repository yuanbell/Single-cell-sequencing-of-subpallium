library(Seurat)
library(ggplot2)
library(data.table)

# (1)#################read data and seurat workflow###################
project <- readRDS('~/results/project.rds')


##select progenitor clusters
project <- subset(project,idents=c('P1','P2','P3','P4','P5','P6'))

##1)split data by GW and seurat workflow by each samples to tsne dimensition reduce
##2)calculate VZ and SVZ markers genes
project<-SplitObject(project,split.by = 'Week')
diffs2<-list()
GW<-c('GW09','GW10','GW11','GW12')
for (i in 1:length(project)) {
  project[[i]] <- NormalizeData(object = project[[i]],verbose = TRUE,scale.factor = 10000)
  project[[i]] <- FindVariableFeatures(project[[i]],selection.method = "vst",
                                       nfeatures = 1000)
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  project[[i]]<- CellCycleScoring(project[[i]], s.features = s.genes, 
                             g2m.features = g2m.genes, set.ident = TRUE)
  project[[i]] <- ScaleData(object = project[[i]], verbose = TRUE,
                            vars.to.regress = c("S.Score", "G2M.Score"))
  project[[i]] <- RunPCA(object = project[[i]], npcs = 10,verbose = T)
  project[[i]] <- FindNeighbors(project[[i]], dims = 1:10,reduction = 'pca')
  project[[i]] <- RunTSNE(object = project[[i]],dims=1:10,reduction ='pca' )
  project[[i]]@active.ident <- actor(as.character(project[[i]]$celltype),
                                     levels =c( 'P3','P2','P5','P4','P1','pre-IN'))
  newidents<-c('VZ','VZ','SVZ','SVZ','SVZ','SVZ')
  names(newidents)<-levels(project[[i]])
  project[[i]]<-RenameIdents(project[[i]],newidents)
  diffs <- FindAllMarkers(object = project[[i]], 
                              logfc.threshold = 0.25, 
                              only.pos = TRUE, test.use="bimod",
                              min.pct = 0.1,assay="RNA")
  diffs2[[i]]<-subset(diffs[grep("^RP[L|S]",diffs$gene, ignore.case = FALSE,invert=TRUE),],
                      subset=avg_logFC>0.25 & pct.1 > 0.5  & (pct.1 > pct.2) & pct.2 < 0.4 & p_val_adj < 0.05)
  diffs2[[i]]<-subset(diffs2,subset=!duplicated(gene))
  write.csv(diffs2[[i]],paste0('~/results/progenitor/',GW[i],'_diffs2.csv'))
}

# (2)##################plot###################
##fig2a 
DimPlot(project$GW09,reduction = 'tsne',cols =c('#006400', '#abf1ab','#f19bf1',
                                                '#b34ef3','#a98de2','#8968a0'),pt.size = 2)
ggsave('~/results/GW09_progenitor_tsne.pdf')
DimPlot(project$GW11,reduction = 'tsne',cols =c('#006400', '#abf1ab','#f19bf1',
                                                '#b34ef3','#a98de2','#8968a0'),pt.size = 2)
ggsave('~/results/GW11_progenitor_tsne.pdf')

##figS3d
DimPlot(project$GW12,reduction = 'tsne',cols =c('#006400', '#abf1ab','#f19bf1',
                                                '#b34ef3','#a98de2','#8968a0'),pt.size = 2)
ggsave('~/results/GW12_progenitor_tsne.pdf')

##figS3a-c
FeaturePlot(project$GW09,features = c('GAD2','DLX2'),reduction = 'tsne')
ggsave('~/results/GW09_progenitor_featureplot.pdf')
FeaturePlot(project$GW11,features = c('GAD2','DLX2'),reduction = 'tsne')
ggsave('~/results/GW11_progenitor_featureplot.pdf')
FeaturePlot(project$GW12,features = c('GAD2','DLX2'),reduction = 'tsne')
ggsave('~/results/GW12_progenitor_featureplot.pdf')

##fig2b-c and figS3d
vinplot_genes<-as.character(rev(read.csv('~/data/pro_violin.csv',header = F)[,1]))

vinplot_pro<-function(GW,vinplot_genes){
  vinplot_cm<-as.matrix(GetAssayData(project[[GW]],slot = 'data')[vinplot_genes,])
  vinplot_data<-as.data.frame(t(vinplot_cm))
  vinplot_data$celltype<-project[[GW]]@active.ident
  vinplot_data<-melt(vinplot_data,id.vars = 'celltype')
  col<-c('#abf1ab','#006400','#b34ef3','#f19bf1',
         '#a98de2','#8968a0')
  p1<-ggplot(vinplot_data,aes(x=celltype,y=value,fill=celltype))+
    geom_violin(scale='width')+facet_wrap(~variable,scales = 'free_y',ncol=1,strip.position = 'right',as.table=F,shrink = F)+
    theme_minimal_vgrid()+
    theme(axis.text.y=element_blank(),strip.text.y = element_text(angle = 0,size=10),axis.ticks.y=element_blank())+
    scale_fill_manual(values = col)+
    guides(fill=FALSE)+
    ylab('Normalizes expression (Log-transformed)')
  p1
  ggsave(paste0('~/results/',GW,'_vinlin_markers.pdf'),height = 8,width = 4)
}

##fig2b
vinplot_pro(GW = 'GW09',vinplot_genes = vinplot_genes)

##fig2c
vinplot_pro(GW = 'GW11',vinplot_genes = vinplot_genes)

##figS3d
vinplot_pro(GW = 'GW12',vinplot_genes = vinplot_genes)

##figS3e
project<-readRDS("~/LAST/project.rds")
project<-subset(project,idents=c('P1','P2','P3','P4','P5','pre-IN'))
VlnPlot(project,features = c('NKX2-1','LHX6','NR2F1','NR2F2','MEIS2','ZFHX3','ETV1'),
                             group.by = 'Week',pt.size = 0,
        cols = c('#a1d99b','#fcc5c0','#9ecae1','#be8abf'))
ggsave('~/results/progenitor_region_genes.pdf')
