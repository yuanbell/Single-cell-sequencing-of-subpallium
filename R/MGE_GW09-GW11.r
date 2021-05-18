library(Seurat)
library(harmony)
library(destiny)
library(slingshot)
library(RColorBrewer)
# (1)###############read data and cluster analysis
##read data
project<-readRDS('~/results/project.rds')

##Get cells annotated as MGE1 and MGE2
MGE_project<-subset(project,ident=c('MGE1','MGE2'))

###Get cell from GW9 and GW11
MGE_project<-subset(MGE_project,cells=row.names(subset(MGE_project@meta.data,
                                                       MGE_project@meta.data$Week=='GW09'|MGE_project@meta.data$Week=='GW11')))

###Clustering analysis
MGE_project<- FindVariableFeatures(object =MGE_project,selection.method = "vst",nfeatures = 1000)
MGE_project <- NormalizeData(object =MGE_project, verbose = TRUE,scale.factor = 10000)
MGE_project <- ScaleData(object = MGE_project, verbose = TRUE,
                         features = MGE_project@assays$RNA@var.features)
MGE_project <- RunPCA(object = MGE_project, npcs =30, 
                      verbose = TRUE,
                      features=VariableFeatures(MGE_project))
MGE_project<-RunHarmony(MGE_project,"Week", plot_convergence = TRUE,dims.use = 1:30
)
MGE_project <- FindNeighbors(MGE_project, reduction = "harmony", dims = 1:30) 
MGE_project <- RunTSNE(object = MGE_project, reduction = "harmony",dims=1:30)
MGE_project<- FindClusters(resolution = 1,MGE_project)

##Perform difference analysis and find marker genes
diffs <- FindAllMarkers(object = MGE_project, logfc.threshold = 0.25, 
                        only.pos = TRUE, min.pct = 0.1,assay="RNA")
diffs$cluster <-as.factor(as.numeric(as.vector(diffs$cluster))+1)
diffs2<-subset(diffs[grep("^RP[L|S]",diffs$gene, ignore.case = FALSE,invert=TRUE),],
               subset=avg_logFC>0.75 & (pct.1 > pct.2) & pct.1> 0.5)
diffs2<-subset(diffs2,subset=!duplicated(gene))

##annoated cluster
new_celltype<-c('ZEB2+/MAF+','POU3F2+/CNTNAP2+','NR2F1+/MEIS2+','LHX8+/NKX2-1+','ANGPT2+/CRABP1+')
names(new_celltype)<-levels(MGE_project)
MGE_project<-RenameIdents(MGE_project,new_celltype)


# (2)###################DiffusionMap and Trajectory analysis ####################

##Diffusion
ct <- as.ExpressionSet(as.data.frame(Embeddings(MGE_project,'harmony')))
ct$celltype<-MGE_project@active.ident
dm <- DiffusionMap(ct,sigma='global',k=5,density_norm = T,n_pcs = NA)

##set colors
cell_cluster<-as.data.frame(MGE_project@active.ident)
colnames(cell_cluster)<-'Cluster'
cell_cluster$Cluster<-as.factor(cell_cluster$Cluster)
colnames(cell_cluster)<-'Cluster'
cols=c(brewer.pal(12,"Paired"),brewer.pal(8,"Dark2"),brewer.pal(9,"Set1"),brewer.pal(8,"Set2"),brewer.pal(12,"Set3"),brewer.pal(8,"Pastel2"),brewer.pal(9,"Pastel1"),brewer.pal(8,"Accent"))
col=unique(cols)[c(1,5,9,13,16)]
cluster<-unique(cell_cluster$Cluster)
for(i in 1:length(cluster)){
  index=which(cell_cluster$Cluster==cluster[i])
  cell_cluster[index,2]<-col[i]
}

##slingshot analysis
##2and3 dimension to visable
sling<-slingshot(data = dm@eigenvectors[,2:3],clusterLabels = MGE_project@active.ident,
                 start.clus='ANGPT2+/CRABP1+',end.clus=c('NR2F1+,MEIS2+','LHX8+,NKX2-1+'))
##1and2 dimension to calculate pseudotime
sling2<-slingshot(data = dm@eigenvectors[,1:2],clusterLabels = MGE_project@active.ident,
                 start.clus='ANGPT2+/CRABP1+',end.clus=c('NR2F1+,MEIS2+','LHX8+,NKX2-1+'))

##this function are reference from https://github.com/IStevant/XX-XY-mouse-gonad-scRNA-seq
rankKeepNA <- function(x) {
  return(
    ifelse(
      is.na(x),
      NA,
      rank(
        x, 
        na.last = TRUE, 
        ties.method="random"
      )
    )
  )
}
get_pseudotime <- function(pseudotime, wthres=wthres, ranked=TRUE){
  pseudoT <- list()
  for(lineage in 1:length(pseudotime@curves))local({
    curve <- pseudotime@curves[[lineage]]
    lambda <- curve$lambda
    weight <- curve$w
    ps <- curve$lambda
    ps[weight < wthres] <- NA
    if (ranked==TRUE){
      ps <- rankKeepNA(ps)
    }
    pseudoT[[lineage]] <<- ps
  })
  df <- t(do.call("rbind",pseudoT))
  colnames(df) <- names(pseudotime@curves)
  return(df)
}

##get pseudotime
pseudotime <- get_pseudotime(sling2, wthres=0.8)
pseudotime_lin <- pseudotime[,"curve1"]
max_pseudotime <- max(pseudotime_lin, na.rm = TRUE)
pseudotime_lin1_percent <- (pseudotime_lin*100)/max_pseudotime
pseudotime_lin <- pseudotime[,"curve2"]
max_pseudotime <- max(pseudotime_lin, na.rm = TRUE)
pseudotime_lin2_percent <- (pseudotime_lin*100)/max_pseudotime
pseudotime[,"curve1"] <- pseudotime_lin1_percent
pseudotime[,"curve2"] <- pseudotime_lin2_percent
colnames(pseudotime)<-c('lin1','lin2')

##calculate pseudotime associated genes
cal_pseu_genes<-function(qvalue,lin){
  expr<-GetAssayData(MGE_project,slot = 'data')[,]
  sample_sheet <- data.frame(cells=colnames(expr), stages=MGE_project$Week, 
                             cellType=MGE_project@active.ident)
  rownames(sample_sheet)<- colnames(expr)
  gene_annotation <- as.data.frame(rownames(expr))
  rownames(gene_annotation)<- rownames(expr)
  colnames(gene_annotation)<- "genes"
  pd <- new("AnnotatedDataFrame", data = sample_sheet)
  fd <- new("AnnotatedDataFrame", data = gene_annotation)
  monocl <- newCellDataSet(
    as(expr_matrix, "sparseMatrix"),
    phenoData = pd,
    featureData = fd,
    lowerDetectionLimit=0.5,
    expressionFamily=negbinomial.size()
  )
  monocle_project <- detectGenes(monocle_project, min_expr = 5)
  monocle_project <- monocle_project[fData(monocle_project)$num_cells_expressed > 10, ]
  monocle_project <- estimateSizeFactors(monocle_project)
  monocle_project <- estimateDispersions(monocle_project)
  
  ##add pseudotime to monocle_project
  lineage_pseudotime_all <- pseudotime[,2]
  lineage <- pseudotime[!is.na(pseudotime[,lin]),2,drop=FALSE]
  monocle_project$Pseudotime <- lineage_pseudotime_all
  monocle_project <- monocle_project[,row.names(lineage)]
  monocle_project <- detectGenes(monocle_project, min_expr = 5)
  monocle_project <- monocle_project[fData(monocle_project)$num_cells_expressed >= 10, ]
  ##calculate pseudotime associated genes
  DE_genes<- differentialGeneTest(monocle_project, 
    fullModelFormulaStr="~sm.ns(Pseudotime, df=3)")
  significant_genes<-subset(DE_genes,DE_genes$qval<0.001)
  return(significant_genes)
}
significant_genes_lin1<-cal_pseu_genes(0.001,'lin1')
significant_genes_lin2<-cal_pseu_genes(0.001,'lin2')
##select TF genes
homo_TF<-read.csv('~/data/homo_TF.csv')
lin1_TF<-intersect(homo_TF$Symbol,signature_genes_lin1$genes)
lin2_TF<-intersect(homo_TF$Symbol,signature_genes_lin2$genes)



## (3)######################plot#########################

##fig4a
DimPlot(MGE_project, reduction = "tsne", pt.size = 1,
        cols = c('#A6CEE3','#FB9A99','#1B9E77','#CAB2D6','#E7298A'))
ggsave('~/results/MGE_celltype_tsne.pdf',width = 6.6,height = 4)

##fig4b
heatmap_genes<-c('ANGPT2','CRABP1','ETV1','ZEB2','MAF','ERBB4','POU3F2','PEG10','CNTNAP2','NR2F1','GSTP1',
                 'MEIS2','LHX8','NKX2-1','ROBO2')
MGE_average<-AverageExpression(MGE_project,return.seurat = T)
cm<-GetAssayData(MGE_average,slot = 'data')[heatmap_genes,]
cm<-as.data.frame(cm)
pheatmap::pheatmap(cm,scale = 'row',cluster_rows = F,cluster_cols = F,border_color = '#3f3f44',
                   color =colorRampPalette(c('#3da4ab','#fcf9ea','#fe8a71'))(100),
                   filename = '~/results/MGE_heatmap.pdf')

##fig4c
pdf('~/results/sling_shot_MGE.pdf')
plot(dm@eigenvectors[,1:2], pch=20,cex=1,col=cell_cluster$V2)
lines(SlingshotDataSet(sling), lwd=2, col='black')
legend('bottomright',legend=unique(cell_cluster$Cluster),col=unique(cell_cluster$V2),pch=20,xjust=1,yjust=1,cex = 0.7)
dev.off()

##fig4g-i
expr<-GetAssayData(MGE_project,slot = 'data')[row.names(MGE_project) %in% c('EBF1','MEIS2','ISL1','IGFBP5','TAC1','ZFHX3','ZNF503',
                                                                          'LHX8','ZIC1','ZIC2','ZIC4','NKX2-1','LHX6'),]
for (i in 1:nrow(expr)) {
  pseudotime_data <- data.frame(
    pseudotime=numeric(),
    lineage=numeric(),
    gene=numeric()
  )
  for (lineages in 1:ncol(as.data.frame(pseudotime))){
    lineage <- pseudotime[,lineages]
    sub_data <-expr[i,]
    
    data <- data.frame(
      pseudotime=lineage,
      lineage=paste("Lineage ",lineages, sep=""),
      gene=sub_data
    )
    
    colnames(data) <- c("pseudotime", "lineage","gene")
    
    pseudotime_data <- rbind(pseudotime_data, data)
  }
  ggplot(pseudotime_data, aes(x=pseudotime, y=gene))+
    geom_smooth(aes(group=lineage, color=lineage, fill=lineage), na.rm = TRUE, method="loess", span=0.75)+
    ylab("Exprision") +
    theme_bw() +
    theme(
      axis.text=element_text(size=16),
      axis.title=element_text(size=16),
      legend.text = element_text(size =16),
      legend.title=element_blank(),
      plot.title = element_text(size=18, face="bold.italic", hjust = 0.5),
      aspect.ratio=0.5,
      legend.position="bottom",
      strip.text.x = element_text(size = 16))+
    ggtitle(row.names(expr)[i])+
    scale_color_manual(values = c("#ffa17f","#654062"))+
    scale_fill_manual(values = c("#ffa17f","#654062"))
  ggsave(filename = paste0('~/results/',row.names(expr)[i],'.pdf'),height = 5,width = 8)
}

##figS8b
MGE_project@reductions$dm<-CreateDimReducObject(dm@eigenvectors,key = 'DC')
FeaturePlot(MGE_project,reduction = 'dm',features = 'MAF',dims = 2:3,combine = T,order = T)+
  xlim(range(df_plot$DC2))+
  ylim(range(df_plot$DC3))+
  scale_color_gradient2(high = '#fe8a71',mid = '#fcf9ea',low = '#3da4ab',
                        midpoint = (max(df_plot[,'MAF'])/3))
