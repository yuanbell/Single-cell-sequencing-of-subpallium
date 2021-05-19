library(Seurat)
library(destiny)
library(data.table)
# (1)#################read data and diffusionmap analysis###################
##read data
project <- readRDS('~/results/project.rds')
MGE_P<-subset(project,ident=c('P1','P2','P3','P4','P5','P6','MGE1','MGE2'))
CGE_P<-subset(project,ident=c('P1','P2','P3','P4','P5','P6','CGE1'))
LGE_P<-subset(project,ident=c('P1','P2','P3','P4','P5','P6','LGE1','LGE2','LGE3'))

##diffusion map dimension reduce
dm_analysis<-function(proejct){
  ct <- as.ExpressionSet(as.data.frame(Embeddings(project,'harmony')))
  ct$celltype<-project@active.ident
  dm <- DiffusionMap(ct,sigma='local',k=1000,density_norm = T,n_pcs = NA)
  return(dm)
}

dm_MGE<-dm_analysis(MGE_P)
dm_CGE<-dm_analysis(CGE_P)
dm_LGE<-dm_analysis(LGE_P)


# (2)###################plot####################
##fig1f
cell_cluster<-as.data.frame(MGE_P@active.ident)
colnames(cell_cluster)<-'Cluster'
col=as.character(read.csv('~/data/MGE_P_cols.csv',header = F)$V1)
for(i in 1:length(levels(cell_cluster$Cluster))){
  index=which(cell_cluster$Cluster==levels(cell_cluster$Cluster)[i])
  cell_cluster[index,2]<-col[i]
}
plot(dm_MGE@eigenvectors[,c(1,2)],col=cell_cluster$V2,pch=20,cex=0.5)
legend('bottomright',legend=unique(cell_cluster$Cluster),
       col=unique(cell_cluster$V2),pch=20,xjust=1,yjust=1,cex = 0.7)

##fig1g
cell_cluster<-as.data.frame(CGE_P@active.ident)
colnames(cell_cluster)<-'Cluster'
col=as.character(read.csv('~/data/CGE_P_cols.csv',header = F)$V1)
for(i in 1:length(levels(cell_cluster$Cluster))){
  index=which(cell_cluster$Cluster==levels(cell_cluster$Cluster)[i])
  cell_cluster[index,2]<-col[i]
}
plot(dm_CGE@eigenvectors[,c(1,2)],col=cell_cluster$V2,pch=20,cex=0.5)
legend('bottomright',legend=unique(cell_cluster$Cluster),
       col=unique(cell_cluster$V2),pch=20,xjust=1,yjust=1,cex = 0.7)

##fig1h
cell_cluster<-as.data.frame(LGE_P@active.ident)
colnames(cell_cluster)<-'Cluster'
col=as.character(read.csv('~/data/LGE_P_cols.csv',header = F)$V1)
for(i in 1:length(levels(cell_cluster$Cluster))){
  index=which(cell_cluster$Cluster==levels(cell_cluster$Cluster)[i])
  cell_cluster[index,2]<-col[i]
}
plot(dm_LGE@eigenvectors[,c(1,2)],col=cell_cluster$V2,pch=20,cex=0.5)
legend('bottomright',legend=unique(cell_cluster$Cluster),
       col=unique(cell_cluster$V2),pch=20,xjust=1,yjust=1,cex = 0.7)
