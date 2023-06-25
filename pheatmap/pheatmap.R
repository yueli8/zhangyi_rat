library(pheatmap)
setwd("~/zhangyi_rat/rat")

mat<-read.table("heatmap_ogrd8_vs_ogdr0_down", head=TRUE, row.names = 1)
x<-as.matrix(mat)
x<-t(x)
pheatmap(log((x+1),2),cellwidth=6, cellheight=10, cluster_cols=T,cluster_rows = T, 
         color=colorRampPalette(c("green", "black", "red"))(100),fontsize=7)


mat<-read.table("heatmap_ogrd10_vs_ogdr0_down", head=TRUE, row.names = 1)
x<-as.matrix(mat)
x<-t(x)
pheatmap(log((x+1),2),cellwidth=10, cellheight=7, cluster_cols=T,cluster_rows = T, 
         color=colorRampPalette(c("green", "black", "red"))(100),fontsize=7)

       
