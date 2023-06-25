library(stringr)
library(sva)
library(data.table)
library(readxl)
library(DESeq2)
library(DESeq)
library(pamr)
library(ggpubr)
library(Seurat)
library(cowplot)
library(scater)
library(scran)
library(BiocParallel)
library(BiocNeighbors)
library(ggplot2)
library(gplots)
library(pca3d)
library(rgl)
library(scatterplot3d)
library(FactoMineR)
library(ggfortify)
library(useful)
library(tidyverse)
library(kableExtra)
library(xfun)
library(psych)
library(limma)
library(calibrate)
library(pheatmap)

#repeat use avergeif command
setwd("~/zhangyi_rat")
DF<-read.table("rat01",header = TRUE)
b<-aggregate(DF[, -c(1:2)], by=list(DF$EntrzID, DF$Name), mean)
write.table(b, file="rat01_averageif.txt")

#计算基因长度
#思路1：计算基因在染色体的起始和结束之差
#思路2：计算每个基因的最长转录本或所有外显子之和
##############################################

#### 方法1 简单把基因在染色体上的起始位置和结束位置之差用作标准化的长度。

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("biomaRt")
library(biomaRt)

#查看基因组参数
mart = useMart('ensembl')
listDatasets(mart)

#count to TPM
library(biomaRt)
library(GenomicFeatures)
library(xlsx)
setwd("~/zhangyi_rat")
ensembl_list<-read.table("tmp11",header=FALSE)
dim(ensembl_list)

rat<-useMart("ensembl",dataset="rnorvegicus_gene_ensembl")
#connection problem,try many time. shutdown the computer,restart computer.
gene_coords=getBM(attributes=c("rgd_symbol","ensembl_gene_id","start_position","end_position"), filters="ensembl_gene_id", values=ensembl_list, mart=rat)
gene_coords$size=gene_coords$end_position - gene_coords$start_position
head(gene_coords)

a<-"tmp01_averageif.txt"
counts01<-data.frame(fread(a),check.names=FALSE,row.names=1)
counts<-as.matrix(counts01)
head(counts)[,1:12]

featureLength01<-read.table("feature_length01",head=TRUE)
featureLength02<-t(featureLength01)
featureLength<-as.vector(featureLength02)

Counts_to_tpm <- function(counts, featureLength) {
  
  # Ensure valid arguments.
  stopifnot(length(featureLength) == nrow(counts))
  
  # Compute effective lengths of features in each library.
  effLen <- featureLength
  
  # Process one column at a time.
  tpm <- do.call(cbind, lapply(1:ncol(counts), function(i) {
    rate = log(counts[,i]) - log(effLen)
    denom = log(sum(exp(rate)))
    exp(rate - denom + log(1e6))
  }))
  
  # Copy the row and column names from the original matrix.
  colnames(tpm) <- colnames(counts)
  rownames(tpm) <- rownames(counts)
  return(tpm)
}

tpm_counts<-Counts_to_tpm(counts,featureLength)
write.table(tpm_counts,"tpm_counts")



