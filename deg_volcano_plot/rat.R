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

#averageif
#setwd("~/zhangyi_rat")
#DF<-read.table("rat01",header = TRUE)
#b<-aggregate(DF[, -c(1:2)], by=list(DF$EntrzID, DF$Name), mean)
#write.table(b, file="rat01_averageif.txt")


#read file
#a<-"tmp01_averageif.txt"
#counts01<-data.frame(fread(a),check.names=FALSE,row.names=1)
#counts<-as.matrix(counts01)
#head(counts)[,1:12]


#pca 
setwd("~/zhangyi_rat/rat")
rat<-read.table("rat_tpm",header=TRUE)
write.csv(rat,"rat.csv")
colnames(rat)

#gse184053_01<-rat[,-1]
pca<-prcomp(t(rat))

plot(pca)
{plot(pca$x[,1],pca$x[,2])
  points(pca$x[1:3,1],pca$x[1:3,2], main="Top 2 PCs",col=1)
  points(pca$x[4:6,1],pca$x[4:6,2], main="Top 2 PCs",col=2)
  points(pca$x[7:9,1],pca$x[7:9,2], main="Top 2 PCs",col=3)
  points(pca$x[10:12,1],pca$x[10:12,2], main="Top 2 PCs",col=4)
}
which(pca$x[,1]< -50)
which(pca$x[,2]< -50)
percentVar <- pca$sdev^2 / sum( pca$sdev^2)


#volcano plot ogdr0_vs_con
cts<-read.table("ogdr0_vs_con",head=TRUE) 
#setup DESeqDataSetFromMatrix
cts_round<-round(cts)
countData<-cts_round
condition <- factor(c(rep("p_p",3),rep("p_h",3)), levels = c("p_p","p_h"))
coldata<-data.frame(row.names=colnames(countData),condition)
condition
coldata
dds<-DESeqDataSetFromMatrix(countData=countData, colData=coldata, design=~condition)
head(dds)
featureData <- data.frame(gene=rownames(cts))
mcols(dds) <- DataFrame(mcols(dds), featureData)
mcols(dds)
keep <- rowSums(counts(dds)) >= 10#keep sum counts >=10
dds <- dds[keep,]
dim(dds)
dds$condition <- factor(dds$condition, levels = c("p_p","p_h"))
dds$condition <- relevel(dds$condition, ref = "p_h")
dds$condition <- droplevels(dds$condition)
dds <- estimateSizeFactors(dds)
dds_norm<-counts(dds,normalized=TRUE)
head(dds_norm)
write.csv(dds_norm,"ogdr0_vs_con_norm.csv")
dds2<-DESeq(dds,fitType = "local")
resultsNames(dds2)
res <- results(dds2)
res <- results(dds2, contrast=c("condition","p_p","p_h"))
resultsNames(dds2)
write.csv(res,"ogdr0_vs_con_deg.csv")

res <- read.csv("ogdr0_vs_con_deg.csv", header=TRUE)
head(res)
with(res, plot(log2FoldChange, -log10(padj), pch=20, main="Volcano plot", xlim=c(-8,12)))
# Add colored points: red if pvalue<0.05, orange of log2FC>1, green if both)
with(subset(res, padj<.05 ), points(log2FoldChange, -log10(padj), pch=20, col="red"))
with(subset(res, abs(log2FoldChange)>1), points(log2FoldChange, -log10(padj), pch=20, col="orange"))
with(subset(res, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(padj), pch=20, col="green"))
#with(subset(res, P.Value<.05 & abs(log2FC)>1), textxy(logFC, -log10(P.Value), labs=Gene, cex=.6))
abline(h=1.3,v=1,lty=3)
abline(v=-1,lty=3)


#volcano plot ogdr10_vs_con
cts<-read.table("ogdr10_vs_con",head=TRUE) 
#setup DESeqDataSetFromMatrix
cts_round<-round(cts)
countData<-cts_round
condition <- factor(c(rep("p_p",3),rep("p_h",3)), levels = c("p_p","p_h"))
coldata<-data.frame(row.names=colnames(countData),condition)
condition
coldata
dds<-DESeqDataSetFromMatrix(countData=countData, colData=coldata, design=~condition)
head(dds)
featureData <- data.frame(gene=rownames(cts))
mcols(dds) <- DataFrame(mcols(dds), featureData)
mcols(dds)
keep <- rowSums(counts(dds)) >= 10#keep sum counts >=10
dds <- dds[keep,]
dim(dds)
dds$condition <- factor(dds$condition, levels = c("p_p","p_h"))
dds$condition <- relevel(dds$condition, ref = "p_h")
dds$condition <- droplevels(dds$condition)
dds <- estimateSizeFactors(dds)
dds_norm<-counts(dds,normalized=TRUE)
head(dds_norm)
write.csv(dds_norm,"ogdr10_vs_con_norm.csv")
dds2<-DESeq(dds,fitType = "local")
resultsNames(dds2)
res <- results(dds2)
res <- results(dds2, contrast=c("condition","p_p","p_h"))
resultsNames(dds2)
write.csv(res,"ogdr10_vs_con_deg.csv")

res <- read.csv("ogdr10_vs_con_deg.csv", header=TRUE)
head(res)
with(res, plot(log2FoldChange, -log10(padj), pch=20, main="Volcano plot", xlim=c(-8,12)))
# Add colored points: red if pvalue<0.05, orange of log2FC>1, green if both)
with(subset(res, padj<.05 ), points(log2FoldChange, -log10(padj), pch=20, col="red"))
with(subset(res, abs(log2FoldChange)>1), points(log2FoldChange, -log10(padj), pch=20, col="orange"))
with(subset(res, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(padj), pch=20, col="green"))
#with(subset(res, P.Value<.05 & abs(log2FC)>1), textxy(logFC, -log10(P.Value), labs=Gene, cex=.6))
abline(h=1.3,v=1,lty=3)
abline(v=-1,lty=3)

#volcano plot ogdr8_vs_con
cts<-read.table("ogdr8_vs_con",head=TRUE) 
#setup DESeqDataSetFromMatrix
cts_round<-round(cts)
countData<-cts_round
condition <- factor(c(rep("p_p",3),rep("p_h",3)), levels = c("p_p","p_h"))
coldata<-data.frame(row.names=colnames(countData),condition)
condition
coldata
dds<-DESeqDataSetFromMatrix(countData=countData, colData=coldata, design=~condition)
head(dds)
featureData <- data.frame(gene=rownames(cts))
mcols(dds) <- DataFrame(mcols(dds), featureData)
mcols(dds)
keep <- rowSums(counts(dds)) >= 10#keep sum counts >=10
dds <- dds[keep,]
dim(dds)
dds$condition <- factor(dds$condition, levels = c("p_p","p_h"))
dds$condition <- relevel(dds$condition, ref = "p_h")
dds$condition <- droplevels(dds$condition)
dds <- estimateSizeFactors(dds)
dds_norm<-counts(dds,normalized=TRUE)
head(dds_norm)
write.csv(dds_norm,"ogdr8_vs_con_norm.csv")
dds2<-DESeq(dds,fitType = "local")
resultsNames(dds2)
res <- results(dds2)
res <- results(dds2, contrast=c("condition","p_p","p_h"))
resultsNames(dds2)
write.csv(res,"ogdr8_vs_con_deg.csv")

res <- read.csv("ogdr8_vs_con_deg.csv", header=TRUE)
head(res)
with(res, plot(log2FoldChange, -log10(padj), pch=20, main="Volcano plot", xlim=c(-8,12)))
# Add colored points: red if pvalue<0.05, orange of log2FC>1, green if both)
with(subset(res, padj<.05 ), points(log2FoldChange, -log10(padj), pch=20, col="red"))
with(subset(res, abs(log2FoldChange)>1), points(log2FoldChange, -log10(padj), pch=20, col="orange"))
with(subset(res, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(padj), pch=20, col="green"))
#with(subset(res, P.Value<.05 & abs(log2FC)>1), textxy(logFC, -log10(P.Value), labs=Gene, cex=.6))
abline(h=1.3,v=1,lty=3)
abline(v=-1,lty=3)


#volcano plot ogdr10_vs_ogdr0
cts<-read.table("ogdr10_vs_ogdr0",head=TRUE) 
#setup DESeqDataSetFromMatrix
cts_round<-round(cts)
countData<-cts_round
condition <- factor(c(rep("p_p",3),rep("p_h",3)), levels = c("p_p","p_h"))
coldata<-data.frame(row.names=colnames(countData),condition)
condition
coldata
dds<-DESeqDataSetFromMatrix(countData=countData, colData=coldata, design=~condition)
head(dds)
featureData <- data.frame(gene=rownames(cts))
mcols(dds) <- DataFrame(mcols(dds), featureData)
mcols(dds)
keep <- rowSums(counts(dds)) >= 10#keep sum counts >=10
dds <- dds[keep,]
dim(dds)
dds$condition <- factor(dds$condition, levels = c("p_p","p_h"))
dds$condition <- relevel(dds$condition, ref = "p_h")
dds$condition <- droplevels(dds$condition)
dds <- estimateSizeFactors(dds)
dds_norm<-counts(dds,normalized=TRUE)
head(dds_norm)
write.csv(dds_norm,"ogdr10_vs_ogdr0_norm.csv")
dds2<-DESeq(dds,fitType = "local")
resultsNames(dds2)
res <- results(dds2)
res <- results(dds2, contrast=c("condition","p_p","p_h"))
resultsNames(dds2)
write.csv(res,"ogdr10_vs_ogdr0_deg.csv")

res <- read.csv("ogdr10_vs_ogdr0_deg.csv", header=TRUE)
head(res)
with(res, plot(log2FoldChange, -log10(padj), pch=20, main="Volcano plot", xlim=c(-8,12)))
# Add colored points: red if pvalue<0.05, orange of log2FC>1, green if both)
with(subset(res, padj<.05 ), points(log2FoldChange, -log10(padj), pch=20, col="red"))
with(subset(res, abs(log2FoldChange)>1), points(log2FoldChange, -log10(padj), pch=20, col="orange"))
with(subset(res, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(padj), pch=20, col="green"))
#with(subset(res, P.Value<.05 & abs(log2FC)>1), textxy(logFC, -log10(P.Value), labs=Gene, cex=.6))
abline(h=1.3,v=1,lty=3)
abline(v=-1,lty=3)


#volcano plot ogdr8_vs_ogdr0
cts<-read.table("ogdr8_vs_ogdr0",head=TRUE) 
#setup DESeqDataSetFromMatrix
cts_round<-round(cts)
countData<-cts_round
condition <- factor(c(rep("p_p",3),rep("p_h",3)), levels = c("p_p","p_h"))
coldata<-data.frame(row.names=colnames(countData),condition)
condition
coldata
dds<-DESeqDataSetFromMatrix(countData=countData, colData=coldata, design=~condition)
head(dds)
featureData <- data.frame(gene=rownames(cts))
mcols(dds) <- DataFrame(mcols(dds), featureData)
mcols(dds)
keep <- rowSums(counts(dds)) >= 10#keep sum counts >=10
dds <- dds[keep,]
dim(dds)
dds$condition <- factor(dds$condition, levels = c("p_p","p_h"))
dds$condition <- relevel(dds$condition, ref = "p_h")
dds$condition <- droplevels(dds$condition)
dds <- estimateSizeFactors(dds)
dds_norm<-counts(dds,normalized=TRUE)
head(dds_norm)
write.csv(dds_norm,"ogdr8_vs_ogdr0_norm.csv")
dds2<-DESeq(dds,fitType = "local")
resultsNames(dds2)
res <- results(dds2)
res <- results(dds2, contrast=c("condition","p_p","p_h"))
resultsNames(dds2)
write.csv(res,"ogdr8_vs_ogdr0_deg.csv")

res <- read.csv("ogdr8_vs_ogdr0_deg.csv", header=TRUE)
head(res)
with(res, plot(log2FoldChange, -log10(padj), pch=20, main="Volcano plot", xlim=c(-8,12)))
# Add colored points: red if pvalue<0.05, orange of log2FC>1, green if both)
with(subset(res, padj<.05 ), points(log2FoldChange, -log10(padj), pch=20, col="red"))
with(subset(res, abs(log2FoldChange)>1), points(log2FoldChange, -log10(padj), pch=20, col="orange"))
with(subset(res, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(padj), pch=20, col="green"))
#with(subset(res, P.Value<.05 & abs(log2FC)>1), textxy(logFC, -log10(P.Value), labs=Gene, cex=.6))
abline(h=1.3,v=1,lty=3)
abline(v=-1,lty=3)

#volcano plot ogdr10_vs_ogdr8
cts<-read.table("ogdr10_vs_ogdr8",head=TRUE) 
#setup DESeqDataSetFromMatrix
cts_round<-round(cts)
countData<-cts_round
condition <- factor(c(rep("p_p",3),rep("p_h",3)), levels = c("p_p","p_h"))
coldata<-data.frame(row.names=colnames(countData),condition)
condition
coldata
dds<-DESeqDataSetFromMatrix(countData=countData, colData=coldata, design=~condition)
head(dds)
featureData <- data.frame(gene=rownames(cts))
mcols(dds) <- DataFrame(mcols(dds), featureData)
mcols(dds)
keep <- rowSums(counts(dds)) >= 10#keep sum counts >=10
dds <- dds[keep,]
dim(dds)
dds$condition <- factor(dds$condition, levels = c("p_p","p_h"))
dds$condition <- relevel(dds$condition, ref = "p_h")
dds$condition <- droplevels(dds$condition)
dds <- estimateSizeFactors(dds)
dds_norm<-counts(dds,normalized=TRUE)
head(dds_norm)
write.csv(dds_norm,"ogdr10_vs_ogdr8_norm.csv")
dds2<-DESeq(dds,fitType = "local")
resultsNames(dds2)
res <- results(dds2)
res <- results(dds2, contrast=c("condition","p_p","p_h"))
resultsNames(dds2)
write.csv(res,"ogdr10_vs_ogdr8_deg.csv")

res <- read.csv("ogdr10_vs_ogdr8_deg.csv", header=TRUE)
head(res)
with(res, plot(log2FoldChange, -log10(padj), pch=20, main="Volcano plot", xlim=c(-8,12)))
# Add colored points: red if pvalue<0.05, orange of log2FC>1, green if both)
with(subset(res, padj<.05 ), points(log2FoldChange, -log10(padj), pch=20, col="red"))
with(subset(res, abs(log2FoldChange)>1), points(log2FoldChange, -log10(padj), pch=20, ol="orange"))
with(subset(res, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(padj), pch=20, col="green"))
#with(subset(res, P.Value<.05 & abs(log2FC)>1), textxy(logFC, -log10(P.Value), labs=Gene, cex=.6))
abline(h=1.3,v=1,lty=3)
abline(v=-1,lty=3)
