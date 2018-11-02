# ref_genome .gz files need gunzip command (.fa .gft .gff files)
# ref_genome establish ref_index with hisat2 build
# raw sequencing data NCBI-geo sra2fastq (prefetch or fastq-dump)
# .fastq quality control with fastqc
# alignment with hisat2 and samtools obtaining .srt.bam file
# read_count calculate using htseq-count
# igv check read_count
# combination read_count of each sample 
# differential expression genes analysis with DESeq2

setwd("C:/mydirectory")

# DEGs
## prepare cts and coldata variables
cts <- read.table("merge_count.txt",header=T)
rownames(cts)=cts[,1]
cts <- cts[, -c(1)]

coldata <- read.table("sample_info.txt",header=T)
coldata$condition
coldata$condition = factor(coldata$condition, c("EV","DNMT3B")) # control and treatment groups
coldata$condition

library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ condition)
dds <- DESeq(dds)
resultsNames(dds)
res <- results(dds, name="condition_DNMT3B_vs_EV") # condition_treat_vs_control
res <- lfcShrink(dds, coef="condition_DNMT3B_vs_EV")
resdf <- as.data.frame(res)
res <- res[order(res$padj),] # oder by padj
resdf <- as.data.frame(res)

write.csv(x=resdf,file="resdf.csv") #export csv file

##PCA
install.packages("ggplot2")
library(ggplot2)
vsd <- vst(dds, blind=FALSE)
head(vsd)
pcaData <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)

percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
ggsave(filename="DESeq2_PCA.pdf")


##MA plot
library(DESeq2)
library(ggplot2)
plotMA(res, ylim=c(-2,2))
ggsave(filename="DESeq2_MA.pdf")



##Heatmap - - count matrix
install.packages("pheatmap")
library("pheatmap")
library(grid)
library("RColorBrewer")
select <- order(rowMeans(counts(dds,normalized=T)), decreasing=T)[1:20]
df <- as.data.frame(colData(dds)[,c("condition")])
df
#ntd <- normTransform(dds)
#ntd
vsd
pheatmap(assay(vsd)[select,],cluster_rows = FALSE, show_rownames = FALSE, cluster_cols = FALSE)


##heatmap -- sample-sample distances
sampleDists <- dist(t(assay(vsd)))
library("RColorBrewer")
library(ggplot2)
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,clustering_distance_rows=sampleDists,clustering_distance_cols=sampleDists,col=colors)
#ggsave(filename="Heat_dis.pdf")

