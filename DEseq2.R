# RNA-seq data analization with DESeq2
setwd("/Users/tingyang/Desktop/UCSD/SP23/BENG285/Team_3_LUAD/P5/1yr")
# Load packages
library("biomaRt")
library("dplyr")
library("tidyverse")
library("DESeq2")
library("ggplot2")
library("ggrepel")
library("EnhancedVolcano")
library("apeglm")

# Load gene expression data (raw GE data: remove recurrent and normal tumors; remove 0 value GE genes)
# Total 511 patients 
cts <- as.matrix(read.csv("cts.txt", sep = "\t", row.names = 1))

# Take a look at our GE dataframe
as_tibble(cts)

coldata <- read.csv("coldata.txt", sep = "\t", header = TRUE, row.names = 1)

# Check if the row names in coldata are the same as the column names of cts data
all(rownames(coldata) %in% colnames(cts))

diff_names <- setdiff(rownames(coldata), colnames(cts))

# Print the difference
print(diff_names)

# Drop the TCGA.44.6144
coldata <- coldata[-which(rownames(coldata) == "TCGA.44.6144"), ]

# Check if order is the same
all(rownames(coldata) == colnames(cts))

# Check if the row names in coldata are the same as the column names of cts data
all(rownames(coldata) %in% colnames(cts))

# Convert the label column to a factor
coldata$Label <- factor(coldata$Label)

# Construct a DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ Label)
dds


# Pre-filtering (remove genes with less than 10 reads)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Genes differential expressed between patients
dds <- DESeq(dds)
res <- results(dds)
res

##############################################################################

# p-values and adjusted p-values
summary(res)

## We can: reorder the results by the smallest p value
resOrdered <- res[order(res$pvalue),]

# See the number of genes with padj < 0.05
sum(res$padj < 0.05, na.rm=TRUE)

top_100 <- head(res[order(res$padj),], 100)
write.table(top_100, "top_100_DEgenes.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

d <- plotCounts(dds, gene=which.min(res$padj), intgroup="Label",returnData=TRUE)

ggplot(d, aes(x=Label, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))

# Log fold change shrinkage for visualization and ranking
resultsNames(dds)
resLFC <- lfcShrink(dds, coef="Label_1_vs_0", type="apeglm")
resLFC

########################
##### Volcano plot #####
########################

fc <- 1
q <- 0.05

keyvals <- rep('grey75', nrow(res))
names(keyvals) <- rep('INSIGNIFICANT', nrow(res))

keyvals[which(abs(res$log2FoldChange) > fc & res$padj > q)] <- 'grey50'
  
keyvals[which(abs(res$log2FoldChange) > -fc & res$padj < q)] <- 'grey25'
  
keyvals[which(res$log2FoldChange < -fc & res$padj < q)] <- 'royalblue'
names(keyvals)[which(res$log2FoldChange < -fc & res$padj < q)] <- 'DOWN'
  
keyvals[which(res$log2FoldChange > fc & res$padj < q)] <- 'red'
names(keyvals)[which(res$log2FoldChange > fc & res$padj < q)] <- 'UP'
    
unique(keyvals)
unique(names(keyvals))
    
topGene <- res[which(res$padj < 10^-8),]

EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'padj',
                selectLab = c(row.names(topGene)),
                xlab = bquote(~Log[2]~ 'fold change'),
                ylab = bquote(~Log[10]~ adjusted~italic(P)),
                title = 'Clinical endpoint DSS 1yr',
                subtitle = 'EnhancedVolcano',
                pCutoff = 10e-8,
                FCcutoff = 1.0,
                pointSize = 3.0,
                colCustom = keyvals,
                labSize = 4.0,
                labCol = 'black',
                labFace = 'bold',
                colAlpha = 0.5,
                legendPosition = 'right',
                legendLabSize = 11,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                xlim = c(-10,10),
                ylim = c(-1,15),
                border = 'full',
                colConnectors = 'grey')


# Log fold change
EnhancedVolcano(resLFC,
                lab = rownames(resLFC),
                x = 'log2FoldChange',
                y = 'pvalue')
