library(DESeq2)
library(tidyverse)
library(ggplot2)  



counts_raw <- read.table("counts.txt", header = TRUE, comment.char = "#", sep = "\t", check.names = FALSE)
counts <- counts_raw[, c(1, 7:ncol(counts_raw))]
# Set gene IDs as rownames
rownames(counts) <- counts$Geneid
counts$Geneid <- NULL
colnames(counts) <- gsub(".bam", "", colnames(counts))
head(counts)

# colnames(counts)[1] <- ""
# rownames(counts) <- counts$Geneid
# head(counts)
# counts <- counts[, -1]

meta <- read.csv("../metadata.csv", row.names = 1)
meta <- meta[colnames(counts), ]  # reorder to match count columns

# rownames(meta) <- meta$SampleID   # SRR4341972, SRR4341973, ...
meta <- meta[meta$Tissue != "Embryo", ]
counts <- counts[, rownames(meta)]  # Filter count matrix too
all(colnames(counts) %in% rownames(meta))  # Should be TRUE

# meta$group <- factor(paste(meta$Tissue, meta$DAA, sep = "_"))
# dds <- DESeqDataSetFromMatrix(countData = counts,
#                               colData = meta,
#                               design = ~ group)

# dds <- DESeqDataSetFromMatrix(countData = counts,
#                               colData = meta,
#                               design = ~ Condition)
# keep<- rowSums(counts(dds)) >=10
# dds<-dds[keep,]
# dds <- DESeq(dds)


meta$Tissue <- factor(meta$Tissue)
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = meta,
                              design = ~ Tissue)
dds <- DESeq(dds)
res<-results(dds)
res
summary(res)
res_0.01<-results(dds,alpha=0.01)

summary(res_0.01)

resLFC <- lfcShrink(dds, coef="Tissue_Seed_vs_PodWall", type="apeglm")
resLFC
summary(resLFC)
# Compare Seed_10 vs PodWall_10
# res1 <- results(dds, contrast = c("Tissue", "Seed_10", "PodWall_10"))

# Save the MA plot
png("MA_plot.png", width=800, height=600)
plotMA(res, main="DESeq2 MA Plot", ylim=c(-5,5))
dev.off()


png("MA_plot_LFC.png", width=800, height=600)
plotMA(resLFC, main="DESeq2 MA Plot", ylim=c(-5,5))
dev.off()
# Convert to data frame and clean
res_df <- as.data.frame(res)
res_df <- na.omit(res_df)
res_df$significant <- res_df$padj < 0.05

# Save volcano plot
png("Volcano_plot.png", width=800, height=600)

with(res_df, plot(log2FoldChange, -log10(padj),
                  pch=20,
                  col=ifelse(significant, "red", "gray"),
                  main="Volcano Plot",
                  xlab="log2(Fold Change)",
                  ylab="-log10(Adjusted P-value)"))

abline(h = -log10(0.05), col = "blue", lty = 2)   # FDR threshold
abline(v = c(-1, 1), col = "black", lty = 3)       # fold change thresholds

dev.off()


# If you have many samples, vst is faster
vsd <- vst(dds, blind = TRUE)  # You can use rlog(dds) instead if dataset is small

# Optional: check a few transformed values
head(assay(vsd)[, 1:5])

# PCA plot and save to PNG
png("PCA_plot.png", width = 800, height = 600)
plotPCA(vsd, intgroup = "Tissue")  # or "group" if you used meta$group
dev.off()















