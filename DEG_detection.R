library(data.table)
library(dplyr)
heatmap.sample <- read.table("heatmap_sample.txt", header=T)
#chr7.gene <- read.table("Chr7_genes.txt", header=T, sep="\t")
chr7.gene <- read.table("Chr7q_genes.txt", header=T, sep="\t")
allsample.allgene <- fread("MLL_l2CPM.txt", header=T, sep="\t")
allsample.gene <- filter(allsample.allgene, !V1 %in% chr7.gene$Gene)
allgene <- allsample.gene$V1
allsample.gene$V1 <- NULL
allsample.gene <- as.matrix(allsample.gene)
rownames(allsample.gene) <- allgene
allsample.gene <- as.data.frame(allsample.gene)
match.sample <- match(heatmap.sample$Sample, colnames(allsample.gene))
match.sample <- na.omit(match.sample)
sample.gene <- allsample.gene[, match.sample]
variance <- c(rep(0,12429))
sample.gene <- t(sample.gene)
for (i in 1:12429) {
  variance[i] <- mad(sample.gene[,i])
  }

sample.gene <- as.data.frame(t(sample.gene))
sample.gene$Var <- variance
sample.gene <- sample.gene[order(sample.gene$Var),]
selected.sample <- colnames(sample.gene)
selected.sample <- as.data.frame(selected.sample)
colnames(selected.sample) <- "ID"
anno.sample <- read.table("del7.anno.txt", header=T, sep="\t")
anno.selected.sample <- left_join(selected.sample, anno.sample)
write.table(anno.selected.sample, "anno.txt", sep="\t", quote=F)

library(limma)
limma.table <- sample.gene
limma.table$Var <- NULL
design <- read.table("anno_20220131.txt", header=T ,sep="\t")
rownames(design) <- design$ID
design$ID <- NULL
fit <- lmFit(limma.table, design=design)
contrast.matrix <- makeContrasts(del7-NK, NK-del7, levels = design)
fit2 <- contrasts.fit(fit,contrast.matrix)
fit2 <- eBayes(fit2)
temp <- topTable(fit2, number=12429)
deg.list <- rownames(temp)

#volcano plot
temp$color <- c(rep("0", 12429))
temp$color[temp$adj.P.Val<0.00001&temp$del7...NK>1] <- 1
temp$color[temp$adj.P.Val<0.00001&temp$del7...NK<(-1)] <- 2
library(ggplot2)
g <- ggplot(temp, aes(x=del7...NK, y=-log10(adj.P.Val), color=color))
g <- g + geom_point(size=3)
g <- g + scale_color_manual(values=c("#00000032", "#8b0000", "#00008b"))
g <- g + scale_y_continuous(limits=c(0,12), expand=c(0,0))
g <- g + labs(x="log2 FC", y="-log FDR")
g <- g + geom_hline(yintercept=5, linetype="dotted")
g <- g + geom_vline(xintercept=1, linetype="dotted") + geom_vline(xintercept=-1, linetype="dotted")
g <- g + theme_classic()
g <- g + theme(legend.position="none")
g