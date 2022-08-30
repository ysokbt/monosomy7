library(tidyverse)
expr.allsample.allgene <- read.table("expr.txt", header=T ,sep="\t")
allsample <- expr.allsample.allgene$ID

del7.expr <- read.table("del7_expr_adjusted.txt", header=T, sep="\t")
ID.del7 <- del7.expr$ID
Other.expr <- filter(expr.allsample.allgene, !ID %in% ID.del7)
ID.other <- Other.expr$ID

del7.expr$ID <- NULL
gene.selected <- colnames(del7.expr)
del7.expr <- read.table("del7_expr_adjusted.txt", header=T, sep="\t")

match.other <- match(gene.selected, colnames(Other.expr))
Other.expr.selected <- Other.expr[,match.other]
Other.expr.selected <- as.matrix(Other.expr.selected)
rownames(Other.expr.selected) <- ID.other
Other.expr.selected <- as.data.frame(Other.expr.selected)

del7.expr$ID <- NULL
rownames(del7.expr) <- ID.del7

expr.allsample <- rbind(del7.expr, Other.expr.selected)
expr.allsample <- as.matrix(expr.allsample)
expr.allsample <- scale(expr.allsample)
expr.allsample <- t(expr.allsample)
expr.allsample <- as.data.frame(expr.allsample)

heatmap.gene <- read.table("selected_genes.txt", header=T, sep="\t")
match.gene <- match(heatmap.gene$Gene, rownames(expr.allsample))
expr.heatmap <- expr.allsample[match.gene,]

expr.heatmap <- read.table("expr.allsample.heatmap.txt", header=T, row.names=1, sep="\t")

anno <- read.table("all.anno.txt", header=T, row.names=1, sep="\t")
ann_color <- list()
ann_color$Chr7.status <- c("red", "grey")
names(ann_color$Chr7.status) <- c("-7/del(7q)", "Other")
ann_color$Chr7.description <- c("Green", "blue", "purple", "grey")
names(ann_color$Chr7.description) <- c("iso", "plus", "CK", "Other")

library(pheatmap)
library(gplots) 
bk=unique(c(seq(-4,-2, length=100),seq(-2,2, length=400),seq(2,4,length=100)))
hmcols=colorpanel(n=length(bk)-1,low="blue",mid="lightyellow",high="red")


pheatmap(expr.heatmap, col=hmcols, breaks=bk, clustering_method="ward.D2", annotation=anno, annotation_colors=ann_color,
         show_colnames=F, show_rownames=F, fontsize_col = 4)

expr.minimal <- expr.heatmap[1:4,]

pheatmap(expr.minimal, col=hmcols, breaks=bk, clustering_method="ward.D2", annotation=anno, annotation_colors=ann_color,
         show_colnames=F, show_rownames=F, fontsize_col = 4)
