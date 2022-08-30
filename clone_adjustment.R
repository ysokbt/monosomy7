library(dplyr)
#Read expression matrix of all samples and all genes in chromosome 7
expr.allsample.allgene <- read.table("expr.txt", header=T ,sep="\t")
allsample <- expr.allsample.allgene$ID

#Exclude low expressed genes
rownames(expr.allsample.allgene) <- expr.allsample.allgene$ID
expr.allsample.allgene$ID <- NULL
allgene.list <- colnames(expr.allsample.allgene)
allgene.list <- as.data.frame(allgene.list)
allgene.list$percentage <- c(rep(1, 694))
for (i in 1:694) {
  allgene.list$percentage[i] <- sum(expr.allsample.allgene[,i]==0)/692*100
}
include.gene.list <- filter(allgene.list, percentage<90)

expr.allsample.allgene <- read.table("expr.txt", header=T ,sep="\t")
match.include <- match(include.gene.list$allgene.list, colnames(expr.allsample.allgene))
expr.allsample.allgene <- expr.allsample.allgene[,match.include]
expr.allsample.allgene$ID <- allsample

#Read slope information
slope <- read.table("coeff_order.txt", header=T ,sep="\t")
match.slope <- match(colnames(expr.allsample.allgene), slope$Gene)
slope.match <- slope[match.slope,]
slope.match <- na.omit(slope.match)

#Order by alphabet
slope.match <- slope.match[order(slope.match$Gene),]

#Generate expression matrix of monosomy7
del7 <- read.table("del7.anno.txt", header=T, sep="\t")
del7.list <- filter(del7, Chr7.status=="-7/del(7q)")
del7.list <- del7.list$ID
del7.list <- as.data.frame(del7.list)
del7.expr <- left_join(del7.list, expr.allsample.allgene, by=c("del7.list"="ID"))

#Selected genes in monosomy7 samples by gene list in slope
match.del7 <- match(slope.match$Gene, colnames(del7.expr))
match.del7 <- na.omit(match.del7)
del7.expr.selected <- del7.expr[, match.del7]
del7.expr.selected <- as.matrix(del7.expr.selected)
rownames(del7.expr.selected) <- del7.list$del7.list
del7.expr.selected <- as.data.frame(del7.expr.selected)

#Adjusted to 100% clone by slope, clonal size, and expression levels (slope*(100-clonal size)+expression levels)
#Read 100-clonal size
clonal.size <- read.table("s_100.txt", header=T, sep="\t")

#Matched by expr matrix
match.clone <- match(rownames(del7.expr.selected), clonal.size$ID)
clonal.size <- clonal.size[match.clone,]

#Generate matrix for adjustment
matrix.adj <- matrix(c(rep(1,48*587)), nrow=48, ncol=587)
rownames(matrix.adj) <- rownames(del7.expr.selected)
colnames(matrix.adj) <- colnames(del7.expr.selected)
for (i in 1:587) {
  matrix.adj[,i] <- slope.match[i,2]*clonal.size$s_100
}

#Adjusted
del7.adj <- as.matrix(del7.expr.selected) + matrix.adj

#Generate expression matrix of normal diploid
NK <- read.table("NK_true_list.txt", header=T ,sep="\t")
NK.list <- NK$id
NK.list <- as.data.frame(NK.list)
NK.expr <- left_join(NK.list, expr.allsample.allgene, by=c("NK.list"="ID"))

#Selected genes in monosomy7 samples by gene list in slope
match.NK <- match(slope.match$Gene, colnames(NK.expr))
match.NK <- na.omit(match.NK)
NK.expr.selected <- NK.expr[, match.NK]
NK.expr.selected <- as.matrix(NK.expr.selected)
rownames(NK.expr.selected) <- NK.list$NK.list
NK.expr.selected <- as.data.frame(NK.expr.selected)

#Generate compare matrix
matrix.comp <- matrix(c(rep(1,8*587)), nrow=587, ncol=8)
rownames(matrix.comp) <- colnames(del7.adj)
colnames(matrix.comp) <- c("25percentile", "30percentile", "35percentile", "40percentile", "45percentile", 
                           "50percentile", "55percentile", "60percentile")

for (i in 1:587) {
  temp <- quantile(NK.expr.selected[,i],0.25)
  matrix.comp[i,1] <- sum(del7.adj[,i]<temp, na.rm=T)/48*100
}

for (i in 1:587) {
  temp <- quantile(NK.expr.selected[,i],0.30)
  matrix.comp[i,2] <- sum(del7.adj[,i]<temp, na.rm=T)/48*100
}

for (i in 1:587) {
  temp <- quantile(NK.expr.selected[,i],0.35)
  matrix.comp[i,3] <- sum(del7.adj[,i]<temp, na.rm=T)/48*100
}

for (i in 1:587) {
  temp <- quantile(NK.expr.selected[,i],0.40)
  matrix.comp[i,4] <- sum(del7.adj[,i]<temp, na.rm=T)/48*100
}

for (i in 1:587) {
  temp <- quantile(NK.expr.selected[,i],0.45)
  matrix.comp[i,5] <- sum(del7.adj[,i]<temp, na.rm=T)/48*100
}

for (i in 1:587) {
  temp <- quantile(NK.expr.selected[,i],0.50)
  matrix.comp[i,6] <- sum(del7.adj[,i]<temp, na.rm=T)/48*100
}

for (i in 1:587) {
  temp <- quantile(NK.expr.selected[,i],0.55)
  matrix.comp[i,7] <- sum(del7.adj[,i]<temp, na.rm=T)/48*100
}

for (i in 1:587) {
  temp <- quantile(NK.expr.selected[,i],0.60)
  matrix.comp[i,8] <- sum(del7.adj[,i]<temp, na.rm=T)/48*100
}


match.slope.adj <- match(rownames(matrix.comp), slope.match$Gene)
slope.match <- slope[match.slope.adj,]

matrix.comp.slope <- cbind(matrix.comp, slope.match$Coefficient)
colnames(matrix.comp.slope) <- c("25percentile", "30percentile", "35percentile", "40percentile", "45percentile", 
                                 "50percentile", "55percentile", "60percentile", "Slope")

#Add position information
position.chr7 <- read.table("genes_chr7.txt", header=T ,sep="\t")
position.chr7$Chromosome <- NULL
position.chr7$Start <- NULL
position.chr7$End <- NULL
matrix.comp.slope <- as.data.frame(matrix.comp.slope)
matrix.comp.slope$Gene <- rownames(matrix.comp.slope)
matrix.comp.slope.position <- left_join(matrix.comp.slope, position.chr7, by=c("Gene"="Gene.name"))

write.table(matrix.comp.slope.position, "matrix.txt", sep="\t", quote=F)
