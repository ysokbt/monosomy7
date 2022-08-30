library(tidyverse)
library(ggplot2)
library(ggpmisc)
library(gridExtra)

#Data transformation
CNV <- read.table("ID_CNV.txt", header=T, sep="\t")
CNV_ID <- CNV$ID

CPM <- read.table("CPM.txt", header=T, sep="\t")
CPM_long <- CPM %>% tidyr::gather(key=Gene, value=CPM, -ID)

Region <- read.table("Region.txt", header=T, sep="\t")
Region_long <- Region %>% tidyr::gather(key=Gene, value=CPM, -ID)

CPM_Region <- left_join(CPM_long, Region_long, by=c("ID", "Gene"))
colnames(CPM_Region) <- c("ID", "Gene", "CPM", "Deletion")

CPM_Region_CNV_all <- left_join(CPM_Region, CNV)

CPM_Region_CNV_all <- na.omit(CPM_Region_CNV_all)
CPM_Region_CNV <- CPM_Region_CNV_all %>% dplyr::filter(Deletion==1)

Gene_unique <- CPM_Region_CNV %>% distinct(Gene)
temp <- c(1:668)
temp <- as.data.frame(temp)
colnames(temp) <- c("Number")
Gene_number <- cbind(Gene_unique, temp)
CPM_Region_CNV_Number <- left_join(CPM_Region_CNV, Gene_number, by=c("Gene"))

#Coefficients detection
coeff <- as.data.frame(matrix(nrow=668, ncol=5))
colnames(coeff) <- c("Gene", "Coefficient", "R-squared", "p.value", "FDR")
for (i in 1:668) {
  temp <- CPM_Region_CNV_Number %>% dplyr::filter(Number==i)
  coeff[i,1] <- temp[1,2]
  temp_2 <- lm(temp$CPM~temp$CNV)
  coeff[i,2] <- temp_2[["coefficients"]][["temp$CNV"]]
  temp_3 <- summary(temp_2)
  coeff[i,3] <- temp_3$r.squared
  temp_4 <- temp_3$fstatistic
  coeff[i,4] <- 1-pf(temp_4["value"],temp_4["numdf"],temp_4["dendf"])
}
coeff$FDR <- p.adjust(coeff$p.value, method="BH")
coeff_order <- order(coeff$Coefficient)
write.table(coeff_order, "coeff_order.txt", sep="\t", quote=F)
