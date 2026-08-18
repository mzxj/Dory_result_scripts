library(ggplot2)
library(grid)
library(gridExtra)
library(dplyr)
library(Matrix)
genes <- read.csv('SourceData/Fetal-Liver_gene.csv')
ctsnml <- read.table('SourceData/GeneExpinCelltype_update.tsv')
colnames(ctsnml) <- genes$x
genexp <- data.frame(ctsnml[, !colSums(ctsnml)==0])
max_in_Hepatocyte <- sapply(genexp, function(x) x[3] == max(x))
genemaxexp <- genexp[, max_in_Hepatocyte]

rgset <- c(1,2,15, 16, 17)
TFset <- data.frame()
for(rg in rgset){
  rgb <- read.table(paste0('SourceData/TFnameUnique/TFname_r', rg, '.txt'))
  rexp <- genemaxexp[, which(colnames(genemaxexp) %in% rgb$V1)]
  TFset <- rbind(TFset, cbind(colnames(rexp), rep(rg, dim(rexp)[2])))
}
TFset <- as.data.frame(TFset)
pairs <- distinct(TFset, regulator=V1, region=V2) # remove duplicates
regulators <- sort(unique(TFset$V1))
i <- match(pairs$regulator, regulators)
j <- match(pairs$region, rgset)
S <- sparseMatrix(i = i, j = j, x = 1,
                  dims = c(length(regulators), length(rgset)),
                  dimnames = list(regulators, rgset))
S_dense <- as.matrix(S)


melted_data <- reshape2::melt(S_dense)
melted_data$Var2 <- factor(melted_data$Var2, levels = c(17, 16, 15, 2,1))
melted_data$Var1 <- melted_data$Var1
pdf('Binding.pdf', height = 2.2, width = 6.4)
ggplot(melted_data, aes(x = Var1, y = Var2, fill = factor(value))) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c("lightgray", "#39593f"),labels = c("No binding", "Binding")) +
  labs(x = "", y = "", fill = "") +
  theme_minimal() +
  theme(axis.title.y = element_text(size = 10),
        axis.text.x = element_text(angle = 45,hjust = 1, size=10),
        axis.text = element_text(size = 10),
        legend.position = 'top',
        legend.margin = margin(t = 0, b = -10, l = 0, r = 0))
dev.off()

