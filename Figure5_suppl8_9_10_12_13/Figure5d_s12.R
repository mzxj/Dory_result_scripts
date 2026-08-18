library(ggplot2)
library(grid)
library(gridExtra)
library(tidyr)
library(dplyr)
library(Matrix)

## considering the match with region bar, plot one by one.


GeneSet <- c("Anxa8", "Fhit", "Oit1", "Mycn", "Myc") ### Oit1 is Fam3d. in the paper figures, the name is Fam3d
########
## gid=1
############
gid=1
gene <- GeneSet[gid]
regulators <- read.table(paste0('SourceData/Regulators_Track/', gene, '_regulators.tsv'), sep = "\t", header = FALSE)
regulators <- paste0(toupper(substring(regulators$V1, 1, 1)), tolower(substring(regulators$V1, 2)))
path <- paste0('SourceData/TFnameUnique/', gene)
rglt_list <- list.files(path)
rglt_rg_set <- c()
for(rgid in 1:length(rglt_list)){
  rg <- sub(".*_(R[0-9]+)\\..*", "\\1", rglt_list[rgid])
  rgfile <- paste0(path, '/', rglt_list[rgid])
  if (file.info(rgfile)$size > 0) {
    rgltset <- read.table(rgfile, header = FALSE)
    regulators_in_rglt_rg <- rgltset$V1[which(rgltset$V1 %in% regulators)]
    rglt_rg_set <- rbind(rglt_rg_set, cbind(regulators_in_rglt_rg, rep(rg, length(regulators_in_rglt_rg))))
  }
}
rglt_rg_set <- as.data.frame(rglt_rg_set)

regions <- paste0("R", 1:length(rglt_list))
pairs <- distinct(rglt_rg_set, regulator=regulators_in_rglt_rg, region=V2) 
i <- match(pairs$regulator, regulators)
j <- match(pairs$region, regions)
S <- sparseMatrix(i = i, j = j, x = 1,
                  dims = c(length(regulators), length(regions)),
                  dimnames = list(regulators, regions))
S_dense <- as.matrix(S)


melted_data <- reshape2::melt(as.matrix(S_dense))

melted_data$Var2 <- factor(melted_data$Var2, levels = paste0("R",1:length(rglt_list)))
melted_data$Var1 <- melted_data$Var1
p <-ggplot(melted_data, aes(x = Var2, y = Var1, fill = factor(value))) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c("lightgray", "orange"),labels = c("No binding", "Binding")) +  
  labs(x = "", y = "", fill = "") +
  theme_minimal() +
  theme(axis.title.y = element_text(size = 10),
        axis.text.x = element_text(angle = 90, size=7),
        axis.text = element_text(size = 10),
        legend.position = 'none',
        legend.margin = margin(t = 0, b = -10, l = 0, r = 0))
pdf(paste0('Binding_', gene, '.pdf'), height = 2.3, width = 6.96)
print(p)
dev.off()
  



########
## gid=2
############
gid=2
gene <- GeneSet[gid]
regulators <- read.table(paste0('SourceData/Regulators_Track/', gene, '_regulators.tsv'), sep = "\t", header = FALSE)
regulators <- paste0(toupper(substring(regulators$V1, 1, 1)), tolower(substring(regulators$V1, 2)))
path <- paste0('SourceData/TFnameUnique/', gene)
rglt_list <- list.files(path)
rglt_rg_set <- c()
for(rgid in 1:length(rglt_list)){
  rg <- sub(".*_(R[0-9]+)\\..*", "\\1", rglt_list[rgid])
  rgfile <- paste0(path, '/', rglt_list[rgid])
  if (file.info(rgfile)$size > 0) {
    rgltset <- read.table(rgfile, header = FALSE)
    regulators_in_rglt_rg <- rgltset$V1[which(rgltset$V1 %in% regulators)]
    rglt_rg_set <- rbind(rglt_rg_set, cbind(regulators_in_rglt_rg, rep(rg, length(regulators_in_rglt_rg))))
  }
}
rglt_rg_set <- as.data.frame(rglt_rg_set)

regions <- paste0("R", 1:length(rglt_list))
pairs <- distinct(rglt_rg_set, regulator=regulators_in_rglt_rg, region=V2) 
i <- match(pairs$regulator, regulators)
j <- match(pairs$region, regions)
S <- sparseMatrix(i = i, j = j, x = 1,
                  dims = c(length(regulators), length(regions)),
                  dimnames = list(regulators, regions))
S_dense <- as.matrix(S)

melted_data <- reshape2::melt(as.matrix(S_dense))

melted_data$Var2 <- factor(melted_data$Var2, levels = paste0("R",1:length(rglt_list)))
melted_data$Var1 <- melted_data$Var1
p <-ggplot(melted_data, aes(x = Var2, y = Var1, fill = factor(value))) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c("lightgray", "orange"),labels = c("No binding", "Binding")) +  
  labs(x = "", y = "", fill = "") +
  theme_minimal() +
  theme(axis.title.y = element_text(size = 10),
        axis.text.x = element_text(angle = 90, size=7),
        axis.text = element_text(size = 10),
        legend.position = 'none',
        legend.margin = margin(t = 0, b = -10, l = 0, r = 0))
pdf(paste0('Binding_', gene, '.pdf'), height = 0.9, width = 6.8)
print(p)
dev.off()
  
  
  
  
  
  
########
## gid=3
############  
gid=3
gene <- GeneSet[gid]
regulators <- read.table(paste0('SourceData/Regulators_Track/', gene, '_regulators.tsv'), sep = "\t", header = FALSE)
regulators <- paste0(toupper(substring(regulators$V1, 1, 1)), tolower(substring(regulators$V1, 2)))
path <- paste0('SourceData/TFnameUnique/', gene)
rglt_list <- list.files(path)
rglt_rg_set <- c()
for(rgid in 1:length(rglt_list)){
  rg <- sub(".*_(R[0-9]+)\\..*", "\\1", rglt_list[rgid])
  rgfile <- paste0(path, '/', rglt_list[rgid])
  if (file.info(rgfile)$size > 0) {
    rgltset <- read.table(rgfile, header = FALSE)
    regulators_in_rglt_rg <- rgltset$V1[which(rgltset$V1 %in% regulators)]
    rglt_rg_set <- rbind(rglt_rg_set, cbind(regulators_in_rglt_rg, rep(rg, length(regulators_in_rglt_rg))))
  }
}
rglt_rg_set <- as.data.frame(rglt_rg_set)

regions <- paste0("R", 1:length(rglt_list))
pairs <- distinct(rglt_rg_set, regulator=regulators_in_rglt_rg, region=V2) 
i <- match(pairs$regulator, regulators)
j <- match(pairs$region, regions)
S <- sparseMatrix(i = i, j = j, x = 1,
                  dims = c(length(regulators), length(regions)),
                  dimnames = list(regulators, regions))
S_dense <- as.matrix(S)


melted_data <- reshape2::melt(as.matrix(S_dense))

melted_data$Var2 <- factor(melted_data$Var2, levels = paste0("R",1:length(rglt_list)))
melted_data$Var1 <- melted_data$Var1
p <-ggplot(melted_data, aes(x = Var2, y = Var1, fill = factor(value))) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c("lightgray", "orange"),labels = c("No binding", "Binding")) +  
  labs(x = "", y = "", fill = "") +
  theme_minimal() +
  theme(axis.title.y = element_text(size = 10),
        axis.text.x = element_text(angle = 90, size=7),
        axis.text = element_text(size = 10),
        legend.position = 'none',
        legend.margin = margin(t = 0, b = -10, l = 0, r = 0))
pdf(paste0('Binding_', gene, '.pdf'), height = 2.6, width = 6.96)
print(p)
dev.off()

  
  
  
  
  
########
## gid=4
############  
gid=4
gene <- GeneSet[gid]
regulators <- read.table(paste0('SourceData/Regulators_Track/', gene, '_regulators.tsv'), sep = "\t", header = FALSE)
regulators <- paste0(toupper(substring(regulators$V1, 1, 1)), tolower(substring(regulators$V1, 2)))
path <- paste0('SourceData/TFnameUnique/', gene)
rglt_list <- list.files(path)
rglt_rg_set <- c()
for(rgid in 1:length(rglt_list)){
  rg <- sub(".*_(R[0-9]+)\\..*", "\\1", rglt_list[rgid])
  rgfile <- paste0(path, '/', rglt_list[rgid])
  if (file.info(rgfile)$size > 0) {
    rgltset <- read.table(rgfile, header = FALSE)
    regulators_in_rglt_rg <- rgltset$V1[which(rgltset$V1 %in% regulators)]
    rglt_rg_set <- rbind(rglt_rg_set, cbind(regulators_in_rglt_rg, rep(rg, length(regulators_in_rglt_rg))))
  }
}
rglt_rg_set <- as.data.frame(rglt_rg_set)

regions <- paste0("R", 1:length(rglt_list))
pairs <- distinct(rglt_rg_set, regulator=regulators_in_rglt_rg, region=V2) 
i <- match(pairs$regulator, regulators)
j <- match(pairs$region, regions)
S <- sparseMatrix(i = i, j = j, x = 1,
                  dims = c(length(regulators), length(regions)),
                  dimnames = list(regulators, regions))
S_dense <- as.matrix(S)


melted_data <- reshape2::melt(as.matrix(S_dense))

melted_data$Var2 <- factor(melted_data$Var2, levels = paste0("R",1:length(rglt_list)))
melted_data$Var1 <- melted_data$Var1
p <-ggplot(melted_data, aes(x = Var2, y = Var1, fill = factor(value))) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c("lightgray", "#FFA50099"),labels = c("No binding", "Binding")) +  
  labs(x = "", y = "", fill = "") +
  theme_minimal() +
  theme(axis.title.y = element_text(size = 10),
        axis.text.x = element_text(angle = 90, size=7),
        axis.text = element_text(size = 10),
        legend.position = 'none',
        panel.grid = element_blank(),
        legend.margin = margin(t = 0, b = -10, l = 0, r = 0))
pdf(paste0('Binding_', gene, '.pdf'), height = 1.75, width = 6.96)
print(p)
dev.off()


########
## gid=5
############ 
gid=5
gene <- GeneSet[gid]
regulators <- read.table(paste0('SourceData/Regulators_Track/', gene, '_regulators.tsv'), sep = "\t", header = FALSE)
regulators <- paste0(toupper(substring(regulators$V1, 1, 1)), tolower(substring(regulators$V1, 2)))
path <- paste0('SourceData/TFnameUnique/', gene)
rglt_list <- list.files(path)
rglt_rg_set <- c()
for(rgid in 1:length(rglt_list)){
  rg <- sub(".*_(R[0-9]+)\\..*", "\\1", rglt_list[rgid])
  rgfile <- paste0(path, '/', rglt_list[rgid])
  if (file.info(rgfile)$size > 0) {
    rgltset <- read.table(rgfile, header = FALSE)
    regulators_in_rglt_rg <- rgltset$V1[which(rgltset$V1 %in% regulators)]
    rglt_rg_set <- rbind(rglt_rg_set, cbind(regulators_in_rglt_rg, rep(rg, length(regulators_in_rglt_rg))))
  }
}
rglt_rg_set <- as.data.frame(rglt_rg_set)

regions <- paste0("R", 1:length(rglt_list))
pairs <- distinct(rglt_rg_set, regulator=regulators_in_rglt_rg, region=V2) 
i <- match(pairs$regulator, regulators)
j <- match(pairs$region, regions)
S <- sparseMatrix(i = i, j = j, x = 1,
                  dims = c(length(regulators), length(regions)),
                  dimnames = list(regulators, regions))
S_dense <- as.matrix(S)


melted_data <- reshape2::melt(as.matrix(S_dense))

melted_data$Var2 <- factor(melted_data$Var2, levels = paste0("R",1:length(rglt_list)))
melted_data$Var1 <- melted_data$Var1
p <-ggplot(melted_data, aes(x = Var2, y = Var1, fill = factor(value))) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c("lightgray", "orange"),labels = c("No binding", "Binding")) +  
  labs(x = "", y = "", fill = "") +
  theme_minimal() +
  theme(axis.title.y = element_text(size = 10),
        axis.text.x = element_text(angle = 90, size=7),
        axis.text = element_text(size = 10),
        legend.position = 'none',
        legend.margin = margin(t = 0, b = -10, l = 0, r = 0))
pdf(paste0('Binding_', gene, '.pdf'), height = 1.2, width = 6.79)
print(p)
dev.off()


  
  
  


