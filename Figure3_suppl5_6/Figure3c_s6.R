#set.seed(42)

############################################
## in our data, permute the ABchange
###########################################
celltype <- data.frame(ID=c(paste0("ct", seq(2,7))),
                       name = c("Endothelial_cell", "Hepatocyte", "Macrophage", "Megakaryocyte", "Erythroblast", "Proerythroblast"))
celltypename <- c()
i <- 2
while(i < 7){
  j <- i + 1
  while(j <= 7){
    celltypename <- append(celltypename, paste0(celltype$name[j-1], ' to ', celltype$name[i-1]))
    j <- j + 1
  }
  i <- i + 1
}

PglPath <- 'SourceData/DiffScore_Matrix/'
ABchangePath<- 'SourceData/ABchange_Results/'
ABfiles <- list.files(ABchangePath, '.tsv')
ctset <- gsub(".tsv", "", ABfiles)

num_permutations <- 1000
library(ggplot2)
library(grid)
library(gridExtra)
plt <- list()
k <- 1
for(ct in ctset){
  dataAB <- read.table(paste0(ABchangePath, ct, '.tsv'), header = TRUE)
  data <- read.table(paste0(PglPath, 'WpCauchypMatrix_', ct, '.tsv'))
  datapfinal <- data[lower.tri(data, diag = FALSE)]
  actual_cor <- cor(dataAB$ABchange, datapfinal)
  perm_cor_values <- replicate(num_permutations, cor(datapfinal, sample(dataAB$ABchange, replace = FALSE)))

  df_perm <- data.frame(perm_cor_values)
  
  plt[[k]] <- ggplot(df_perm, aes(x = perm_cor_values)) +
    geom_histogram(binwidth = 0.001, fill = "darkgray", color = "#4d4d4d") +  
    geom_vline(xintercept = actual_cor, color = "#d62728", linetype = "dashed", size = 0.5) +  
    labs(title = paste0(celltypename[k], "\nPermutation Test for Correlation"),
         x = "Correlation Coefficient",
         y = "Frequency") +
    xlim(-1, 1) +  # Set x-axis limits
    theme_minimal()+
    theme(axis.title = element_text(size = 8),
          axis.text = element_text(size = 8),
          plot.title = element_text(size = 8))
  pdf(paste0('permuteABchange_', ct,'.pdf'), height = 1.5, width =2.5)
  print(plt[[k]])
  dev.off()
  k <- k + 1
}






