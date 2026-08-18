celltype <- data.frame(ID=c(paste0("ct", seq(2,7))),
                       name = c("Endothelial_cell", "Hepatocyte", "Macrophage", "Megakaryocyte", "Erythroblast", "Proerythroblast"))
celltypename <- c()
i <- 2
while(i < 7){
  j <- i + 1
  while(j <= 7){
    #celltypename <- append(celltypename, paste0(celltype$name[i-1], ' VS ', celltype$name[j-1]))
    celltypename <- append(celltypename, paste0(celltype$name[j-1], ' to ', celltype$name[i-1]))
    j <- j + 1
  }
  i <- i + 1
}

library(grid)
library(gridExtra)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
PglPath <- 'SourceData/DiffScore_Matrix/'
ABchangePath<- 'SourceData/ABchange_Results/'
ABfiles <- list.files(ABchangePath, '.tsv')
ctset <- gsub(".tsv", "", ABfiles)


thrd <- 3
plt <- list()
k <- 1
for(ct in ctset){
#ct<- "ct4VSct7"
#celltypename <- "Proerythroblast to Macrophage"
  dataAB <- read.table(paste0(ABchangePath, ct, '.tsv'), header = TRUE)
  data <- read.table(paste0(PglPath, 'WpCauchypMatrix_', ct, '.tsv'))
  datapfinal <- data[lower.tri(data, diag = FALSE)]
  dataplot <- data.frame(ABchange = dataAB$ABchange, 
                         Diffscore = datapfinal)
  colorblue <- densCols(dataAB$ABchange, datapfinal)
  model <- lm(Diffscore ~ ABchange, data = dataplot)
  dataplot$predicted_y <- predict(model)
  dataplot$residuals <- residuals(model)
  residual_sd <- sd(dataplot$residuals) # Determine the standard deviation of residuals
  
  results <- data.frame(index = numeric(), coef_change = numeric(), r2_change = numeric())
  for (i in 1:dim(dataplot)[1]) {
    temp_data <- dataplot[-i, ]
    temp_model <- lm(Diffscore ~ ABchange, data = temp_data)
    coef_change <- as.numeric(abs(coef(temp_model)) - abs(coef(model) ))[2]
    r2_change <- summary(temp_model)$r.squared - summary(model)$r.squared 
    results <- rbind(results, data.frame(index = i, coef_change = coef_change, r2_change = r2_change))
  }
  results1 <- cbind(results, dataplot$residuals)

  colorred <- which(abs(dataplot$residuals) > thrd * residual_sd & results1$r2_change > 0)
  colorblue[colorred] <- "#faa943" #"#E41A1C"  # "#874f8d"  # "#EE6A33" #"#E41A1C"
  
  dataplot <- data.frame(ABchange = dataAB$ABchange, 
                         Diffscore = datapfinal,
                         PointDense = colorblue)
  cor_result <- cor.test(dataplot$ABchange, dataplot$Diffscore)
  r_value <- round(cor_result$estimate, 3)
  p_value <- signif(cor_result$p.value, 3)
  

  plt[[k]] <- ggplot(dataplot, aes(x=ABchange, y=Diffscore))+
    geom_point(color=I(colorblue),  size=1)+
    theme_bw()+
    theme(legend.position = "none",
          axis.title = element_text(size = 7),
          axis.text = element_text(size = 7),
          plot.title = element_text(size = 7))+ # top right bottom left
    geom_smooth(method = "lm", se = TRUE, color = "#377EB8", fill="darkgray", size=0.5) +
    stat_cor(method = "pearson", label.x.npc = 0.3, label.y.npc = 0.95, size= 2)+
    labs(title = paste0(celltypename[k]))
  pdf(paste0('changecol_', ct, '.pdf'), height = 1.9, width = 1.9)    
  print(plt[[k]])
  dev.off()
  
  k <- k + 1
}

