TFlist <- c("Ctnnb1", "Nr1h3", "Nrf1", "Srebf1", "Smad4", "Hnf4a")

TFexp <- read.table('SourceData/TFexp.tsv', sep = "\t", header = TRUE)
marker <- TFexp[, which(colnames(TFexp) %in% TFlist)]

library(ggplot2)
for(i in 1:dim(marker)[2]){
  data <- data.frame(Celltype = rownames(marker), Expression = marker[,i])
  data$Celltype <- factor(data$Celltype, levels = c("Hepatocyte", "Erythroblast", "Macrophage", "Megakaryocyte", "Proerythroblast"))
  plt <- ggplot(data, aes(x=Celltype, y=Expression, fill = Celltype))+
    geom_bar(stat = "identity")+
    scale_fill_manual(values=c("#39593f", "#617c52",  "#889f64", "#b0c277", "#d7e589"))+
    theme_bw()+
    theme(legend.position = "none",
          panel.grid.major = element_blank(),  # Remove major gridlines
          panel.grid.minor = element_blank(),
          axis.title.y = element_text(size = 10),
          axis.text.x = element_text(angle = 60,vjust = 1, hjust = 1, size=10),
          axis.text = element_text(size = 10),
          plot.margin = margin(t = 1, r = 1, b = 0, l = 1))+
    labs(x="")+
    scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
    expand_limits(y = max(data$Expression) + 0.4) +
    annotate("text", x = length(data$Celltype), y = max(data$Expression)+0.4, label = colnames(marker)[i], hjust = 1, vjust = 1, size=3)
  pdf(paste0('Exp_',colnames(marker)[i], '.pdf'),  height = 2.4, width = 1.2)
  print(plt)
  dev.off()
}
