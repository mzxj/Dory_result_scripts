library(ggplot2)
library(grid)
library(gridExtra)
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
  TFset <- rbind(TFset, cbind(colnames(rexp), rep(paste0("Region", rg), dim(rexp)[2])))
}



TFuniq <- unique(TFset$V1)
TFuniqExp <- genexp[, which(colnames(genexp) %in% TFuniq)]

plt <- list()
for(i in 1:length(TFuniq)){
  data <- data.frame(Celltype = rownames(TFuniqExp), Expression = TFuniqExp[,i])
  data$Celltype <- factor(data$Celltype, levels = c("Hepatocyte", "Erythroblast", "Macrophage", "Megakaryocyte", "Proerythroblast"))
  plt[[i]] <- ggplot(data, aes(x=Celltype, y=Expression, fill = Celltype))+
    geom_bar(stat = "identity")+
    scale_fill_manual(values=c("#39593f", "#617c52",  "#889f64", "#b0c277", "#d7e589"))+
    theme_bw()+
    theme(legend.position = "none",
          panel.grid.major = element_blank(),  # Remove major gridlines
          panel.grid.minor = element_blank(),
          axis.title.y = element_text(size = 8),
          axis.text.x = element_text(angle = 60,vjust = 1, hjust = 1, size=8),
          axis.text = element_text(size = 10),
          plot.margin = margin(t = 1, r = 5, b = 0, l = 5))+ 
    labs(x="")+
    scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
    annotate("text", x = length(data$Celltype), y = max(data$Expression), label = (TFuniq[i]), hjust = 1, vjust = 1, size=3)
}


pdf(paste0('Exp1.pdf'),  height = 2, width = 7)
grid.arrange(grobs=plt[1:5],ncol = 5, nrow = 1)
dev.off()
pdf(paste0('Exp2.pdf'),  height = 2, width = 7)
grid.arrange(grobs=plt[6:10],ncol = 5, nrow = 1)
dev.off()
pdf(paste0('Exp3.pdf'),  height = 2, width = 7)
grid.arrange(grobs=plt[11:15],ncol = 5, nrow = 1)
dev.off()
pdf(paste0('Exp4.pdf'),  height = 2, width = 7)
grid.arrange(grobs=plt[16:20],ncol = 5, nrow = 1)
dev.off()
pdf(paste0('Exp5.pdf'),  height = 2, width = 7)
grid.arrange(grobs=plt[21:25],ncol = 5, nrow = 1)
dev.off()
pdf(paste0('Exp6.pdf'),  height = 2, width = 7)
grid.arrange(grobs=plt[26],ncol = 5, nrow = 1)
dev.off()


