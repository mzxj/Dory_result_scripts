# gene/TF expression
library(dplyr)
library(ggplot2)
library(ggrepel)
convert_to_vector <- function(region_string) {
  # Remove the "region:" prefix and split by comma
  region_numbers <- unlist(strsplit(sub("region:", "", region_string), ","))
  # Convert the split strings to numeric and create a vector
  region_vector <- as.numeric(region_numbers)
  return(region_vector)
}

geneset0 <- read.table('SourceData/GeneTSSandPromoterRegion_update.tsv', header = TRUE)
geneset1 <- geneset0[which(geneset0$Gene %in%c("Anxa8", "Fhit", "Oit1", "Mycn", "Myc")), ]
geneset1$label <- geneset1$Gene
geneset1$label[which(geneset1$label == "Oit1")] <- "Fam3d"


RNA <- read.csv('SourceData/Download_Supplementary5_bulkRNAseq.csv')
for(g in 1:dim(geneset1)[1]){
  pmt_rg <- convert_to_vector(geneset1[["Promoter2kInRegion"]][g])
  gene <- geneset1$Gene[g]
  tflist <- read.table(paste0("SourceData/DiffRegion_p1e.05_DiffTF_log2(2)/", gene, "_R", pmt_rg, "_DR_DTF.txt"))
  ## extract rna
  datause <- RNA[which(RNA$Gene %in% tflist$V1), ]
  colorset <- rep("gray", dim(RNA)[1])
  colorset[which(RNA$Gene %in% tflist$V1)] <- "orange"
  dataplot <- data.frame(gene = RNA$Gene,
                         logFC = RNA$log2FoldChange..LUAD.AdenomaG.,
                         negLogPadj = -log10(RNA$padj),
                         colorset = colorset)
  p <- ggplot(dataplot, aes(x = logFC, y= negLogPadj))+
    geom_point(data=subset(dataplot, colorset == "gray"), color="gray", alpha=0.8, size=0.7)+
    geom_point(data=subset(dataplot, colorset == "orange"), color="orange", size=0.7)+
    geom_vline(xintercept = log2(2), linetype = "dashed", color="#025259", linewidth = 0.2)+
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color="#025259", linewidth = 0.2)+
    geom_text_repel(data=subset(dataplot, colorset == "orange"), aes(label = gene), color="darkorange", vjust = -0.5, size=4)+
    theme_bw()+
    xlim(-5, 15)+
    labs(x="log2FoldChange", y="-log10(p.adj)", title = geneset1$label[g])+
    theme(panel.grid.major = element_blank(),  
          panel.grid.minor = element_blank())+
    annotate("rect", xmin=log2(2), xmax = Inf, ymin = -log10(0.05), ymax = Inf, alpha = 0.1, fill = "orange")   # #f29325 background on the left
  pdf(paste0('RegulatorExp_', geneset1$label[g], '.pdf'), width = 3.5, height = 4)
  print(p)
  dev.off()  
}



