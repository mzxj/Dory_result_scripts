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
thrd <- 5
## plot mean or median
DS <- read.table('SourceData/CPD_CID_DiffScore.txt', sep = "\t", header = TRUE)
geneset <- geneset0[which(geneset0$GeneName %in% c("Anxa8", "Myc", "Mycn", "Fhit", "Oit1")), ] #since the DiffRegion Num are different, so plot them one by one.
geneset$label <- geneset$Gene
geneset$label[which(geneset$label == "Oit1")] <- "Fam3d"
for(k in c(2,3)){
  gene <- geneset$Gene[k]
  genename <- geneset$GeneName[k]
  geneid <- geneset$GeneID[k]
  pmt_rg <- convert_to_vector(geneset[["Promoter2kInRegion"]][k])
  S3DisTrace <- read.table(paste0('SourceData/DistanceS3_DistanceDF/S3PromoterDistance_', genename, '_R', pmt_rg, '.tsv'))
  S3Dis <- as.data.frame(apply(S3DisTrace, 2, median, na.rm = TRUE))
  S2DisTrace <- read.table(paste0('SourceData/DistanceS2_DistanceDF/S2PromoterDistance_', genename, '_R', pmt_rg, '.tsv'))
  S2Dis <- as.data.frame(apply(S2DisTrace, 2, median, na.rm = TRUE))
  
  Difmat <- read.table(paste0("SourceData/Matrix/WpCauchypMatrix_State3Gene", geneid, "VSState2Gene", geneid,".tsv"))
  Dif <- t(Difmat[pmt_rg, ])
  
  dataplot <- data.frame(S3Dis=S3Dis[-pmt_rg, ],
                         S2Dis=S2Dis[-pmt_rg, ],
                         diff=Dif[-pmt_rg, ],
                         id=c(seq(1, pmt_rg-1), seq(pmt_rg+1, length(Dif))))
  colorset <- rep("black", dim(dataplot)[1])
  colorset[which(dataplot$diff < - thrd)] <- "#d60b0e"
  ylimin <- floor(min(dataplot$S3Dis, na.rm = TRUE) * 10) / 10
  dataplot[which(dataplot$diff < - thrd),]
  
  p <- ggplot(dataplot)+
    geom_segment(aes(x = diff, y = S2Dis, xend = diff, yend = S3Dis),
                 arrow = arrow(length = unit(0.05, "inches")), 
                 size = 0.3, color = colorset)+
    theme_bw()+
    xlim(min(dataplot$diff) - 1 , max(dataplot$diff) +0.5)+
    ylim(ylimin, NA)+
    theme(panel.grid.major = element_blank(),  # Remove major gridlines
          panel.grid.minor = element_blank(),
          axis.title.y = element_text(hjust = 0.1))+
    labs(x = "DiffScore", y= 'Median NML spatial distance', title = geneset$label[k])+
    annotate("rect", xmin=-Inf, xmax = 0, ymin = -Inf, ymax = Inf, alpha = 0.1, fill = "#d60b0e") +  
    annotate("rect", xmin=0, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.1, fill = "#0f70bf") +  
    geom_vline(xintercept = -5, linetype = "dashed", color = "#d60b0e", alpha=0.8) +
    annotate("text", x = -5.8, y = ylimin + 0.03, label = "Differential", alpha=0.8, color = "#d60b0e", size = 2.5)
  
  pdf(paste0(geneset$label[k], '_R', pmt_rg, '.pdf'), width = 3.5, height = 2.5)
  print(p)
  dev.off()
}


for(k in c(1,4,5)){
  gene <- geneset$Gene[k]
  genename <- geneset$GeneName[k]
  geneid <- geneset$GeneID[k]
  pmt_rg <- convert_to_vector(geneset[["Promoter2kInRegion"]][k])
  S3DisTrace <- read.table(paste0('SourceData/DistanceS3_DistanceDF/S3PromoterDistance_', genename, '_R', pmt_rg, '.tsv'))
  S3Dis <- as.data.frame(apply(S3DisTrace, 2, median, na.rm = TRUE))
  S2DisTrace <- read.table(paste0('SourceData/DistanceS2_DistanceDF/S2PromoterDistance_', genename, '_R', pmt_rg, '.tsv'))
  S2Dis <- as.data.frame(apply(S2DisTrace, 2, median, na.rm = TRUE))
  
  Difmat <- read.table(paste0("SourceData/Matrix/WpCauchypMatrix_State3Gene", geneid, "VSState2Gene", geneid,".tsv"))
  Dif <- t(Difmat[pmt_rg, ])
  
  dataplot <- data.frame(S3Dis=S3Dis[-pmt_rg, ],
                         S2Dis=S2Dis[-pmt_rg, ],
                         diff=Dif[-pmt_rg, ],
                         id=c(seq(1, pmt_rg-1), seq(pmt_rg+1, length(Dif))))
  colorset <- rep("black", dim(dataplot)[1])
  colorset[which(dataplot$diff < - thrd)] <- "#d60b0e"
  ylimin <- floor(min(dataplot$S3Dis, na.rm = TRUE) * 10) / 10
  
  dataplot[which(dataplot$diff < - thrd),]
  
  p <- ggplot(dataplot)+
    geom_segment(aes(x = diff, y = S2Dis, xend = diff, yend = S3Dis),
                 arrow = arrow(length = unit(0.05, "inches")), # Arrow customization
                 size = 0.3, color = colorset)+
    theme_bw()+
    xlim(min(dataplot$diff) - 1 , max(dataplot$diff) +0.5)+
    ylim(ylimin, NA)+
    theme(panel.grid.major = element_blank(),  # Remove major gridlines
          panel.grid.minor = element_blank(),
          axis.title.y = element_text(hjust = 0.1))+
    labs(x = "DiffScore", y= 'Median NML spatial distance', title = genename)+
    annotate("rect", xmin=-Inf, xmax = 0, ymin = -Inf, ymax = Inf, alpha = 0.1, fill = "#d60b0e") +  # Purple background on the left
    annotate("rect", xmin=0, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.1, fill = "#0f70bf") +  # Purple background on the left
    geom_vline(xintercept = -5, linetype = "dashed", color = "#d60b0e", alpha=0.8) +
    annotate("text", x = -6.2, y = ylimin + 0.03, label = "Differential", alpha=0.8, color = "#d60b0e", size = 2.5)
  
  pdf(paste0( geneset$label[k], '_R', pmt_rg, '.pdf'), width = 3.5, height = 2.5)
  print(p)
  dev.off()
}
