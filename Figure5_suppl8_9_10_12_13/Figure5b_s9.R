if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("plotgardenerData")
library("TxDb.Mmusculus.UCSC.mm10.knownGene")
library("org.Mm.eg.db")
library("plotgardener")
library(ggplot2)

convert_to_vector <- function(region_string) {
  # Remove the "region:" prefix and split by comma
  region_numbers <- unlist(strsplit(sub("region:", "", region_string), ","))
  # Convert the split strings to numeric and create a vector
  region_vector <- as.numeric(region_numbers)
  return(region_vector)
}


thrd <- 5
geneset0 <- read.table('SourceData/GeneTSSandPromoterRegion_update.tsv', header = TRUE)
geneset <- geneset0[which(geneset0$GeneName %in% c("Anxa8", "Fhit", "Oit1", "Mycn", "Myc")), ]
# oit1 is fam3d
geneset$label <- geneset$Gene
geneset$label[which(geneset$label == "Oit1")] <- "Fam3d"
for(k in 1:dim(geneset)[1]){
  gene <- geneset$Gene[k]
  genename <- geneset$GeneName[k]
  geneid <- geneset$GeneID[k]
  gloci<-read.table(paste0("SourceData/GeneRegions/Regions_", gene, ".txt"),header=F,as.is=T)
  gloci$V1 <- paste0("chr", gloci$V1)
  
  ## Input barplot
  input_data_NML <- read.table(paste0("SourceData/Matrix/WpCauchypMatrix_State3Gene", geneid, "VSState2Gene", geneid,".tsv"))
  pmt_rg <- convert_to_vector(geneset[["Promoter2kInRegion"]][k])
  data_NML <- data.frame(gloci[,1:4], score=t(input_data_NML[pmt_rg, ]))
  colnames(data_NML) <- c("chrom", "start", "end", "regionID", "score")
  data_NML$end <- data_NML$end -1
  data_range <- pgParams(range = c(min(data_NML$score, na.rm = TRUE), max(data_NML$score, na.rm = TRUE)),
                         assembly = "mm10")
  data_arc <- data.frame(chrom1 = character(), start1 = numeric(), end1 = numeric(),
                         chrom2 = character(), start2 = numeric(), end2 = numeric())
  slcrg <- which(data_NML$score < -thrd)
  for(rg in slcrg){
    data_arc <- rbind(data_arc, data.frame(
      chrom1 = data_NML$chrom[pmt_rg], start1 = (data_NML$start[pmt_rg] + data_NML$end[pmt_rg])/2, end1 = (data_NML$start[pmt_rg] + data_NML$end[pmt_rg])/2,
      chrom2 = data_NML$chrom[rg], start2 = (data_NML$start[rg] +data_NML$end[rg])/2, end2 = (data_NML$start[rg] +data_NML$end[rg])/2
    ))
  }
  data_arc$length <- (data_arc$start2 - data_arc$start1)/1000
  data_arc$h <- data_arc$length/max(data_arc$length)
  data_scale <- data.frame(chrom1 = data_NML$chrom[1], start1 = data_NML$start[2], end1 = data_NML$start[2] + 500,
                           chrom2 = data_NML$chrom[1], start2 = data_NML$end[2] - 500, end2 = data_NML$end[2] )
  ## gene region
  paramsbig <- pgParams(
    chrom = unique(data_NML$chrom),
    chromstart = data_NML$start[1] , chromend = data_NML$end[dim(gloci)[1]] ,
    assembly = "mm10", width = 6, x = 0.25
  )

  pdf(paste0('Promoter2k_', genename, '.pdf'), width = 6.5, height = 1.23)
  pageCreate(width = 6.5, height = 1.23, default.units = "inches", showGuides = FALSE)
  ## Plot data bar signal
  fillrect <- rep("lightgray", dim(gloci)[1])
  fillrect[slcrg] <- "#fbe6e5" 
  fillrect[pmt_rg] <- "darkgray"
  plotRect(x = seq(0.25+6/40/2, 6.25-6/40/2, 6/40), y = 0.1,height = 0.1, width = 6/40, 
            default.units = "inches", fill = fillrect)
  plotText(label = c(1, seq(5, 40, 5)),  fontcolor = "black", fontsize = 10,
           x = c(0.25+6/40/2 , seq(0.25+6/40/2 + 6/40*4 , 6.25-6/40/2 , 6/40*5)), y = 0.25) # 0.02 is to move label a little bit right not at the rectangle left frame.
  plotPairsArches(
    data = data_arc, params = paramsbig,
    linecolor = "black", lwd = 0.5, fill = 'black',
    flip = TRUE, curvature = 5,
    archHeight = "h", alpha = 1,
    x = 0.25, y = 0.15, height = 0.4,
    just = c("left", "top"),
    default.units = "inches")
  plotPairs(
    data = data_scale, params = paramsbig, fill = "black", lwd=1,
    y = 0.4, height = 0.1, just = c("left", "top"), default.units = "inches"
  )
  plotText(
    label = paste0(round((data_NML$end[1] - data_NML$start[1])/1000),"kb"),  x = 0.25+6/40*3 +6/40/1.5, y = 0.45, fontsize = 10
  )
  genesPlot <- plotGenes(
    params = paramsbig,
    geneHighlights = data.frame(
      "gene" = c(geneset$label[k]),
      "color" = c("#018f52")),
    geneBackground = "#3d5c6f",
    y = 0.5, height = 0.5,
    just = c("left", "top"), default.units = "inches")
  annoGenomeLabel(params = paramsbig,
    plot = genesPlot,  y = 1.05, scale = "bp",
    just = c("left", "top"))
  dev.off()
  
  data_NML$regionID <- factor(data_NML$regionID, levels = paste0("Region", seq(1, dim(data_NML)[1], 1)))
  pdf(paste0('Promoter2k_', genename, '_bar.pdf'), width = 6.58, height = 1.07)
  ## ggplot bar plot
  p <- ggplot(data_NML, aes(x = regionID, y = score, fill = score > 0)) +
    geom_bar(stat = "identity", color = "black") +
    scale_fill_manual(values = c("#d60b0e", "#0f70bf")) +
    geom_hline(yintercept = -thrd, color = "darkgray", linetype = "dashed") +  
    labs(x = NULL,y = "DiffScore",fill = NULL) +
    theme_classic()+
    theme(
      panel.grid.major = element_blank(),     
      panel.grid.minor = element_blank(),     
      panel.background = element_blank(),     
      axis.line.x = element_blank(),          
      axis.text.x = element_blank(),          
      axis.ticks.x = element_blank(),         
      axis.title.x = element_blank(),
      axis.text.y = element_text(size=10),
      legend.position = "none", 
    )
  print(p)
  dev.off()
}



