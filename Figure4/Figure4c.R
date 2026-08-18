##################################################################################
# !!! note: replace the HNF4A ChIP-seq track with the genome browser track
##################################################################################

library(GenomicRanges)
library(plotgardener)
library(rtracklayer)

cCRE <- data.table::fread('SourceData/mm10-chr19-44193676-44293675-strain_B6NTacB6NCrl_Lap_liver_tissue_embryo_14.5_days_ENCDO469XMJ.bed', header = TRUE)
cCRE_use <- cCRE[, c(1,2,3,5)]
colnames(cCRE_use) <- c("chr", "start", "end", "classification")
data <- cCRE_use[which(cCRE_use$classification %in% c("dELS", "pELS")), ]
bed_ELS <- GRanges(seqnames = data$chr,
                   ranges = IRanges(start = data$start, end = data$end),
                   name = data$classification)
bigwig_file <- "SourceData/GSM2406340_HNF4A_ChIP-seq_ZT10_HNF4A_floxed_rep2_chr19_mm10.bw"
bw_data <- import(bigwig_file)
signal_df <- as.data.frame(bw_data)


pdf(paste0('enhancer_HNF4A_update.pdf'), width = 6.8, height = 2)
pageCreate(width = 6, height = 1.6, default.units = "inches", showGuides = FALSE)
plotText(label = 'Hnf4a\nChIP-seq',  fontcolor = "#8a4b21", fontsize = 10,just = c("right", "top"),
         x = 0.22, y = 0.5)
plotSignal(
  data = bw_data, binCap = FALSE, chrom = "chr19", chromstart = 44193676, chromend = 44293675,
  score = score(bw_data), assembly = "mm10", x = 0.3, y = 0.1, orientation = "h", width = 5.5, height = 0.8,
  linecolor = "#8a4b21", just = c("left", "top")
)
plotText(label = 'ELS',  fontcolor = "#3d4a77", fontsize = 10,just = c("right", "top"), #"#eb8d56" "#d6af6e"
         x = 0.22, y = 0.98)
plotRanges(
  data = bed_ELS,
  chrom = "chr19",   # Modify chromosome based on your BED file
  chromstart = 44193676, chromend = 44293675,
  x = 0.3, y = 0.95,   # Position on plot
  width = 5.5, height = 0.15, collapse = TRUE,
  fill = "#3d4a77", alpha = 0.9
)

plotRect(x = seq(0.3+5.5/20/2, 0.3+5.5/20/2 + 5.5/20*4, 5.5/20), y = 1.2,height = 0.1, width = 5.5/20, 
         default.units = "inches", fill = c("#fbe6e5", "#fbe6e5",rep("lightgray", 3) )) ##7bc4c5  #c5d1d2
plotRect(x = seq(0.3+5.5/20/2 + 5.5/20*6, 5.8-5.5/20/2, 5.5/20), y = 1.2,height = 0.1, width = 5.5/20, 
         default.units = "inches", fill = c(rep("lightgray", 9), rep("#fbe6e5", 3), "lightgray", "darkgray" ))
plotText(label = c(1, 3),  fontcolor = "black", fontsize = 10,
         x = c(0.3+5.5/20/2 , 0.3+5.5/20/2 + 5.5/20*2), y = 1.32) # 0.02 is to move label a little bit right not at the rectangle left frame.
plotText(label = c(seq(6, 19, 3)),  fontcolor = "black", fontsize = 10,
         x = c(seq(0.3+5.5/20/2 + 5.5/20*6 , 5.8-5.5/20/2, 5.5/20*3)), y = 1.32)

plotGenomeLabel(
  chrom = "chr19",
  chromstart = 44193676, chromend = 44293675,
  assembly = "mm10",
  x = 0.3, y = 1.45, length = 5.5,
  default.units = "inches"
)
dev.off()



