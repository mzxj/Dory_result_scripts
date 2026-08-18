###############################
##!! note: there is a small issue for the gene annotation in this library for 'Mir3084−2'
################################

library("TxDb.Mmusculus.UCSC.mm10.knownGene")
library("org.Mm.eg.db")
library("plotgardener")
library(ggplot2)
gloci<-read.csv("SourceData/region.csv",header=F,as.is=T)
colnames(gloci) <- c("chr", "start", "end", "id")
# Input barplot
input_data<- read.table(paste0("SourceData/WpCauchypMatrix_5kbHepatocyteVS5kbNonHepatocyte.tsv"))
pmt_rg <- 19
data_use <- data.frame(gloci[,1:4], score=t(input_data[pmt_rg, ]))
colnames(data_use) <- c("chrom", "start", "end", "regionID", "score")
data_range <- pgParams(range = c(min(data_use$score, na.rm = TRUE), max(data_use$score, na.rm = TRUE)),
                       assembly = "mm10")
thrd <- -log10(0.05)
## region arches
data_arc <- data.frame(chrom1 = character(), start1 = numeric(), end1 = numeric(),
                       chrom2 = character(), start2 = numeric(), end2 = numeric())
slcrg <- which(data_use$score < -thrd)
for(rg in slcrg){
  data_arc <- rbind(data_arc, data.frame(
    chrom1 = data_use$chrom[pmt_rg], start1 = (data_use$start[pmt_rg] + data_use$end[pmt_rg])/2, end1 = (data_use$start[pmt_rg] + data_use$end[pmt_rg])/2,
    chrom2 = data_use$chrom[rg], start2 = (data_use$start[rg] +data_use$end[rg])/2, end2 = (data_use$start[rg] +data_use$end[rg])/2
  ))
}
data_arc$length <- abs(data_arc$start2 - data_arc$start1)/1000
data_arc$h <- data_arc$length/max(data_arc$length)
data_scale <- data.frame(chrom1 = data_use$chrom[1], start1 = data_use$start[1]+100, end1 = data_use$start[1] + 200,
                         chrom2 = data_use$chrom[1], start2 = data_use$end[1] , end2 = data_use$end[1] +100 )
## gene region
paramsbig <- pgParams(
  chrom = unique(data_use$chrom),
  chromstart = data_use$start[1] , chromend = data_use$end[dim(gloci)[1]] + 10000, ## adding two region legth at the end
  assembly = "mm10", width = 3, x = 0.15
)
pdf(paste0('Promoter2k_Scd2.pdf'), width = 3.5, height = 1.6)
pageCreate(width = 3.5, height = 1.6, default.units = "inches", showGuides = FALSE)
## Plot data bar signal
plotRect(x = seq(0.15+3/22/2, 0.15+3/22/2 + 3/22*4, 3/22), y = 0.3,height = 0.1, width = 3/22, 
         default.units = "inches", fill = c("#fbe6e5", "#fbe6e5",rep("lightgray", 3) ))
plotRect(x = seq(0.15+3/22/2 + 3/22*6, 3.15-3/22/2 - 3/22*2, 3/22), y = 0.3,height = 0.1, width = 3/22, 
         default.units = "inches", fill = c(rep("lightgray", 9), rep("#fbe6e5", 3), "lightgray", "darkgray" ))

plotText(label = c(1, 3),  fontcolor = "black", fontsize = 10,
         x = c(0.15+3/22/2 , 0.15+3/22/2 + 3/22*2), y = 0.45) # 0.02 is to move label a little bit right not at the rectangle left frame.
plotText(label = c(seq(6, 19, 3)),  fontcolor = "black", fontsize = 10,
         x = c(seq(0.15+3/22/2 + 3/22*6 , 3.15-3/22/2 -3/22*2, 3/22*3)), y = 0.45)
plotPairsArches(
  data = data_arc[3:5, ], params = paramsbig,
  linecolor = "black", lwd = 0.5, fill = 'black',
  flip = FALSE, curvature = 5,
  archHeight = "h", alpha = 1,
  x = 0.15, y = 0.25, height = 0.08,
  just = c("left", "bottom"),
  default.units = "inches")
plotPairsArches(
  data = data_arc[1:2, ], params = paramsbig,
  linecolor = "black", lwd = 0.5, fill = 'black',
  flip = TRUE, curvature = 50,
  archHeight = "h", alpha = 1,
  x = 0.15, y = 0.35, height = 0.4,
  just = c("left", "top"),
  default.units = "inches")
plotPairs(
  data = data_scale, params = paramsbig, fill = "black", lwd=1,
  y = 0.63, height = 0.1, just = c("left", "top"), default.units = "inches"
)
plotText(
  label = paste0(round((data_use$end[1] - data_use$start[1])/1000),"kb"),  x = 0.15+3/22*1 + 3/22, y = 0.67, fontsize = 10
)
genesPlot <- plotGenes(
  params = paramsbig,
  geneHighlights = data.frame(
    #"gene" = c(genename), c(geneset$label[k]),
    "gene" = "Scd2",
    "color" = c("#018f52")),
  geneBackground = "#3d5c6f",
  y = 0.8, height = 0.5,
  just = c("left", "top"), default.units = "inches")
## Annotate genome label
annoGenomeLabel(params = paramsbig,
                plot = genesPlot,  y = 1.35, scale = "bp",
                just = c("left", "top"))
dev.off()




#######
## bar plot
################
data_use$regionID <- factor(data_use$regionID, levels = seq(1, 19, 1))
data_use1 <- data.frame(regionID = c(1:20),
                        score = c(data_use$score[1:5], NA, data_use$score[6:19]))
data_use1$regionID <- factor(data_use1$regionID, levels = seq(1, 20, 1))
pdf(paste0('Promoter2k_Scd2_bar.pdf'), width = 3.3, height = 1)
## ggplot bar plot
p <- ggplot(data_use1, aes(x = regionID, y = score, fill = score > 0)) +
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
    axis.title.y = element_text(size=10),
    legend.position = "none", 
  )
print(p)
dev.off()
