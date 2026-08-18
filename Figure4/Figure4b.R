library(ggplot2)
library(ggrepel)
pmt_rg <- 19
DisTrace <- read.table(paste0('SourceData/PromoterRow/NonHepatocyte_r', pmt_rg,'.tsv'))
NonHepDis <- as.data.frame(apply(DisTrace, 2, median, na.rm = TRUE))

DisTrace <- read.table(paste0('SourceData/PromoterRow/Hepatocyte_r', pmt_rg,'.tsv'))
HepDis <- as.data.frame(apply(DisTrace, 2, median, na.rm = TRUE))

Difmat <- read.table(paste0("SourceData/WpCauchypMatrix_5kbHepatocyteVS5kbNonHepatocyte.tsv"))
Dif <- t(Difmat[pmt_rg, ])

dataplot <- data.frame(NonHepDis=NonHepDis[1:18, ],
                       HepDis=HepDis[1:18, ],
                       diff=Dif[1:18, ],
                       id=c(1, 2, rep("", 12), 15, 16, 17, ""))
ylimin <- 0.1


pdf('Medianthrd0.05.pdf', width = 3.2, height = 2.2)
ggplot(dataplot, aes(x=diff, y=distance))+
  geom_segment(aes(x = diff, y = NonHepDis, xend = diff, yend = HepDis),
    arrow = arrow(length = unit(0.05, "inches")), # Arrow customization
    size = 0.3, color = c("#d60b0e", "#d60b0e", rep("black", 12), rep("#d60b0e", 3), "black"    ) 
  )+
  theme_bw()+
  ylim(ylimin, NA)+
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  labs(x = "DiffScore", y= 'Median spatial distance')+
  annotate("rect", xmin=-Inf, xmax = 0, ymin = -Inf, ymax = Inf, alpha = 0.1, fill = "#d60b0e") +  
  annotate("rect", xmin=0, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.1, fill = "#0f70bf") +  
  geom_vline(xintercept = log10(0.05), linetype = "dashed", color = "#d60b0e", alpha=0.8) +
  annotate("text", x = log10(0.05)-1, y = ylimin+0.01, label = "Differential", alpha=0.8, color = "#d60b0e", size = 4)#+
dev.off()

