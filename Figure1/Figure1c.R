##################
##########bar plot
####################
dt <- read.table('SourceData/Rep4_WpCauchyPGreatLess_rep4ct3VSrep4ct4.tsv', sep = "\t", header = TRUE)
for(i in 1:dim(dt)[1]){
  if(!is.na(dt$pless[i])){
    a <- c(-log10(dt$pgreat[i]), log10(dt$pless[i]))
    dt$pf[i] <- a[which.max(abs(a))]
    dt$pfs[i] <- which.max(abs(a))
    dt$pfi[i] <- a[which.min(abs(a))]
    dt$pfis[i] <- which.min(abs(a))
  }else{
    dt$pf[i] <- NA
    dt$pfs[i] <- NA
    dt$pfi[i] <- NA
    dt$pfis[i] <- NA
  }
}
dt$pfs <- ifelse(dt$pfs == 1, "great", "less")
dt$pfis <- ifelse(dt$pfis == 1, "great", "less")
dt1 <- data.frame(p_final=dt$pf, p_final_sign=dt$pfs, p_other=dt$pfi, p_other_sign=dt$pfis)
dt2 <- dt1[order(dt1$p_final, decreasing = TRUE),]
dt2$ID <- paste0(seq(1, dim(dt)[1]))

data_o <- reshape2::melt(dt2, id=c("ID", "p_final_sign", "p_other_sign"))
data <- data.frame(ID=data_o$ID, variable=data_o$variable, value=data_o$value, value_sign=ifelse(data_o$variable=="p_final", data_o$p_final_sign, data_o$p_other_sign))
data$ID <- factor(data$ID, levels = paste0( seq(dim(dt)[1], 1)))

dt3 <- dt2[, c(1,5)]
dataline <- dt3
dataline$ID <- factor(dataline$ID, levels = paste0( seq(dim(dt)[1], 1)))
n005g <- length(which(dt$pgreat < 10^(-4)))
n005l <- length(which(dt$pless < 10^(-4)))
library(ggplot2)
pdf(paste0('BarPlot_addLine_revisedAnnot1.pdf'), width = 2, height = 2.56)
ggplot(data,aes(y=ID))+ 
  geom_bar(aes(x=value, color = value), stat = "identity")+
  scale_color_gradient2(low="#d60b0e", mid="white", high="#0f70bf", midpoint = 0, na.value = "grey", limits=c(-abs(max(max(data$value, na.rm = TRUE), min(data$value, na.rm = TRUE))), abs(max(max(data$value, na.rm = TRUE), min(data$value, na.rm = TRUE)))) )+
  geom_vline(xintercept = -log10(10^(-4)), col="grey", linetype = "dashed")+
  geom_vline(xintercept = log10(10^(-4)), col="grey", linetype = "dashed")+
  geom_line(data=dataline, aes(x=p_final,y=ID, group = 1))+
  theme(panel.background = element_rect(fill = "#f2f2f2", 
                                        colour = "#f2f2f2",
                                        size = 0.1, linetype = "solid"),
        panel.grid.major = element_line(size = 0.1, linetype = 'solid',
                                        colour = "#f2f2f2"), 
        panel.grid.minor = element_line(size = 0.1, linetype = 'solid',
                                        colour = "#f2f2f2"),
        legend.position = "none")+
  labs(fill="", y="Region pairs ID", x="DiffScore")+
  xlim(-8, 8)+
  annotate('text', y=dt2$ID[floor(1225*1/4)], x=-log10(10^(-6.5)), label=paste0("DRP\nGreater direction"), size=3, color="black", angle=270)+
  annotate('text', y=dt2$ID[floor(1225*3/4)], x=log10(10^(-6.5)), label=paste0("DRP\nSmaller direction"), size=3, color="black", angle=90)+
  scale_y_discrete(breaks = c(seq(1, 1225, by = 250), 1225), labels=c(seq(1, 1225, by=250), 1225))
dev.off()
