
#########################
#heatmap
##########################
library(ggplot2)
dt <- read.table('SourceData/Rep4_WpCauchypMatrix_rep4ct3VSrep4ct4.tsv', sep = "\t")
colnames(dt) <- paste0("region", seq(1:dim(dt)[1]))
rownames(dt) <- paste0("region", seq(1:dim(dt)[1]))
data <- reshape2::melt(dplyr::mutate(as.data.frame(dt), index=row.names(dt)), id="index")
colnames(data) <- c("T1", "T2", "value")
data$T1 <- factor(data$T1, levels=paste0("region", seq(1,dim(dt)[1],1)))
data$T2 <- factor(data$T2, levels=paste0("region", seq(dim(dt)[1],1,-1)))
png('Fig1_diff_Dory.png', height = 1700, width = 2210, res=600)
ggplot(data,aes(x=T1,y=T2,fill=value))+ 
  scale_fill_gradient2(low="#d60b0e", mid="white", high="#0f70bf", midpoint = 0, na.value = "grey", limits=c(-max(abs(data$value), na.rm = TRUE), max(abs(data$value), na.rm = TRUE)))+
  geom_raster()+
  labs(fill="DiffScore",
       x="Region ID",
       y="Region ID")+
  scale_x_discrete(breaks = paste0("region", seq(5, 50, by = 5)), labels=seq(5, 50, by=5))+
  scale_y_discrete(breaks = paste0("region", seq(5, 50, by = 5)), labels=seq(5, 50, by=5))
dev.off()



