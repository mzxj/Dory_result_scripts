
#######################################################
###### bar plot. show the trace number in erythroblast and proerythroblast
########################################################
library(ggplot2)
library(grid)
library(gridExtra)
data <- data.frame(
  Category = c('Proerythroblast', 'Erythroblast'),  
  Value = c(1490, 2615) ## number from the data statistics 
)
data$Category <- factor(data$Category, levels = c( 'Proerythroblast', 'Erythroblast'))

pdf(paste0('f2a_tracenum.pdf'), width = 1.2, height = 2.1)
ggplot(data, aes(x = Category, y = Value, fill = Category)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c(Proerythroblast="#bcd1bc",  Erythroblast="#f3dac0")) +
  labs(x = NULL, y = "Trace number", fill=NULL) +
  theme_bw()+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size=10),
        text=element_text(size=10))
dev.off()

