#######################################################
###### bar plot. 
########################################################
library(ggplot2)
library(grid)
library(gridExtra)
# Create a data frame for the bar plot
data <- data.frame(
  Category = c('Up', 'Down'),  
  Value = c(21, 10) # number from our results
)
data$Category <- factor(data$Category, levels = c('Up', 'Down'))
p <-ggplot(data, aes(x = Category, y = Value, fill = Category)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#F28E7A",  "#76B7B2")) +
  labs(x = NULL, y = "Gene count", fill=NULL, title = "Shorter direction") +
  theme_classic()+
  theme(legend.position = "none",
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10),
        plot.title = element_text(size = 10))
pdf(paste0('genenum.pdf'), width = 1.6, height = 1.5)
plot(p)
dev.off()


#################
## background
###################
# Create a data frame for the bar plot
data <- data.frame(
  Category = c('Up', 'Down'), 
  Value = c(94, 105)
)
data$Category <- factor(data$Category, levels = c('Up', 'Down'))
p <-ggplot(data, aes(x = Category, y = Value, fill = Category)) +
  geom_bar(stat = "identity") +
  #scale_fill_manual(values = c("#F28E7A",  "#76B7B2")) +
  scale_fill_manual(values = c("#F5B3A2",  "#B0D4D2")) +
  scale_y_continuous(position = "right") +
  labs(x = NULL, y = "Gene count", fill=NULL, title = "Non-shorter direction") +
  theme_classic()+
  theme(legend.position = "none",
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10),
        plot.title = element_text(size = 10))
pdf(paste0('genenum_bg.pdf'), width = 1.6, height = 1.5 )
plot(p)
dev.off()

