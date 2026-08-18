library(ggplot2)
library(gridExtra)
library(grid)

#########################
##the whole sample
#########################
lists <- "Rep4_WpCauchyPGreatLess_rep4ct6VSrep4ct7.tsv"
lists1 <- sub("WpCauchyPGreatLess_", "", lists)
objs <- as.character(t(data.frame(strsplit(lists1, "[.]"))[1,]))
for (i in 1:length(lists)) {
  assign(objs[i],read.table(paste0("SourceData/",lists[i]), sep="\t", header = TRUE)) ##assign 函数 Assign a Value to a Name
}

objslog <- paste0("obj", 1:length(objs), "log")
for(k in 1:length(objslog)){
  data <- c()
  for(i in 1:dim(get(objs[k]))[1]){
    a1 <- c(-log10(get(objs[k])$pgreat[i]), log10(get(objs[k])$pless[i]))
    if(!is.na(get(objs[k])$pgreat[i])){
      a2 <-  a1[which.max(abs(a1))]
    }else{
      a2 <- NA
    }
    data <- rbind(data, cbind(-log10(get(objs[k])$pgreat[i]), log10(get(objs[k])$pless[i]), a2))
  }
  colnames(data) <- c("pgreat", "pless", "pfinal")
  assign(objslog[k], as.data.frame(data))
}

pfinalMat <- c()
for(k in 1:length(objslog)){
  pfinalMat <- cbind(pfinalMat, get(objslog[k])$pfinal)
}

pfinalMat_all <- pfinalMat


###########################################
### sub
###########################################
celltype <- data.frame(ID=c(paste0("ct", seq(6,7))),
                       name = c( "Erythroblast", "Proerythroblast"))

subs <- paste0("sub", seq(100, 1000, 100))
pcor_allsubs <- c()
for(sub in subs){
  cti <- 6
  ctj <- 7
  for(k in 1:20){
    data <- read.table(paste0("SourceData/", sub, "/Rep4_WpCauchyPGreatLess_rep4ct", cti, sub, "random", k, "VSrep4ct", ctj, sub, "random", k, ".tsv"), header = TRUE)
    dtin <- c()
    for(i in 1:dim(data)[1]){
      a1 <- c(-log10(data$pgreat[i]), log10(data$pless[i]))
      if(!is.na(data$pgreat[i])){
        a2 <-  a1[which.max(abs(a1))]
      }else{
        a2 <- NA
      }
      dtin <- rbind(dtin, a2)
    }
    pfinalMat_sub <- dtin
    pcor_pearson <- cor(as.numeric(pfinalMat_sub), as.numeric(pfinalMat_all), method = "pearson", use = "complete.obs")
    cordf<- data.frame(Subs = sub,
                       Random = k,
                       cor_pearson = pcor_pearson)
    pcor_allsubs <- rbind(pcor_allsubs, cordf)
  }
}

dataplot <- reshape2::melt(pcor_allsubs, id=c("Subs", "Random"))
colnames(dataplot) <- c("Subs", "Random", "Method", "Number")

dataplot$Subs <- sub("sub", "", dataplot$Subs)
dataplot$Subs <- factor(dataplot$Subs, levels =  c(seq(100,1000,100)))

pdf('cor_boxplot_pearson.pdf',width = 2.9, height = 2.3)
ggplot(dataplot, aes(x=Subs, y=Number))+
  geom_boxplot(outlier.size = 0.5, color = "#3d5c6f")+
  labs(x="Number of sub traces",  
       y="Pearson\ncorrelation coefficient", 
       title ="PCC(sub traces, the whole traces)")+ 
  ylim(0, 1)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=10),
        axis.text.y = element_text(size=10),
        axis.title = element_text(size = 10),
        text=element_text(size=10),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),
        legend.position = "none",
        plot.title = element_text(size = 10))
dev.off()
