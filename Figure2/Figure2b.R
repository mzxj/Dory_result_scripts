library(ggplot2)
library(gridExtra)
library(grid)
library(ggpubr)
library(RColorBrewer)
#########################
##the whole sample
#########################
lists <- "Rep4_WpCauchyPGreatLess_rep4ct6VSrep4ct7.tsv"
lists1 <- sub("WpCauchyPGreatLess_", "", lists)
objs <- as.character(t(data.frame(strsplit(lists1, "[.]"))[1,]))
for (i in 1:length(lists)) {
  assign(objs[i],read.table(paste0("SourceData/",lists[i]), sep="\t", header = TRUE)) 
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


##################################
##### sub sample
#####################################
celltype <- data.frame(ID=c(paste0("ct", seq(6,7))), name = c( "Erythroblast", "Proerythroblast"))

suball <- seq(100, 1000, 100)
cti <- 6
ctj <- 7

correcords <- c()
sub <- 500
correcord <- c()
k <- 2
data0 <- read.table(paste0("SourceData/Rep4_WpCauchyPGreatLess_rep4ct", cti,"sub", sub, "random", k, "VSrep4ct", ctj,"sub", sub, "random", k, ".tsv"), header = TRUE)

data <- c()
for(i in 1:dim(data0)[1]){
  a1 <- c(-log10(data0$pgreat[i]), log10(data0$pless[i]))
  if(!is.na(data0$pgreat[i])){
    a2 <-  a1[which.max(abs(a1))]
  }else{
    a2 <- NA
  }
  data <- rbind(data, cbind(-log10(data0$pgreat[i]), log10(data0$pless[i]), a2))
}
colnames(data) <- c("pgreat", "pless", "pfinal")
data <- as.data.frame(data)
pfinalMat <- data$pfinal
pfinalMat_sub <- pfinalMat

data <- data.frame(sub358 = pfinalMat_sub,
                   whole = pfinalMat_all[,1],
                   PointDense = densCols(pfinalMat_sub, pfinalMat_all[,1]))
cor_result <- cor.test(data$sub358, data$whole)
r_value <- round(cor_result$estimate, 3)
p_value <- signif(cor_result$p.value, 3)

plt <- ggplot(data, aes(x=sub358, y=whole, color=I(PointDense)))+
  geom_point(size=0.5)+
  labs(x=paste0("DiffScore (", sub, " traces)"), 
       y=paste0("DiffScore (the whole traces)"))+
  theme_bw()+
  theme(axis.title = element_text(size = 10),
        panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10))+
  annotate("text", x = -6, y = 13, label = paste0("R = ", round(cor_result$estimate,3), "\np < 2.2e-16"), size = 3, hjust=0)

pdf(paste0(sub,'case', k, '.pdf'), width = 2.2, height = 2.2)
print(plt)
dev.off()


