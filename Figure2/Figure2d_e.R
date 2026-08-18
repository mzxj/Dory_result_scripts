##############################################
#reproducibility score for each subsample
############################################
library(ggplot2)
celltype <- data.frame(ID=c(paste0("ct", seq(6,7))),
                       name = c("Erythroblast", "Proerythroblast"))

top <- 50
objs1 <- seq(100, 1000, 100) ##sub id
objs2 <- paste0("random", seq(1,20,1))

cti <- 6
ctj <- 7

datasubs <- c()
for(obj1 in objs1){
  sigid_p <- list()
  sigid_n <- list()
  for(obj2 in objs2){
    data_o <- read.table(paste0('SourceData/sub', obj1, '/Rep4_WpCauchyPGreatLess_rep4ct', cti, 'sub', obj1, obj2, 'VSrep4ct', ctj,'sub', obj1, obj2, '.tsv'), header = TRUE)
    data <- data.frame(k=seq(1, dim(data_o)[1], 1), -log10(data_o))
    for(i in 1:dim(data)[1]){
      a1 <- c(data$pgreat[i], -data$pless[i])
      if(!is.na(data$pgreat[i])){
        data$pfinal[i] <-  a1[which.max(abs(a1))]
      }else{
        data$pfinal[i] <- NA
      }
    }
    sigid_p[[paste0(obj2, '_sig_p')]] <- order(data$pfinal, decreasing = TRUE)[1:top]
    sigid_n[[paste0(obj2, '_sig_n')]] <- order(data$pfinal, decreasing = FALSE)[1:top]
  }
  
  rep_p <- data.frame(celltypeid=paste0('ct',cti, 'VSct', ctj),
                      celltype = paste0(celltype$name[1], 'VS', celltype$name[2]),
                      subsample = obj1)
  rep_n <- data.frame(celltypeid=paste0('ct',cti, 'VSct', ctj),
                      celltype = paste0(celltype$name[1], 'VS', celltype$name[2]),
                      subsample = obj1)
  rep <- data.frame(celltypeid=paste0('ct',cti, 'VSct', ctj),
                    celltype = paste0(celltype$name[1], 'VS', celltype$name[2]),
                    subsample = obj1)
  for(k in 1:length(objs2)){
    p <- length(intersect(unlist(sigid_p[[k]]), as.vector(unique(unlist(sigid_p[-k])))))
    n <- length(intersect(unlist(sigid_n[[k]]), as.vector(unique(unlist(sigid_n[-k])))))
    rep_p <- cbind(rep_p, p/length(unlist(sigid_p[[k]])))
    rep_n <- cbind(rep_n, n/length(unlist(sigid_n[[k]])))
    rep <- cbind(rep, (p + n)/(length(sigid_p[[k]])+length(sigid_n[[k]])))
  }
  colnames(rep) <- c('celltypeid', 'celltype', 'subsample', paste0("reproducibility_random", seq(1:20)))
  colnames(rep_p) <- c('celltypeid', 'celltype', 'subsample', paste0("reproducibility_random", seq(1:20)))
  colnames(rep_n) <- c( 'celltypeid', 'celltype', 'subsample', paste0("reproducibility_random", seq(1:20)))
  
  dataout_both <- data.frame(direction=c(rep("Greater", dim(rep_p)[1]),
                                         rep("Less", dim(rep_n)[1])),
                             rbind(rep_p, rep_n))
  
  dataall <- rbind(data.frame(direction = rep("Both", dim(rep)[1]), rep),
                   dataout_both)
  
  datasubs <- rbind(datasubs, dataall)
}


directions <- data.frame(obj= c("Greater", "Less"),
                         #name = c("longer", "shorter"),
                         name = c("Greater", "Less"),
                         values = c("#0f70bf",  "#d60b0e"),
                         cols = c("#5aadd0", "#af7aa1"))

#######
##all sub traces
#################
labelname <- c("Greater", "Smaller")
for(i in 1:dim(directions)[1]){
  dataplot0 <- datasubs[which(datasubs$direction == directions[i, 1] & datasubs$subsample %in% seq(100, 1000, 100)), ]
  dataall_plot <- dataplot0 [, c(-1, -2, -3)]
  dataallplot <- reshape2::melt(dataall_plot, id=c( 'subsample'))
  dataallplot$subsample <- factor(dataallplot$subsample, levels = seq(100, 1000, 100))
  
  pdf(paste0('reproducibility_allsubtraces_', labelname[i], '.pdf'),width = 3.2, height = 2.3)
  plt<- ggplot(dataallplot, aes(x=subsample, y=value))+
    geom_boxplot(outliers = FALSE, alpha=0.5, size=0.2, fill=directions[i, 3])+
    geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.5, stroke=0.5, fill=directions[i, 3]) + 
    theme_bw()+
    labs(x=paste0("Trace count"), y="Reproducibility rate",title=paste0("Top 50/1225 region pairs (",labelname[i], ")" ))+
    theme(legend.position = 'none',
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=10),
          axis.text.y = element_text(size=10),
          axis.title = element_text(size = 10),
          text=element_text(size=10),
          panel.grid.major = element_blank(),  
          panel.grid.minor = element_blank(),
          plot.title = element_text(size = 10))+
    scale_y_continuous(limits = c(0, 1))#+
  print(plt)
  dev.off()
}


