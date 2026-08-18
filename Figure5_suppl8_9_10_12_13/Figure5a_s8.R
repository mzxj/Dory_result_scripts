
#########################################################
# check the statistical significance
#############################################################
library(ggplot2)
library(gridExtra)
library(grid)
library(emmeans)
library(car)
library(dplyr)
library(patchwork)

geneset0 <- read.table('SourceData/GeneTSSandPromoterRegion_update.tsv', header = TRUE)
geneset0$label <- geneset0$Gene
geneset0$label[which(geneset0$label == "Oit1")] <- "Fam3d"

geneset_CPD <- geneset0[which(geneset0$CandidateProgressionDrivers.CandidateTumorSuppressors.CandidateTumorInitiation == 'CPD' | geneset0$GeneName == 'Myc'), ]

geneset <- geneset_CPD
state_pvalue <- c()
posthoc <- c()
for(k in 1:dim(geneset)[1]){
  gene <- geneset$Gene[k]
  genename <- geneset$GeneName[k]
  geneid <- geneset$GeneID[k]
  
  temp_env <- new.env() 
  load(paste0('SourceData/DistanceCompareS2_Background/', genename, '.RData'), envir = temp_env)
  state2 <- temp_env[[genename]] 
  
  load(paste0('SourceData/DistanceCompare_Background/', genename, '.RData'), envir = temp_env)
  state3 <- temp_env[[genename]]
  
  data_list <- list()
  for (i in seq_along(state2)) {
    temp_df <- data.frame(dis = state2[[i]], regionID = i-1, state = "state2")
    data_list[[length(data_list) + 1]] <- temp_df
  }
  for (i in seq_along(state3)) {
    temp_df <- data.frame(dis = state3[[i]], regionID = i-1, state = "state3")
    data_list[[length(data_list) + 1]] <- temp_df
  }
  
  final_data <- do.call(rbind, data_list) 
  final_data1 <- final_data[!is.na(final_data$dis), ] 
  final_data1$regionID <- as.integer(final_data1$regionID)
  
  
  df <- final_data1
  df$state  <- factor(df$state, levels = c("state2", "state3"))
  df$regionID <- factor(df$regionID)   
  fit <- aov(dis ~ state * regionID, data = df)
  a <- as.data.frame(Anova(fit, type = 3) )
  state_pvalue <- rbind(state_pvalue, cbind(genename, a$`Pr(>F)`[2]))
  emm_by_region <- emmeans(fit, ~ state | regionID)
  b <- pairs(emm_by_region, adjust = "BH")
  contr <- test(pairs(emm_by_region, adjust = "none"), side = ">") |> as.data.frame()
  tab <- as.data.frame(contr) |>
    transmute(
      regionID,
      estimate,          # LUAD - Adenoma (because of factor order)
      SE, df, t = t.ratio, p = p.value,
      p_adj = p.adjust(p, method = "BH"), neglog10 = -log10(p_adj),
      dir = ifelse(estimate < 0, "LUAD > Adenoma", "LUAD < Adenoma"),
      star = cut(p_adj, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf), labels = c("***","**","*"," "))
    )
  tab$neglog10_cap <- pmin(tab$neglog10, 10)

  p <- ggplot(final_data1, aes(x = factor(regionID), y = dis, color = state)) +
    geom_boxplot(outliers = FALSE) +
    labs(x = NULL, y = "Spatial distance", color = NULL, title = genename) +
    theme_classic()+
    scale_color_manual(values=c( '#449945', '#ea7827'),
                       labels=c('Adenoma', 'LUAD'))+
    theme(legend.position = 'none', #legend.position = c(0.1, 0.9),
          axis.title = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          axis.text.x = element_blank(),
          plot.margin = margin(2, 5, 0, 5),
          plot.title = element_text(hjust = 0.5))+
    scale_x_discrete(breaks = c(1, seq(5, length(state2), by = 5)), labels=c(1, seq(5, length(state2), by=5)))
  
  
  P <- ggplot(tab, aes(x = regionID, y = 1, fill = neglog10_cap)) +
    geom_tile(color = "grey70",height = 0.1) +
    geom_text(aes(label = star), size = 1.8, angle=90, vjust= 0.8) +
    scale_fill_gradient(
      name = expression(-log[10]("FDR p")),
      low = "white", high = "purple", limits = c(0, 10)  # grayscale to avoid red/blue
    ) +
    scale_y_continuous(NULL, breaks = NULL) +
    labs(x = 'Genomic distance', y = NULL) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size=10),
      legend.position = "none",
      plot.margin = margin(0, 5, 2, 5),
      axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 1)
    )+
    scale_x_discrete(breaks = c(1, seq(5, length(state2), by = 5)), labels=c(1, seq(5, length(state2), by=5)))
  
  
  combined <- p / P + plot_layout(heights = c(7, 1))
  
  pdf(paste0('posthoc_', geneset$GeneName[k],'.pdf'), height = 1.8, width = 3.5)
  print(combined)
  dev.off()
  
}




write.table(state_pvalue, 'ANOVA_pvalue.tsv', quote=FALSE, sep="\t", row.names = FALSE)




