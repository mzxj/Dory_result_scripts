
################################################
### box plot for ABchange score distribution 
#################################################

library(ggplot2)
library(dplyr)
ABchangePath<- 'SourceData/ABInteractionChange_Results/'
ABfiles <- list.files(ABchangePath, '.tsv')
ctset <- gsub(".tsv", "", ABfiles)

df <- c()
for(ct in ctset){
  dataAB <- read.table(paste0(ABchangePath, ct, '.tsv'), header = TRUE)
  ABchange_p <- dataAB[which(dataAB$ABchange > 0), ]
  ABchange_n <- dataAB[which(dataAB$ABchange < 0),]
  p_AAAB <- ABchange_p[which(ABchange_p$AB_ct1 %in% c('AA') & ABchange_p$AB_ct2 %in% c('AB', 'BA')), 'ABchange']
  p_BBAB <- ABchange_p[which(ABchange_p$AB_ct1 %in% c('BB') & ABchange_p$AB_ct2 %in% c('AB', 'BA')), 'ABchange']
  p_sAAwAA <- ABchange_p[which(ABchange_p$AB_ct1 %in% c('AA') & ABchange_p$AB_ct2 %in% c('AA')), 'ABchange']
  p_sAAwBB <- ABchange_p[which(ABchange_p$AB_ct1 %in% c('AA') & ABchange_p$AB_ct2 %in% c('BB')), 'ABchange']
  p_sBBwBB <- ABchange_p[which(ABchange_p$AB_ct1 %in% c('BB') & ABchange_p$AB_ct2 %in% c('BB')), 'ABchange']
  p_sBBwAA <- ABchange_p[which(ABchange_p$AB_ct1 %in% c('BB') & ABchange_p$AB_ct2 %in% c('AA')), 'ABchange']
  p_ABAB <- ABchange_p[which(ABchange_p$AB_ct1 %in% c('AB', 'BA') & ABchange_p$AB_ct2 %in% c('AB', 'BA')), 'ABchange']
  
  n_ABAA <- ABchange_n[which(ABchange_n$AB_ct1 %in% c('AB', 'BA') & ABchange_n$AB_ct2 %in% c('AA')), 'ABchange']
  n_ABBB <- ABchange_n[which(ABchange_n$AB_ct1 %in% c('AB', 'BA') & ABchange_n$AB_ct2 %in% c('BB')), 'ABchange']
  n_wAAsAA <- ABchange_n[which(ABchange_n$AB_ct1 %in% c('AA') & ABchange_n$AB_ct2 %in% c('AA')), 'ABchange']
  n_wBBsAA <- ABchange_n[which(ABchange_n$AB_ct1 %in% c('BB') & ABchange_n$AB_ct2 %in% c('AA')), 'ABchange']
  n_wBBsBB <- ABchange_n[which(ABchange_n$AB_ct1 %in% c('BB') & ABchange_n$AB_ct2 %in% c('BB')), 'ABchange']
  n_wAAsBB <- ABchange_n[which(ABchange_n$AB_ct1 %in% c('AA') & ABchange_n$AB_ct2 %in% c('BB')), 'ABchange']
  n_ABAB <- ABchange_n[which(ABchange_n$AB_ct1 %in% c('AB', 'BA') & ABchange_n$AB_ct2 %in% c('AB', 'BA')), 'ABchange']
  label_0 <- c(rep('p_AAAB', length(p_AAAB)), rep('p_BBAB', length(p_BBAB)), rep('p_sAAwAA', length(p_sAAwAA)), rep('p_sAAwBB', length(p_sAAwBB)), rep('p_sBBwBB', length(p_sBBwBB)), rep('p_sBBwAA', length(p_sBBwAA)), rep('p_ABAB', length(p_ABAB)), 
             rep('n_ABAA', length(n_ABAA)), rep('n_ABBB', length(n_ABBB)), rep('n_wAAsAA', length(n_wAAsAA)), rep('n_wBBsAA', length(n_wBBsAA)), rep('n_wBBsBB', length(n_wBBsBB)), rep('n_wAAsBB', length(n_wAAsBB)), rep('n_ABAB', length(n_ABAB)))
  label <- c(rep('AA<-AB', length(p_AAAB)), rep('BB<-AB', length(p_BBAB)), rep('sAA<-wAA', length(p_sAAwAA)), rep('sAA<-wBB', length(p_sAAwBB)), rep('sBB<-wBB', length(p_sBBwBB)), rep('sBB<-wAA', length(p_sBBwAA)), rep('wAB<-sAB', length(p_ABAB)), 
             rep('AB<-AA', length(n_ABAA)), rep('AB<-BB', length(n_ABBB)), rep('wAA<-sAA', length(n_wAAsAA)), rep('wBB<-sAA', length(n_wBBsAA)), rep('wBB<-sBB', length(n_wBBsBB)), rep('wAA<-sBB', length(n_wAAsBB)), rep('sAB<-wAB', length(n_ABAB)))
  value <- c(p_AAAB, p_BBAB, p_sAAwAA, p_sAAwBB, p_sBBwBB, p_sBBwAA, p_ABAB,
             n_ABAA, n_ABBB, n_wAAsAA, n_wBBsAA, n_wBBsBB, n_wAAsBB, n_ABAB)
  df <- rbind(df, cbind(label_0, label, value))
}

df <- as.data.frame(df)
df$label <- as.character(df$label)
df$value <- as.numeric(df$value)

lvl <- df %>%
  arrange(label_0) %>%          
  distinct(label) %>%           
  pull(label)

df <- df %>%
  mutate(
    label_f = factor(label, levels = lvl),
    x_pos   = as.integer(label_f)
  ) %>%
  add_count(label_f, name = "n")  # per-group counts

label_expr <- parse(text = gsub("<-", "%<-%", levels(df$label_f)))

pdf('ABchangeScore_p_n_box2.pdf', height = 3, width = 5)
ggplot(df, aes(x = x_pos, y = value, group = label_f)) +
  geom_boxplot(aes(fill = n), width = 0.8, outlier.size = 0.1, size = 0.3, color = "black") +
  scale_x_continuous(breaks = sort(unique(df$x_pos)),
                     labels = label_expr) +  #levels(df$label_f)
  scale_fill_viridis_c(name = "Count (per group)") +
  labs(x = "", y = "ABchange score") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )+
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray")
dev.off()

