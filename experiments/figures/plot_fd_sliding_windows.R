## plotting fd sliding windows

library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggridges)
library(viridis)
library(reshape2)

spp_list <- c("lsta.lmic",
              "lmic.lsta",
              "lang.lmar",
              "lmar.lang",
              "lmic.lmar",
              "lmar.lmic",
              "lmic.lang",
              "lang.lmic",
              "lsta.lmar",
              "lmar.lsta",
              "lsta.lang",
              "lang.lsta")

for(s in spp_list) {
  fd <- read_csv(paste0("C:/Users/jrick/Dropbox/i/projects/lates_admixture/results/fd_sliding_windows/latesWGS_july2020.beagle.w500.",s))
  fd$fd[fd$D < 0] <- 0
  fd <- fd %>% mutate(spp=s)
  
  assign(paste0(s,"_fd"),fd)
  
  cutoff1 <- quantile(fd$fd[fd$fd >= 0 & fd$fd <= 1],
                      probs=c(0.975),na.rm=T)
  cutoff2 <- quantile(fd$fd[fd$fd >= 0 & fd$fd <= 1],
                      probs=c(0.995),na.rm=T)
  
  plot <- ggplot(fd, aes(x=mid, y=fd)) +
    geom_line() +
    facet_wrap(~scaffold, scales="free_x") +
    ylim(0,1) +
    theme_pubr() +
    theme(axis.title=element_blank(),
          axis.text.x=element_text(size=12),
          axis.text.y=element_blank(),
          legend.position="none") +
    geom_hline(yintercept=cutoff1,lty=2,col="red")+
    geom_hline(yintercept=cutoff2,lty=2,col="blue")
  assign(paste0(s,"_plot"),plot)
  #print(plot)
  
  props <- data.frame(spp=s,
             cutoff_95=length(fd$fd[fd$fd > cutoff1])/length(fd$fd),
             cutoff_99=length(fd$fd[fd$fd > cutoff2])/length(fd$fd))
  assign(paste0(s,"_props"),props)
}

lang.lmar_plot
lmar.lang_plot
lmic.lsta_plot
lsta.lmic_plot

all_props <- bind_rows(lsta.lmic_props,
                       lmic.lsta_props,
                       lang.lmar_props,
                       lmar.lang_props,
                       lmic.lang_props,
                       lang.lmic_props,
                       lsta.lmar_props,
                       lmar.lsta_props,
                       lsta.lang_props,
                       lang.lsta_props)

ggplot(all_props, aes(x=spp)) +
  geom_boxplot(aes(y=cutoff_95)) +
  theme_pubr() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=14),
        legend.position="none") +
  ylab("proportion of genome admixed with 95% threshold")

########

all_fd <- bind_rows(lsta.lmic_fd,
                    lmic.lsta_fd,
                    lang.lmar_fd,
                    lmar.lang_fd,
                    lmic.lang_fd,
                    lang.lmic_fd,
                    lsta.lmar_fd,
                    lmar.lsta_fd,
                    lsta.lang_fd,
                    lang.lsta_fd,
                    lmar.lmic_fd,
                    lmic.lmar_fd
                    )

ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

pdf("fd_sliding_windows_byChrom_fdM_w500.pdf", width=11, height=5)
for (chr in unique(all_fd$scaffold)){
  cutoff <- quantile(all_fd$fdM[all_fd$scaffold == chr],probs=c(0.975),na.rm=T)
    
  plot <- all_fd %>%
    filter(scaffold == chr) %>%
    drop_na() %>%
    mutate(color=ifelse(fdM > cutoff, spp, "ns")) %>%
    ggplot(aes(x=mid, y=fdM)) +
    geom_point(aes(color=spp),alpha=0.5) +
    #ylim(0,1) +
    theme_pubr() +
    theme(axis.title=element_blank(),
          axis.text.x=element_text(size=12),
          axis.text.y=element_text(size=12),
          legend.position="right") +
    #scale_color_manual(name="Spp pair",values=c(ggplotColours(12),"gray80")) +
    ggtitle(chr) #+
    #facet_wrap(~spp, ncol=3)
  print(plot)
}
dev.off()

all_fd %>%
  ggplot(aes(x=fd,y=spp,fill=..x..)) +
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01)+
  theme_ridges() +
  scale_fill_viridis(name="fd", option="C") +
  theme(axis.title.x = element_blank(),
        legend.position = "none",
        axis.text = element_text(size=10))  +
  xlim(0,0.5)+
  facet_wrap(~scaffold, scales="free_x",
             strip.position="right")

length(fd$fd[fd$fd > cutoff1])/length(fd$fd)
length(fd$fd[fd$fd > cutoff2])/length(fd$fd)

## looking at correlation between species
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}

cormat <- all_fd %>%
  pivot_wider(names_from=spp,values_from=c(ABBA,BABA,D,fd,fdM)) %>%
  select(fd_lsta.lmic:fd_lang.lsta) %>%
  cor(.,use="pairwise.complete.obs") %>%
  reorder_cormat() %>%
  get_upper_tri() %>%
  melt(.,na.rm=T)
  

ggheatmap <- ggplot(cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()

ggheatmap + 
  #geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.6, 0.7),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))

corrplot2 <- function(data,
                      method = "pearson",
                      sig.level = 0.05,
                      order = "original",
                      diag = FALSE,
                      type = "upper",
                      tl.srt = 90,
                      number.font = 1,
                      number.cex = 1,
                      mar = c(0, 0, 0, 0)) {
  library(corrplot)
  data_incomplete <- data
  #data <- data[complete.cases(data), ]
  mat <- cor(data, method = method, use="pairwise.complete.obs")
  cor.mtest <- function(mat, method) {
    mat <- as.matrix(mat)
    n <- ncol(mat)
    p.mat <- matrix(NA, n, n)
    diag(p.mat) <- 0
    for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
        tmp <- cor.test(mat[, i], mat[, j], method = method)
        p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
      }
    }
    colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
    p.mat
  }
  p.mat <- cor.mtest(data, method = method)
  col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
  corrplot(mat,
           method = "color", col = col(200), number.font = number.font,
           mar = mar, number.cex = number.cex,
           type = type, order = order,
           addCoef.col = "black", # add correlation coefficient
           tl.col = "black", tl.srt = tl.srt, # rotation of text labels
           # combine with significance level
           p.mat = p.mat, sig.level = sig.level, insig = "blank",
           # hide correlation coefficients on the diagonal
           diag = diag
  )
}

df <- all_fd %>%
  pivot_wider(names_from=spp,values_from=c(ABBA,BABA,D,fd,fdM)) %>%
  select(fdM_lsta.lmic:fdM_lang.lsta)

corrplot2(
  data = df,
  method = "pearson",
  sig.level = 0.05,
  order = "original",
  diag = FALSE,
  type = "upper",
  tl.srt = 75
)
