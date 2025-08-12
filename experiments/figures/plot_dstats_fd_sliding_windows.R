## plotting window f_d statistics from dsuite, sept 2021

library(tidyverse)
library(ggpubr)
library(ggstance)
library(ggdist)

## import windowed dstats per chromosome
dir <- "~/Dropbox/i/projects/lates_admixture/results/dstats_feb2023/"
files <- list.files(dir, pattern="*100_20.txt",full.names=TRUE)

read_file <- function(file) {
  trio <- gsub(pattern=".*/(L[a-z]+_L[a-z]+_L[a-z]+)_localFstats.*","\\1",x=file,fixed=FALSE)
  dat <- read_table2(file) %>%
    mutate(f_d_corr = case_when(D < 0 ~ 0,
                                TRUE ~ f_d),
           trio = trio,
           chr_num = gsub(pattern="ASB_LG([0-9]+)","\\1",chr,fixed=FALSE))
  return(dat)
}

dstats <- files %>%
  map_dfr(~read_file(.)) %>%
  mutate(chr_num = case_when(chr_num == "16_LG22" ~ "16",
                             chr_num == "7_1" ~ "7.1",
                             chr_num == "7_2" ~ "7.2",
                             TRUE ~ chr_num)) %>%
  mutate(zscore = (f_d_corr - mean(f_d_corr, na.rm=T)) / sd(f_d_corr, na.rm=T),
         pval = pnorm(-zscore)) 

dstat_sum <- dstats %>%
  group_by(trio,chr) %>%
  summarize(mean = mean(f_d_corr, na.rm = TRUE),
            median = median(f_d_corr, na.rm = TRUE),
            thr95 = quantile(f_d_corr,0.95, na.rm = TRUE),
            thr99 = quantile(f_d_corr,0.99, na.rm = TRUE),
            thr995 = quantile(f_d_corr,0.995, na.rm = TRUE),
            max = max(f_d_corr, na.rm = TRUE))
dstats %>%
  summarize(mean=mean(f_d_corr, na.rm = TRUE),
            median = median(f_d_corr, na.rm = TRUE),
            thr95 = quantile(f_d_corr,0.95, na.rm = TRUE),
            thr99 = quantile(f_d_corr,0.99, na.rm = TRUE),
            thr995 = quantile(f_d_corr,0.995, na.rm = TRUE),
            max = max(f_d_corr, na.rm = TRUE))

dstats %>%
  filter(D > 0) %>%
  ggplot(aes(x=f_d_corr,fill=chr)) +
    geom_density(alpha=0.5) +
    facet_wrap(~trio)

dstats %>%
  filter(chr == "ASB_LG14") %>%
  ggplot(aes(x=(windowStart+windowEnd)/2,y=f_d_corr,col=trio)) +
    geom_point(alpha=0.5, size=1.5) +
    theme_custom() +
    facet_wrap(~trio, scales="free_x") +
    ylim(0,1) +
  xlab("Position") +
  ylab("Corrected fd Statistic") +
    geom_hline(yintercept=quantile(dstats$f_d_corr,0.99, na.rm = TRUE),lty=2) +
    geom_hline(yintercept=quantile(dstats$f_d_corr,0.95, na.rm = TRUE),lty=3) +
    geom_hline(yintercept=quantile(dstats$f_d_corr,0.995, na.rm = TRUE),lty=1)

## attempt at manhattan plot
data_cum <- dstats %>%
  mutate(chr_num = as.numeric(chr_num)) %>%
  arrange(desc(chr_num)) %>% 
  group_by(chr_num,chr) %>%
  summarize(max_bp = max(windowEnd)) %>%
  ungroup() %>%
  mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>%
  dplyr::select(chr_num, chr, bp_add)
dstats_cum <- dstats %>%
  inner_join(data_cum, by="chr") %>%
  mutate(windowMid_cum = (windowStart + windowEnd)/2 + bp_add)

axis_set <- dstats_cum %>%
  group_by(chr) %>%
  summarize(center = mean(windowMid_cum)) %>% 
  mutate(chr=gsub("ASB_","",chr))

colors <- c("#264653","#2a9d8f","#e9c46a","#f4a261","#e76f51")

## plotting f_d (with negative values set to zero)
fig3d <- dstats_cum %>% 
  mutate(chr=factor(gsub("ASB_","",chr),
                    levels=c( "LG1","LG2","LG3","LG4","LG5","LG6",
                              "LG7_1","LG7_2","LG8","LG9","LG10","LG11","LG12",
                              "LG13","LG14","LG15","LG16_LG22","LG17","LG18","LG19",
                              "LG20","LG21","LG23","LG24"))) %>%
  mutate(signif = case_when(pval < 0.01 ~ "p0.01",
                            pval >= 0.01 & pval < 0.05 ~ "p0.05",
                            TRUE ~ "ns")) %>%
  ggplot(aes(x=windowMid_cum, y=f_d_corr,col=chr)) +
    geom_point(data= . %>% filter(signif == "ns"),alpha=0.2, size=0.6) +
    geom_point(data= . %>% filter(signif == "p0.05"),alpha=0.6, size=0.6) +
    geom_point(data= . %>% filter(signif == "p0.01"),alpha=0.6, size=0.6) +
#    geom_hline(yintercept=quantile(dstats$f_d_corr,0.90),lty=2, col="gray50") +
    geom_hline(yintercept=-qnorm(0.01) * sd(dstats_cum$f_d_corr,na.rm=T) + mean(dstats_cum$f_d_corr,na.rm=T),lty=2, col="gray20") +
#    geom_hline(yintercept=quantile(dstats$f_d_corr,0.95),lty=4, col="gray50") +
    scale_x_continuous(label = axis_set$chr, breaks = axis_set$center, 
                       limits = c(-10,max(dstats_cum$windowMid_cum)+10),
                       expand = c(0.01,0.01)) +
    scale_color_manual(values = rep(colors, unique(length(axis_set$chr)))) +
    #scale_color_viridis_d() +
    scale_y_continuous(limits=c(0,1),expand=c(0.01,0.01)) +
    xlab("") +
    ylab("Corrected f_d Statistic") +
    theme_custom() +
    theme(legend.position="none",
          axis.text.x = element_text(angle=60, size=16, vjust=1, hjust=1)) +
    facet_wrap(~trio, ncol=1)

## plotting f_dM
dstats_cum %>%
  ggplot(aes(x=windowMid_cum, y=f_dM,col=as.factor(chr))) +
  geom_point(alpha=0.3, size=0.6) +
  #geom_hline(yintercept=quantile(dstats$f_dM,0.99),lty=2, col="gray50") +
  #geom_hline(yintercept=quantile(dstats$f_dM,0.995),lty=3, col="gray50") +
  #geom_hline(yintercept=quantile(dstats$f_dM,0.95),lty=4, col="gray50") +
  scale_x_continuous(label = axis_set$chr, breaks = axis_set$center, 
                     limits = c(-10,max(dstats_cum$windowMid_cum)+10),
                     expand = c(0,0)) +
  scale_color_manual(values = rep(colors, unique(length(axis_set$chr)))) +
  #ylim(0,1) +
  ylab("f_dM Statistic") +
  theme_custom() +
  theme(legend.position="none",
        axis.text.x = element_text(angle=90, size=8, vjust=1)) +
  facet_wrap(~trio, ncol=1)

## extracting locations with significant signal of admixture
## and writing locations to a bed file
dstats %>%
  #filter(f_d_corr < quantile(f_d_corr,0.95,na.rm=T)) %>%
  filter(zscore < 2.33) %>%
  select(chr,windowStart,windowEnd) %>%
  write_delim(file="~/Dropbox/i/projects/lates_admixture/results/dstats_feb2023/low_fd_windows_95quantile.bed",delim="\t",col_names=FALSE)

dstats %>%
  #filter(f_d_corr >= quantile(f_d_corr,0.90)) %>%
  filter(zscore >= 2.33) %>%
  select(chr,windowStart,windowEnd) %>%
  write_delim(file="~/Dropbox/i/projects/lates_admixture/results/dstats_feb2023/high_fd_windows_90quantile.bed",delim="\t",col_names=FALSE)

dstats %>%
  # mutate(set = case_when(f_d_corr >= quantile(f_d_corr,0.90) ~ "omit",
  #                        f_d_corr < quantile(f_d_corr,0.90) ~ "keep",
  #                        TRUE ~ "other")) %>%
  mutate(set = case_when(zscore >= 2.33 ~ "omit",
                         zscore < 2.33 ~ "keep",
                         TRUE ~ "other")) %>%
  select(chr,windowStart,windowEnd,set) %>%
  write_delim(file="~/Dropbox/i/projects/lates_admixture/results/dstats_feb2023/fd_windows_90quantile.bed",delim="\t",col_names=FALSE)

## quantifying significant windows per comparison
dstats %>%
  mutate(window_length=windowEnd-windowStart,
         signif = case_when(pval >= 0.05 ~ "no",
                            pval < 0.05 & pval >= 0.01 ~ "p0.05",
                            pval < 0.01 ~ "p0.01")) %>%
  group_by(trio,signif,chr) %>%
  summarize(n_window = n(),
            total_length = sum(window_length)) %>%
  group_by(trio,chr) %>%
  mutate(prop_window = n_window/sum(n_window),
         prop_length = total_length/sum(total_length)) %>%
  filter(signif == "p0.01") %>%
  ggplot(aes(x=trio,y=prop_length)) +
  geom_point(aes(col=chr)) +
  geom_boxplot(aes(),size=0.8,alpha=0.5) +
  theme_custom() +
  theme(panel.grid.major=element_line(linetype="dashed"))

dstats_signif <- dstats %>%
  mutate(window_length=windowEnd-windowStart,
         signif = case_when(pval >= 0.05 ~ "no",
                            pval < 0.05 & pval >= 0.01 ~ "p0.05",
                            pval < 0.01 ~ "p0.01")) %>%
  group_by(chr,trio) %>%
  arrange(windowStart) %>%
  mutate(same = case_when(signif == lag(signif) ~ "same",
                          signif != lag(signif) ~ "diff",
                          TRUE ~ "diff"),
         window_group = 1)

# identifying consecutive windows
dstats_signif_blocks <- dstats_signif %>%
  ungroup() %>%
  filter(chr == 1)
for (t in unique(dstats_signif$trio)){
  for (c in unique(dstats_signif$chr)) {
    subset <- dstats_signif[dstats_signif$chr == c & dstats_signif$trio == t,]
    for (i in 2:nrow(subset)) {
      if (subset$same[i] == "same") {
        subset$window_group[i] <- subset$window_group[i-1]
      } else {
        subset$window_group[i] <- subset$window_group[i-1] + 1
      }
    }
    dstats_signif_blocks <- dstats_signif_blocks %>%
      add_row(subset)
  }
}

dstats_signif_sum <- dstats_signif_blocks %>%
  group_by(chr,trio,window_group,signif) %>%
  summarize(window_min = min(windowStart),
            window_max = max(windowEnd),
            window_total_length = window_max - window_min,
            mean_fd_corr = mean(f_d_corr),
            median_fd_corr = median(f_d_corr),
            mean_pval = mean(pval))

dstats_signif_sum %>%
  filter(signif == "p0.01") %>% 
  ggplot() +
  geom_jitter(aes(x=trio,y=window_total_length), alpha=0.3, width=0.05, height=1) +
  ggdist::stat_halfeye(aes(x=trio,y=window_total_length, col=trio), 
                       slab_alpha=0.5, point_size = 3, 
                       point_interval="mean_qi", position = position_nudge(x=0.1),
                       adjust=0.8) +
  theme_custom() +
  theme(panel.grid.major = element_line(linetype="dashed"),
        axis.title.y = element_blank(),
        legend.position = "none") +
  ylab("Admixture Window Length (bp)") +
  coord_flip() +
  geom_text(data=summarize(group_by(filter(dstats_signif_sum,signif == "p0.01"),trio),
                           mean = mean(window_total_length), max = max(window_total_length)),
            aes(y=100000,x=trio,label=paste("mean = ",scales::comma(round(mean,1)),"bp\nmax = ",scales::comma(round(max,0)),"bp"),color=trio),
            nudge_x=0.4, family = "Open Sans")

window_length_table <- dstats_signif_sum %>%
  filter(signif == "p0.01") %>%
  group_by(trio,chr) %>%
  summarize(mean_window_length = mean(window_total_length),
            total_window_length = sum(window_total_length),
            mean_fd_corr = mean(mean_fd_corr)) %>%
  arrange(trio,chr) %>%
  kableExtra::kable(format = "simple",digits = 3,align="l",caption="table caption here")

dstats_signif_sum %>%
  group_by(trio,signif,chr) %>%
  summarize(mean_window_length = mean(window_total_length),
            total_window_length = sum(window_total_length),
            mean_fd_corr = mean(mean_fd_corr)) %>%
  group_by(trio,chr) %>%
  mutate(prop_total_length = total_window_length/sum(total_window_length)) %>%
  filter(signif == "p0.01") %>%
  ggplot(aes(x=trio,y=prop_total_length)) +
  geom_point(aes(col=chr,fill=chr),size=4,alpha=0.6,pch=19) +
  geom_text(data=. %>% filter(chr %in% c("ASB_LG14","ASB_LG1","ASB_LG20","ASB_LG18","ASB_LG3")),
            aes(x=trio,y=prop_total_length,label=chr),
            position=position_nudge(x=0.1), 
            hjust=0, size=3) +
  #geom_point(aes(col=chr),size=3,alpha=0.6) +
  stat_summary(size=1.2,fun.data="mean_cl_boot",
               position=position_nudge(x=-0.2)) +
  theme_custom() +
  ylab("Proportion of Sequenced Genome Admixed") + 
  theme(axis.text.x = element_text(angle=60,hjust=1),
        axis.title.x = element_blank()) +
  coord_cartesian(xlim=c(1,6.5),clip="off")

## what about grouping by P2 taxon? because positive fd = sharing among P2+P3
dstats_signif_sum %>%
  group_by(trio,signif,chr) %>%
  summarize(mean_window_length = mean(window_total_length),
            total_window_length = sum(window_total_length),
            mean_fd_corr = mean(mean_fd_corr)) %>%
  mutate(p1 = gsub("(L[a-z]+)_(L[a-z]+)_(L[a-z]+)","\\1",trio),
         p2 = gsub("(L[a-z]+)_(L[a-z]+)_(L[a-z]+)","\\2",trio),
         p3 = gsub("(L[a-z]+)_(L[a-z]+)_(L[a-z]+)","\\3",trio)) %>%
  group_by(trio,chr) %>%
  mutate(prop_total_length = total_window_length/sum(total_window_length)) %>%
  group_by(p2,chr,signif) %>%
  summarize(mean_mean_window_length = mean(mean_window_length),
         mean_total_window_length = mean(total_window_length),
         mean_prop_total_length = mean(prop_total_length)) %>%
  filter(signif == "p0.01") %>%
  ggplot(aes(x=p2,y=mean_prop_total_length)) +
  geom_point(aes(col=chr,fill=chr),size=12,alpha=0.6,pch="_") +
  geom_text(data=. %>% filter(chr %in% c("ASB_LG14","ASB_LG1","ASB_LG20","ASB_LG18","ASB_LG3")),
            aes(x=p2,y=mean_prop_total_length,label=chr),
            position=position_nudge(x=0.1), 
            hjust=0, size=3) +
  #geom_point(aes(col=chr),size=3,alpha=0.6) +
  stat_summary(size=1.2,fun.data="mean_cl_boot",
               position=position_nudge(x=-0.2)) +
  theme_custom() +
  ylab("Proportion of Sequenced Genome Admixed") + 
  theme(axis.text.x = element_text(angle=60,hjust=1),
        axis.title.x = element_blank()) +
  coord_cartesian(xlim=c(1,3.5),clip="off")
  

#######################################
##### OLD CODE BELOW
#######################################
## identifying significant peaks
dstats_cum %>%
  filter(chr == "ASB_LG14") %>%
  ggplot(aes(x=windowStart, y=f_d_corr, color=trio)) +
  geom_point(alpha=0.5,size=0.6) +
  geom_hline(yintercept=quantile(dstats$f_d_corr,0.995),lty=4, col="gray50") +
  geom_hline(yintercept=quantile(dstats$f_d_corr,0.95),lty=4, col="gray50") +
  #ylim(0,1) +
  #xlim(4.5e6,5.5e6) +
  xlab("") +
  ylab("Corrected f_d Statistic") +
  theme_custom() +
  theme(legend.position="none",
        axis.text.x = element_text(angle=90, size=8, vjust=1))

dstats_cum %>%
  filter(chr == "ASB_LG14") %>%
  filter(f_d_corr > quantile(dstats$f_d_corr, 0.995) &
           windowStart < 5.5e6 & windowStart > 4.5e6) %>%
  group_by(trio) %>%
  summarize(min=min(windowStart),
            max=max(windowEnd),
            mean_fd=mean(f_d_corr))

summaryMatrix_fd <- dstats_signif_sum %>%
  group_by(trio,signif,chr) %>%
  summarize(mean_window_length = mean(window_total_length),
            total_window_length = sum(window_total_length),
            mean_fd_corr = mean(mean_fd_corr)) %>%
  group_by(trio,chr) %>%
  mutate(prop_total_length = total_window_length/sum(total_window_length)) %>%
  filter(signif == "p0.01") %>%
  mutate(tax1 = gsub("(L[a-z]+)_(L[a-z]+)_(L[a-z]+)","\\1",trio),
         tax2 = gsub("(L[a-z]+)_(L[a-z]+)_(L[a-z]+)","\\2",trio),
         tax3 = gsub("(L[a-z]+)_(L[a-z]+)_(L[a-z]+)","\\3",trio)) %>%
  group_by(tax2,tax3) %>%
  summarize(mean_prop_total_length = mean(prop_total_length)) %>%
  ggplot(aes(x=factor(tax2,levels=c("Lsta","Lmic","Lang","Lmar","Lnil")),
             y=factor(tax3,levels=c("Lsta","Lmic","Lang","Lmar","Lnil")),
             fill = mean_prop_total_length))+
  geom_tile(color = "white")+
  geom_tile(aes(x=factor(tax3,levels=c("Lsta","Lmic","Lang","Lmar","Lnil")),
                y=factor(tax2,levels=c("Lsta","Lmic","Lang","Lmar","Lnil")), 
                fill = mean_prop_total_length), color = "white") +
  # scale_fill_gradient2(low = "white", high = "red", mid = "yellow",na.value = "grey50",
  #                      midpoint = max(introgressionSummary$mean_prop)/2, limit = c(0,max(introgressionSummary$mean_prop)), 
  #                      name="Average introgression fraction") +
  scale_fill_viridis_c(na.value="gray20",option="B",begin=0.3,direction=-1,
                       name="Average introgression proportion") +
  geom_abline(slope = 1, intercept=0)+
  geom_vline(xintercept=seq(1.5,5+0.5,1), alpha=0.3)+
  geom_hline(yintercept=seq(1.5,5+0.5,1), alpha=0.3)+
  labs(x="", y="")+
  geom_text(aes(label=round(mean_prop_total_length,5)),family="Open Sans Light",size=8) +
  geom_text(aes(x=factor(tax3,levels=c("Lsta","Lmic","Lang","Lmar","Lnil")), 
                y=factor(tax2,levels=c("Lsta","Lmic","Lang","Lmar","Lnil")),
                label=round(mean_prop_total_length,5)),family="Open Sans Light",size=8) +
  #scale_x_discrete(position = "top") +
  theme_custom() +
  theme(panel.grid = element_blank(),
        legend.position = "top",
        legend.text = element_text(angle=60,hjust=1,size=rel(1.5)),
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'),
        axis.text = element_text(size=rel(2))) + #change legend key width +
  scale_x_discrete(drop=FALSE, labels=c("Lsta","Lmic","Lang","Lmar","Lnil")) +
  scale_y_discrete(drop=FALSE, labels=c("Lsta","Lmic","Lang","Lmar","Lnil"))
