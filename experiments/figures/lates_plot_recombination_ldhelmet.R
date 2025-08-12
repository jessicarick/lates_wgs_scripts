## plotting recombination rates estimated from ldhelmet
## j.rick, sept 2021

read_recomb <- function(x) {
  chr <- gsub(".*_(ASB_LG.*)\\.run[0-9]\\.post\\.txt","\\1",x)
  recomb <- read_table2(paste0(dir,x),
                        comment="#",skip=3,
                        col_names=c("left_snp", "right_snp", "mean", "p0.025", "p0.500", "p0.0975")) %>%
    add_column(scaffold=chr,.before=1) 
  recomb <- recomb %>%
    mutate(mid=(left_snp+right_snp)/2)
  return(recomb)
}

dir <- "~/Dropbox/i/projects/lates_admixture/results/ldhelmet_may2023/"
files_r1 <- list.files(dir,pattern="output.*run1")
files_r2 <- list.files(dir,pattern="output.*run2")
files_r3 <- list.files(dir,pattern="output.*run3")

recomb_r1 <- files_r1 %>%
  map_df(~read_recomb(.))
recomb_r2 <- files_r2 %>%
  map_df(~read_recomb(.))
recomb_r3 <- files_r3 %>%
  map_df(~read_recomb(.))

recomb_r1 %>%
  ggplot(aes(x=mid,y=mean)) +
  geom_point(alpha=0.5) +
  facet_wrap(~scaffold,ncol=4,scales="free_x")
recomb_r2 %>%
  ggplot(aes(x=mid,y=mean)) +
  geom_point(alpha=0.5) +
  facet_wrap(~scaffold,ncol=4,scales="free_x")

recomb_all <- recomb_r1 %>%
  full_join(recomb_r2,by=c("scaffold", "left_snp", "right_snp", "mid"),suffix=c(".r1",".r2")) %>%
  full_join(recomb_r3,by=c("scaffold", "left_snp", "right_snp", "mid"),suffix=c("",".r3"))

## how well do estimates from using different block sizes match up?
recomb_all %>%
  ggplot(aes(x=mean.r1,y=mean.r2)) +
  geom_point(alpha=0.5) +
  geom_point(aes(x=mean.r1,y=mean.r2),col="darkturquoise",alpha=0.5) +
  #geom_hex(bins=200) +
  theme_custom()
for (chr in unique(recomb_all$scaffold)){
  p1 <- recomb_all %>%
    filter(scaffold == chr) %>%
    ggplot(aes(x=mid)) +
    # geom_point(aes(y=mean.r1),alpha=0.5) +
    # geom_point(aes(y=mean.r2),col="violetred4",alpha=0.5) +
    geom_point(aes(y=mean),col="darkturquoise",alpha=0.5) +
    #stat_summary_bin(aes(y=mean.r1),bins=500,geom="area",fill="gray50",alpha=0.5) +
    #stat_summary_bin(aes(y=mean.r2),bins=500,geom="area",fill="violetred4",alpha=0.5) +
    stat_summary_bin(aes(y=-mean),bins=500,geom="area",fill="darkturquoise",alpha=0.5) +
    theme_custom()
  print(p1)
}

recomb_sum <- recomb_all %>%
  summarize(mean_r1 = mean(mean.r1),
            max_r1 = max(mean.r1),
            mean_r2 = mean(mean.r2),
            max_r2 = max(mean.r2),
            mean_r3 = mean(mean),
            max_r3 = max(mean),
            overall_mean = mean(c(mean.r1,mean.r2,mean),na.rm=TRUE),
            q0.025 = quantile(c(mean.r1,mean.r2,mean),0.025,na.rm=TRUE),
            q0.975 = quantile(c(mean.r1,mean.r2,mean),0.975,na.rm=TRUE))
recomb_sum %>%
  ggplot(aes(x=scaffold)) +
  geom_point(aes(y=mean_r1),pch=19) +
  geom_point(aes(y=mean_r2),pch=3) +
  geom_point(aes(y=mean_r3),pch=6) +
  geom_point(aes(y=overall_mean,col=scaffold),size=4,alpha=0.5) +
  #geom_linerange(aes(ymin=q0.025,ymax=q0.975,x=scaffold)) +
  theme_custom() +
  theme(legend.position="none")
  
recomb <- recomb_r1
colnames(recomb)[colnames(recomb) == "scaffold"] <- "chr"
recomb %>%
  #group_by()
  ggplot(aes(x=mid,y=mean*2)) +
  #geom_point(color="gray80") +
  stat_summary_bin(bins=300,fill="gray50",geom="area") +
  theme_custom() +
  ylab("Mean Estimated Recomb Rate") +
  stat_summary_bin(bins=300,data=dstats[dstats$chr %in% unique(recomb$chr),], aes(x=windowStart,y=-f_d_corr,fill=chr),
                   inherit.aes=FALSE,geom="area") +
  geom_hline(yintercept=0) +
  facet_wrap(~chr, ncol=4, scales="free") +
  scale_fill_manual(values = rep(colors, length(unique(recomb$chr)))) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(n.breaks=8, labels=abs) +
  theme(legend.position = "none")

recomb %>%
  group_by(chr) %>%
  summarize(mean_mean = mean(mean),
            windows = n(),
            median = median(mean),
            q0.025 = quantile(mean,0.025),
            q0.975 = quantile(mean,0.975),
            min = min(mean),
            max = max(mean),
            sd = sd(mean,na.rm=TRUE)) %>%
  arrange(mean_mean) %>%
  ggplot() +
  #geom_point(data=dstat_sum,aes(x=chr,y=mean),col="turquoise") +
  geom_boxplot(data=dstat_sum,aes(x=chr,y=mean),alpha=0.5) +
  #geom_point(aes(x=chr,y=mean_mean),size=3,col=colors[3]) +
  #geom_line(aes(x=chr,y=mean*6,group="none"),lty=2,col=colors[3]) +
  geom_pointrange(aes(x=chr,y=mean_mean*2,ymin=q0.025*2, ymax=q0.975*2),fatten=5,col=colors[2],position=position_nudge(x=0.2)) +
  theme_custom() +
  theme(axis.text.x = element_text(angle=60,hjust=1),
        axis.title.x = element_blank()) +
  scale_y_continuous(name = "estimated recombination rate",
                     sec.axis=sec_axis(trans=~./2,name="corrected f_d statistic"))
  
rel_pos_recomb <- recomb %>%
  group_by(chr) %>%
  mutate(max = max(right_snp)) %>%
  ungroup() %>%
  mutate(rel_mid = mid/max) 
rel_pos_dstat <- dstats %>%
  group_by(chr) %>%
  mutate(max = max(windowEnd)) %>%
  ungroup() %>%
  mutate(rel_mid = ((windowStart+windowEnd)/2)/max)

## binned recombination rate vs f_d plot
## panel B for figure 3
fig3b <- ggplot() +
  #geom_point(aes(col=chr),alpha=0.3) +
  #stat_summary_bin(aes(col=chr),alpha=0.5, bins=12) +
  stat_summary_bin(data=rel_pos_recomb, aes(x=rel_mid,y=mean), geom="line",lty=2,alpha=0.5, binwidth=0.05, fun="mean") +
  stat_summary_bin(data=rel_pos_recomb, aes(x=rel_mid,y=mean), binwidth=0.05,size=1.5, fun.data="mean_cl_normal") +
  stat_summary_bin(data=rel_pos_dstat,aes(x=rel_mid,y=f_d_corr/2), geom="line",lty=2,alpha=0.5,binwidth=0.05,col=colors[5], fun="mean") +
  stat_summary_bin(data=rel_pos_dstat,aes(x=rel_mid,y=f_d_corr/2), binwidth=0.05,size=1.5,col=colors[5], fun.data="mean_cl_normal") +
  theme_custom() +
  theme(legend.position="top",
        axis.title.y = element_text(color = "black"),
        axis.title.y.right = element_text(color = colors[5]),
        axis.text.y.right = element_text(colour = colors[5])) +
  xlab("relative chromosomal position") +
  ylab("estimated mean recombination rate (1/bp)") +
  scale_y_continuous(sec.axis=sec_axis(trans=~.*2,name="corrected f_d statistic"))

## combine with fd and twisst results
## requires "dstats" dataframe from plot_dstats_fd_sliding_windows.R
## and "combined_weights" dataframe from lates_wgs_twisst_sept2021.R
## note! will need to fix chr names to match
all_binned <- tibble()  
for (i in unique(recomb$chr)){
  max_bp <- max(max(recomb %>% select(right_snp)),
                max(dstats %>% filter(chr == i) %>% select(windowEnd)))
  binned_dstats <- dstats %>%
    filter(chr == i) %>%
    mutate(windowMid = (windowStart + windowEnd)/2) %>%
    group_by(ints = cut(windowMid, breaks = seq(0, max_bp, 5000))) %>%
    summarise(mean = mean(f_d_corr),
              median = median(f_d_corr),
              n=n())
  binned_recomb <- recomb %>%
    filter(chr == i) %>%
    group_by(ints = cut(mid, breaks = seq(0, max_bp, 5000))) %>%
    summarise(mean = mean(mean),
              median = median(mean),
              n=n())
  binned_topo13 <- combined_weights %>%
    filter(scaffold == i & topo == "topo13") %>%
    group_by(ints = cut(mid, breaks = seq(0, max_bp, 5000))) %>%
    summarise(mean = mean(weighting),
              median = median(weighting),
              n=n())
  binned_combined <- binned_dstats %>%
    full_join(binned_recomb,by="ints",suffix=c(".fd",".recomb")) %>%
    full_join(binned_topo13,by="ints",suffix=c("",".twisst")) %>%
    mutate(chr = i)
  all_binned <- rbind(all_binned,binned_combined)
}

## panel C for figure 3
fig3c <- all_binned %>%
  rename(mean.twisst = mean, median.twisst = median, n.twisst = n) %>% drop_na() %>%
  filter(mean.fd > 0 & mean.recomb > 0 & mean.recomb < 1) %>% 
  ggplot(aes(x=mean.recomb,y=mean.fd)) +
  geom_point(aes(col=mean.recomb),color=colors[2],alpha=0.2,size=3) +
  #geom_hex(bins=100) +
  #geom_smooth(method="gam") +
  geom_smooth(method="lm",color="black",lty=2, se = FALSE,lwd=3) +
  #scale_color_viridis_c(direction = -1,end=0.9,option="plasma",values=c(0,0.0001,0.001,0.01,0.1,0.2,0.5,1.0))+
  #scale_color_gradient(low="gray10",high="black") +
  theme_custom() +
  theme(legend.position = "none") +
  xlab("log Estimated Recombination Rate (1/bp)") +
  ylab("log Estimated f_d Statistic") +
  scale_x_continuous(trans="log10") +
  scale_y_continuous(trans="log10") +
  annotation_logticks(sides="bl")

all_binned %>%
  filter(mean.fd > 0 & mean.recomb > 0 & mean.recomb < 1) %>% 
  mutate(log_recomb = log10(mean.recomb),log_fd=log10(mean.fd)) %>%
  ggscatterhist("log_recomb","log_fd",log="xy",alpha=0.1,margin.plot="density") 

all_binned %>%
  #filter(mean.x > 0 & mean.y > 0) %>%
  ggplot(aes(x=mean.recomb, y=mean.fd) ) +
  geom_hex(bins = 50) +
  #scale_fill_continuous(type = "magma") +
  scale_fill_distiller(palette= "Spectral", direction=-1) +
  theme_custom() +
  #scale_x_continuous(trans="log10") +
  #scale_y_continuous(trans="log10") +
  #geom_smooth(method="lm",col="black",lty=2) +
  geom_smooth(method="gam",formula = y ~ s(log(x)),col="black",lty=2) +
  xlim(0,1) + ylim(0,0.5) +
  xlab("Estimated Recombination Rate (1/bp)") +
  ylab("Estimated f_d Statistic") +
  ggpmisc::stat_poly_eq(formula = y ~ log(x),
                        eq.with.lhs = "italic(y)~`=`~",
                        aes(label = paste(..eq.label..,..rr.label..,..p.value.label..,sep = "*plain(\",\")~")),
                        parse = TRUE, label.x = 0.9, label.y=0.9, size=5, hjust=1)

all_binned %>%
  filter(chr == "ASB_LG14") %>%
  ggplot() +
  geom_point(aes(x=ints,y=mean.fd), color="lightblue", alpha=0.5) +
  geom_point(aes(x=ints,y=-median.recomb), color="darkred", alpha=0.5) +
  theme_custom() + facet_wrap(~chr, scales="free") +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank())
           
library(patchwork)
(fig3a + fig3b + fig3c) / fig3d
