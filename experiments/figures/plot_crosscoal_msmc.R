## script for plotting msmc/msmc-im cross-coalescence
library(tidyverse)
library(ggrepel)
library(ggpubr)
source("theme_custom.R")

colors <- c(rev(c("#976ab0","#de5f7d","#f7a76a","#ffd166","#83d483","#078879","#3a63a4",
                           "darkblue")),"gray")

#### MSMC-IM #####
dir <- "~/Dropbox/i/projects/lates_admixture/results/msmc_may2023"
files <- list.files(dir, pattern="estimates.txt")[-1]

read_crosscoal <- function(x) {
  ind_pair <- str_split(gsub("msmc_([A-Z]+.*)_09[0-9]+21.*","\\1",x),"\\.",simplify=TRUE)
  dat <- read_table2(paste0(dir,x)) %>%
    add_column(ind1 = ind_pair[1],ind2 = ind_pair[2],ind_pair=paste(ind_pair[1],ind_pair[2],sep="_"), .before=1)
  est_split <- read_lines(paste0(dir,gsub("estimates","fittingdetails",x)),n_max=1) %>%
    gsub(".*\\[(.*)\\].*","\\1",.) %>%
    gsub("\\s","",.) %>%
    str_split(",", simplify=TRUE)
  dat <- dat %>%
    mutate(est_split_0.25 = est_split[1],
           est_split_0.5 = est_split[2],
           est_split_0.75 = est_split[3])
  return(dat)
}

spp_ind <- tibble(ind=c("CEW16_022",    "CEW17_120",    "CEW17_127",    "CEW17_151",    "JAR18_139370", "JAR18_139523",
                        "JAR18_139895", "JAR18_139902", "JAR18_139933", "JAR18_140014", "JAR18_140015", "JAR18_140016", 
                        "JAR18_140024","AMNH216682","SRR3140997","CEW19_134","CEW19_135"),
                  spp=c("Lmic","Lang","Lsta","Lsta",
                        "Lmic","Lang","Lmar","Lang","Lmar","Lmar",
                        "Lmic","Lmic","Lsta","Lnil","Lcal","Lmic","Lang"))

# ibd_length <- ibd %>%
#   left_join(spp_ind,by=c("ind1"="ind")) %>%
#   left_join(spp_ind,by=c("ind2"="ind"),suffix=c("1","2")) %>%
#   mutate(length=end-start,
#          ind_pair=paste0(ind1,"-",ind2),
#          spp_pair=paste0(spp1,"-",spp2)) %>%
#   group_by(ind_pair,chrom,spp_pair) %>%
#   summarise(total_length=sum(length))

all_crosscoal <- files %>%
  map_df(~read_crosscoal(.))

all_dat <- all_crosscoal %>%
  left_join(spp_ind,by=c("ind1"="ind")) %>%
  left_join(spp_ind,by=c("ind2"="ind"),suffix=c("1","2")) %>%
  mutate(spp_pair = paste(spp1,spp2,sep="-")) %>%
  mutate(spp_pair2 = case_when(spp_pair == "Lmic-Lmar" ~ "Lmar-Lmic",
                               spp_pair == "Lmic-Lang" ~ "Lang-Lmic",
                               spp_pair == "Lmar-Lang" ~ "Lang-Lmar",
                               TRUE ~ spp_pair)) %>%
  mutate(class = case_when(spp1 == spp2 ~ "intra_spp",
                           spp1 == "Lnil" ~ "lnil_LT",
                           (spp1 == "Lmar" | spp2 == "Lmar") & (spp1 == "Lmic" | spp2 == "Lmic") ~ "Lmar-Lmic",
                           (spp1 == "Lmar" | spp2 == "Lmar") & (spp1 == "Lang" | spp2 == "Lang") ~ "Lang-Lmar",
                           (spp1 == "Lmar" | spp2 == "Lmar") & (spp1 == "Lsta" | spp2 == "Lsta") ~ "Lmar-Lsta",
                           (spp1 == "Lang" | spp2 == "Lang") & (spp1 == "Lmic" | spp2 == "Lmic") ~ "Lang-Lmic",
                           (spp1 == "Lsta" | spp2 == "Lsta") & (spp1 == "Lmic" | spp2 == "Lmic") ~ "Lmic-Lsta",
                           (spp1 == "Lang" | spp2 == "Lang") & (spp1 == "Lsta" | spp2 == "Lsta") ~ "Lang-Lsta",
                           #spp1 != "Lnil" & spp1 != spp2 ~ "intra_LT",
                           TRUE ~ "other"))
lnil_div <- all_dat %>% 
  filter(left_time_boundary > 0) %>%
  filter(class == "lnil_LT") %>%
  ggplot(aes(x=left_time_boundary*3,y=M,col=spp_pair)) +
  geom_line(aes(group=ind_pair)) +
  scale_x_continuous(trans="log10") +
  theme_custom() +
  geom_hline(yintercept=0.5,lty=3) +
  geom_text(data=filter(all_dat,left_time_boundary == 0 & spp1 == "Lnil"),inherit.aes=FALSE,
            aes(x=as.numeric(est_split_0.5)*3,y=0.5,label=round(as.numeric(est_split_0.5)*3,0),col=spp_pair),
            angle=85, nudge_x = -0.04,nudge_y = 0.05) +
  ylab("M(t)") +
  xlab("Time (years before present)") +
  theme(legend.position=c(0.9,0.2))
lnil_div

intra_spp_div <- all_dat %>%
  filter(class == "intra_spp") %>%
  ggplot(aes(x=left_time_boundary*3,y=M,col=spp_pair)) +
  geom_line(aes(group=ind_pair)) +
  scale_x_continuous(trans="log10") +
  theme_custom() +
  geom_hline(yintercept=0.5,lty=3) +
  ylab("M(t)") +
  xlab("Time (years before present)") +
  theme(legend.position=c(0.9,0.2))
intra_spp_div

intra_lt_div <- all_dat %>%
  filter(class == "intra_LT") %>%
  ggplot(aes(x=left_time_boundary*3,y=M,col=spp_pair)) +
  geom_line(aes(group=ind_pair)) +
  scale_x_continuous(trans="log10") +
  theme_custom() +
  geom_hline(yintercept=0.5,lty=3) +
  ylab("M(t)") +
  xlab("Time (years before present)") +
  theme(legend.position=c(0.9,0.2))
intra_lt_div

mean_dat <- all_dat %>%
  filter(as.numeric(est_split_0.5) > 1000) %>%
  group_by(class) %>%
  summarize(mean=mean(as.numeric(est_split_0.5),na.rm=T),
            median=median(as.numeric(est_split_0.5)),
            mean_0.25 = mean(as.numeric(est_split_0.25),na.rm=T),
            mean_0.75 = mean(as.numeric(est_split_0.75),na.rm=T))

all_div <- all_dat %>%
  filter(class != "intra_spp") %>%
  ggplot(aes(x=left_time_boundary*3,y=M,col=class)) +
  geom_line(aes(group=ind_pair),alpha=0.8) +
  scale_x_continuous(trans="log10",limits=c(1e4,1.3e6)) +
  theme_custom() +
  theme(strip.text = element_text(size=14),
        legend.position=c(0.2,0.2)) +
  geom_hline(yintercept=0.5,lty=3) +
  #geom_text(data=filter(mean_dat,class!="intra_spp"),aes(x=1e7,y=0.5,label=round(mean,0)),nudge_y=0.025) +
  #geom_text(data=filter(mean_dat,class!="intra_spp"),aes(x=1e7,y=0.5,label=paste0("(",round(mean_0.25,0),"-",round(mean_0.75,0),")")),nudge_y=-0.025,size=3) +
  #facet_wrap(~class,ncol=3) +
  ylab("M(t)") +
  xlab("Time (years before present)") +
  annotation_logticks(sides="b")
all_div

all_dat %>%
  filter(class != "intra_spp") %>%
  ggplot(aes(x=left_time_boundary*3,y=m,col=class)) +
  geom_step(aes(group=ind_pair),alpha=0.8) +
  scale_x_continuous(trans="log10",limits=c(1e4,1e7)) +
  theme_custom() +
  theme(strip.text = element_text(size=14),
        legend.position="none") +
  #geom_hline(yintercept=0.5,lty=3) +
  #geom_text(data=filter(mean_dat,class!="intra_spp"),aes(x=1e7,y=0.5,label=round(mean,0)),nudge_y=0.025) +
  #geom_text(data=filter(mean_dat,class!="intra_spp"),aes(x=1e7,y=0.5,label=paste0("(",round(mean_0.25,0),"-",round(mean_0.75,0),")")),nudge_y=-0.025,size=3) +
  #facet_wrap(~class,ncol=3) +
  ylab("m(t)") +
  xlab("Time (years before present)") +
  coord_cartesian() +
  annotation_logticks(sides="b") +
  facet_wrap(~spp_pair,scales="free_y")

## attempting to use spline like d.marques
eval_spline <- function (x,p=0.5) {
  dat <- data.frame(x=x$left_time_boundary,y=x$M)
  if (dat$y[1] > 0.2) {
    return(NA)
  } else {
    colnames(dat) <- c("x","y")
    spline <- smooth.spline(dat,spar=0.3)
    m <- uniroot(function(x) predict(spline, x, deriv = 0)$y - p,
                 interval = c(10, 1e6))$root
    return(m)
  }
}

roots <- data.frame(ind_pair = unique(all_dat$ind_pair),
                    m50 = numeric(length(unique(all_dat$ind_pair))),
                    m99 = numeric(length(unique(all_dat$ind_pair))))
for (i in unique(all_dat$ind_pair)) {
  m50 = eval_spline(all_dat[all_dat$ind_pair == i,],0.5)
  m99 = eval_spline(all_dat[all_dat$ind_pair == i,],0.99)
  roots[roots$ind_pair == i,] <- c(i,m50,m99)
}

all_dat %>%
  left_join(roots,by="ind_pair") %>%
  group_by(ind_pair,spp_pair2,class) %>%
  summarize(est_split_0.5 = mean(as.numeric(est_split_0.5)),
            m50 = mean(as.numeric(m50)),
            m99 = mean(as.numeric(m99))) %>%
  ggplot() +
  geom_boxplot(aes(x=class,y=m50,fill=class),alpha=0.5) +
  geom_jitter(aes(x=class,y=m50),alpha=0.5,width=0.1,height=1) +
  theme_custom() +
  theme(legend.position="none")

#######
#### MSMC CROSSCOAL ####
#######

dir <- "~/Dropbox/i/projects/lates_admixture/results/msmc_may2023/"
files_crosscoal <- list.files(dir,"crosscoal.final.txt")[-1]
read_crosscoal2 <- function(x) {
  ind_pair <- str_split(gsub(".*msmc_latesWGS_2.bamlist_([A-Z]+[0-9]+.*)_050[0-9]23.*","\\1",x),"\\.",simplify=TRUE)
  dat <- read_table2(paste0(dir,x)) %>%
    add_column(ind1 = ind_pair[1],ind2 = ind_pair[2],ind_pair=paste(ind_pair[1],ind_pair[2],sep="_"), .before=1) %>%
    add_column(p_par = gsub(".*_06[0-9]+22_([0-9]+seg)\\.crosscoal.*","\\1",x))
  return(dat)
}

dat_msmc <- files_crosscoal %>%
  map_df(~read_crosscoal2(.)) %>%
  filter(ind1 != "CEW17_120" & ind2 != "CEW17_120") %>% 
  mutate(crosscoal = 2*lambda_01/(lambda_00+lambda_11)) %>%
  left_join(spp_ind,by=c("ind1"="ind")) %>%
  left_join(spp_ind,by=c("ind2"="ind"),suffix=c("1","2")) %>%
  mutate(spp_pair = paste(spp1,spp2,sep="-")) %>%
  mutate(spp_pair2 = case_when(spp_pair == "Lmic-Lmar" ~ "Lmar-Lmic",
                               spp_pair == "Lmic-Lang" ~ "Lang-Lmic",
                               spp_pair == "Lmar-Lang" ~ "Lang-Lmar",
                               spp_pair == "Lsta-Lang" ~ "Lang-Lsta",
                               spp_pair == "Lsta-Lmic" ~ "Lmic-Lsta",
                               spp_pair == "Lsta-Lmar" ~ "Lmar-Lsta",
                               TRUE ~ spp_pair)) %>%
  mutate(class = case_when(spp1 == spp2 ~ "intra_spp",
                           spp1 == "Lnil" ~ "lnil_LT",
                           spp2 == "Lcal" & spp1 == "Lnil" ~ "lnil_lcal",
                           spp2 == "Lcal" & spp1 != "Lnil" ~ "lcal_LT",
                           (spp1 == "Lmar" | spp2 == "Lmar") & (spp1 == "Lmic" | spp2 == "Lmic") ~ "Lmar-Lmic",
                           (spp1 == "Lmar" | spp2 == "Lmar") & (spp1 == "Lang" | spp2 == "Lang") ~ "Lang-Lmar",
                           (spp1 == "Lmar" | spp2 == "Lmar") & (spp1 == "Lsta" | spp2 == "Lsta") ~ "Lmar-Lsta",
                           (spp1 == "Lang" | spp2 == "Lang") & (spp1 == "Lmic" | spp2 == "Lmic") ~ "Lang-Lmic",
                           (spp1 == "Lsta" | spp2 == "Lsta") & (spp1 == "Lmic" | spp2 == "Lmic") ~ "Lmic-Lsta",
                           (spp1 == "Lang" | spp2 == "Lang") & (spp1 == "Lsta" | spp2 == "Lsta") ~ "Lang-Lsta",
                           #spp1 != "Lnil" & spp1 != spp2 ~ "intra_LT",
                           TRUE ~ "other")) %>%
  mutate(class = factor(class,levels=c("Lang-Lmar","Lang-Lmic","Lang-Lsta","Lmar-Lmic","Lmar-Lsta",
                                       "Lmic-Lsta","lnil_LT","lcal_LT","lnil_lcal","intra_spp")))
mu <- 3.5e-9

cc_timing <- dat_msmc %>%
  #filter(p_par == "30seg") %>%
  mutate(right_time_boundary=as.numeric(right_time_boundary),
         left_time_boundary=as.numeric(left_time_boundary),
         left_years = (left_time_boundary/mu)*3,
         right_years = (right_time_boundary/mu)*3,
         mid_years = (left_years + right_years) / 2) %>%
  group_by(ind_pair,p_par) %>%
  filter(crosscoal >= 0.5 & left_years > 1e3) %>%
  slice(1) 
  
plot_cc <- dat_msmc %>%
  #filter(p_par == "30seg") %>%
  mutate(left_years = (left_time_boundary/mu)*3) %>%
  filter(class != "intra_spp") %>%
  ggplot(aes(x=left_years,y=crosscoal,col=class)) +
  #annotate("rect",xmin=5.8e5,xmax=7.8e5,ymin=-2,ymax=2,fill="gray80",alpha=0.5) +
  #annotate("rect",xmin=9.9e5,xmax=1.35e6,ymin=-2,ymax=2,fill="gray80",alpha=0.5) +
  geom_step(aes(group=interaction(ind_pair,p_par)),alpha=0.8) +
  theme_custom() +
  scale_x_continuous(trans="log10") +
  theme(legend.position="right") + 
  coord_cartesian(ylim=c(0,2),clip="on",xlim=c(1000,1e7)) +
  geom_hline(yintercept = 0.5, lty=2) +
  annotation_logticks(sides="b",outside=FALSE,scaled=TRUE) +
  ylab("Cross-Coalescence Rate") +
  xlab("Time (years before present)") +
  #theme(legend.position=c(0.2,0.8)) +
  #scale_color_manual(values=colors) +
  guides(color = guide_legend(override.aes = list(size=2))) +
  #geom_point(data=cc_timing,aes(x=left_years,y=0.5,col=class),alpha=0.4,size=4) +
  facet_wrap(~class) +
  geom_vline(xintercept=1e5,col="lightgray")

plot_ne <- dat_msmc %>%
  #filter(p_par == "30seg") %>%
  filter(ind1 != "CEW17_120" & ind2 != "CEW17_120") %>% 
  mutate(left_years = (left_time_boundary/mu)*3) %>%
  #filter(class != "intra_spp") %>%
  ggplot(aes(x=left_years,y=(1/lambda_00)/(2*mu))) +
  annotate("rect",xmin=5.8e5,xmax=7.8e5,ymin=-2,ymax=2,fill="gray80",alpha=0.5) +
  annotate("rect",xmin=9.9e5,xmax=1.35e6,ymin=-2,ymax=2,fill="gray80",alpha=0.5) +
  geom_step(aes(group=ind_pair,col=ind1),alpha=0.8) +
  geom_step(aes(y=(1/lambda_11)/(2*mu),group=ind_pair,col=ind2),alpha=0.8) +
  theme_custom() +
  scale_x_continuous(trans="log10") +
  theme(legend.position="right") + 
  coord_cartesian(ylim=c(0,5e5),clip="on",xlim=c(1e4,5e6)) +
  annotation_logticks(sides="b",outside=FALSE,scaled=TRUE) +
  ylab("Ne") +
  xlab("Time (years before present)") 

ggarrange(plot_cc,all_div,nrow=1)

dat_msmc %>%
  mutate(left_years = (as.numeric(left_time_boundary)/mu)*3,
         right_years = (as.numeric(right_time_boundary)/mu)*3) %>%
  group_by(spp_pair2) %>%
  filter(crosscoal >= 0.5 & left_years > 1e5) %>%
  slice(1) %>% 
  group_by(spp_pair2,class) %>%
  summarize(mean_left_yr = mean(left_years),
            mean_right_yr = mean(right_years)) %>% 
  group_by(class) %>%
  summarize(mean_mean_left_yr = mean(mean_left_yr),
            mean_mean_right_yr = mean(mean_right_yr))

order <- c("Lang-Lmar", "Lmic-Lsta",
           "Lang-Lmic","Lmar-Lmic","Lang-Lsta","Lmar-Lsta",
           "Lnil-Lang","Lnil-Lmar","Lnil-Lmic","Lnil-Lsta")
plot_cc_div <- dat_msmc %>%
  #filter(p_par == "30seg" & class != "intra_spp") %>%
  filter(class != "intra_spp" & spp2 != "Lcal") %>% 
  mutate(left_years = (as.numeric(left_time_boundary)/mu)*3,
         right_years = (as.numeric(right_time_boundary)/mu)*3) %>%
  group_by(ind_pair) %>%
  filter(crosscoal >= 0.5 & left_years > 1e5) %>%
  mutate(spp_pair2 = factor(spp_pair2, levels=order)) %>%
  slice(1) %>% # select first range where crosscoal reaches 0.5
  ggplot(aes(x=spp_pair2)) +
    #geom_segment(aes(x=spp_pair2,xend=spp_pair2,
    #                  y=left_years/1e6,yend=right_years/1e6), size=1.5,col="gray50") +
    geom_point(aes(x=spp_pair2, y=(left_years+right_years)/2/1e6,
                   col=class,fill=class),
               size=3, alpha=0.7) +
    stat_summary(aes(x=spp_pair2,y=(left_years+right_years)/2/1e6,
                     group=spp_pair2),
                 shape=17,size=4,geom="point") +
    theme_custom() + 
    theme(legend.position="none",
          panel.grid.major = element_line(linetype="dashed")) +
    ylab(expression(Cross-coalescence~T[50]~(Mya))) +
    xlab("Species Pair") +
    scale_color_manual(values=colors) +
    coord_flip() + 
    geom_rect(data=adm_ests_ci %>% filter(!grepl("TADM",param)), 
              aes(ymax=ucl_95_yr/1e6,ymin=lcl_95_yr/1e6, 
                                    xmin=c(0.5,1.5,2.5,6.5),xmax=c(1.5,2.5,6.5,10.5)),
              alpha=0.1,inherit.aes=FALSE) +
    geom_segment(data=adm_ests_ci %>% filter(!grepl("TADM",param)), 
                 aes(y=median_est_yr/1e6,yend=median_est_yr/1e6,
                     x=c(0.5,1.5,2.5,6.5),xend=c(1.5,2.5,6.5,10.5)),
                 lwd=2,col="darkgray",alpha=0.5)

ggarrange(plot_cc,plot_cc_div,nrow=1,
          widths = c(2,1))
      
