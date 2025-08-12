## plotting parametric bootstrap estimates ----
library(tidyverse)
library(ggridges)
library(viridis)
library(ggdist)
library(patchwork)
library(ggpubr)
source("theme_custom.R")

file <- "~/Dropbox/i/projects/lates_admixture/results/mods_recompute_lhoods_latesWGS_2_082922/bestmod_noMarAng_noStaMar_noMicAng_noNilSta_bs_nonparam_200runs.bestlhoods"
file_alt <- "~/Dropbox/i/projects/lates_admixture/results/mods_recompute_lhoods_newSFS_dp30_091321/fbranch_model_new_altFossil.bootstraps_200runs.bestlhoods"

bs <- read_table(file) %>%
  rownames_to_column(var="bs_num") %>%
  mutate(diffLhood = MaxObsLhood - MaxEstLhood)
  
bs_long <- bs %>%
  pivot_longer(cols=POP0:MaxObsLhood,names_to="param",values_to="est")
bs_CI <- bs_long %>%
  #filter(diffLhood < quantile(diffLhood,0.9)) %>%
  group_by(param) %>%
  summarize(mean = mean(est),
            median = median(est),
            sd = sd(est),
            lcl = mean - 1.96*sd,
            ucl = mean + 1.96*sd,
            q0.025 = quantile(est,probs=0.025),
            q0.975 = quantile(est,probs=0.975),
            lcl_IQR = median - 1.5*(q0.975-q0.025),
            ucl_IQR = median + 1.5*(q0.975-q0.025),
            q0.05 = quantile(est,probs=0.05),
            q0.95 = quantile(est,probs=0.95)
            ) %>%
  mutate(linesty = case_when(grepl("DIV",param) ~ 2,
                             TRUE ~ 3))

bs_long %>%
  ggplot(aes(x=diffLhood,y=est)) +
  geom_jitter(aes(col=param),alpha=0.5) +
  facet_wrap(~param,scales="free_y")+
  #theme_custom() +
  theme(legend.position="none") +
  geom_hline(data=bs_CI,aes(yintercept=q0.05),lty=3) +
  geom_hline(data=bs_CI,aes(yintercept=q0.95),lty=3) +
  geom_hline(data=bs_CI,aes(yintercept=median)) +
  geom_hline(data=bs_CI,aes(yintercept=mean),lty=2)

bs_CI %>%
  filter(mean > 100) %>% # remove mig rates
  ggplot() +
  geom_pointrange(aes(x=param,y=mean,ymin=q0.025,ymax=q0.975)) +
  theme_custom() 

adm_dat <- bs %>%
  select(contains(c("bs_num","diffLhood","DIV","TADM"))) %>%
  #filter(diffLhood < 3000) %>%
  pivot_longer(cols=contains(c("DIV","TADM")),names_to="param",values_to="est") %>%
  filter(param != "TADM13") 

param.adm.plot <- adm_dat %>%
  ggplot(aes(x=reorder(param, est),y=est*3,fill=reorder(param, est))) +
  #geom_density_ridges_gradient(alpha=0.5) +
  ggdist::stat_halfeye(adjust=2,justification=-0.1,point_color=NA,trim=FALSE,.width=0,width=1.5,alpha=0.7) +
  # geom_point(data=filter(adm_dat,diffLhood == min(diffLhood)),aes(y=est*3,x=reorder(param, est)),
  #            size=4,position=position_nudge(x=0.15),pch=4) +
  geom_point(data=filter(bs_CI,param %in% unique(adm_dat$param)),aes(y=median*3,x=param),inherit.aes=FALSE,
             size=4,position=position_nudge(x=0.15),pch=4) +
  #geom_hline(data=filter(adm_dat,diffLhood == min(diffLhood)), aes(yintercept=est*3),lty=3) +
  geom_segment(data=filter(bs_CI,param %in% unique(adm_dat$param)),aes(y=median*3,yend=median*3,x=0,xend=param,lty=as.factor(linesty)),inherit.aes=FALSE) +
  scale_linetype_manual(values=c(2,3)) +
  geom_point(alpha=0.5,size=1,position=position_jitter(seed=1,width=0.05)) +
  #geom_text(data=filter(bs_CI,param %in% unique(adm_dat$param)),aes(x=0,y=median*3,label=round(median*3/1000000,2)),hjust=1,col="gray50",nudge_x=-0.2,inherit.aes=FALSE) +
  #coord_flip() +
  #theme_custom() +
  theme(legend.position="none",
        axis.text.y = element_blank(),
        axis.text.x = element_text(angle=60,hjust=1),
        plot.margin = margin(0.5, 0.5, 0.5, 1, "cm"),
        panel.border = element_blank(),
        axis.ticks = element_blank()) +
  ylab("") +
  xlab("") +
  coord_cartesian(clip="off",xlim=c(0,11),expand=FALSE) +
  theme_custom()

mig_dat <- bs %>%
  select(contains(c("bs_num","diffLhood","MIG"))) %>%
  #filter(diffLhood < 3000) %>%
  pivot_longer(cols=contains(c("MIG")),names_to="param",values_to="est") 

param.mig.plot <- mig_dat %>%
  ggplot(aes(x=reorder(param, est),y=est,fill=reorder(param, est))) +
  #geom_density_ridges_gradient(alpha=0.5) +
  ggdist::stat_halfeye(adjust=1.5,justification=-0.1,point_color=NA,trim=TRUE,.width=0,width=1.5,alpha=0.7) +
  geom_point(data=filter(mig_dat,diffLhood == min(diffLhood)),aes(y=est,x=reorder(param, est)),
             size=4,position=position_nudge(x=0.15),pch=4) +
  #geom_hline(data=filter(adm_dat,diffLhood == min(diffLhood)), aes(yintercept=est*3),lty=3) +
  #geom_segment(data=filter(mig_dat,diffLhood == min(diffLhood)),aes(y=est,yend=est,x=0,xend=reorder(param, est)),lty=3) +
  geom_point(alpha=0.5,size=1,position=position_jitter(seed=1,width=0.05)) +
  #theme_custom() +
  theme(legend.position="none",
        #axis.text.y = element_blank(),
        #axis.text.x = element_text(angle=60,hjust=1),
        plot.margin = margin(0.5, 0.5, 0.5, 1, "cm")) +
  #scale_y_continuous(trans="log") +
  ylab("") +
  xlab("") 

mig_dat %>%
  pivot_wider(id_cols=bs_num:diffLhood, names_from=param, values_from=est) %>%
  mutate(RAT14 = MIG14 / MIG41,
         RAT15 = MIG15 / MIG51,
         RAT12 = MIG12 / MIG21,
         RAT23 = MIG23 / MIG32,
         RAT34 = MIG34 / MIG43,
         RAT25 = MIG25 / MIG52) %>%
  select(bs_num:diffLhood,RAT14:RAT25) %>% 
  pivot_longer(cols=RAT14:RAT25, names_to = "param", values_to = "est") %>%
  filter(est < 100) %>%
  ggplot(aes(x=est,y=param)) +
  stat_summary()

## plotting nonparametric bootstrap estimates ----
np_file <- "~/Dropbox/i/projects/lates_admixture/results/mods_recompute_lhoods_newSFS_dp30_091321/fbranch_model_new_bs_nonparam.bestlhoods"
np_bs <- read_table2(np_file) %>%
  rownames_to_column(var="bs_num") %>%
  mutate(diffLhood = MaxObsLhood - MaxEstLhood)

np_bs_long <- np_bs %>%
  pivot_longer(cols=POP0:MaxObsLhood,names_to="param",values_to="est")
np_bs_CI <- np_bs_long %>%
  #filter(diffLhood < quantile(diffLhood,0.9)) %>%
  group_by(param) %>%
  summarize(mean = mean(est),
            median = median(est),
            sd = sd(est),
            lcl = mean - 1.96*sd,
            ucl = mean + 1.96*sd,
            q0.025 = quantile(est,probs=0.025),
            q0.975 = quantile(est,probs=0.975),
            lcl_IQR = median - 1.5*(q0.975-q0.025),
            ucl_IQR = median + 1.5*(q0.975-q0.025),
            q0.05 = quantile(est,probs=0.05),
            q0.95 = quantile(est,probs=0.95)
  ) 

np_bs_long %>%
  ggplot(aes(x=diffLhood,y=est)) +
  geom_jitter(aes(col=param),alpha=0.5) +
  facet_wrap(~param,scales="free_y")+
  theme_custom() +
  theme(legend.position="none") +
  geom_hline(data=np_bs_CI,aes(yintercept=q0.05),lty=3) +
  geom_hline(data=np_bs_CI,aes(yintercept=q0.95),lty=3) +
  geom_hline(data=np_bs_CI,aes(yintercept=median)) +
  geom_hline(data=np_bs_CI,aes(yintercept=mean),lty=2)

np_bs_CI %>%
  filter(mean > 10) %>% # remove mig rates
  ggplot() +
  geom_pointrange(aes(x=param,y=mean,ymin=q0.025,ymax=q0.975)) +
  theme_custom()

np_adm_dat <- np_bs %>%
  select(contains(c("bs_num","diffLhood","DIV","TADM"))) %>%
  #filter(diffLhood < 200) %>%
  pivot_longer(cols=contains(c("DIV","TADM")),names_to="param",values_to="est") 

nparam.adm.plot <- np_adm_dat %>%
  ggplot(aes(x=reorder(param, est),y=est*3,fill=reorder(param, est))) +
  #geom_density_ridges_gradient(alpha=0.5) +
  ggdist::stat_halfeye(adjust=2,justification=-0.1,point_color=NA,trim=FALSE,.width=0,width=1.5,alpha=0.7) +
  # geom_point(data=filter(adm_dat,diffLhood == min(diffLhood)),aes(y=est*3,x=reorder(param, est)),
  #            size=4,position=position_nudge(x=0.15),pch=4) +
  geom_point(data=filter(np_bs_CI,param %in% unique(np_adm_dat$param)),aes(y=median*3,x=param),inherit.aes=FALSE,
             size=4,position=position_nudge(x=0.15),pch=4) +
  #geom_hline(data=filter(adm_dat,diffLhood == min(diffLhood)), aes(yintercept=est*3),lty=3) +
  geom_segment(data=filter(np_bs_CI,param %in% unique(np_adm_dat$param)),aes(y=median*3,yend=median*3,x=0,xend=param),inherit.aes=FALSE,lty=3) +
  geom_point(alpha=0.5,size=1,position=position_jitter(seed=1,width=0.05)) +
  geom_text(data=filter(np_bs_CI,param %in% unique(np_adm_dat$param)),aes(x=0,y=median*3,label=round(median*3/1000000,2)),hjust=1,col="gray50",nudge_x=-0.2,inherit.aes=FALSE) +
  #coord_flip() +
  theme_custom() +
  theme(legend.position="none",
        axis.text.y = element_blank(),
        axis.text.x = element_text(angle=60,hjust=1),
        plot.margin = margin(0.5, 0.5, 0.5, 1, "cm")) +
  ylab("") +
  xlab("") +
  coord_cartesian(clip="off",xlim=c(0,12),expand=FALSE)

np_mig_dat <- np_bs %>%
  select(contains(c("bs_num","diffLhood","MIG"))) %>%
  #filter(diffLhood < 3000) %>%
  pivot_longer(cols=contains(c("MIG")),names_to="param",values_to="est") 

nparam.mig.plot <- np_mig_dat %>%
  ggplot(aes(x=reorder(param, est),y=est,fill=reorder(param, est))) +
  #geom_density_ridges_gradient(alpha=0.5) +
  ggdist::stat_halfeye(adjust=1.5,justification=-0.1,point_color=NA,trim=TRUE,.width=0,width=1.5,alpha=0.7) +
  geom_point(data=filter(np_mig_dat,diffLhood == min(diffLhood)),aes(y=est,x=reorder(param, est)),
             size=4,position=position_nudge(x=0.15),pch=4) +
  #geom_hline(data=filter(adm_dat,diffLhood == min(diffLhood)), aes(yintercept=est*3),lty=3) +
  #geom_segment(data=filter(mig_dat,diffLhood == min(diffLhood)),aes(y=est,yend=est,x=0,xend=reorder(param, est)),lty=3) +
  geom_point(alpha=0.5,size=1,position=position_jitter(seed=1,width=0.05)) +
  theme_custom() +
  theme(legend.position="none",
        #axis.text.y = element_blank(),
        #axis.text.x = element_text(angle=60,hjust=1),
        plot.margin = margin(0.5, 0.5, 0.5, 1, "cm")) +
  #scale_y_continuous(trans="log") +
  ylab("") +
  xlab("") 

library(patchwork)
ggarrange(param.adm.plot, nparam.adm.plot)

### combined halfeye plot (wish me luck!)
both.adm.plot <- np_adm_dat %>%
  ggplot(aes(x=reorder(param, est),y=est*3/1e6,fill=reorder(param, est))) +
  #geom_density_ridges_gradient(alpha=0.5) +
  ggdist::stat_halfeye(adjust=2,point_color=NA,trim=TRUE,.width=0, width=1.5,alpha=0.9,
                       side="right",justification=-0.01,scale=0.7) +
  ggdist::stat_halfeye(data=adm_dat,aes(x=reorder(param, est),y=est*3/1e6,fill=reorder(param, est)), 
                       adjust=2,point_color=NA,trim=TRUE,.width=0,width=1.5,alpha=0.5,
                       side="left",scale=0.7) +
  # geom_point(data=filter(adm_dat,diffLhood == min(diffLhood)),aes(y=est*3,x=reorder(param, est)),
  #            size=4,position=position_nudge(x=0.15),pch=4) +
  # geom_point(data=filter(np_bs_CI,param %in% unique(np_adm_dat$param)),aes(y=median*3,x=param),inherit.aes=FALSE,
  #            size=4,position=position_nudge(x=0),pch=4) +
  # geom_point(data=filter(bs_CI,param %in% unique(adm_dat$param)),aes(y=median*3,x=param),inherit.aes=FALSE,
  #            size=4,position=position_nudge(x=-0),pch=2,color="gray50") +
  stat_summary(data=adm_dat,aes(x=reorder(param, est),y=est*3/1e6),
               fun="mean",fun.min=function(z) {quantile(z,0.025)},fun.max=function(z) {quantile(z,0.975)},col="gray50") +
  stat_summary(aes(x=reorder(param, est),y=est*3/1e6),
               fun="mean",fun.min=function(z) {quantile(z,0.025)},fun.max=function(z) {quantile(z,0.975)},col="black") +
  #geom_hline(data=filter(adm_dat,diffLhood == min(diffLhood)), aes(yintercept=est*3),lty=3) +
  #geom_segment(data=filter(np_bs_CI,param %in% unique(np_adm_dat$param)),aes(y=median*3,yend=median*3,x=0,xend=param),inherit.aes=FALSE,lty=3) +
  #geom_point(alpha=0.5,size=1,position=position_jitter(seed=1,width=0.05)) +
  geom_text(data=filter(np_bs_CI,param %in% unique(np_adm_dat$param)),
            aes(x=param,y=median*3/1e6,label=round(median*3/1000000,2)),
            hjust=0,col="black",nudge_x=0.1,inherit.aes=FALSE) +
  geom_text(data=filter(bs_CI,param %in% unique(adm_dat$param)),
            aes(x=param,y=median*3/1e6,label=round(median*3/1000000,2)),
            hjust=0,col="gray50",nudge_x=-0.4,inherit.aes=FALSE) +
  annotate(x=-0.2,hjust=1,y=max(np_adm_dat$est*3)/1e6,geom="text",label="Time\n(Mya)",
           family="Open Sans",size=rel(4)) +
  #coord_flip() +
  theme_custom() +
  theme(legend.position="none",
        #axis.text.y = element_blank(),
        axis.text.x = element_text(angle=60,hjust=1),
        plot.margin = margin(0.5, 0.5, 0.5, 1.5, "cm"),
        panel.grid.major = element_line(linetype="dashed")) +
  ylab("") +
  xlab("") +
  coord_cartesian(clip="off",xlim=c(0,12),ylim=c(-100000/1e6,max(np_adm_dat$est*3.1)/1e6),expand=FALSE)
both.adm.plot

#####------------######
### Alt fossil calibration
######-----------#######

bs_alt <- read_table(file_alt) %>%
  rownames_to_column(var="bs_num") %>%
  mutate(diffLhood = MaxObsLhood - MaxEstLhood)


bs_alt_long <- bs_alt %>%
  pivot_longer(cols=POP0:MaxObsLhood,names_to="param",values_to="est")
bs_alt_CI <- bs_alt_long %>%
  #filter(diffLhood < quantile(diffLhood,0.9)) %>%
  group_by(param) %>%
  summarize(mean = mean(est),
            median = median(est),
            sd = sd(est),
            lcl = mean - 1.96*sd,
            ucl = mean + 1.96*sd,
            q0.025 = quantile(est,probs=0.025),
            q0.975 = quantile(est,probs=0.975),
            lcl_IQR = median - 1.5*(q0.975-q0.025),
            ucl_IQR = median + 1.5*(q0.975-q0.025),
            q0.05 = quantile(est,probs=0.05),
            q0.95 = quantile(est,probs=0.95)
  ) %>%
  mutate(linesty = case_when(grepl("DIV",param) ~ 2,
                             TRUE ~ 3))

bs_alt_long %>%
  ggplot(aes(x=diffLhood,y=est)) +
  geom_jitter(aes(col=param),alpha=0.5) +
  facet_wrap(~param,scales="free_y")+
  theme_custom() +
  theme(legend.position="none") +
  geom_hline(data=bs_alt_CI,aes(yintercept=q0.05),lty=3) +
  geom_hline(data=bs_alt_CI,aes(yintercept=q0.95),lty=3) +
  geom_hline(data=bs_alt_CI,aes(yintercept=median)) +
  geom_hline(data=bs_alt_CI,aes(yintercept=mean),lty=2)
  
# plot normal bootstraps and alt bootstraps together
bs_alt_CI %>%
  add_column(model="fbranch_new_altFossil") %>%
  add_row(bs_CI %>% add_column(model="fbranch_new")) %>%
  filter(mean > 100) %>% # remove mig rates
  ggplot() +
  geom_pointrange(aes(x=param,y=mean,ymin=q0.025,ymax=q0.975,col=model),
                  position = position_dodge(width = 0.25),
                  size=1) +
  theme_custom() 

adm_dat_alt <- bs_alt %>%
  select(contains(c("bs_num","diffLhood","DIV","TADM"))) %>%
  #filter(diffLhood < 3000) %>%
  pivot_longer(cols=contains(c("DIV","TADM")),names_to="param",values_to="est") 

param.adm.plot_alt <- adm_dat_alt %>%
  ggplot(aes(x=reorder(param, est),y=est*3,fill=reorder(param, est))) +
  #geom_density_ridges_gradient(alpha=0.5) +
  ggdist::stat_halfeye(adjust=2,justification=-0.1,point_color=NA,trim=FALSE,.width=0,width=1.5,alpha=0.7) +
  # geom_point(data=filter(adm_dat,diffLhood == min(diffLhood)),aes(y=est*3,x=reorder(param, est)),
  #            size=4,position=position_nudge(x=0.15),pch=4) +
  geom_point(data=filter(bs_alt_CI,param %in% unique(adm_dat_alt$param)),aes(y=median*3,x=param),inherit.aes=FALSE,
             size=4,position=position_nudge(x=0.15),pch=4) +
  #geom_hline(data=filter(adm_dat,diffLhood == min(diffLhood)), aes(yintercept=est*3),lty=3) +
  geom_segment(data=filter(bs_alt_CI,param %in% unique(adm_dat_alt$param)),aes(y=median*3,yend=median*3,x=0,xend=param,lty=as.factor(linesty)),inherit.aes=FALSE) +
  scale_linetype_manual(values=c(2,3)) +
  geom_point(alpha=0.5,size=1,position=position_jitter(seed=1,width=0.05)) +
  #geom_text(data=filter(bs_CI,param %in% unique(adm_dat$param)),aes(x=0,y=median*3,label=round(median*3/1000000,2)),hjust=1,col="gray50",nudge_x=-0.2,inherit.aes=FALSE) +
  #coord_flip() +
  theme_custom() +
  theme(legend.position="none",
        axis.text.y = element_blank(),
        axis.text.x = element_text(angle=60,hjust=1),
        plot.margin = margin(0.5, 0.5, 0.5, 1, "cm"),
        panel.border = element_blank(),
        axis.ticks = element_blank()) +
  ylab("") +
  xlab("") +
  coord_cartesian(clip="off",xlim=c(0,11),expand=FALSE)

mig_dat_alt <- bs_alt %>%
  select(contains(c("bs_num","diffLhood","MIG"))) %>%
  #filter(diffLhood < 3000) %>%
  pivot_longer(cols=contains(c("MIG")),names_to="param",values_to="est") 

param.mig.plot_alt <- mig_dat_alt %>%
  ggplot(aes(x=reorder(param, est),y=est,fill=reorder(param, est))) +
  #geom_density_ridges_gradient(alpha=0.5) +
  ggdist::stat_halfeye(adjust=1.5,justification=-0.1,point_color=NA,trim=TRUE,.width=0,width=1.5,alpha=0.7) +
  geom_point(data=filter(mig_dat_alt,diffLhood == min(diffLhood)),aes(y=est,x=reorder(param, est)),
             size=4,position=position_nudge(x=0.15),pch=4) +
  #geom_hline(data=filter(adm_dat,diffLhood == min(diffLhood)), aes(yintercept=est*3),lty=3) +
  #geom_segment(data=filter(mig_dat,diffLhood == min(diffLhood)),aes(y=est,yend=est,x=0,xend=reorder(param, est)),lty=3) +
  geom_point(alpha=0.5,size=1,position=position_jitter(seed=1,width=0.05)) +
  theme_custom() +
  theme(legend.position="none",
        #axis.text.y = element_blank(),
        #axis.text.x = element_text(angle=60,hjust=1),
        plot.margin = margin(0.5, 0.5, 0.5, 1, "cm")) +
  #scale_y_continuous(trans="log") +
  ylab("") +
  xlab("") 

### combined halfeye plot
both.adm.plot <- adm_dat_alt %>%
  ggplot(aes(x=param,y=est*3/1e6,fill=param)) +
  #geom_density_ridges_gradient(alpha=0.5) +
  ggdist::stat_halfeye(adjust=1.5,point_color=NA,trim=TRUE,.width=0, width=1.5,alpha=0.9,
                       side="right",justification=-0.01,scale=0.6) +
  ggdist::stat_halfeye(data=adm_dat,aes(x=param,y=est*3/1e6,fill=param), 
                       adjust=1.5,point_color=NA,trim=TRUE,.width=0,width=1.5,alpha=0.5,
                       side="left",scale=0.6) +
  # geom_point(data=filter(adm_dat,diffLhood == min(diffLhood)),aes(y=est*3,x=param),
  #            size=4,position=position_nudge(x=0.15),pch=4) +
  # geom_point(data=filter(np_bs_CI,param %in% unique(np_adm_dat$param)),aes(y=median*3,x=param),inherit.aes=FALSE,
  #            size=4,position=position_nudge(x=0),pch=4) +
  # geom_point(data=filter(bs_CI,param %in% unique(adm_dat$param)),aes(y=median*3,x=param),inherit.aes=FALSE,
  #            size=4,position=position_nudge(x=-0),pch=2,color="gray50") +
  stat_summary(data=adm_dat,aes(x=param,y=est*3/1e6),
               fun="median",fun.min=function(z) {quantile(z,0.025)},fun.max=function(z) {quantile(z,0.975)},col="gray50",
               position=position_nudge(x=-0.05)) +
  stat_summary(aes(x=param,y=est*3/1e6),
               fun="median",fun.min=function(z) {quantile(z,0.025)},fun.max=function(z) {quantile(z,0.975)},col="black",
               position=position_nudge(x=0.07)) +
  #geom_hline(data=filter(adm_dat,diffLhood == min(diffLhood)), aes(yintercept=est*3),lty=3) +
  #geom_segment(data=filter(np_bs_CI,param %in% unique(np_adm_dat$param)),aes(y=median*3,yend=median*3,x=0,xend=param),inherit.aes=FALSE,lty=3) +
  #geom_point(alpha=0.5,size=1,position=position_jitter(seed=1,width=0.05)) +
  geom_text(data=filter(bs_alt_CI,param %in% unique(adm_dat_alt$param)),
            aes(x=param,y=median*3/1e6,label=round(median*3/1000000,2)),
            hjust=0,col="black",nudge_x=0.12,inherit.aes=FALSE) +
  geom_text(data=filter(bs_CI,param %in% unique(adm_dat$param)),
            aes(x=param,y=median*3/1e6,label=round(median*3/1000000,2)),
            hjust=0,col="gray50",nudge_x=-0.5,inherit.aes=FALSE) +
  annotate(x=-0.2,hjust=1,y=max(adm_dat_alt$est*3)/1e6,geom="text",label="Time\n(Mya)",
           family="Open Sans",size=rel(4)) +
  #coord_flip() +
  theme_custom() +
  theme(legend.position="none",
        #axis.text.y = element_blank(),
        axis.text.x = element_text(angle=60,hjust=1),
        plot.margin = margin(0.5, 0.5, 0.5, 1.5, "cm"),
        panel.grid.major = element_line(linetype="dashed")) +
  ylab("") +
  xlab("") +
  coord_cartesian(clip="off",xlim=c(0,12),ylim=c(-100000/1e6,max(adm_dat$est*3.1)/1e6),expand=FALSE) +
  geom_point(data=altFossil_events %>% group_by(param) %>% slice(1),
             aes(x=param,y=time*3/1e6),
             pch=4,size=3,inherit.aes = FALSE, position=position_nudge(x=0.07))  +
  geom_point(data=fbranch_events %>% group_by(param) %>% slice(1),
             aes(x=param,y=time*3/1e6),
             pch=4,size=3,inherit.aes = FALSE, position=position_nudge(x=-0.05),
             color="gray50")
both.adm.plot


adm_ests_ci <- adm_dat %>%
  mutate(est_yr = 3*est) %>%
  group_by(param) %>%
  summarize(median_est = median(est),
            mean_est = mean(est),
            lcl_95 = quantile(est,0.025),
            ucl_95 = quantile(est,0.975),
            median_est_yr = median(est_yr),
            mean_est_yr = mean(est_yr),
            lcl_95_yr = quantile(est_yr,0.025),
            ucl_95_yr = quantile(est_yr,0.975))

mig_ests_ci <- mig_dat %>%
  #mutate(est_yr = 3*est) %>%
  group_by(param) %>%
  summarize(median_est = median(est),
            mean_est = mean(est),
            lcl_95 = quantile(est,0.025),
            ucl_95 = quantile(est,0.975))
