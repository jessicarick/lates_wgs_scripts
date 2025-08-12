## script for plotting fastsimcoal likelihood distributions
## J. Rick, Jan 2021

library(tidyverse)
library(ggridges)
library(viridis)
library(ggpubr)
library(ggdist)

read_and_add <- function(flnm) {
  read_csv(flnm, col_names="lhood") %>%
    mutate(model = gsub(".*mods_recompute_lhoods_twopop_[a-z]+\\.[a-z]+/(twopop.*)\\.lhoods","\\1",flnm),
           spp = gsub(".*mods_recompute_lhoods_twopop_([a-z]+\\.[a-z]+)/twopop.*\\.lhoods","\\1",flnm))
}

folder <- paste0("~/Dropbox/i/projects/lates_admixture/results/mods_recompute_lhoods_twopop")

lhoods <- list.files(folder,full.names=TRUE,recursive=TRUE,pattern="\\.lhoods") %>% 
  map_df(~read_and_add(.)) %>%
#  mutate(model_lab = case_when(model == "noadmix" ~ "No geneflow",
#                               model == "lnil-lt_admix" ~ "Single shared admixture Lnil->LT",
#                               model == "lnil-all_admix_lt_admix" ~ "Lnil-LT admixture + LT-LT admixture",
#                               model == "lnil-all_admix_asymm" ~ "Bidirectional Lnil-LT admixture",
#                               model == "lnil-all_admix" ~ "Unidirectional Lnil-LT admixture")) %>%
  mutate(model = as.factor(model),
         spp = as.factor(spp))

#new_folder <- paste0("C:/Users/jrick/Dropbox/i/projects/lates_admixture/results/mods_recompute_lhoods_sixpop2")

# lhoods_new <- list.files(new_folder,full.names=TRUE,recursive=TRUE,pattern="\\.lhoods") %>% 
#   map_df(~read_and_add(.)) %>%
#   mutate(model = as.factor(model))

#maxEstLhoods <- read_csv("~/Dropbox/i/projects/lates_admixture/results/mods_recompute_lhoods_sixpop/maxEstLhoods.csv")
lhoodMax <- lhoods %>%
  filter(!grepl("gf",model)) %>%
  group_by(spp) %>%
  arrange(-lhood) %>%
  slice(1)
lhoodMaxMedian <- lhoods %>%
  filter(!grepl("gf",model)) %>%
  group_by(model,spp) %>%
  summarize(median=median(lhood)) %>%
  group_by(spp) %>%
  arrange(-median) %>%
  slice(1)
lhoodMaxMean <- lhoods %>%
  filter(!grepl("gf",model)) %>%
  group_by(model,spp) %>%
  summarize(mean=mean(lhood)) %>%
  group_by(spp) %>%
  arrange(-mean) %>%
  slice(1)
density <- lhoods %>%
    #filter(spp %in% c("lang.lnil","lmar.lnil","lmic.lnil","lsta.lnil")) %>%
    filter(!grepl("gf",model)) %>%
    ggplot(aes(x=-lhood,y=model, fill=..x..)) +
    geom_density_ridges_gradient(scale = 2, rel_min_height = 0.01) +
    theme_ridges() +
    scale_fill_viridis(name="Likelihood", option="C") +
    facet_wrap(~spp,scales="free_x") +
    theme(axis.title = element_blank(),
          legend.position = "none",
          axis.text.x=element_text(size=12),
          panel.grid.major.x = element_blank()) 
density_plot_lnil <- lhoods %>%
  filter(grepl("lnil",spp)) %>%
  filter(!grepl("gf",model) & !grepl("noadmix",model)) %>%
  ggplot(aes(x=-lhood,y=model)) +
  geom_density_ridges(scale = 1.2, rel_min_height=0.001) +
  #geom_point() +
  theme_custom() +
  scale_fill_viridis(name="- log Likelihood", option="C") +
  facet_grid(~spp,scales="free_x") +
  theme(axis.title.y = element_blank(),
        legend.position = "none",
        axis.text.x=element_text(size=10),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(),
        plot.margin = margin(1,1,10,10),
        strip.text = element_text(size=rel(1.5))) +
  xlab("- Likelihood") +
  geom_vline(data=filter(lhoodMax,grepl("lnil",spp)), 
             aes(xintercept=-lhood),lty=2) #+
  geom_text(data=filter(lhoodMaxMedian,grepl("lnil",spp)), 
             aes(x=-median,y=model,label=paste0("best: ",model)),hjust=0,family="Open Sans Light") 
density_plot_lt <- lhoods %>%
  filter(!grepl("lnil",spp)) %>%
  filter(!grepl("gf",model) & !grepl("noadmix",model)) %>%
  ggplot(aes(x=-lhood,y=model)) +
  geom_density_ridges(scale = 1.2, rel_min_height=0.001) +
  #geom_point() +
  theme_custom() +
  scale_fill_viridis(name="- log Likelihood", option="C") +
  facet_grid(~spp,scales="free_x") +
  theme(axis.title.y = element_blank(),
        legend.position = "none",
        axis.text.x=element_text(size=10),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(),
        plot.margin = margin(1,1,10,10),
        strip.text = element_text(size=rel(1.5))) +
  xlab("- Likelihood") +
  scale_x_log10() +
  geom_vline(data=filter(lhoodMax,!grepl("lnil",spp)), 
             aes(xintercept=-lhood),lty=2) #+
  geom_text(data=filter(filter(lhoodMaxMedian,!grepl("gf",model)),!grepl("lnil",spp)), 
            aes(x=-median,y=model,label=paste0("best: ",model)),hjust=0,family="Open Sans Light")
  
#density
ggarrange(density_plot_lnil,density_plot_lt,ncol=1)

## adding AIC

aic <- read_table("~/Dropbox/i/projects/lates_admixture/results/mods_recompute_lhoods_twopop/mods_twopops.AIC") %>%
  mutate(AIC = as.numeric(AIC), deltaL = as.numeric(deltaL)) %>%
  mutate(mod2=factor(model, levels = levels(lhoods$model))) %>%
  drop_na() %>%
  #filter(!grepl("noMono",model)) %>%
  filter(!grepl("gf",model)) %>%
  #filter(model %in% lhoods$model) %>%
  group_by(spp) %>%
  mutate(deltaAIC=AIC25-min(AIC25,na.rm=T),
         best=case_when(deltaAIC > 10 ~ "FALSE",
                        TRUE ~ "TRUE")) %>%
  ungroup(spp)

aic_plot1 <- ggplot(aic, aes(x=mod2,y=deltaAIC,col=best)) +
  geom_point(size=3) +
  geom_hline(yintercept=0,lty=2) +
  coord_flip() +
  theme_custom() +
  theme(axis.title=element_blank(),
        axis.text.x=element_text(size=12),
        axis.text.y=element_text(),
        legend.position="none",
        panel.grid.major=element_line(linetype="dashed")) +
  #geom_hline(yintercept=min(aic$AIC200), lty=2) +
  facet_grid(~spp,scales="free_x")

aic_plot_lnil <- aic %>%
  filter(spp %in% c("lang.lnil","lmar.lnil","lmic.lnil","lsta.lnil")) %>%
  ggplot(aes(x=mod2,y=deltaAIC,col=best)) +
  geom_point(size=3) +
  coord_flip() +
  theme_custom() +
  theme(axis.title.y=element_blank(),
        axis.text.x = element_text(size=10),
        legend.position="none",
        panel.grid.major=element_line(linetype="dashed"),
        plot.margin = margin(1,1,10,10),
        strip.text = element_text(size=rel(1.5))) +
  geom_hline(yintercept=0, lty=2) +
  facet_grid(~spp,scales="free_x")
aic_plot_lt <- aic %>%
  filter(!(spp %in% c("lang.lnil","lmar.lnil","lmic.lnil","lsta.lnil"))) %>%
  ggplot(aes(x=mod2,y=deltaAIC,col=best)) +
  geom_point(size=3) +
  coord_flip() +
  theme_custom() +
  theme(axis.title.y=element_blank(),
        axis.text.x = element_text(size=10),
        legend.position="none",
        panel.grid.major=element_line(linetype="dashed"),
        plot.margin = margin(1,1,10,10),
        strip.text = element_text(size=rel(1.5))) +
  geom_hline(yintercept=0, lty=2) +
  facet_grid(~spp,scales="free_x") 

ggarrange(aic_plot_lnil,aic_plot_lt,ncol=1)
ggarrange(density_plot_lnil,aic_plot_lnil,ncol=1)
ggarrange(density_plot_lt,aic_plot_lt,ncol=1)

lhoods %>%
  ggplot(aes(x=reorder(model, lhood, FUN=mean),y=lhood,fill=model)) +
  geom_boxplot() +
  theme_custom() +
  theme(axis.title = element_blank(),
        legend.position = "none",
        axis.text.x=element_text(size=12),
        panel.grid.major.x = element_blank()) +
  coord_flip() +
  facet_wrap(~spp,scale="free_x")

############### old code below ##################
boxes <- ggboxplot(lhoods, x="model",y="lhood",fill="model",order=levels(reorder(lhoods$model,lhoods$lhood,FUN=median)),
          legend="none") +
  #geom_hline(yintercept=median(lhoods$lhood[lhoods$model == "lnil-all_admix_lt_admix_asymm2_sixpop"]), lty=2, color="gray20", alpha=0.75) +
  #theme(axis.text.x=element_text(size=12,angle=0)) +
  coord_flip() +
  theme_custom() +
  theme(legend.position = "none",
        legend.text = element_text(size=14),
        legend.title = element_blank(),
        axis.title.y = element_blank()) +
  geom_boxplot(data=lhoods,
               aes(x=model,y=lhood),
               outlier.color="gray50") +
  geom_point(data=maxEstLhoods,
             aes(x=model,y=MaxEstLhood),pch=8,cex=3)

aic_plot3 <- ggplot(aic, aes(x=reorder(mod2,-AIC50),y=deltaAIC)) +
  geom_point(size=5) +
  theme_custom() +
  theme(axis.title.y=element_blank(),
        axis.text.x=element_text(),
        axis.text.y=element_text(),
        legend.position="none") +
  scale_x_discrete(position="top") +
  coord_flip()

ggarrange(boxes,aic_plot3,nrow=1)

#####

calc_overlap <- function(mod1,mod2) {
  lhoods1 <- lhoods[lhoods$model == mod1,]
  lhoods2 <- lhoods[lhoods$model == mod2,]
  maxmod <- if (max(lhoods1) > max(lhoods2)) {
    "mod1"
  } else if (max(lhoods1) < max(lhoods2)) {
    "mod2"
  }
  
}

# plot of 95% intervals for all models
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

p1.lt <- lhoods %>%
  filter(model != "twopop_3admix_uni2")  %>%
  filter(!(spp %in% c("lang.lnil","lmar.lnil","lmic.lnil","lsta.lnil"))) %>%
  group_by(spp,model) %>%
  summarise(lowerCI = quantile(lhood,0.05),
         upperCI = quantile(lhood,0.95),
         mean = mean(lhood),
         median = median(lhood),
         min = min(lhood),
         max = max(lhood),
         mode = getmode(lhood)) %>%
  mutate(model2 = reorder(model, median)) %>%
  #left_join(aic) %>%
  ggplot(aes(col=model)) +
  geom_segment(aes(x=model,xend=model,y=lowerCI,yend=upperCI),
              lwd=2,alpha=0.7) +
#  geom_point(aes(x=reorder(model,-AIC50),y=mean),size=3,col="black",
#             position = position_nudge(x=0.4),shape=4) +
  theme_custom() +
  theme(legend.position = "none",
        axis.text.x = element_text()) +
  ylab("Model Likelihood") +
  coord_flip() +
  facet_wrap(~spp,scales="free_x",ncol=1,strip.position="right")

p2.lt <- p1.lt + ggdist::stat_halfeye(data=lhoods %>%
                                  filter(model != "twopop_3admix_uni2" &
                                           !(spp %in% c("lang.lnil","lmar.lnil","lmic.lnil","lsta.lnil"))),
                                aes(x=model,y=lhood),
    adjust=0.5, width=0.7, justification=-0.3,.width=0,point_color=NA, alpha=0.5) +
  geom_point(data=lhoods %>%
               filter(model != "twopop_3admix_uni2" & 
                        !(spp %in% c("lang.lnil","lmar.lnil","lmic.lnil","lsta.lnil"))),
             aes(x=model,y=lhood), # circles
             size=1.3,alpha=0.3,position=position_jitter(seed=1,width=0.1)) +
  geom_boxplot(data=lhoods %>%
                 filter(model != "twopop_3admix_uni2" &
                          !(spp %in% c("lang.lnil","lmar.lnil","lmic.lnil","lsta.lnil"))),
               aes(x=model,y=lhood),
               width=0.15, outlier.shape=NA, col="black", fill=NA)  +
  theme(strip.text=element_text(size=rel(1.5)),
        axis.text = element_text(size=8))
  #geom_hline(aes(yintercept=max(lowerCI),group=spp),lty=2)
  # geom_point(data=lhoods,aes(x=model,y=lhood), # tick marks
  #            shape=95, size=10, alpha=0.3) +

ggarrange(p2.lt,aic_plot_lt)

p1.lnil <- lhoods %>%
  filter(model != "twopop_3admix_uni2")  %>%
  filter((spp %in% c("lang.lnil","lmar.lnil","lmic.lnil","lsta.lnil"))) %>%
  group_by(spp,model) %>%
  summarise(lowerCI = quantile(lhood,0.05),
            upperCI = quantile(lhood,0.95),
            mean = mean(lhood),
            median = median(lhood),
            min = min(lhood),
            max = max(lhood),
            mode = getmode(lhood)) %>%
  mutate(model2 = reorder(model, median)) %>%
  #left_join(aic) %>%
  ggplot(aes(col=model)) +
  geom_segment(aes(x=model,xend=model,y=lowerCI,yend=upperCI),
               lwd=2,alpha=0.7) +
  #  geom_point(aes(x=reorder(model,-AIC50),y=mean),size=3,col="black",
  #             position = position_nudge(x=0.4),shape=4) +
  theme_custom() +
  theme(legend.position = "none",
        axis.text.x = element_text()) +
  ylab("Model Likelihood") +
  coord_flip() +
  facet_wrap(~spp,scales="free_x",ncol=1,strip.position="right")

p2.lnil <- p1.lt + ggdist::stat_halfeye(data=lhoods %>%
                                        filter(model != "twopop_3admix_uni2" &
                                                 (spp %in% c("lang.lnil","lmar.lnil","lmic.lnil","lsta.lnil"))),
                                      aes(x=model,y=lhood),
                                      adjust=0.5, width=0.7, justification=-0.3,.width=0,point_color=NA, alpha=0.5) +
  geom_point(data=lhoods %>%
               filter(model != "twopop_3admix_uni2" & 
                        (spp %in% c("lang.lnil","lmar.lnil","lmic.lnil","lsta.lnil"))),
             aes(x=model,y=lhood), # circles
             size=1.3,alpha=0.3,position=position_jitter(seed=1,width=0.1)) +
  geom_boxplot(data=lhoods %>%
                 filter(model != "twopop_3admix_uni2" &
                          (spp %in% c("lang.lnil","lmar.lnil","lmic.lnil","lsta.lnil"))),
               aes(x=model,y=lhood),
               width=0.15, outlier.shape=NA, col="black", fill=NA)  +
  theme(strip.text=element_text(size=rel(1.5)),
        axis.text = element_text(size=8))
#geom_hline(aes(yintercept=max(lowerCI),group=spp),lty=2)
# geom_point(data=lhoods,aes(x=model,y=lhood), # tick marks
#            shape=95, size=10, alpha=0.3) +
    

ggarrange(p2.lnil,aic_plot_lnil)

# top_mods <- lhoods %>%
#   group_by(model) %>%
#   summarise(lowerCI = quantile(lhood,0.05),
#             upperCI = quantile(lhood,0.95),
#             mean = mean(lhood),
#             median = median(lhood),
#             min = min(lhood),
#             max = max(lhood),
#             mode = getmode(lhood)) %>%
#   mutate(model2 = reorder(model, median, order=TRUE)) %>%
#   filter(upperCI > max(lowerCI)) %>%
#   select(model)
# 
# ##########33
# read_and_add <- function(flnm) {
#   model <- str_split_fixed(str_split_fixed(basename(flnm),"\\.",2)[1,1],"_maxL",2)[1,1]
#   par <- scan(flnm, character(0), sep = "\n", strip.white = TRUE)
#   data <- parse_fsc_par(par = par, start = 2)
#   as_tibble(data$hist_event_info$processed_hist_events) %>%
#     mutate(migrants=as.numeric(migrants),
#            model=model,
#            time_yr=time*3,
#            param=length(par))
# }
# 
# folder <- paste0("~/Dropbox/i/projects/lates_admixture/results/mods_recompute_lhoods_newSFS_dp30_091321")
# 
# parfiles <- list.files(folder,full.names=TRUE,recursive=TRUE,pattern="_maxL\\.par") %>% 
#   map_dfr(~read_and_add(.)) %>%
#   filter(!grepl("noMono",model)) %>%
#   filter(!grepl("alt",model)) %>%
#   filter(!grepl("noFossil",model))
# 
# sum_stats <- parfiles %>%
#   mutate(category = case_when(source == 2 & sink == 1 & migrants == 1 ~ "Lnil-LT",
#                               source == 4 & sink == 2 & migrants == 1 ~ "LT-LT",
#                               source == 3 & sink == 2 & migrants == 1 ~ "Lsta-Lmic",
#                               source == 5 & sink == 4 & migrants == 1 ~ "Lmar-Lang",
#                               TRUE ~ "other")) %>%
#   filter(model %in% top_mods$model) %>%
#   group_by(category) %>%
#   summarise(lowerCI = quantile(time_yr,0.05),
#             upperCI = quantile(time_yr,0.95),
#             mean = mean(time_yr),
#             median = median(time_yr),
#             min = min(time_yr),
#             max = max(time_yr),
#             mode = getmode(time_yr))
# 
# parfiles %>%
#   mutate(category = case_when(source == 2 & sink == 1 & migrants == 1 ~ "Lnil-LT",
#                    source == 4 & sink == 2 & migrants == 1 ~ "LT-LT",
#                    source == 3 & sink == 2 & migrants == 1 ~ "Lsta-Lmic",
#                    source == 5 & sink == 4 & migrants == 1 ~ "Lmar-Lang",
#                    TRUE ~ "other")) %>%
#   filter(category != "other") %>%
#   filter(model %in% top_mods$model) %>%
#   #filter(grepl("bestmod",model) | grepl("fbranch",model) | model %in% c("lnil-all_admix_asymm_sixpop_noShared_allsites","lnil-all_admix_asymm_sixpop_allsites","lnil-shared_admix_allsites","lnil-all_admix_asymm_sixpop_allsites_newSFS")) %>%
#   #ggdensity(x="time_yr",fill="category",col="category",add="mean",rug=TRUE)
#   ggplot(aes(x=time_yr)) +
#     geom_density(aes(fill=category), alpha=0.5, trim=FALSE, adjust=0.9)  +
#     geom_errorbarh(inherit.aes=FALSE, data=sum_stats[sum_stats$category != "other",],
#                    aes(xmin=lowerCI,xmax=upperCI,y=-0.00000005,color=category),height=0.00000005,size=2,show.legend=FALSE,alpha=0.5) +
#     geom_point(inherit.aes=FALSE, data=sum_stats[sum_stats$category != "other",],
#                aes(x=median,y=-0.00000005,color=category),size=3) +
#     geom_text_repel(inherit.aes=FALSE, data=sum_stats[sum_stats$category != "other",],
#              aes(x=median,y=-3.5*0.00000005,color=category,label=formatC(median,format="g",digits=2)),size=3,direction="y",ylim=c(-5*0.00000005,-0.00000005)) +
#     theme_custom() +
#     theme(legend.position = "none") +
#     xlim(0,2500000) +
#     ylab("Density") +
#     xlab("Time (years before present)") 
# 
# lhood_aic <- function(lhood,mod,folder) {
#   est <- readLines(paste0(folder,mod,".est"))
#   k <- (grep("RULES",est))-(grep("//all Ns are",est)+1)
#   aic<-2*k-2*(lhood/log10(exp(1)))
#   return(aic)
# }
# #aic <- lhood_aic(lhoods[1,1],lhoods[1,2],folder="~/Dropbox/i/projects/lates_admixture/results/mods_recompute_lhoods_newSFS_dp30/")
# lhoods_aic <- lhoods %>%
#   mutate(aic = apply(lhoods, 1, function(x) lhood_aic(as.numeric(x[1]),x[2],folder="~/Dropbox/i/projects/lates_admixture/results/mods_recompute_lhoods_newSFS_dp30/"))) %>%
#   mutate(delta_aic = max(aic)-aic)
# lhoods_aic$aic <- as.numeric(lhoods_aic$aic)
# 
# p1_aic <- lhoods_aic %>%
#   group_by(model) %>%
#   summarise(lowerCI = quantile(aic,0.05),
#             upperCI = quantile(aic,0.95),
#             mean = mean(aic),
#             median = median(aic),
#             min = min(aic),
#             max = max(aic),
#             mode = getmode(aic)) %>%
#   mutate(model2 = reorder(model, median, order=TRUE)) %>%
#   left_join(aic) %>%
#   ggplot(aes(col=model)) +
#   # geom_segment(aes(x=model,xend=model,y=lowerCI,yend=upperCI),
#   #              lwd=2,alpha=0.7) +
#   geom_point(aes(x=reorder(mod2,-deltaAIC),y=mean),size=3,col="black",
#              position = position_nudge(x=0.4),shape=4) +
#   theme_custom() +
#   theme(legend.position = "none",
#         axis.text.x = element_text()) +
#   ylab("Model AIC") +
#   coord_flip()
# 
# p2_aic <- p1_aic + ggdist::stat_halfeye(data=lhoods_aic,aes(x=model,y=aic),
#                                 adjust=0.5, width=0.7, justification=-0.3,.width=0,point_color=NA, alpha=0.5) +
#   geom_point(data=lhoods_aic,aes(x=model,y=aic), # circles
#              size=1.3,alpha=0.3,position=position_jitter(seed=1,width=0.1)) +
#   geom_boxplot(data=lhoods_aic,aes(x=model,y=aic),
#                width=0.15, outlier.shape=NA, col="black", fill=NA) 
# 
# 
# ggarrange(p2_aic,aic_plot3,nrow=1)
