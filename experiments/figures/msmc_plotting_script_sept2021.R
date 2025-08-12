### MSMC Plotting Script
## J. Rick, Sept 2021

library(readr)
library(tidyverse)

pop.cols <- c("#E76F51","#5CCDAF","#FFC000","#AAA7BC","#858585","#0A0A0A")

## plotting the original msmc output
dir <- "~/Dropbox/i/projects/lates_admixture/results/msmc_may2023"
finals <- list.files(path=dir,pattern="msmc_output\\.msmc_050823.*\\.final\\.txt",full.names=TRUE)
msmc_output <- list()
msmc_output_df <- data.frame()

for (f in finals) {
  ind <- gsub(".*msmc_050823_([A-Z]+[0-9]+.*)\\.final\\.txt","\\1",f)
  log <- gsub("\\.final\\.txt","\\.log",f)
  mut_rate <- as.numeric(scan(log,skip=3,nlines=1,what=character())[2])
  recomb_rate <- as.numeric(scan(log,skip=4,nlines=1,what=character())[2])
  df <- read_table2(paste0(f)) %>%
    mutate(ind=ind,
           mu=mut_rate/100000,
           #mu=3.5e-9,
           recomb=recomb_rate,
           spp=case_when(ind %in% c("CEW16_022","JAR18_139370","CEW19_134") ~ "lmic",
                         ind %in% c("CEW17_127","CEW17_151","JAR18_140024") ~ "lsta",
                         ind %in% c("JAR18_139523","JAR18_139902","CEW19_135","CEW17_120") ~ "lang",
                         ind %in% c("JAR18_139895","JAR18_139933","JAR18_140014") ~ "lmar",
                         ind == "SRR3140997" ~ "lcal",
                         ind == "AMNH216682" ~ "lnil"))
  # plot(((df$left_,time_boundary/(mu.mar))*gen.mar),(1/df$lambda)/(2*(mu.mar)), 
  #      col="red", type="s", xlab="Years before Present", 
  #      ylab="Effective Population Size",lwd=2,ylim=c(0,100000),main=ind)
  msmc_output[[ind]] <- df
  msmc_output_df <- rbind(msmc_output_df,df)
}

gen.mic <- 4 # 4 years, per fishbase
gen.mar <- 3.5 # 3-4 years, per fishbase
gen.sta <- 2.5 # average 28 months: http://www.fao.org/3/S7047E/S7047E08.htm
gen.ang <- 4 # 4 years, per fishbase
gen.nil <- 4 # age at maturity 2yr for males, 4 for females: https://doi.org/10.1007/BF00005874
gen.cal <- 3 # age at maturity 3-4yr for males, 6-8yr for females: http://www.fao.org/fishery/culturedspecies/Lates_calcarifer/en

ne_plot <- msmc_output_df %>%
  mutate(right_time_boundary = as.numeric(right_time_boundary)) %>%
  mutate(spp = factor(spp, levels = c("lmar", "lang", "lsta", "lmic", "lnil", "lcal"))) %>%
  mutate(right_time_corrected =  case_when(spp == "lmic" ~ (right_time_boundary/mu)*gen.mic,
                                       spp == "lmar" ~ (right_time_boundary/mu)*gen.mar,
                                       spp == "lsta" ~ (right_time_boundary/mu)*gen.sta,
                                       spp == "lang" ~ (right_time_boundary/mu)*gen.ang,
                                       spp == "lnil" ~ (right_time_boundary/mu)*gen.nil,
                                       spp == "lcal" ~ (right_time_boundary/mu)*gen.cal),
         pop_corrected=(1/lambda)/(2*mu)) %>%
  #filter(right_time_corrected > 0 & pop_corrected > 0) %>%
  ggplot(aes(x=right_time_corrected,y=pop_corrected,col=ind,lty=spp,group=ind)) +
  #annotate(geom="rect",xmin=77000,xmax=731000,ymin=3e3,ymax=3e6,fill="gray50",alpha=0.2) +
  #geom_point() +
  #geom_vline(xintercept=c(1000000,2000000),lty=2,col="gray80") +
  #geom_rect(xmin=1000000/1e6,xmax=2000000/1e6,ymin=-1000,ymax=400000,fill="gray90",alpha=0.1,col="gray80",lty=2) +
  geom_step(linewidth=1.2,alpha=0.6,lty=1) +
  #geom_smooth(aes(group=spp),span=0.25,method="loess",alpha=0.2) +
  coord_cartesian(ylim=c(1e4,1e7),xlim=c(10000,5e7)) +
  #coord_trans(x="log10",xlim=c(20000/1e6,2700000/1e6),ylim=c(1000,300000),expand=FALSE) +
  #scale_x_continuous(trans='log10',limits=c(100,4000000)) +
  #coord_cartesian(xlim=c(20000,2700000),ylim=c(1000,300000)) +
  theme_custom() + 
  scale_x_log10(labels = scales::comma,expand=c(0,0)) +
  scale_y_log10(expand=c(0,0)) +
  annotation_logticks(sides="bl") +
  theme(axis.text.x = element_text(angle=0),
        legend.position = "right",
        legend.text = element_text(size=rel(1.5)),
        panel.grid.major.x = element_line(linetype="dashed")) +
  #scale_color_manual(values=pop.cols,labels=c("L. mariae","L. angustifrons", "L. stappersii", "L. microlepis",
  #                                            "L. niloticus", "L. calcarifer")) +
  guides(color = guide_legend(override.aes = list(size=2))) +
  ylab("Effective Population Size (Ne)") +
  xlab("") + 
  #annotate(geom="rect",xmin=1000,xmax=40000,ymin=1e4,ymax=1e7,fill="white",alpha=0.7)+
  #annotate(geom="rect",xmin=2100000,xmax=5000000,ymin=1e4,ymax=1e7,fill="white",alpha=0.7) +
  geom_vline(xintercept=1000, linetype="dashed", col="gray80") +
  geom_vline(xintercept=10000, linetype="dashed", col="gray80") +
  geom_vline(xintercept=1000000, linetype="dashed", col="gray80") 
  facet_wrap(~spp)
  
mut_plot <- msmc_output_df %>% 
  group_by(ind,spp) %>%
  summarize(mu=mean(mu)) %>%
  ggplot(aes(x=spp,y=mu,col=spp,fill=spp)) + 
  geom_point(size=4,alpha=0.7,pch=19) + 
  geom_hline(yintercept=3.5e-9, lty=3) + # this is the estimate i've used in fastsimcoal
  theme_custom() +
  ylab("Estimated Mutation Rate") +
  xlab("") +
  scale_color_manual(values=pop.cols,labels=c("L. mariae","L. angustifrons", "L. stappersii", "L. microlepis",
                                              "L. niloticus", "L. calcarifer")) +
  theme(legend.position = "none")

recomb_plot <- msmc_output_df %>% 
  ggplot(aes(x=spp,y=recomb,col=spp)) + 
  geom_point(size=3) + 
  #geom_hline(yintercept=3.5e-9, lty=3) + 
  theme_custom()

msmc_output_df$recomb / (msmc_output_df$mu*100000) # recomb is 0.25*mu 

ggarrange(ne_plot,mut_plot,nrow=1,widths=c(3,1))

########################           
## species combined runs
########################

## plotting the original msmc output
dir <- "~/Dropbox/i/projects/lates_admixture/results/msmc_091721_whatshap_beagle/"
finals_spp <- list.files(path=dir,pattern="msmc_output.l[a-z].*091721_30seg.final.txt")
msmc_output_spp <- list()
msmc_output_spp_df <- data.frame()

for (f in finals_spp) {
  spp <- gsub(".*msmc_([a-z]+)_091721_30seg.*","\\1",f)
  log <- list.files(path=dir,pattern=paste0(spp,"_091721_30seg.log"),full.names = TRUE)
  mut_rate <- as.numeric(scan(log,skip=3,nlines=1,what=character())[2])
  recomb_rate <- as.numeric(scan(log,skip=4,nlines=1,what=character())[2])
  df <- read_table2(paste0(dir,f)) %>%
    mutate(spp=spp,
           mu=mut_rate/100000,
           recomb=recomb_rate
           # spp=case_when(ind %in% c("CEW16_022","JAR18_139370","CEW19_134") ~ "lmic",
           #               ind %in% c("CEW17_127","CEW17_151","JAR18_140024") ~ "lsta",
           #               ind %in% c("JAR18_139523","JAR18_139902","CEW19_135") ~ "lang",
           #               ind %in% c("JAR18_139895","JAR18_139933","JAR18_140014") ~ "lmar",
           #               ind == "SRR3140997" ~ "lcal",
           #               ind == "AMNH216682" ~ "lnil")
           )
  # plot(((df$left_,time_boundary/(mu.mar))*gen.mar),(1/df$lambda)/(2*(mu.mar)), 
  #      col="red", type="s", xlab="Years before Present", 
  #      ylab="Effective Population Size",lwd=2,ylim=c(0,100000),main=ind)
  msmc_output_spp[[ind]] <- df
  msmc_output_spp_df <- rbind(msmc_output_spp_df,df)
}

gen.mic <- 4 # 4 years, per fishbase
gen.mar <- 3.5 # 3-4 years, per fishbase
gen.sta <- 2.3 # average 28 months: http://www.fao.org/3/S7047E/S7047E08.htm
gen.ang <- 4 # 4 years, per fishbase
gen.nil <- 4.5 # age at maturity 2yr for males, 4 for females: https://doi.org/10.1007/BF00005874
gen.cal <- 2.5 # age at maturity 3-4yr for males, 6-8yr for females: http://www.fao.org/fishery/culturedspecies/Lates_calcarifer/en

msmc_output_spp_df %>%
  mutate(left_time_corrected=case_when(spp == "lmic" ~ (left_time_boundary/mu)*gen.mic,
                                       spp == "lmar" ~ (left_time_boundary/mu)*gen.mar,
                                       spp == "lsta" ~ (left_time_boundary/mu)*gen.sta,
                                       spp == "lang" ~ (left_time_boundary/mu)*gen.ang,
                                       spp == "lnil" ~ (left_time_boundary/mu)*gen.nil,
                                       spp == "lcal" ~ (left_time_boundary/mu)*gen.cal),
         pop_corrected=(1/lambda)/(2*mu)) %>%
  ggplot(aes(x=left_time_corrected,y=pop_corrected,col=spp,group=spp)) +
  #geom_point() +
  geom_step(size=1.2,alpha=0.7,lty=1) +
  #geom_smooth(aes(group=spp),span=0.5,method="loess") +
  scale_x_continuous(trans='log10',limits=c(10000,5000000)) +
  coord_cartesian(xlim=c(80000,2000000),ylim=c(1000,500000)) +
  theme_custom()



##-------------------------------------------##
## plotting the bootstrap replicates
list.bootstraps <- NA

for (i in seq(1,50,1)) {
  name <- NA
  name <- read_delim(paste0("D:/mt_moran/MSMC_091817/msmc_output.091817.CEW16_022.091817.bootstrap_",i,".final.txt"),"\t", escape_double = FALSE, trim_ws = TRUE)
  name[,5:10] <- NA
  colnames(name)<- c("time_index","left_time_boundary","right_time_boundary","lambda","GenerationsMin","GenerationsMax","YearsMin","YearsMax","scaledPopSize","PopSize")
  
  name$GenerationsMin <- (name$left_time_boundary)/mu.mar  ## dividing by mutation rate to get generations
  name$YearsMin <- (name$GenerationsMin)/gen.mar  ## dividing by generation time (from fishbase) to get years
  
  name$GenerationsMax <- (name$right_time_boundary)/mu.mar  ## dividing by mutation rate to get generations
  name$YearsMax <- (name$GenerationsMax)/gen.mar  ## dividing by generation time (from fishbase) to get years
  
  name$scaledPopSize <- 1/name$lambda
  name$PopSize <- name$scaledPopSize/(2*mu.mar)
  
  assign(paste0("msmc_bootstrap.",i),name)
  
}

list.bootstraps <- lapply(ls(patt='^msmc_bootstrap.'),get)


plot((msmc_output$left_time_boundary/((mu.mar)*gen.mar)),((1/msmc_output$lambda)/(2*(mu.mar))), col="red", type="s", lwd=2, log="x",ylim=c(0,100000))

colors <- rainbow(50, s=0.5)
for (i in 1:50) {
  yearsmin <- as.data.frame(list.bootstraps[i])$YearsMin
  popsize <- as.data.frame(list.bootstraps[i])$PopSize
  points((yearsmin),(popsize),col=colors[i], type="s")
}

