## Attempt to characterize differences in MSMC trajectories
## J. Rick, 7 June 2022

library(tidyverse)
library(here)
library(BINCOR)
library(TSdist)
?bin_cor

mu <- 3.5e-9

## loading msmc results
here::i_am("scripts/lates_admixture_analyses/msmc_trajectory_comparison.R")
dir <- here("results","msmc_091721_whatshap_beagle")
finals <- list.files(path=dir,pattern="30seg\\.final\\.txt",full.names=TRUE)
msmc_output <- list()
msmc_output_df <- data.frame()

for (f in finals) {
  ind <- gsub(".*msmc_([A-Z]+[0-9]+.*)_[0-9]+_30seg.*","\\1",f) # pulling out individual name
  #log <- gsub("\\.final\\.txt","\\.log",f)
  #mut_rate <- as.numeric(scan(log,skip=3,nlines=1,what=character())[2]) # reads inferred mu estimate
  recomb_rate <- as.numeric(scan(log,skip=4,nlines=1,what=character())[2])
  df <- read_table2(paste0(f)) %>%
    mutate(ind=ind,
           #mu=mut_rate/100000,
           mu=3.5e-9,
           recomb=recomb_rate,
           spp=case_when(ind %in% c("CEW16_022","JAR18_139370","CEW19_134") ~ "lmic",
                         ind %in% c("CEW17_127","CEW17_151","JAR18_140024") ~ "lsta",
                         ind %in% c("JAR18_139523","JAR18_139902","CEW19_135") ~ "lang",
                         ind %in% c("JAR18_139895","JAR18_139933","JAR18_140014") ~ "lmar",
                         ind == "SRR3140997" ~ "lcal",
                         ind == "AMNH216682" ~ "lnil"))
  # plot(((df$left_,time_boundary/(mu.mar))*gen.mar),(1/df$lambda)/(2*(mu.mar)), 
  #      col="red", type="s", xlab="Years before Present", 
  #      ylab="Effective Population Size",lwd=2,ylim=c(0,100000),main=ind)
  msmc_output[[ind]] <- df
  msmc_output_df <- rbind(msmc_output_df,df)
}

######-----------######
## now, moving to calculating distances among trajectories

dist.all.df <- tibble(ind1=character(),
                      ind2=character(),
                      dist_acf=numeric(),
                      dist_cor=numeric(),
                      dist_ccor=numeric())
for (i in 1:length(unique(msmc_output_df$ind))) {
  for (j in 1:length(unique(msmc_output_df$ind))) {
      ind1 <- unique(msmc_output_df$ind)[i]
      ind2 <- unique(msmc_output_df$ind)[j]
      
      # trajectory for ind1
      msmc_dat1 <- msmc_output_df %>%
        mutate(right_time_boundary = as.numeric(right_time_boundary)) %>%
        mutate(spp = factor(spp, levels = c("lmar", "lang", "lsta", "lmic", "lnil", "lcal"))) %>%
        mutate(right_time_corrected =  case_when(spp == "lmic" ~ (right_time_boundary/mu)*gen.mic,
                                                 spp == "lmar" ~ (right_time_boundary/mu)*gen.mar,
                                                 spp == "lsta" ~ (right_time_boundary/mu)*gen.sta,
                                                 spp == "lang" ~ (right_time_boundary/mu)*gen.ang,
                                                 spp == "lnil" ~ (right_time_boundary/mu)*gen.nil,
                                                 spp == "lcal" ~ (right_time_boundary/mu)*gen.cal),
               pop_corrected=(1/lambda)/(2*mu)) %>%
        filter(right_time_corrected > 0 & pop_corrected > 0) %>%
        filter(ind == ind1 & pop_corrected < 1e7 & right_time_corrected < 3e6) %>%
        select(right_time_corrected,pop_corrected) %>%
        as.data.frame(.)
      
      # trajectory for ind2
      msmc_dat2 <- msmc_output_df %>%
        mutate(right_time_boundary = as.numeric(right_time_boundary)) %>%
        mutate(spp = factor(spp, levels = c("lmar", "lang", "lsta", "lmic", "lnil", "lcal"))) %>%
        mutate(right_time_corrected =  case_when(spp == "lmic" ~ (right_time_boundary/mu)*gen.mic,
                                                 spp == "lmar" ~ (right_time_boundary/mu)*gen.mar,
                                                 spp == "lsta" ~ (right_time_boundary/mu)*gen.sta,
                                                 spp == "lang" ~ (right_time_boundary/mu)*gen.ang,
                                                 spp == "lnil" ~ (right_time_boundary/mu)*gen.nil,
                                                 spp == "lcal" ~ (right_time_boundary/mu)*gen.cal),
               pop_corrected=(1/lambda)/(2*mu)) %>%
        filter(right_time_corrected > 0 & pop_corrected > 0) %>%
        filter(ind == ind2 & pop_corrected < 1e7 & right_time_corrected < 3e6) %>%
        select(right_time_corrected,pop_corrected) %>%
        as.data.frame(.)
      
      # interpolate data points for a fixed interval from 20000 to 1e6
      # should probably play with the choices for start/stop and spacing
      interp_dat1 <- approx(msmc_dat1, method = "linear", xout=seq(20000,1e6,by=20000),
                            rule = 1, f = 0, ties = mean, na.rm = TRUE)
      iterp_dat1_df <- data.frame(x=interp_dat1$x,y=interp_dat1$y)
      interp_dat2 <- approx(msmc_dat2, method = "linear", xout=seq(20000,1e6,by=20000),
                            rule = 1, f = 0, ties = mean, na.rm = TRUE)
      iterp_dat2_df <- data.frame(x=interp_dat2$x,y=interp_dat2$y)
      
      # three different distance metrics for time series
      # from TSdist package
      dist.acf <- ACFDistance(iterp_dat1_df$y,iterp_dat2_df$y)
      dist.cor <- CorDistance(iterp_dat1_df$y,iterp_dat2_df$y)
      dist.ccor <- CCorDistance(iterp_dat1_df$y,iterp_dat2_df$y,lag.max=10)
      
      dist.all <- tibble(ind1=ind1,
                            ind2=ind2,
                            dist_acf=dist.acf,
                            dist_cor=dist.cor,
                            dist_ccor=dist.ccor)
      dist.all.df <- dist.all.df %>%
        add_row(dist.all)
    }
}

# plot correlation among distance metrics
ggplot(data=dist.all.df) +
  geom_point(aes(x=dist_acf,y=dist_cor,col=ind1))

# create distance matrix
dist.all.mat <- dist.all.df %>%
  select(ind1,ind2,dist_cor) %>%
  pivot_wider(id_cols=ind1,names_from=ind2,values_from=dist_cor) %>%
  column_to_rownames("ind1") %>%
  as.matrix()

# pcoa on distance matrix, and plotting the result
dist.all.pcoa <- ape::pcoa(dist.all.mat)
pcoa.df <- as.tibble(dist.all.pcoa$vectors) %>%
  add_column(ind=rownames(dist.all.pcoa$vectors),.before=1) %>%
  mutate(spp=case_when(ind %in% c("CEW16_022","JAR18_139370","CEW19_134") ~ "lmic",
                       ind %in% c("CEW17_127","CEW17_151","JAR18_140024") ~ "lsta",
                       ind %in% c("JAR18_139523","JAR18_139902","CEW19_135") ~ "lang",
                       ind %in% c("JAR18_139895","JAR18_139933","JAR18_140014") ~ "lmar",
                       ind == "SRR3140997" ~ "lcal",
                       ind == "AMNH216682" ~ "lnil"))
ggplot(data=pcoa.df, aes(x=Axis.1,y=Axis.2,col=ind)) +
  geom_point(aes(shape=spp),cex=3)

#######################################
####### unused code below #############
#######################################

dat1 <- dat_msmc %>%
  filter(p_par == "30seg") %>%
  mutate(left_years = (left_time_boundary/mu)*3,
         ne = (1/lambda_00)/(2*mu)) %>%
  filter(ind_pair == "AMNH216682_CEW16_022" & left_years > 20000 & left_years < 1e7) %>%
  select(left_years,ne) %>%
  as.data.frame(.)
dat2 <- dat_msmc %>%
  filter(p_par == "30seg") %>%
  mutate(left_years = (left_time_boundary/mu)*3,
         ne = (1/lambda_00)/(2*mu)) %>%
  filter(ind_pair == "CEW16_022_CEW17_127" & left_years > 20000 & left_years < 1e7) %>%
  select(left_years,ne) %>%
  as.data.frame(.)

test.bincor <- bin_cor(dat1,dat2,FLAGTAU=3,"test.bincor")

binned_dat <- test.bincor$Binned_time_series %>%
  as.data.frame(.) %>%
  drop_na()
cor_ts(binned_dat[,1:2],binned_dat[,c(1,3)],KoCM="spearman")
ccf_ts(binned_dat[,1:2],binned_dat[,c(1,3)])

ggplot(dat1) +
  geom_point(aes(x=left_years,y=ne)) +
  geom_point(data=dat2,aes(x=left_years,y=ne), col="turquoise") +
  geom_point(data=as.data.frame(test.bincor$Binned_time_series),aes(x=tau.t.mean,y=mean.ts1),pch=2,col="red") +
  geom_point(data=as.data.frame(test.bincor$Binned_time_series),aes(x=tau.t.mean,y=mean.ts2),pch=4,col="green")
test.bincor$Binned_time_series
