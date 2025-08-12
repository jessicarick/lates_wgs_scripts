## plotting and comparing msfs for fastsimcoal
library(tidyverse)
library(here)
library(ggpubr)

here::i_am("scripts/lates_admixture_analyses/sfs_comparison.R")

old_sfs <- read_table(here("data","latesWGS_july2020_allsites.nocontam_dp5_miss2_q30_MSFS.obs"),skip=2,
                      col_names=FALSE)
 old_sfs <- read_table(here("data","latesWGS_2.all_miss2_q30_MSFS.obs"),skip=2,
                       col_names=FALSE)
new_sfs <- read_table(here("data","new.latesWGS_2.summed_MSFS.obs"),skip=2,
                      col_names=FALSE)

both_sfs <- tibble(old=t(old_sfs),new=t(new_sfs)[,1])

both_sfs %>%
  filter(old < 1e6) %>%
  ggplot(aes(x=old,y=new)) +
  geom_point()

p1 <- both_sfs %>%
  mutate(sfs_entry=seq(1,nrow(both_sfs))) %>%
  filter(old < 1e6 & sfs_entry < 5001) %>%
  ggplot(aes(x=sfs_entry)) +
  geom_point(aes(y=old),alpha=0.5) +
  geom_point(aes(y=new),col="blue",alpha=0.5) +
  theme_bw() +
  coord_cartesian(xlim=c(-1,5001),expand=FALSE) +
  ylab("value")

p2 <- both_sfs %>%
  mutate(sfs_entry=seq(1,nrow(both_sfs))) %>%
  filter(old < 1e6 & sfs_entry > 5000 & sfs_entry < 10001) %>%
  ggplot(aes(x=sfs_entry)) +
  geom_point(aes(y=old),alpha=0.5) +
  geom_point(aes(y=new),col="blue",alpha=0.5) +
  theme_bw() +
  coord_cartesian(xlim=c(4999,10001),expand=FALSE) +
  ylab("value")

p3 <- both_sfs %>%
  mutate(sfs_entry=seq(1,nrow(both_sfs))) %>%
  filter(old < 1e6 & sfs_entry > 10000 & sfs_entry < 15001) %>%
  ggplot(aes(x=sfs_entry)) +
  geom_point(aes(y=old),alpha=0.5) +
  geom_point(aes(y=new),col="blue",alpha=0.5) +
  theme_bw() +
  coord_cartesian(xlim=c(9999,15001),expand=FALSE) +
  ylab("value")

p4 <- both_sfs %>%
  mutate(sfs_entry=seq(1,nrow(both_sfs))) %>%
  filter(old < 1e6 & sfs_entry > 15000 ) %>%
  ggplot(aes(x=sfs_entry)) +
  geom_point(aes(y=old),alpha=0.5) +
  geom_point(aes(y=new),col="blue",alpha=0.5) +
  theme_bw() +
  coord_cartesian(xlim=c(14999,21610),expand=FALSE) +
  ylab("value")

ggarrange(p1,p2,p3,p4,
          ncol=1)
