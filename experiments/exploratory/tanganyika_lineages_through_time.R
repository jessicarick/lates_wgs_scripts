## attempting to make a cool plot
## with lineages through time by clade
## nov 2021

library(tidyverse)
library(ggpubr)
source("theme_custom.R")

dat <- read_csv("../../data/Tanganyika_Divergence_Data.csv")

dat %>%
  filter(Cumulative_Lineages > 0) %>%
  ggplot(aes(x=Time_End,y=Cumulative_Lineages,col=Clade)) +
  geom_point() +
  geom_smooth(aes(col=Clade)) +
  theme_custom()

new_years <- seq(10,0,by=-0.1)
new_dates <- data.frame(Time_Start=new_years[-101],
                  Time_End=new_years[-1])

all_dates <- tibble(Time_Start = numeric(),
                    Time_End = numeric(),
                    Clade = character(),
                    Group = character(),
                    Cumulative_Lineages = numeric())
for (c in unique(dat$Clade)) {
  clade_dates <- new_dates %>%
    add_column(Clade = c,
               Group = dat$Group[dat$Clade == c][1],
               Cumulative_Lineages = numeric(nrow(new_dates)))
  for(i in 1:nrow(clade_dates)) {
    start <- clade_dates[i,1]
    end <- clade_dates[i,2]
    
    seg_data <- dat %>%
      arrange(Time_Start) %>%
      filter(Clade == c & Time_Start >= start) %>%
      slice(1)
    clade_dates$Cumulative_Lineages[i] <- seg_data$Cumulative_Lineages
  }
  all_dates <- all_dates %>%
    add_row(clade_dates) %>%
}

summary(all_dates)

all_dates %>%
  group_by(Group,Time_Start,Time_End) %>%
  summarize(sum_Cumulative_Lineages = sum(Cumulative_Lineages)) %>%
  group_by(Group) %>%
  mutate(Prop_Lineages = sum_Cumulative_Lineages/max(sum_Cumulative_Lineages),
         Group2 = factor(Group, levels=rev(order))) %>%
  ggplot(aes(x=Time_Start,y=Prop_Lineages,col=Group2)) +
  geom_line(alpha=0.9,lwd=2) +
  # stat_smooth(
  #   geom = 'area', method = 'loess', span = 0.1,
  #   alpha = 0.5) +
  theme_custom() +
  scale_x_reverse() +
  geom_vline(xintercept=1.68,alpha=0.1, size=5) +
  xlab("Million Years Before Present") +
  ylab("Proportion of Extant Lineages") +
  theme(axis.line.x = element_line(size=2),
        axis.line.y = element_blank(),
        panel.border = element_blank(),
        axis.ticks.length.x = unit(0.5,"cm")) +
  scale_color_brewer(palette = "Spectral",aesthetics = c("fill","color")) 

library(streamgraph)
order <- all_dates %>% 
  filter(Time_End == 0) %>% 
  group_by(Group) %>%
  summarize(sum = sum(Cumulative_Lineages)) %>%
  arrange(sum) %>% 
  pull(Group)

all_dates %>%
  mutate(Time_Start = -1*Time_Start,
         Time_End = -1*Time_End,
         Group2 = factor(Group, levels=rev(order))) %>%
  filter(Group != "Cichlids") %>%
  group_by(Group2,Time_Start,Time_End) %>%
  summarize(sum_Cumulative_Lineages = sum(Cumulative_Lineages)) %>%
  streamgraph(key="Group2",value="sum_Cumulative_Lineages",
             date="Time_End",scale="continuous",offset="zero",
             interactive=FALSE) %>%
  sg_fill_brewer("YlGnBu") %>%
  sg_colors(axis_color="black")

all_dates %>%
  mutate(Time_Start = -1*Time_Start,
         Time_End = -1*Time_End,
         Group2 = factor(Group, levels=rev(order))) %>%
  group_by(Group2,Time_Start,Time_End) %>%
  summarize(sum_Cumulative_Lineages = sum(Cumulative_Lineages)) %>%
  ggplot(aes(x=Time_End,y=sum_Cumulative_Lineages,col=Group2,fill=Group2)) +
  #geom_point() +
  geom_smooth() +
  theme_custom() +
  scale_color_brewer(palette="YlGnBu") +
  scale_fill_brewer(palette="YlGnBu")

