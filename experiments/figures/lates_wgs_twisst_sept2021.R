## twisst plotting for WGS
source("../plot_twisst.R")
library(ggpubr)
library(chromoMap)

group <- "latesWGS"

# weights files
dir <- "~/Dropbox/i/projects/lates_admixture/results/twisst_apr2023/"
files <- list.files(path=dir, pattern=".weights.csv.gz")

window_data_file <- "../../results/twisst_apr2023/latesWGS_2.bamlist.all_miss0.5_maf0.01.recode.beagle.raxml.w100.data.tsv"

weights.all <- read.table(paste0(dir,files[1]), header = T)
#normalise rows so weights sum to 1
weights.all <- weights.all / apply(weights.all, 1, sum)
#retrieve the names of the topologies
topoNames = names(weights.all)

window_data = read.table(window_data_file, header = T)

#exclude any rows where data is missing
good_rows = which(is.na(apply(weights.all, 1, sum)) == F)
weights.all <- weights.all[good_rows, ]
window_data = window_data[good_rows, ]

lnil.within <- paste0("topo",c(16,17,18,28,29,30,31,32,33,43,44,45,46,47,48,53,54,58,59,60,61,62,63,19,22,25,34,37,40,49,52,55,64,67,70,20,21,23,24,26,27,35,36,38,39,41,42,50,51,56,57,65,66,68,69,71,72,73,74,75,76:105))
sp.topos <- paste0("topo",c(1,4,9,12,13))

read_and_combine <- function(file) {
  chr <- gsub(".*(ASB_LG[0-9]+).*","\\1",file)
  weights <- read_table2(paste0(dir,file),col_names=TRUE,comment="#")
  weights <- as_tibble(weights/apply(weights,1,sum))  %>%
    rownames_to_column(var = "window") %>%
    drop_na()
  window_data <- read_table2(paste0(dir,gsub("\\.weights\\.csv\\.gz","\\.data\\.tsv",file)), col_names = TRUE) %>%
    rownames_to_column(var = "window")
  
  #exclude any rows where data is missing
  #good_rows <- which(is.na(apply(weights[,-1], 1, sum)) == F)
  #weights <- weights[good_rows, ]
  #window_data <- window_data[good_rows, ]
  
  all_dat <- window_data %>%
    right_join(weights, by="window") %>%
    pivot_longer(cols=starts_with("topo"), names_to="topo",values_to="weighting")
  return(all_dat)
}

combined_weights <- files %>%
  map_df(~read_and_combine(.))
summary(combined_weights)

## which topos are most common?
best_topos <- combined_weights %>%
  group_by(scaffold,window) %>%
  summarise(best_topo=topo[which.max(weighting)])

best_topos %>%
  ggplot(aes(x=best_topo)) +
  geom_histogram(stat="count",bins=105) +
  theme_custom() +
  theme(axis.text.x=element_text(angle=90))
  
dom_topos <- best_topos %>%
  group_by(best_topo) %>%
  tally() %>%
  arrange(desc(n))

# combined_weights %>%
#   filter(topo %in% lnil.within) %>%
#   mutate(chr = scaffold) %>%
#   ggplot() +
#   #geom_point(alpha=0.5, size=0.6) +
#   stat_summary_bin(aes(x=mid,y=weighting,fill=topo),
#                    fun="mean",geom="area",
#                    inherit.aes=FALSE,position="stack",bins=50) +
#   stat_summary_bin(data=dstats,aes(x=windowStart,y=f_d_corr/2),
#                    inherit.aes=FALSE,bins=50,fun.data="mean_se") +
#   stat_summary_bin(data=dstats,aes(x=windowStart,y=f_d_corr/2),
#                    inherit.aes=FALSE,bins=50,fun="mean",geom="line",lty=2) +
#   #geom_line() +
#   scale_fill_manual(values = rep(colors, unique(length(lnil.within)))) +
#   theme_custom() +
#   theme(legend.position="none") +
#   scale_y_continuous(name="Topology Weighting",
#                      sec.axis=sec_axis(trans=~.*2,name="Corrected f_d Statistic")) +
#   facet_wrap(~chr, scales="free", ncol=4)

## plotting only topo 13
topo13 <- combined_weights %>%
  filter(topo == "topo13") %>%
  mutate(chr_num=gsub(pattern="ASB_LG([0-9]+)","\\1",scaffold,fixed=FALSE)) %>%
  mutate(chr_num = case_when(chr_num == "16_LG22" ~ "16",
                             chr_num == "7_1" ~ "7.1",
                             chr_num == "7_2" ~ "7.2",
                             TRUE ~ chr_num))

data_cum <- topo13 %>%
  mutate(chr_num = as.numeric(chr_num)) %>%
  arrange(chr_num) %>% 
  group_by(scaffold) %>%
  summarize(max_bp = max(end)) %>%
  mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>%
  select(scaffold, bp_add)

topo13_cum <- topo13 %>%
  inner_join(data_cum, by="scaffold") %>%
  mutate(mid_cum = (mid + bp_add))

axis_set <- topo13_cum %>%
  group_by(scaffold) %>%
  summarize(center = mean(mid_cum))

colors <- c("#264653","#2a9d8f","#e9c46a","#f4a261","#e76f51")

topo13_cum %>%
  ggplot(aes(x=mid_cum, y=weighting,col=as.factor(scaffold))) +
  stat_summary_bin(alpha=0.5, size=0.6, fun="median", geom="point",binwidth=5000) +
  #geom_hline(yintercept=quantile(topo13$weighting,0.995),lty=3, col="gray50") +
  scale_x_continuous(label = axis_set$scaffold, breaks = axis_set$center, 
                     limits = c(-10,max(topo13_cum$mid_cum)+10),
                     expand = c(0.01,0.01)) +
  scale_color_manual(values = rep(colors, unique(length(axis_set$scaffold)))) +
  scale_y_continuous(expand=c(0.01,0.01)) +
  xlab("") +
  ylab("TWISST Weighting") +
  theme_custom() +
  theme(legend.position="none",
        axis.text.x = element_text(angle=60, size=16, vjust=1, hjust=1))

## topo13 smoothed plot
#use loess to smooth weights.
chr_cum <- window_data %>%
  mutate(chr_num = gsub("ASB_LG","",scaffold),
         chr_num = case_when(chr_num == "16_LG22" ~ "16",
                             chr_num == "7_1" ~ "7.1",
                             chr_num == "7_2" ~ "7.2",
                             TRUE ~ chr_num)) %>%
  arrange(chr_num) %>% 
  group_by(scaffold) %>%
  summarize(max_bp = max(end)) %>%
  mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>%
  select(scaffold, bp_add)
windows_cum <- window_data %>%
  inner_join(chr_cum, by="scaffold") %>%
  mutate(mid_cum = (mid + bp_add))


all_cum_smooth = smooth.weights(
  window_positions = windows_cum$mid_cum,
  weights_dataframe = weights.all,
  span = 0.1,
  window_sites = windows_cum$sites,
  na.action(na.omit)
)

plot.weights(
  weights_dataframe = all_cum_smooth[-c(278479,278480),] %>% select(topo13),
  positions = windows_cum$mid_cum[-c(278479,278480)],
  line_cols = topo_cols,
  lwd = 1,
  fill_cols = scales::alpha(topo_cols,0.5),
  stacked = FALSE,
  xlab = "",
  ylab = ""
)

topo13_smooth <- topo13_cum %>%
  filter(scaffold=="ASB_LG1") %>%
  ggplot(aes(x=mid, y=weighting, col=as.factor(scaffold))) +
  stat_summary_bin(alpha=0.5, size=0.6, fun="median", geom="point",binwidth=500) +
  #geom_hline(yintercept=quantile(topo13$weighting,0.995),lty=3, col="gray50") +
  geom_smooth(method="loess", se=FALSE, span=0.05, aes(group=as.factor(scaffold)),
              color="blue") +
  # scale_x_continuous(label = axis_set$scaffold, breaks = axis_set$center, 
  #                    limits = c(-10,max(topo13_cum$mid_cum)+10),
  #                    expand = c(0.01,0.01)) +
  scale_color_manual(values = rep(colors, unique(length(axis_set$scaffold)))) +
  scale_y_continuous(expand=c(0.01,0.01)) +
  xlab("") +
  ylab("TWISST Weighting") +
  theme_custom() +
  theme(legend.position="none",
        axis.text.x = element_text(angle=60, size=16, vjust=1, hjust=1))


##########################
## compared to f_d_corr
## new analysis 092221
#######################

all_combined_sp <- tibble()
for (chrom in unique(combined_weights$scaffold)){
  max_bp <- max(max(combined_weights %>% filter(scaffold == chrom) %>% select(end)),
                max(dstats %>% filter(chr == chrom) %>% select(windowEnd)))
  binned_weights <- combined_weights %>%
    filter(scaffold == chrom) %>%
    filter(topo == "topo13") %>%
    group_by(window, scaffold, start, end, mid, sites) %>%
    summarize(sum = sum(weighting)) %>%
    ungroup() %>%
    group_by(ints = cut(mid, breaks = seq(0, max_bp, 5000))) %>%
    summarise(mean = mean(sum),
              median = median(sum),
              n=n())
  binned_dstats <- dstats %>%
    filter(chr == chrom) %>%
    mutate(windowMid = (windowStart + windowEnd)/2) %>%
    group_by(ints = cut(windowMid, breaks = seq(0, max_bp, 5000))) %>%
    summarise(mean = mean(f_d_corr),
              median = median(f_d_corr),
              n=n())
  combined <- binned_weights %>%
    full_join(binned_dstats, by="ints") %>%
    mutate(scaffold = chrom)

  if (length(all_combined_sp) == 0) {
    all_combined_sp <- combined
  } else {
    all_combined_sp <- all_combined_sp %>%
      add_row(combined)
  }
}

p1 <- all_combined_sp %>%
  drop_na() %>%
  #filter(mean.y > 0) %>%
  ggplot(aes(x=mean.x,y=mean.y)) +
     geom_point(col=colors_reord3[1],alpha=0.4,size=2) +
     #stat_summary_bin(bins=40) +
     #geom_smooth(method="gam", span=2) +
     xlab("mean topology weighting") +
     ylab("mean corrected f_d stat") +
     theme_custom() +
     xlim(0,1) +
     ylim(0,1)
     #geom_hline(yintercept=quantile(all_combined$mean.y,0.95,na.rm=T), lty=2, alpha=0.8) +
     #geom_hline(yintercept=quantile(all_combined$mean.y,0.99,na.rm=T), lty=3, alpha=0.8)
     # ggpmisc::stat_poly_eq(formula = y ~ x,
     #                       eq.with.lhs = "italic(y)~`=`~",
     #                       aes(label = paste(..rr.label..,..p.value.label..,sep = "*plain(\",\")~")),
     #                       parse = TRUE, label.x = 0.9, label.y=0.9, size=5, hjust=1)

# now, introgression topologies
all_combined_introg <- tibble()
for (chrom in unique(combined_weights$scaffold)){
  max_bp <- max(max(combined_weights %>% filter(scaffold == chrom) %>% select(end)),
                max(dstats %>% filter(chr == chrom) %>% select(windowEnd)))
  binned_weights <- combined_weights %>%
    filter(scaffold == chrom) %>%
    filter(topo %in% lnil.within) %>%
    group_by(window, scaffold, start, end, mid, sites) %>%
    summarize(sum = sum(weighting)) %>%
    ungroup() %>%
    group_by(ints = cut(mid, breaks = seq(0, max_bp, 5000))) %>%
    summarise(mean = mean(sum),
              median = median(sum),
              n=n())
  binned_dstats <- dstats %>%
    filter(chr == chrom) %>%
    mutate(windowMid = (windowStart + windowEnd)/2) %>%
    group_by(ints = cut(windowMid, breaks = seq(0, max_bp, 5000))) %>%
    summarise(mean = mean(f_d_corr),
              median = median(f_d_corr),
              n=n())
  combined <- binned_weights %>%
    full_join(binned_dstats, by="ints") %>%
    mutate(scaffold = chrom)
  
  if (length(all_combined_introg) == 0) {
    all_combined_introg <- combined
  } else {
    all_combined_introg <- all_combined_introg %>%
      add_row(combined)
  }
}

p2 <- all_combined_introg %>%
  drop_na() %>%
  #filter(mean.y > 0) %>%
  ggplot(aes(x=mean.x,y=mean.y)) +
  geom_point(col=colors_reord3[3],alpha=0.4) +
  #stat_summary_bin(bins=40) +
  #geom_smooth(method="gam", span=2) +
  xlab("mean topology weighting") +
  ylab("mean corrected f_d stat") +
  theme_custom() +
  xlim(0,1) +
  ylim(0,1) 

## combined plot
ggarrange(p1,p2,common.legend=TRUE,legend="right",
          labels="AUTO",font.label=list(size=30,family="Open Sans"),
          hjust=-3.5,vjust=2)

## both on same plot
fig3a <- all_combined_sp %>%
  drop_na() %>%
  ggplot() +
  geom_point(aes(x=median.x,y=mean.y+0.005),col="black",alpha=0.6,size=1.2) +
  geom_point(data=all_combined_introg,aes(x=median.x,y=-mean.y-0.005),col=colors_reord3[4],alpha=0.6,size=1.2, inherit.aes=TRUE) +
  geom_hline(yintercept=0,lwd=0.8) +
  #stat_summary_bin(bins=40) +
  #geom_smooth(method="gam", span=2) +
  xlab("mean topology weighting") +
  ylab("mean corrected f_d stat") +
  theme_custom() +
  scale_x_continuous(limits=c(0,1),expand=c(0.01,0.01)) 
  ylim(-0.8,0.8)

###############

all_combined %>%
  ggplot(aes(x=n.x,y=n.y,col=scaffold)) +
  geom_point()

ggplot() +
  geom_point(data=dstats[dstats$chr == "ASB_LG14",],aes(x=(windowStart + windowEnd)/2,y=f_d_corr)) +
  geom_point(data=combined_weights[combined_weights$scaffold == "ASB_LG14" & combined_weights$topo == "topo13",],aes(x=mid,y=weighting),col="turquoise")
ggplot() +
  stat_summary_bin(data=dstats[dstats$chr == "ASB_LG14",],aes(x=(windowStart + windowEnd)/2,y=f_d_corr),breaks=seq(0,max(dstats$windowEnd[dstats$chr == "ASB_LG14"]),1000),fun="mean") +
  stat_summary_bin(data=combined_weights[combined_weights$scaffold == "ASB_LG14" & combined_weights$topo == "topo13",],aes(x=mid,y=weighting),col="turquoise",breaks=seq(0,max(dstats$windowEnd[dstats$chr == "ASB_LG14"]),1000),fun="mean") +
  theme_custom()


# par(mar = c(0, 0, 0, 0), mfrow = c(4, 6))
# for (chrom in scaffolds) {
#   weights.chrom <- weights[window_data$scaffold == chrom, ]
#   window_data.chrom <- window_data[window_data$scaffold == chrom, ]
#   
#   ################ plot smoothed data ######
#   #use loess to smooth weights.
#   weights_smooth = smooth.weights(
#     window_positions = combined_weights$mid,
#     weights_dataframe = weights.chrom,
#     span = 0.3,
#     window_sites = window_data.chrom$sites,
#     na.action(na.omit)
#   )
#   #plot smoothed data with polygons stacked
#   #png(file = paste0("lates01_twisst.smooth.stacked",chrom,".png"))
#   #par(mar = c(4,4,1,1))
#   plot.weights(
#     weights_dataframe = weights_smooth,
#     positions = window_data.chrom$mid,
#     line_cols = topo_cols,
#     lwd = 1,
#     fill_cols = scales::alpha(topo_cols,0.5),
#     stacked = TRUE,
#     xlab = "",
#     ylab = ""
#   )
#   abline(h = 0.5, lty = 2)
#   text(x = median(window_data.chrom$mid),
#        y = 0.9,
#        paste(chrom, " smoothed"))
# }
# dev.off()

################ plot topologies #########

# library(ape)
# 
# topos = read.tree(
#   file = "../../results/twisst_052021/latesWGS_july2020_allsites.nocontam.ASB_LG10.beagle.raxml.w500SNP.topologies.trees"
#   )
# 
# #unrooted trees
# 
# pdf(
#   file = paste0(group, ".topos.w500.all.unrooted.pdf"),
#   width = 5,
#   height = 3
# )
# par(mfrow = c(4, 3),
#     mar = c(2, 2, 2, 2),
#     xpd = NA)
# for (n in 1:length(topos)) {
#   plot.phylo(
#     topos[[n]],
#     type = "unrooted",
#     #edge.color = topo_cols[n],
#     edge.width = 5,
#     cex = 1,
#     rotate.tree = 90,
#     adj = .5,
#     label.offset = .5
#   )
#   mtext(side = 3, text = paste0("topo", n))
# }
# dev.off()
# 
# 
# #root trees
# for (i in 1:length(topos)) {
#   topos[[i]] <- root(topos[[i]], "Lcal", resolve.root = T)
# }
# 
# pdf(
#   file = paste0(group, ".topos.all.w500.rooted.pdf"),
#   width = 8,
#   height = 3
# )
# par(mfrow = c(1, 3),
#     mar = c(2, 2, 2, 2),
#     xpd = NA)
# for (n in 1:length(topos)) {
#   plot.phylo(
#     topos[[n]],
#     type = "clad",
#     #edge.color = topo_cols[n],
#     edge.width = 5,
#     label.offset = 0.3,
#     cex = 1
#   )
#   mtext(side = 3, text = paste0("topo", n))
# }
# dev.off()


#weights_locs <- left_join(weights,window_data,by="window")
weights_byTopo <- combined_weights %>%
  group_by(scaffold,topo) %>%
  summarise(sum=sum(weighting),
            sum_greater0.3 = sum(weighting[weighting > 0.3]))

#dom_byTopo <- apply(weights_locs[,-c(1:5)],2,function(x){sum(x > 0.3, na.rm=T)})
barplot(weights_byTopo$sum_greater0.3, las=2, cex.axis=0.5)

par(mar=c(5,4,1,1))
pdf("latesWGS_twisst_plots_w10000.pdf", width=11, height=4)
plot1 <- weights_byTopo %>%
  ggplot(aes(x=topo, y=sum))+
  geom_col() +
  theme_pubr() +
  theme(axis.text.x=element_text(angle=90, size=8)) +
  ggtitle("all chrom")
print(plot1)

for (scaf in unique(weights_byTopo$scaffold)){
  plot <- weights_byTopo %>%
    filter(scaffold==scaf) %>%
    ggplot(aes(x=topo, y=sum))+
    geom_col() +
    theme_pubr() +
    theme(axis.text.x=element_text(angle=90, size=8)) +
    ggtitle(scaf)
  print(plot)
}
dev.off()

weights_topo1 <- weights_locs[weights_locs$topo1 > 0.5,]
weights_topo2 <- weights_locs[weights_locs$topo2 > 0.5,]
weights_topo3 <- weights_locs[weights_locs$topo3 > 0.5,]

sort(apply(weights_locs,2,function(x){sum(x > 0.5, na.rm=T)}))

## categorize topologies by general features
weights_byTopo_cat <- weights_byTopo %>%
  mutate(topo.category = case_when(
    topo %in% paste0("topo",c(1,4)) ~ "lsta.lmic",  
    topo %in% paste0("topo",c(9,12)) ~ "lmar.lang", 
    topo %in% paste0("topo",c(2,7)) ~ "lsta.lmar", 
    topo %in% paste0("topo",c(3,10)) ~ "lmic.lmar", 
    topo %in% paste0("topo",c(5,8)) ~ "lsta.lang", 
    topo %in% paste0("topo",c(6,11)) ~ "lmic.lang", 
    topo %in% paste0("topo",c(13)) ~ "two.clust.lsta.lmic", 
    topo %in% paste0("topo",c(14)) ~ "two.clust.lsta.lmar", 
    topo %in% paste0("topo",c(15)) ~ "two.clust.lmic.lmar",
    topo %in% paste0("topo",c(16,17,18,28,29,30,31,32,33,43,44,45,46,47,48,53,54,58,59,60,61,62,63)) ~ "one.outside",
    topo %in% paste0("topo",c(19,22,25,34,37,40,49,52,55,64,67,70)) ~ "two.outside",
    topo %in% paste0("topo",c(20,21,23,24,26,27,35,3,38,39,41,42,50,51,56,57,65,66,68,69,71,72,73,74,75)) ~ "three.outside", 
    topo %in% paste0("topo",c(76:105)) ~ "sister.to.lnil",
    topo %in% TRUE ~ "UNK"
  )) %>% 
  group_by(scaffold,topo.category) %>%
  summarise(total_category_weight=mean(sum))

weights_byTopo_cat %>%
  ggplot(aes(x=scaffold,y=total_category_weight,col=topo.category)) +
  geom_point() 

## attempting plotting on chromosomes
# using chromMap
library(chromoMap)
best_topos <- combined_weights %>%
  group_by(scaffold,window,start,end,mid,sites) %>%
  summarise(best_topo=topo[which.max(weighting)],
            max = max(weighting))

weights_best_topos <- best_topos %>%
  mutate(topo.category = case_when(
    best_topo %in% paste0("topo",c(1,4)) ~ "lsta.lmic",  
    best_topo %in% paste0("topo",c(9,12)) ~ "lmar.lang", 
    best_topo %in% paste0("topo",c(2,7)) ~ "lsta.lmar", 
    best_topo %in% paste0("topo",c(3,10)) ~ "lmic.lmar", 
    best_topo %in% paste0("topo",c(5,8)) ~ "lsta.lang", 
    best_topo %in% paste0("topo",c(6,11)) ~ "lmic.lang", 
    best_topo %in% paste0("topo",c(13)) ~ "two.clust.lsta.lmic", 
    best_topo %in% paste0("topo",c(14)) ~ "two.clust.lsta.lmar", 
    best_topo %in% paste0("topo",c(15)) ~ "two.clust.lmic.lmar",
    best_topo %in% paste0("topo",c(16,17,18,28,29,30,31,32,33,43,44,45,46,47,48,53,54,58,59,60,61,62,63)) ~ "one.outside",
    best_topo %in% paste0("topo",c(19,22,25,34,37,40,49,52,55,64,67,70)) ~ "two.outside",
    best_topo %in% paste0("topo",c(20,21,23,24,26,27,35,36,38,39,41,42,50,51,56,57,65,66,68,69,71,72,73,74,75)) ~ "three.outside", 
    best_topo %in% paste0("topo",c(76:105)) ~ "sister.to.lnil",
    best_topo %in% TRUE ~ "UNK"
    ),
    topo.category.general = case_when(
      best_topo %in% paste0("topo",c(1,4,12)) ~ "alternate",  
      best_topo %in% paste0("topo",c(9)) ~ "alternate", 
      best_topo %in% paste0("topo",c(2,7)) ~ "other.LT.together", 
      best_topo %in% paste0("topo",c(3,10)) ~ "other.LT.together", 
      best_topo %in% paste0("topo",c(5,8)) ~ "other.LT.together", 
      best_topo %in% paste0("topo",c(6,11)) ~ "other.LT.together", 
      best_topo %in% paste0("topo",c(13)) ~ "raxml.tree", 
      best_topo %in% paste0("topo",c(14)) ~ "other.LT.together", 
      best_topo %in% paste0("topo",c(15)) ~ "other.LT.together",
      best_topo %in% paste0("topo",c(16,17,18,28,29,30,31,32,33,43,44,45,46,47,48,53,54,58,59,60,61,62,63)) ~ "lnil.within.LT",
      best_topo %in% paste0("topo",c(19,22,25,34,37,40,49,52,55,64,67,70)) ~ "lnil.within.LT",
      best_topo %in% paste0("topo",c(20,21,23,24,26,27,35,36,38,39,41,42,50,51,56,57,65,66,68,69,71,72,73,74,75)) ~ "lnil.within.LT", 
      best_topo %in% paste0("topo",c(76:105)) ~ "lnil.within.LT",
      best_topo %in% TRUE ~ "UNK"
    ),
    element_name = paste0(scaffold,"_",start)) %>%
  relocate(element_name)
weights_best_topos$topo.category.general <- factor(weights_best_topos$topo.category.general,
                                                      levels=c("raxml.tree","alternate",
                                                               "other.LT.together","lnil.within.LT"))

chrom_data <- read.csv("Lcal_chrom_lengths.csv",skip=1,header=FALSE)
anno_data <- weights_best_topos %>%
  ungroup() %>%
  select(element_name,scaffold,start,end,topo.category.general) %>% 
  as.data.frame(.)
anno_data$topo.category.general <- factor(anno_data$topo.category.general,
                                          levels=c("raxml.tree","alternate",
                                                   "other.LT.together","lnil.within.LT"))

#colors <- c("#007da3","#008db8","#d1927b","#a44a3f","#72977e","#bedada")
#colors <- wesanderson::wes_palette("FantasticFox1")
colors_reord <- c("#E1D565","#94B46F","#57A3B8","#7B469A")
#colors2 <- c("gray60","#e07a5f","#3d405b","#81b29a","#f2cc8f")
#colors3 <- c("#f25f5c","#ffe066",rep("gray80",2),"#70c1b3")

colors_reord2 <- c("#C6893E","#E1D565","#94B46F","gray80")
colors_reord3 <- c("#29667a","#00A8B8","#5BD7B4","#d6d6d6")
scales::show_col(colors_reord3)

chromoMap(list(chrom_data),
          list(anno_data),
          data_based_color_map=T,
          data_type="categorical",
          data_colors=list(rev(colors_reord3)),
          left_margin=150,
          top_margin=25,
          chr_color=c("gray80"),
          chr_width=18,
          chr_length=10,
          ch_gap=4,
          chr_curve=3,
          text_font_size=16,
          grid_text_size=12,
          canvas_height=1200,
          canvas_width=1200,
          legend = T, 
          lg_x = 150,
          lg_y = 300)

# plot matching histogram of categories
plot_best1 <- weights_best_topos %>% 
  mutate(topo_num=as.integer(gsub("topo([0-9]+)","\\1",best_topo))) %>%
  #filter(max > 0.5) %>%
  ggplot(aes(x=topo_num,fill=topo.category.general)) + 
  geom_histogram(bins=105,col="gray20",alpha=0.9,stat="count") + 
  # theme_bw(base_size=12, base_family="Open Sans Light") %+replace% 
  # theme(
  #   panel.background  = element_blank(),
  #   plot.background = element_rect(fill="white", color=NA), 
  #   legend.background = element_rect(fill="transparent", colour=NA),
  #   legend.key = element_rect(fill="transparent", colour=NA),
  #   panel.grid.major = element_blank(),
  #   panel.grid.minor = element_blank(),
  #   axis.title = element_text(size=18, family="Open Sans"),
  #   axis.text = element_text(size=14),
  #   legend.text = element_text(size=14),
  #   legend.title = element_blank()
  # ) + 
  theme_custom() +
  scale_fill_manual(values=(colors_reord3)) +
  theme(axis.title.x=element_blank())
plot_best2 <- weights_best_topos %>% 
  mutate(topo_num=as.integer(gsub("topo([0-9]+)","\\1",best_topo))) %>%
  filter(max > 0.8) %>%
  ggplot(aes(x=topo_num,fill=topo.category.general)) + 
  geom_histogram(bins=105,col="gray20",alpha=0.9,aes(y=..count../sum(..count..))) + 
  # theme_bw(base_size=12, base_family="Open Sans Light") %+replace% 
  # theme(
  #   panel.background  = element_blank(),
  #   plot.background = element_rect(fill="white", color=NA), 
  #   legend.background = element_rect(fill="transparent", colour=NA),
  #   legend.key = element_rect(fill="transparent", colour=NA),
  #   panel.grid.major = element_blank(),
  #   panel.grid.minor = element_blank(),
  #   axis.title = element_text(size=18, family="Open Sans"),
  #   axis.text = element_text(size=14),
  #   legend.text = element_text(size=14),
  #   legend.title = element_blank()
  # ) + 
  theme_custom() +
  scale_fill_manual(values=(colors_reord3)) +
  xlab("Topology Number") +
  ylab("Proportion of Windows")

ggarrange(plot_best1,plot_best2,ncol=1,labels="AUTO",common.legend=TRUE,
          font.label = list(size=24,family="Open Sans"))

weights_best_topos %>% 
  mutate(topo_num=as.integer(gsub("topo([0-9]+)","\\1",best_topo))) %>%
  filter(max > 0.8) %>%
  group_by(topo.category.general) %>%
  summarize(count = n(),
            windowsize = sum(end-start)) %>%
  mutate(prop = windowsize/sum(windowsize),
         prop_n = count/sum(count))

## adding dstats to chromomap? ## DIDN'T WORK!!
bins <- anno_data %>%
  group_by(scaffold) %>%
  arrange(start) %>%
  ungroup() %>%
  arrange(scaffold) %>%
  select(element_name,scaffold,start,end)
dstats_agg <- do.call(data.frame,aggregate(f_d_corr~chr, dstats, 
                                           FUN=function(x) c(Count=length(x), Mean=mean(x))))

for (i in unique(dstats$chr)[-1]) {
  dstats_chr <- dstats[dstats$chr == i,]
  dstats_chr$interval = cut(dstats_chr$windowStart, breaks=bins$start[bins$scaffold == i], right = FALSE, dig.lab=10)
  dstats2 <- dstats2 %>%
    add_row(dstats_chr)
}

dstat_anno <- dstats2 %>%
  rename(scaffold = chr) %>%
  group_by(scaffold,interval) %>%
  summarize(mean_fd = mean(f_d_corr, na.rm=TRUE),
            sd_fd = sd(f_d_corr)) %>%
  mutate(interval = as.character(interval)) %>%
  mutate(start = substring(gsub("([0-9+]),[0-9]+\\)","\\1",interval),2),
         end = gsub(".*,([0-9]+)\\)","\\1",interval),
         element_name = paste(scaffold,start,sep="_")) %>%
  select(scaffold,element_name,mean_fd) 

new_anno <- anno_data %>%
  left_join(dstat_anno,by=c("scaffold","element_name")) %>%
  relocate("mean_fd",.before = "topo.category.general")

chromoMap(list(chrom_data),
          list(new_anno),
          data_based_color_map=T,
          data_type="numeric",
          plots = "scatter",
          data_colors=list(rev(colors_reord3)),
          plot_filter = list(c("col","byCategory")),
          left_margin=150,
          top_margin=25,
          chr_color=c("gray80"),
          chr_width=18,
          chr_length=10,
          ch_gap=4,
          chr_curve=3,
          text_font_size=16,
          grid_text_size=12,
          canvas_height=1200,
          canvas_width=1200,
          legend = T, 
          lg_x = 150,
          lg_y = 300)