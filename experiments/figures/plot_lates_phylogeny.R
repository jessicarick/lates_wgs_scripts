## phylogeny plotting

library(treeio)
library(ape)
library(phytools)
library(tidyverse)
library(ggtree)
library(here)

here::i_am("scripts/lates_admixture_analyses/plot_lates_phylogeny.R")
tree_text <- "(CEW19_135:1,(JAR18_139523:1,(JAR18_139902:1,((JAR18_139933:1,(JAR18_139895:1,JAR18_140014:1)36.19:0.043819594037420403)99.85:6.064573143486736,((SRR3140997:5,AMNH216682:5)99.21:4.434688082193132,((JAR18_140024:1,(CEW17_127:1,CEW17_151:1)41.35:0.12816839166757135)99.95:7.231605440742371,(CEW16_022:1,(CEW19_134:1,JAR18_139370:1)40.92:0.12079011839649637)99.85:6.07635922173381)57.39:0.44760154911504013)74.13:0.9466574664170476)99.84:6.014996633060428)35.03:0.02581696451676173):1);"
tree <- read.tree(here("results","RAxML_bestTree.latesWGS_2.bamlist.all_miss0.5_maf0.01"))
# tree <- read.tree(text=tree_text[1])
# tree$tip.label <- gsub(".*aln_([A-Z]+[0-9]+.*)\\.sorted\\.bam","\\1",tree$tip.label)

fishinfo <- read.csv('~/Dropbox/i/projects/lates_admixture/data/lates_wgs_metadata.csv',
                     header=TRUE,stringsAsFactors = FALSE)
  
treeinfo <- left_join(as_tibble_col(tree$tip.label, column_name="SampleID"),
                      fishinfo, by="SampleID")

plot(tree)
edgelabels()
nodelabels()
tiplabels()

tree2 <- phytools::reroot(tree, interactive = TRUE)
tree2$root.edge <- 0.01
tree3 <- drop.tip(tree2,"SRR3140997")
tree3$root.edge <- 0.1
plot.phylo(tree3,show.tip.label=TRUE)
edgelabels()
nodelabels()

nodes <- c(26,24,21,18)
#lates_colors <- c("#AAA7BC","#FFC000","#E76F51","#858585","#0A0A0A")
lates_colors <- c("#E76F51","#5CCDAF","#FFC000","#AAA7BC","#858585","#0A0A0A")

# tree layout, with bars to designate clades
p <- ggtree(tree3, layout="fan", open.angle=180, size=1.2, ladderize=FALSE) 
p + #geom_tippoint(aes(color=Location, shape=entropy_dnaID), size=4, alpha=0.8) +
  geom_rootedge(size=1.2) +
  theme_tree() +
  #geom_hilight(node=nodes[1], fill="#e76f51", alpha=0.4) +
  geom_cladelabel(node=nodes[1], color="#e76f51", fontsize=0,
                  label="", align=TRUE, offset=0,
                  offset.text=0, barsize=2) +
  #geom_hilight(node=nodes[2], fill="#2a9d8f",alpha=0.4) +
  geom_cladelabel(node=nodes[2], color="#2a9d8f", fontsize=0,
                  label="", align=TRUE, offset=0,
                  offset.text=0, barsize=2) +
  #geom_hilight(node=nodes[3], fill="#e9c46a",alpha=0.4) +
  geom_cladelabel(node=nodes[3], color="#e9c46a", fontsize=0,
                  label="", align=TRUE, offset=0,
                  offset.text=0, barsize=2) +
  #geom_hilight(node=nodes[4], fill="#70567d",alpha=0.4) +
  geom_cladelabel(node=nodes[4], color="#70567d", fontsize=0,
                  label="", align=TRUE, offset=0,
                  offset.text=0, barsize=2) +
  geom_cladelabel(node=16, color="#000000", fontsize=0,
                  label="", align=FALSE, offset=0.003,
                  offset.text=0, barsize=2)

p <- ggtree(tree3, size=1.2, ladderize=TRUE) 
p %<+% treeinfo + 
  geom_tippoint(size=4, alpha=0.8, 
                aes(col=factor(Species,levels=c("Lates calcarifer","Lates niloticus","Lates microlepis",
                                                "Lates stappersii","Lates angustifrons","Lates mariae"))),
                position=position_nudge(x=0.2)) +
  geom_rootedge(size=1.2) +
  theme_tree() +
  geom_nodelab(nudge_y=0,nudge_x=0.4) +
  scale_color_manual(values=rev(lates_colors)[-1])+
  theme(legend.position="none") 
# circular tree layout, with node colors by location
p2 <- ggtree(tree3, ladderize=FALSE) 
p2


## bootstraps
bs.trees <- read.tree("../../results/RAxML_bootstrap.latesWGS_july2020_rawvariants_miss1.0_maf0.01.out")

class(bs.trees) <- "multiPhylo"
densityTree(bs.trees)
phangorn::densiTree(bs.trees, type="cladogram",alpha=0.2)


