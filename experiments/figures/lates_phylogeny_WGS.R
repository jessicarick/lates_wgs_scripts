## phylogeny plotting

#library(treeio)
library(ape)
library(phytools)
library(tidyverse)
library(ggtree)

tree <- read.tree("../../results/RAxML_bipartitions.latesWGS_july2020_noContam_miss0.5_maf0.01")
#tree$tip.label <- gsub(".*aln_([A-Z]+[0-9]+.*)\\.sorted\\.bam","\\1",tree$tip.label)

fishinfo <- read.csv('../../../lates_admixture/data/lates_combined_all_info_lib.csv',
                     header=TRUE,stringsAsFactors = FALSE)
treeinfo <- left_join(as_tibble_col(tree$tip.label, column_name="Moran_FishID"),
                      fishinfo, by="Moran_FishID")
treeinfo$Location <- as.factor(treeinfo$Location)

plotTree(tree)
edgelabels()
nodelabels()
tiplabels()

tree2 <- phytools::reroot(tree, interactive = TRUE)
tree2$root.edge <- 0.01
tree3 <- drop.tip(tree2,"SRR3140997")
tree3$root.edge <- 0.0005
plot.phylo(tree3,show.tip.label=FALSE,root.edge=TRUE,type="phylogram",label.offset=0.0001,
           edge.width=2)
edgelabels()
nodelabels(tree3$node.label)

nodes <- c(24,22,19,17)
lates_colors <- c("#AAA7BC","#FFC000","#E76F51","#858585","#0A0A0A")


# tree layout, with bars to designate clades
p <- ggtree(tree3) 
p + #geom_tippoint(aes(color=Location, shape=entropy_dnaID), size=4, alpha=0.8) +
  #geom_rootedge(size=1.2) +
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

# circular tree layout, with node colors by location
p2 <- ggtree(tree2, ladderize=FALSE) 
p2


## bootstraps
bs.trees <- read.tree("../../results/RAxML_bootstrap.latesWGS_july2020_rawvariants_miss1.0_maf0.01.out")

class(bs.trees) <- "multiPhylo"
densityTree(bs.trees)
phangorn::densiTree(bs.trees, type="cladogram",alpha=0.2)
