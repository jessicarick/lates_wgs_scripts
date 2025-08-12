library(phytools)
library(ape)
library(phangorn)
library(tidyverse)
library(here)

here::i_am("scripts/lates_admixture_analyses/plot_lates_phylogeny.R")
gtrees <- read.tree(here("results","RAxML_bipartitions.latesWGS_2.bamlist.miss0.5_maf0.01.all_genetrees.tre"))
gtrees.root <- root.multiPhylo(gtrees,"SRR3140997")

gtrees.thin <- gtrees.root[seq(1,length(gtrees.root),by=10)]
gtrees.thin <- gtrees.root[1:1000]
gtrees.thin.noCal <- lapply(gtrees.thin,drop.tip,tip=c("SRR3140997"))
class(gtrees.thin.noCal) <- "multiPhylo"

lang <- c("JAR18_139902","JAR18_139523","CEW19_135")
lmar <- c("JAR18_139933","JAR18_140014","JAR18_139895")
lmic <- c("CEW16_022","JAR18_139370","CEW19_134")
lsta <- c("CEW17_151","CEW17_127","JAR18_140024")
lnil <- "AMNH216682"

gtrees.thin.mono <- sapply(gtrees.thin.noCal, function(x) is.monophyletic(x,tips=c(lang,lmar,lmic,lsta)))
gtrees.thin.lmar.lang <- sapply(gtrees.thin.noCal, function(x) is.monophyletic(x,tips=c(lang,lmar)))
gtrees.thin.lmic.lsta <- sapply(gtrees.thin.noCal, function(x) is.monophyletic(x,tips=c(lmic,lsta)))

class <- case_when(gtrees.thin.lmar.lang & gtrees.thin.lmic.lsta ~ "raxml",
                   gtrees.thin.mono & gtrees.thin.lmar.lang ~ "one_clade",
                   gtrees.thin.mono & gtrees.thin.lmic.lsta ~ "one_clade",
                   gtrees.thin.mono & !(gtrees.thin.lmar.lang | gtrees.thin.lmic.lsta) ~ "mono",
                   TRUE ~ "lnil_within")
class <- as.factor(class)

cols <- viridis::magma(5)
densiTree(gtrees.thin.noCal,alpha=0.2,type="cladogram",scaleX=TRUE,width=1.5,
          jitter = list(amount=1,random=TRUE),tip.color="blue",
          consensus=c("JAR18_139902","JAR18_139523","CEW19_135",
                      "JAR18_139933","JAR18_140014","JAR18_139895",
                      "CEW17_151","CEW17_127","JAR18_140024",
                      "CEW16_022","JAR18_139370","CEW19_134",
                      "AMNH216682"
                      ),
          col="gray20", direction="leftwards", scale.bar=TRUE)
 
# densiTree(gtrees.thin.noCal,alpha=0.1,type="cladogram",scaleX=TRUE,width=1.5,
#           jitter = list(amount=0,random=TRUE),tip.color="blue",
#           consensus=c("JAR18_139902_A","JAR18_139902_B","JAR18_139523_A","JAR18_139523_B","CEW19_135_A","CEW19_135_B",
#                       "JAR18_139933_A","JAR18_139933_B","JAR18_140014_A","JAR18_140014_B","JAR18_139895_A","JAR18_139895_B",
#                       "CEW17_151_A","CEW17_151_B","CEW17_127_A","CEW17_127_B","JAR18_140024_A","JAR18_140024_B",
#                       "CEW16_022_A","CEW16_022_B","JAR18_139370_A","JAR18_139370_B","CEW19_134_A","CEW19_134_B",
#                       "AMNH216682_A","AMNH216682_B"
#                       ),
#           col="gray20", direction="leftwards", scale.bar=FALSE)
