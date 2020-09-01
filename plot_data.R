##script for visualizing phiX deoptimization data associated with R01
##relies on
####R_generate_deopt_seqs.R
####R_run_mfold.R

#libraries
{
  library(stringr)
  library(ggplot2)
  library(gridExtra)
  library(Biostrings)
  library(seqinr)
  data(caitab)
  library(cowplot)
  library(corrplot)
  library(vhica)
  library(RSQLite)
  library(Biobase)
  library(data.table)
  library(plyr)
  library(lme4)
  library(Hmisc)
  library(plyr)
  library(ggsci)
}

##load up WT phiX data
phix.vals <- read.csv("data/phix_seq_vals.csv")

##plot distributions from random codon changes
per25 <- read.csv("data/G_rand_0.25per_vals.csv", header = T, row.names = 1)
per50 <- read.csv("data/G_rand_0.5per_vals.csv", header = T, row.names = 1)
per75 <- read.csv("data/G_rand_0.75per_vals.csv", header = T, row.names = 1)
per100 <- read.csv("data/G_rand_1per_vals.csv", header = T, row.names = 1)
max <- read.csv("data/G_rand_max_vals.csv", header = T, row.names = 1)
max$gene <- str_replace(max$gene,"_1_","_max_")
