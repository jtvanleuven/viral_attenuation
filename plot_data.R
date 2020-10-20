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
  library(reshape2)
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



#look at new sequences
diffcount <- read.csv(file="results/diffcounts.csv")
diffcount.l <- melt(diffcount)
ggplot(data=diffcount.l, aes(x=variable, y=X, fill=value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme(text = element_text(size=8))
ggsave(filename = "plots/AAdiffCount.pdf", width = 6, height= 20, units = "in")

tab.g <- read.csv("results/G_new_diffs.csv", header = T, stringsAsFactors = F)
names(tab.g) <- c(1:176, "id")
plot.g <- melt(tab.g, id.vars = c("id"))
ggplot(data=plot.g, aes(x=variable, y=id, fill=value)) +
  geom_tile() +
  scale_fill_gradient(low = "steelblue", high = "white") +
  theme(text = element_text(size=10), legend.position = "none", axis.text.x=element_text(angle=-45)) +
  xlab("") +
  ylab("")
ggsave(filename = "plots/G_diffs.pdf", width = 16, height= 10, units = "in")

tab.f <- read.csv("results/F_new_diffs.csv", header = T, stringsAsFactors = F)
names(tab.f) <- c(1:428, "id")
plot.f <- melt(tab.f, id.vars = c("id"))
#plot.f <- plot.f[!str_detect(plot.f$id, "highest"),]
ggplot(data=plot.f, aes(x=variable, y=id, fill=value)) +
  geom_tile() +
  scale_fill_gradient(low = "#31a354", high = "white") +
  theme(text = element_text(size=10), legend.position = "none", axis.text.x=element_text(angle=-45)) +
  xlab("") +
  ylab("")
ggsave(filename = "plots/F_diffs.pdf", width = 16, height= 10, units = "in")

tab.h <- read.csv("results/H_new_diffs.csv", header = T, stringsAsFactors = F)
names(tab.h) <- c(1:329, "id")
plot.h <- melt(tab.h, id.vars = c("id"))
ggplot(data=plot.h, aes(x=variable, y=id, fill=value)) +
  geom_tile() +
  scale_fill_gradient(low = "#756bb1", high = "white") +
  theme(text = element_text(size=10), legend.position = "none", axis.text.x=element_text(angle=-45)) +
  xlab("") +
  ylab("")
ggsave(filename = "plots/H_diffs.pdf", width = 16, height= 10, units = "in")
