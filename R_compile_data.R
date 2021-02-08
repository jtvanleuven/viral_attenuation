#libraries
{
library(seqinr)
library(stringr)
library(reshape2)
library(ggplot2)
library(dplyr)
library(readxl)
library(ggsci)
library(cowplot)
library(ggthemes)
library(tidyr)
}

lu2 <- read_xlsx("~/Dropbox/R01 phage research/Data/Lu/Constructs_Transcript Stability_Worst.xlsx", sheet = 5, skip = 3) 
lu2 <- lu2[which(!is.na(lu2[,1])),]
LuComments <- rbind(lu2[,c(1,38)])
LuComments <- LuComments[!is.na(LuComments[,2]),]
lu2 <- lu2[,1:36]  ##drop comments and average fit
lu.l <- melt(lu2, id.vars = c('...1', 'sample #'))
lu.l <- lu.l[which(!str_detect(lu.l$value, '[a-zA-Z]')),]  ##drop cells with comments and stuff
names(lu.l) <- c('construct', 'sample', 'date', 'fitness')
lu.l$date <- openxlsx::convertToDate(as.character(lu.l$date))
lu.l$fitness <- round(as.numeric(lu.l$fitness), digits=2)

emma <- read_xlsx("~/Dropbox/R01 phage research/Data/Emma/Constructs_Shuffles_Randoms.xlsx", sheet = 4, skip = 7)
emma.s <- na.omit(emma[,c("Construct", "Date", "flask #", "fitness")])

emma2 <- read_xlsx("~/Dropbox/R01 phage research/Data/Emma/Constructs_Shuffles_Randoms.xlsx", sheet = 5) 
emma2 <- emma2[which(!is.na(emma2$Construct)),]
EmmaComments <- emma2[,c('Construct', 'Photo?', 'times assembled', 'notes', 'not enough/too many plqs on plate to count')]
emma2 <- emma2[,1:40]  ##drop comments and average fit
emma.l <- melt(emma2, id.vars = c('Sample #', 'Construct'))
emma.l <- emma.l[which(!str_detect(emma.l$value, '[a-zA-Z]')),]  ##drop cells with comments and stuff
names(emma.l) <- c( 'sample', 'construct', 'date', 'fitness')
emma.l$date <- openxlsx::convertToDate(as.character(emma.l$date))
emma.l$fitness <- round(as.numeric(emma.l$fitness), digits=2)
emma.l[emma.l$construct=='wt Sanger2mut (10^3)',]$construct <- 'wt Sanger2mut'

fails <- read_excel("~/Dropbox/R01 phage research/Construct Failures.xlsx", sheet = 1, col_names = F)
fail.s <- data.frame(fails[,1])
fail.s <- fail.s[which(!str_detect(fail.s$...1, '\\.')),]  

lu.cnt <- nrow(lu.l[str_detect(lu.l$construct, "wt"),])
emma.cnt <- nrow(emma.l[str_detect(emma.l$construct, "wt"),])

wt <- rbind(emma.l[str_detect(emma.l$construct, "wt"),], lu.l[str_detect(lu.l$construct, "wt"),])
wt$sci <- c(rep("emma",emma.cnt), rep("lu", lu.cnt))
wt$construct <- "wt"
ggplot(wt, aes(x=sci, y=fitness, color=sci, fill=sci)) +
  geom_violin(alpha=0.2, show.legend = FALSE) +
  geom_jitter(size=3, width = 0.1, show.legend = FALSE) +
  scale_color_aaas() +
  scale_fill_aaas() +
  theme_classic()
#ggsave(filename = "plots/wt_fitness.pdf", width = 5.5, height= 5.5, units = "in")
#ggsave(filename = "plots/wt_fitness.png", width = 5.5, height= 5.5, units = "in")
  

emma.s <- emma.l
##more formatting
emma.s$sci <- rep("emma",nrow(emma.s))
emma.s$gene <- str_split(emma.s$construct, " ", simplify = T)[,1]
emma.s$type <- str_split(emma.s$construct, " ", simplify = T)[,2]
emma.s$treat <- str_split(emma.s$construct, " ", simplify = T)[,3]
emma.s$rep <- str_split(emma.s$construct, " ", simplify = T)[,5]
emma.s[emma.s$gene=="wt",]$treat <- "0%"
emma.s$treat_numeric <- str_replace(emma.s$treat, "%", "")
emma.s$treat_numeric <- str_replace(emma.s$treat_numeric, "max", "125")
emma.s$treat_numeric <- as.numeric(emma.s$treat_numeric)/100
emma.s[emma.s$rep=="",]$rep <- "1"
emma.s$type <- str_replace(emma.s$type, "Sanger2mut", "wt")

lu.s <- lu.l
lu.s$sci <- rep("lu", nrow(lu.s))
lu.s$gene <- str_replace(str_extract(lu.s$construct, "[FGH]_"), "_", "")
lu.s[str_detect(lu.s$construct, "wt"),]$gene <- "wt"
lu.s$type <- "NA"
lu.s[str_detect(lu.s$construct, "shuf"),]$type <- "shuffle"
lu.s[str_detect(lu.s$construct, "worst"),]$type <- "worst"
lu.s[str_detect(lu.s$construct, "Hi"),]$type <- "shuffle_Hi"
lu.s[str_detect(lu.s$construct, "Lo"),]$type <- "shuffle_Lo"
lu.s[str_detect(lu.s$construct, "wt"),]$type <- "wt"
lu.s$treat <- str_extract(lu.s$construct, "[\\d\\.]{2,}")
lu.s[lu.s$gene=="wt",]$treat <- "0"
lu.s[str_detect(lu.s$type, "shuffle"),]$treat <- "0.5"
lu.s$treat <- paste(as.character(as.numeric(lu.s$treat)*100), "%", sep="")
lu.s$rep <- str_replace(str_extract(lu.s$construct, "rep\\d"), "rep", "")
lu.s[str_detect(lu.s$construct, "stab"),]$rep <- str_split(lu.s[str_detect(lu.s$construct, "stab"),]$construct,"_",simplify = T)[,3]
lu.s$rep <- str_replace(lu.s$rep,"Lo","")
lu.s$rep <- str_replace(lu.s$rep,"Hi","")
lu.s[lu.s$gene=="wt",]$rep <- "1"
lu.s$treat_numeric <- as.numeric(str_replace(lu.s$treat, "%", ""))/100

fail.s <- str_replace(fail.s, 'transcript stability_Lo1', 'shuffle-Lo_50%_rep1')
fail.s <- str_replace(fail.s, 'max', '1.25%')
fail.dat <- data.frame(Construct=fail.s,
                       date=NA,
                       sample=1,
                       fitness=0, 
                       sci='lu', 
                       gene=str_split(fail.s, '_', simplify = T)[,1],
                       type=str_split(fail.s, '_', simplify = T)[,2],
                       treat=str_split(fail.s, '_', simplify = T)[,3],
                       rep=str_split(fail.s, '_', simplify = T)[,4],
                       treat_numeric=str_split(fail.s, '_', simplify = T)[,3])
fail.dat$type <- str_replace(fail.dat$type,'-','_')
names(fail.dat) <- c('construct', 'date', 'sample', 'fitness', 'sci', 'gene', 'type', 'treat', 'rep', 'treat_numeric')
fail.dat$treat_numeric <- as.numeric(str_extract(fail.dat$treat_numeric, '\\d+'))/100
fail.dat[fail.dat$treat_numeric == 0.01,]$treat_numeric <- 1.25

emma.s <- emma.s[,names(lu.s)]
fail.dat <- fail.dat[,names(lu.s)]
fit.all <- rbind(emma.s, lu.s, fail.dat)
fit.all <- fit.all[!fit.all$fitness < 0,]
fit.s <- fit.all[,c("fitness", "date", "sci", "gene", "type", "treat", "rep")]
fit.s$date <- as.character(fit.s$date)

#write.csv(fit.all, 'results/fitness_all.csv', quote = F, row.names = F)


#explore data before modeling
ggplot(fit.all, aes(x=treat_numeric, y=fitness, color=type, shape=gene)) +
  geom_jitter(width=0.05, size=2, alpha=0.75) +
  theme_classic() +
  scale_colour_colorblind() +
  scale_fill_colorblind()
#ggsave(filename = "plots/all_oneplot.pdf", width = 8, height= 6, units = "in")
#ggsave(filename = "plots/all_oneplot.png", width = 8, height= 6, units = "in")

fit.s.nostab <- fit.s[!str_detect(fit.s$type, "Lo"),]
fit.s.nostab <- fit.s.nostab[!str_detect(fit.s.nostab$type, "Hi"),]
#fit.s.nostab[fit.s.nostab$type == "wt",]$type <- NA


size=3
p.shuf <- ggplot(fit.s.nostab[fit.s.nostab$type %in% c("wt","shuffle"),], aes(x=treat, y=fitness, shape=gene, color=gene, group=gene))+
  geom_jitter(position = position_jitterdodge(dodge.width=0.4, jitter.width=0.1, jitter.height = 0), size=size, alpha=0.75) +
  #facet_grid(type ~., margins = T, scales = "fixed", drop = T) +
  scale_x_discrete(limits=c("0%", "25%", "50%", "75%", "100%", "max")) +
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank()) +
  ggtitle("Shuffle recoding")
p.rand <- ggplot(fit.s.nostab[fit.s.nostab$type %in% c("wt","random"),], aes(x=treat, y=fitness, shape=gene, color=gene, group=gene))+
  geom_jitter(position = position_jitterdodge(dodge.width=0.4, jitter.width=0.1, jitter.height = 0), size=size, alpha=0.75) +
  #facet_grid(type ~., margins = T, scales = "fixed", drop = T) +
  scale_x_discrete(limits=c("0%", "25%", "50%", "75%", "100%", "max")) +
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank()) +
  ggtitle("Random recoding")
p.worst <- ggplot(fit.s.nostab[fit.s.nostab$type %in% c("wt","worst"),], aes(x=treat, y=fitness, shape=gene, color=gene, group=gene))+
  geom_jitter(position = position_jitterdodge(dodge.width=0.4, jitter.width=0.1, jitter.height = 0), size=size, alpha=0.75) +
  #facet_grid(type ~., margins = T, scales = "fixed", drop = T) +
  scale_x_discrete(limits=c("0%", "25%", "50%", "75%", "100%", "max")) +
  theme_bw() +
  ylim(0,25) +
  ggtitle("Worst codon")
plot_grid( p.shuf, p.rand, p.worst, ncol = 1)


fit.s.stab <- fit.s[fit.s$type %in% c("wt", "shuffle", "shuffle_Hi", "shuffle_Lo"),]
fit.s.stab <- fit.s.stab[fit.s.stab$treat %in% c("0%","50%"),]
p.shuff <- ggplot(fit.s.stab, aes(x=type, y=fitness, shape=gene, color=gene, group=gene))+
  geom_jitter(position = position_jitterdodge(dodge.width=0.4, jitter.width=0.1, jitter.height = 0), size=size, alpha=0.75) +
  #facet_grid(type ~., margins = T, scales = "fixed", drop = T) +
  scale_x_discrete(limits=c("wt", "shuffle_Lo", "shuffle", "shuffle_Hi")) +
  theme_bw() +
  ylim(0,23) +
  ggtitle("Stability")

plot_grid(p.shuf, p.rand, p.worst, p.shuff, ncol = 1)
#ggsave(filename = "plots/all.pdf", width = 8.5, height= 11, units = "in")
#ggsave(filename = "plots/all.png", width = 8.5, height= 11, units = "in")

#figure out how to look at assumptions for model
ggplot(fit.s[!fit.s$gene=="wt",], aes(fitness, fill=treat)) +
  geom_histogram() +
  facet_grid(treat ~ ., margins = T, scales = "free") +
  theme_bw()

with(fit.s[!fit.s$gene=="wt",], tapply(fitness, treat, function(x) {
  sprintf("M (SD) = %1.2f (%1.2f)", mean(x), sd(x))
}))
#means are greater than errors

library(MASS)
mod.1 <- lm(fitness ~ treat + type + gene, data=fit.s)
mod.2 <- lm(fitness ~ treat + type + gene + date + sci, data=fit.s)
AIC(mod.1)
AIC(mod.2)
logLik(mod.1)
logLik(mod.2)
anova(mod.1, mod.2)

fit.s$mod.fit <- mod.2$fitted.values
plot(fit.s$fitness, fit.s$mod.fit)

fit.summary <- fit.all[which(!is.na(fit.all$date)),] %>%   ##dropping constructs that couldn't be made 
  group_by(gene, type, treat_numeric) %>%
  summarise(mean = mean(fitness),
            n = n(),
            sd = sd(fitness),
            se = sd/sqrt(n))
#write.csv(fit.summary, 'results/fit_summary.csv', quote = F, row.names = F)
ggplot(fit.summary, aes(x=mean, y=se, fill=gene, color=gene, shape=gene)) +
  geom_point(alpha=0.75) +
  theme_bw() +
  scale_color_jama() +
  xlab('Mean fitness') +
  ylab('Fitness SE')
#ggsave(filename = "plots/errorByfitness.pdf", width = 4.5, height= 4, units = "in")
#ggsave(filename = "plots/errorByfitness.png", width = 4.5, height= 4, units = "in")

add <- fit.all[which(is.na(fit.all$date)),
               c('gene', 'type', 'treat_numeric', 'fitness')]
add$n <- 1 
add$sd <- 0  
add$se <- 0


ggplot() +
  geom_pointrange(data=fit.summary, 
                  aes(x=treat_numeric, 
                      y=mean, 
                      color=type, 
                      shape=gene, 
                      ymin=mean-sd, 
                      ymax=mean+sd), 
                  alpha=0.75, 
                  position=position_jitter(width=0.07)) +
  geom_point(data=add, aes(x=treat_numeric, y=fitness, color=type, shape=gene),
             alpha=0.75, 
             position=position_jitter(width=0.07),
             size=2) +
  theme_classic() +
  scale_colour_colorblind() +
  scale_fill_colorblind() +
  scale_shape_manual(values = c(15,16,17,18)) +
  xlab("Proportion of codons changed") +
  ylab('Fitness (doublings/hr)')
#ggsave(filename = "plots/all_oneplot.pdf", width = 6, height= 4.5, units = "in")
#ggsave(filename = "plots/all_oneplot.png", width = 6, height= 4.5, units = "in")

ggplot(data = df1, aes(x, y, color = Group)) +
  geom_pointrange(aes(ymin = lb, ymax = ub), 
                  position=position_jitter(width=0.5), 
                  linetype='dotted') +
  theme_bw()
