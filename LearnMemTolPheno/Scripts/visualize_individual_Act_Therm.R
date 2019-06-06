library(ggplot2)
library(cowplot)
library(tidyverse)

source(file="../../Functions/ggplot_theme.R")

rawdLM <- read.table("../ProcessedData/LearnMem_processed.txt",sep="\t", header=TRUE, stringsAsFactors = FALSE)
rawdLM$patRILF <- as.factor(rawdLM$patRIL)

rawdLM <- rawdLM %>% mutate(LearnNorm = quant.norm(Learning),
                            MemoryNorm = quant.norm(Memory),
                            ActNorm = quant.norm(Activity_score))

mod.la <- lm(Learning ~ Activity_score, data=rawdLM)
summary(mod.la)

mod.ma <- lm(Memory ~ Activity_score, data=rawdLM)
summary(mod.ma)

mod.lt <- lm(Learning ~ incapacitation, data=rawdLM)
summary(mod.lt)

mod.mt <- lm(Memory ~ incapacitation, data=rawdLM)
summary(mod.mt)


la <- ggplot(rawdLM, aes(x=Activity_score, y=Learning)) +
  geom_point(size=0.25, alpha=1/5) +
  annotate("text", x = 245, y = 1.1, size=3,label = "italic(R) ^ 2 == 0.07", parse=TRUE) +
  theme(legend.position="none") +
  xlab("Activity Score")+
  my_theme

ma <- ggplot(rawdLM, aes(x=Activity_score, y=Memory)) +
  geom_point(size=0.25, alpha=1/5) +
  annotate("text", x = 245, size=3, y = 1.1, label = "italic(R) ^ 2 == 0.02", parse=TRUE) +
  theme(legend.position="none") +
  xlab("Activity Score")+
  my_theme

lt <- ggplot(rawdLM, aes(x=incapacitation, y=Learning)) +
  geom_point(size=0.25, alpha=1/5) +
  annotate("text", x = 525, size=3, y = 1.1, label = "italic(R) ^ 2 == 0.002", parse=TRUE) +
  theme(legend.position="none") +
  xlab("Thermal Tolerance")+
  my_theme

mt <- ggplot(rawdLM, aes(x=incapacitation, y=Memory)) +
  geom_point(size=0.25, alpha=1/5) +
  theme(legend.position="none") +
  annotate("text", x = 525, y = 1.1, label = "italic(R) ^ 2 == 0.005", parse=TRUE, size=3) +
  xlab("Thermal Tolerance") +
  xlim(c(0,600)) +
  my_theme 
  
pall <- plot_grid(la, ma, lt, mt, labels=c('a.','b.','c.','d.'), label_size = 10, nrow=2)

ggsave(pall, file="../Plots/Pheno_ind_all.pdf", width=6.5, height=6.5)
