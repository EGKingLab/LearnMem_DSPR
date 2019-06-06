library(DSPRqtl)
library(nlme)
library(ggplot2)
library(cowplot)
library(tidyverse)

source(file="../../Functions/ggplot_theme.R")
source("../../Functions/mappingfunctions.R")
source("../../Functions/h2_functs.R")

load(file= "../ProcessedData/sig_ths.rda")

line.ww <- 0.75

rawdLM <- read.table("../../LearnMemTolPheno/ProcessedData/LearnMem_processed.txt",sep="\t", header=TRUE, stringsAsFactors = FALSE)
rawdLM$patRILF <- as.factor(rawdLM$patRIL)
  
rawdLM <- rawdLM %>% mutate(LearnNorm = quant.norm(Learning),
                            MemoryNorm = quant.norm(Memory),
                            ActNorm = quant.norm(Activity_score))

#mod <- lme(Learning ~ incapacitation, random = ~1|patRILF, data=rawdLM)
mod <- lm(Learning ~ incapacitation, data=rawdLM)
summary(mod)
rawdLM$LearnAC <- resid(mod)

#mod <- lme(Memory ~ incapacitation, random = ~1|patRILF, data=rawdLM)
mod <- lm(Memory ~ incapacitation, data=rawdLM)
summary(mod)
rawdLM$MemAC <- resid(mod)


L_MDATA <- rawdLM %>% 
  group_by(patRIL) %>% 
  summarise(Learning_Mean = mean(Learning),
            Memory_Mean = mean(Memory),
            LearnAC_Mean =mean(LearnAC),
            MemAC_Mean =mean(MemAC),
            Act_Mean = mean(Activity_score))
            
L_MDATA$LearnPowerTrans <- L_MDATA$Learning_Mean^3
L_MDATA$MemACTrans <- quant.norm(L_MDATA$MemAC_Mean)

All_3_Traits_Data <- as.data.frame(L_MDATA)

#Fetch founder hyplotype probabilty 
myGenos <- DSPRgenos(design='inbredA', All_3_Traits_Data, id.col='patRIL')

#Look at the possitions. These are in 5.x coordinates
positions <- myGenos$positions

#Put column names of phenotypes here
pp<-myGenos[['phenotype']][,c('LearnPowerTrans','Memory_Mean','LearnAC_Mean', 'MemAC_Mean', 'Act_Mean','MemACTrans')]


Null.Mods<-logLik.multi(lm(as.matrix(pp)~1))

Null.Ls<-unlist(Null.Mods)/log(10)

gg <- myGenos$genolist
Phen.Mods<-lapply(gg, function(mat) logLik.multi(lm(as.matrix(pp) ~ mat[,1] + mat[,2]+ mat[,3]+ mat[,4]+ mat[,5]+ mat[,6]+ mat[,7])))

Phen.Ls<-lapply(Phen.Mods, function(xx) (xx/log(10)) - Null.Ls)

all.LL <- array(dim=c(length(Phen.Ls),dim(Phen.Ls[[1]])[1]))

for(zz in seq(along=Phen.Ls)) all.LL[zz,] <- Phen.Ls[[zz]]
colnames(all.LL)<-colnames(pp)

obs.LL <- array(dim=c(length(Phen.Ls),dim(Phen.Ls[[1]])[1]))

for(zz in seq(along=Phen.Ls)) obs.LL[zz,] <- Phen.Ls[[zz]]
colnames(obs.LL)<-colnames(pp)

obs.LL <- cbind(obs.LL, positions)

save(obs.LL, file="../ProcessedData/ThermCorr_LODscores.rda")


Nrils <- nrow(myGenos$phenotype)

#set value for theresholds 
fdr <- min(fdr.out$fdr$threshold[fdr.out$fdr$fdr<=0.05])
fwer <- fdr.out$fwer

pfdr <- getP(fdr,Nrils,7,Nrils-8)
pfwer <- getP(fwer,Nrils,7,Nrils-8)


peak.i <- apply(obs.LL[,1:4],2,function(x) peakInfo(x,cM=poslist[,c('chr','Gpos')] ,th=5, tol.dist=3))

load(file="../ProcessedData/Peaks_wCIs.rda")

obs.LM_long <-gather(obs.LL, key='Trait', value='LOD', LearnPowerTrans, LearnAC_Mean)
obs.LM_long$Pval <- getP(obs.LM_long$LOD,Nrils,7,Nrils-8)
top <- max(obs.LM_long$Pval)

p1 <- ggplot(obs.LM_long, aes(x=Gaxis, y=Pval, color=Trait)) +
  geom_rect(xmin=66.3,ymin=-10,xmax=120.3,ymax=top+10,fill='aliceblue',color='aliceblue')+
  geom_rect(xmin=174,ymin=-10,xmax=221,ymax=top+10,fill='aliceblue',color='aliceblue')+
  geom_segment(x=-50, xend=275,y=pfdr,yend=pfdr, color='grey30', size=line.ww) +
  annotate("text",label= "FDR", x=295, y=pfdr, color='grey30',size=3) +
  annotate("text",label= "FWER", x=295, y=pfwer, color='grey30',size=3) +
  geom_segment(x=-50, xend=275,y=pfwer,yend=pfwer, color='grey30', size=line.ww, lty=3) +
  geom_line(alpha=0.6, size=line.ww) +
  scale_x_continuous(breaks=c(0,66.3,120,174,221,277,33.15,98.3,145,205,249), 
                     labels=c(0,"66  0",54,"108  0",47,103,'\nX','\n2L','\n2R','\n3L','\n3R'),limits = c(-0.01,max(obs.LM_long$Gaxis)+25))+
  
  theme(axis.ticks.x=element_blank())+
  ylab(expression("-log"[10]*"(P-value)"))+
  xlab("Position (cM)") +
  scale_color_manual(name='Phenotype',breaks=c('LearnPowerTrans','LearnAC_Mean'), labels=c('Learning','Learning (corrected)'), values=c('red','black'))+
  my_theme

#ggsave(p1, filename = "../Plots/Learn_correct_QTL.pdf", height=3,width=7)

obs.LM_long <-gather(obs.LL, key='Trait', value='LOD', Memory_Mean, MemAC_Mean)
obs.LM_long$Pval <- getP(obs.LM_long$LOD,Nrils,7,Nrils-8)
top <- max(obs.LM_long$Pval)

p2 <- ggplot(obs.LM_long, aes(x=Gaxis, y=Pval, color=Trait)) +
  geom_rect(xmin=66.3,ymin=-10,xmax=120.3,ymax=top+10,fill='aliceblue',color='aliceblue')+
  geom_rect(xmin=174,ymin=-10,xmax=221,ymax=top+10,fill='aliceblue',color='aliceblue')+
  geom_segment(x=-50, xend=275,y=pfdr,yend=pfdr, color='grey30', size=line.ww) +
  annotate("text",label= "FDR", x=295, y=pfdr, color='grey30',size=3) +
  annotate("text",label= "FWER", x=295, y=pfwer, color='grey30',size=3) +
  geom_segment(x=-50, xend=275,y=pfwer,yend=pfwer, color='grey30', size=line.ww, lty=3) +
  geom_line(alpha=0.6, size=line.ww) +
  scale_x_continuous(breaks=c(0,66.3,120,174,221,277,33.15,98.3,145,205,249), 
                     labels=c(0,"66  0",54,"108  0",47,103,'\nX','\n2L','\n2R','\n3L','\n3R'),limits = c(-0.01,max(obs.LM_long$Gaxis)+25))+
  
  theme(axis.ticks.x=element_blank())+
  ylab(expression("-log"[10]*"(P-value)"))+
  xlab("Position (cM)") +
  scale_color_manual(name='Phenotype',breaks=c('Memory_Mean','MemAC_Mean'), labels=c('Memory','Memory (corrected)'),values=c('red','black'))+
  my_theme

#ggsave(p1, filename = "../Plots/Memory_correct_QTL.pdf", height=3,width=7)

pall <- plot_grid(p1,p2, labels=c("a.","b."), nrow=2, label_size = 10)
ggsave(pall, filename = "../Plots/LM_Tcorrect_QTL.pdf", height=6.5,width=6.5)

ind <- 8119
obs.LL[(ind-25):(ind+25),]
fdr.out
