library(DSPRqtl)
library(tidyverse)
data("positionlist_wgenetic")
load(file = "../ProcessedData/DSPR_5_6_parsed.rda")

source("../../Functions/mappingfunctions.R")

load(file="../ProcessedData/sig_ths.rda")#fdr.out
load(file="../ProcessedData/Lodscores_3traits.rda")#obs.LL



#th.2 <- fdr.out$fwer
th.l <- fdr.out$fdr$threshold[min(which(fdr.out$fdr$fdr<=0.05))]



#find peaks
peak.i <- apply(obs.LL[,1:3],2,function(x) peakInfo(x,cM=poslist[,c('chr','Gpos')] ,th=th.l, tol.dist=3))



ci.peak <- vector(mode="list", length=3)
names(ci.peak)<-names(peak.i)

for(kk in 1:3){
qq <- poslist[,1:3]
qq$LOD <- obs.LL[,kk]
qq$chr_b <- qq$chr
pp.set <- peak.i[[kk]]
pp.set$lp <- NA
pp.set$up <- NA
pp.set$lpchr <- NA
pp.set$upchr <- NA
pp.set$lg <- NA
pp.set$ug <- NA
pp.set$ulod <- NA
pp.set$llod <- NA
pp.set$chrR6 <- NA
pp.set$PposR6 <- NA
pp.set$lpR6 <- NA
pp.set$upR6 <- NA

pp.set$Ppos <- as.integer(pp.set$Ppos)



for(jj in 1:nrow(pp.set))
{
ff<-findCI(pp.set$chr[jj],pp.set$Ppos[jj], qtldat=qq, method='BCI')
pp.set$up[jj] <- ff$Ppos[2]
pp.set$lp[jj] <- ff$Ppos[1]
pp.set$upchr[jj] <- ff$chr_b[2]
pp.set$lpchr[jj] <- ff$chr_b[1]
pp.set$ug[jj] <- ff$Gpos[2]
pp.set$lg[jj] <- ff$Gpos[1]
pp.set$ulod[jj] <- ff$LOD[2]
pp.set$llod[jj] <- ff$LOD[1] 
pp.set$chrR6[jj] <- coord.table[which(coord.table$R5chr==pp.set$chr[jj] & coord.table$R5pos==as.integer(pp.set$Ppos[jj])),'R6chr']
pp.set$PposR6[jj] <- as.numeric(coord.table[which(coord.table$R5chr==pp.set$chr[jj] & coord.table$R5pos==pp.set$Ppos[jj]),'R6pos'])

pp.set$lpR6[jj] <- as.numeric(coord.table[which(coord.table$R5chr==pp.set$lpchr[jj] & coord.table$R5pos==ff$Ppos[1]),'R6pos'])
pp.set$upR6[jj] <- as.numeric(coord.table[which(coord.table$R5chr==pp.set$upchr[jj] & coord.table$R5pos==ff$Ppos[2]),'R6pos'])

} 

pp.set<-pp.set[order(pp.set$chr,pp.set$Ppos),]

#remove overlapping peaks
pp.set$chrN <- pp.set$chr
pp.set$chrN[pp.set$chrN=='2L']<-'2'
pp.set$chrN[pp.set$chrN=='2R']<-'2'
pp.set$chrN[pp.set$chrN=='3L']<-'3'
pp.set$chrN[pp.set$chrN=='3R']<-'3'
pp.set.final <- pp.set[0,]

for(pp in 1:nrow(pp.set))
{
  f.peak <- pp.set[pp,]
  f.rest <- pp.set[-pp,]
  ww <- which(f.rest$chrN==f.peak$chrN & (f.rest$lg <= f.peak$ug & f.rest$ug >= f.peak$lg) & f.rest$LL > f.peak$LL)
  if(length(ww)==0)
  {
    pp.set.final <- rbind(pp.set.final, f.peak)
  }
}
#(Learn_gene_list_sub$startp <= foc.peak$upR6 & Learn_gene_list_sub$stopp >= foc.peak$lpR6)))
ci.peak[[kk]]<-pp.set.final
}

#get hap means
load(file="../ProcessedData/myGenos.rda")
positions <- myGenos$positions
pp<-myGenos[['phenotype']][,c('LearnPowerTrans','Memory_Mean','Tolsqrtvariable')]
ci.peak.new <- vector(mode="list", length=3)
names(ci.peak.new)<-names(peak.i)

for(kk in 1:3){
  pp.set <- ci.peak[[kk]]
  pp.set[,19:42]<-NA
  colnames(pp.set)[19:42]<- c('A1','A2','A3','A4','A5','A6','A7','A8',
                              'A1se','A2se','A3se','A4se','A5se','A6se','A7se','A8se',
                              'A1N','A2N','A3N','A4N','A5N','A6N','A7N','A8N')
  pheno <- pp[,kk]
  for(jj in 1:nrow(pp.set))
  {
    foc.p <- pp.set[jj,]
    ww <- which(positions$chr==foc.p$chr & positions$Ppos==foc.p$Ppos)
    gg <- myGenos$genolist[[ww]]
    ccs <- summary(lm(pheno ~ gg -1))
    pp.set[jj,19:26] <- ccs$coefficients[,1]
    pp.set[jj,27:34] <- ccs$coefficients[,2]
    pp.set[jj,35:42] <- colSums(gg)
  }
 ci.peak.new[[kk]]<-pp.set 
}

ci.peak <- ci.peak.new
save(ci.peak,file="../ProcessedData/Peaks_wCIs.rda")


load(file="../ProcessedData/Peaks_wCIs.rda")
learn.dist <- ci.peak[[1]]$upR6-ci.peak[[1]]$lpR6
learn.dist <- learn.dist[-8]
mem.dist <- ci.peak[[2]]$upR6-ci.peak[[2]]$lpR6
mean(c(learn.dist, mem.dist))

learn.dist <- ci.peak[[1]]$ug-ci.peak[[1]]$lg
learn.dist <- learn.dist[-8]
mem.dist <- ci.peak[[2]]$ug-ci.peak[[2]]$lg
mean(c(learn.dist, mem.dist))

load(file="../../LearnMemTolPheno/ProcessedData/L_MDATA.rda")
NN <- nrow(L_MDATA)

r2.learn <- LOD_R2(ci.peak[[1]]$LL, NN)
r2.mem <- LOD_R2(ci.peak[[2]]$LL, NN)
mean(c(r2.learn,r2.mem))
max(c(r2.learn,r2.mem))
min(c(r2.learn,r2.mem))

