slideThermo <- function(xx, tt)
{
  st <- 301
  steps <- 50
  width <- 600
  #tol <- 1
  tol <- 5
  set <- 20
  disp <- 20
  dtol <-10
  
  while((set > tol) & (disp > dtol))
  {
    if(st<=(length(xx)-width)){
      set <- var(xx[st:(st+width)])
      disp <- max(abs(diff(xx[st:(st+width)],lag=100)))
      #cat(subR$time[st],"\t",subR$time[st+300],"\n")
      st <- st+50
    }else{
      set<-0
      st<-length(xx)+50
    }
  }
  
  return(tt[(st-50)])  
}

slideRec <- function(st, xx, tt)
{
  steps <- 50
  width <- 600
  
  iis<-seq(st,5701,by=50)
  vvs <- numeric(length(iis))
  for(ii in 1:length(iis))
  {
    st <- iis[ii]
    vvs[ii] <- var(xx[st:(st+width)])
  }
  
  return(max(vvs,na.rm=TRUE))  
}


#set as data frame
LL_Mem <-data.frame('patRIL'= numeric(length=0) , 'chamber'= numeric(length=0), 
                      'diff.test.A'= numeric(length=0), 
                      'diff.test.B'= numeric(length=0), 
                      'diff.train.A'= numeric(length=0), 
                      'diff.train.B'= numeric(length=0), 
                      'LearnAct'= numeric(length=0),
                      'ActStop'= numeric(length=0),
                      'ActRvar'= numeric(length=0),
                      'group'=character(length=0),
                      'file'=character(length=0), stringsAsFactors = FALSE)

basef <- "/home/pwilliams/DSPR/RawData/Learn_Mem_Raw/"
folds <-list.files(basef)

fcheck <- data.frame("file"=character(length=0),"rows"=numeric(length=0), stringsAsFactors=FALSE)


for(kk in folds)
{
  fiset<-list.files(file.path(basef,kk,fsep=""),pattern=".asc$")
  RR <- substr(kk,1,5)
  for(gg in fiset)
  {
    therm.set<-read.table(file.path(basef,kk,gg),sep="\t", header=TRUE, skip=29)
    fcheck <- rbind(fcheck,data.frame("file"=file.path(basef,kk,gg),"rows"=nrow(therm.set), stringsAsFactors=FALSE))
    all.dat <-data.frame('patRIL'= rep(RR,16) , 'chamber'= numeric(length=16), 
                         'diff.test.A'=numeric(length=16), 
                         'diff.test.B'=numeric(length=16), 
                         'diff.train.A'=numeric(length=16),
                         'diff.train.B'=numeric(length=16), 
                         'LearnAct'=numeric(length=16),
                         'ActStop'= numeric(length=16),
                         'ActRvar'= numeric(length=16),
                         'group'=character(length=16),
                         'file'=character(length=16),stringsAsFactors = FALSE)
    chambs<-seq(2,16)
    allR<-therm.set[,1:10]
    all.dat[1,'chamber'] <- allR$chamber[1]
    gg.1 <- substring(gg, 6)
    
    all.dat[1,'group'] <- as.numeric(gsub("[^\\d]+", "", gg.1,perl=TRUE))
    
    all.dat[1,'file'] <- gg
    #first set incap done above
    #target temp for accl - time 0 -300 target
    all.dat[1,'diff.test.A'] <-mean(abs(therm.set[1:301, 7]-therm.set[1:301, 8]))
    all.dat[1,'diff.test.B'] <-mean(abs(therm.set[1:301, 9]-therm.set[1:301, 10]))
    all.dat[1,'diff.train.A'] <-mean(abs(therm.set[302:5701, 7]-therm.set[302:5701, 8]))
    all.dat[1,'diff.train.B'] <-mean(abs(therm.set[302:5701, 9]-therm.set[302:5701, 10]))
    all.dat[1,'LearnAct'] <-var(allR$pos)
    incap <- slideThermo(xx=allR$pos, tt = allR$time)
    all.dat[1,'ActStop']<-incap
    if(incap > 569999){
      all.dat[1,'ActRvar'] <- NA 
    }else{
      all.dat[1,'ActRvar'] <- slideRec(which(allR$time==incap),xx=allR$pos, tt = allR$time)
    }
    
    counter<-2
    for(ii in chambs)
    {
      ss<-counter+9
      ee<-ss+8
      subR<-therm.set[,c(1,ss:ee)]
      colnames(subR)<-colnames(allR)
      all.dat[ii,'chamber'] <- subR$chamber[1]
      gg.1 <- substring(gg, 6)
      all.dat[ii,'group'] <- as.numeric(gsub("[^\\d]+", "", gg.1,perl=TRUE))
      all.dat[ii,'file'] <- gg
      all.dat[ii,'diff.test.A'] <-mean(abs(subR[1:301, 7]-subR[1:301, 8]),na.rm=TRUE)
      all.dat[ii,'diff.test.B'] <-mean(abs(subR[1:301, 9]-subR[1:301, 10]),na.rm=TRUE)
      all.dat[ii,'diff.train.A'] <-mean(abs(subR[302:5701, 7]-subR[302:5701, 8]),na.rm=TRUE)
      all.dat[ii,'diff.train.B'] <-mean(abs(subR[302:5701, 9]-subR[302:5701, 10]),na.rm=TRUE)
      all.dat[ii,'LearnAct'] <-var(subR$pos)
      incap <- slideThermo(xx=subR$pos, tt = subR$time)
      all.dat[ii,'ActStop']<-incap
      if(incap > 569999){
        all.dat[ii,'ActRvar'] <- NA 
      }else{
        all.dat[ii,'ActRvar'] <- slideRec(which(subR$time==incap),xx=subR$pos, tt = subR$time)
      }
      
      
      counter<-counter+9
      
    } #ii
    LL_Mem <- rbind(LL_Mem, all.dat)
    #cat(min(allR$pos),"\t", max(allR$pos),"\n")
    cat(gg, "\n")
  } #gg
}#kk

#save the R object
LearnMem <- LL_Mem
save(LearnMem, file="../ProcessedData/Learn_raw.rda")


hist(LL_Mem$diff.train.A)
hist(LL_Mem$diff.train.B)
hist(LL_Mem$diff.test.A)
hist(LL_Mem$diff.test.B)

ww <- which(LL_Mem$diff.train.A > 1|LL_Mem$diff.train.B>1)

LL_Mem.good<-LL_Mem[-ww,]

nums <- table(LL_Mem.good$patRIL)

#12127

