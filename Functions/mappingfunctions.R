logLik.multi <- function(object, REML = FALSE, ...)
{
  
  all.val<-numeric(length=0)
  for(j in 1:ncol(object$residuals))
  {
    res <- object$residuals[,j] # not resid(object) because of NA methods
    p <- object$rank
    N <- length(res) 
    if(is.null(w <- object$weights)) {
      w <- rep.int(1, N)
    } else {
      ## this is OK as both resids and weights are for the cases used
      excl <- w == 0      # eliminating zero weights
      if (any(excl)) {
        res <- res[!excl]
        N <- length(res)
        w <- w[!excl]
      }
    }
    N0 <- N
    if(REML) N <- N - p
    val <- .5* (sum(log(w)) - N * (log(2 * pi) + 1 - log(N) +
                                     log(sum(w*res^2))))
    if(REML) val <- val - sum(log(abs(diag(object$qr$qr)[1L:p])))
    attr(val, "nall") <- N0 # NB, still omits zero weights
    attr(val, "nobs") <- N
    attr(val, "df") <- p + 1
    class(val) <- "logLik"
    all.val<-rbind(all.val,val)
  }#for close
  all.val
}

pfind<-function(pp, cM, th, tol.dist)
{
  
  d1<-c(diff(pp),0)
  d2<- c(0,-d1[-length(d1)])
  
  init.p<- poslist[which(pp > th & d1 <= 0 & d2 <= 0),]
  init.p$LL<-pp[which(pp > th & d1 <= 0 & d2 <= 0)]
  init.p<-init.p[order(-init.p$LL),]
  #if none - return 0
  cc<-1
  for(kk in 1:nrow(init.p))
  {
    if(cc<=nrow(init.p))
    {
      ind.e<-which(init.p[,'chr']==init.p[cc,'chr'] & abs(init.p[,'Gpos'] - init.p[cc, 'Gpos']) < tol.dist)
      ind.e<-ind.e[-which(ind.e==cc)]
      if(length(ind.e)>0)
      {
        init.p<-init.p[-ind.e,,drop=FALSE]
      }
      cc<-cc+1
    }
  }
  return(nrow(init.p))
  
}

getP<-function(LOD,n,df,ddf)
{
  ff<-(10^((2/n)*LOD)-1)*((ddf)/df)
  pp<- -(pf(ff,df,ddf,lower.tail=FALSE,log=TRUE)/log(10))
  #pp<- pf(ff,df,ddf,lower.tail=FALSE)
  return(pp)
}

peakInfo<-function(pp, cM, th, tol.dist)
{
  
  d1<-c(diff(pp),0)
  d2<- c(0,-d1[-length(d1)])
  
  init.p<- poslist[which(pp > th & d1 <= 0 & d2 <= 0),]
  init.p$LL<-pp[which(pp > th & d1 <= 0 & d2 <= 0)]
  init.p<-init.p[order(-init.p$LL),]
  #if none - return 0
  cc<-1
  for(kk in 1:nrow(init.p))
  {
    if(cc<=nrow(init.p))
    {
      ind.e<-which(init.p[,'chr']==init.p[cc,'chr'] & abs(init.p[,'Gpos'] - init.p[cc, 'Gpos']) < tol.dist)
      ind.e<-ind.e[-which(ind.e==cc)]
      if(length(ind.e)>0)
      {
        init.p<-init.p[-ind.e,,drop=FALSE]
      }
      cc<-cc+1
    }
  }
  return(init.p)
  
}


Conv_5_6<-function(ccc,ppp,coord.table)
  
{
  return(coord.table[which(coord.table$R5chr==ccc & coord.table$R5pos==ppp),c('R6chr','R6pos')])
}


LOD_R2<-function(LOD,n)
{
  return(1-(10^((-2/n)*LOD)))
  
}
