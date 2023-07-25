#### SCR functions

# Distance between traps and activity centers
# Faster than using loops
e2dist <- function (x, y)
{
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}



simSCR0 <-
  function(N=100,K=20,alpha0=-2.5,sigma=.5, discard0=TRUE,array3d=FALSE,rnd=2013){
    set.seed(rnd)
    
    # make trapping grid. Normally you would provide a 2-dimensional matrix
    #   of trap coordinates and read it in like this:
    #    tralocs<-read.csv("traplocs.csv") or similar
    traplocs<- cbind(sort(rep(1:5,5)),rep(1:5,5))
    Dmat<-e2dist(traplocs,traplocs)
    ntraps<-nrow(traplocs)
    plot(traplocs)
    
    # define state-space of point process. (i.e., where animals live).
    # Here "delta" just adds
    # a fixed buffer to the outer extent of the traps.
    buffer<-2
    Xl<-min(traplocs[,1] - buffer)
    Xu<-max(traplocs[,1] + buffer)
    Yl<-min(traplocs[,2] - buffer)
    Yu<-max(traplocs[,2] + buffer)
    
    sx<-runif(N,Xl,Xu)
    sy<-runif(N,Yl,Yu)
    S<-cbind(sx,sy)
    
    # how far is each individual from each trap?
    D<- e2dist(S,traplocs)
    
    #alpha0<- -2.5
    #sigma<- 0.5
    alpha1<- 1/(2*sigma*sigma)
    
    #cloglog.probcap<- alpha0  - alpha1*D*D
    #probcap<- 1-exp(-exp(cloglog.probcap))
    # this is logit model here:
    #probcap<- expit(-2.5 - alpha1*D)
    
    probcap<-plogis(alpha0)*exp(-alpha1*D*D)
    # now generate the encounters of every individual in every trap
    Y<-matrix(NA,nrow=N,ncol=ntraps)
    for(i in 1:nrow(Y)){
      Y[i,]<-rbinom(ntraps,K,probcap[i,])
    }
    
    ## NOTE
    ## Y is a matrix of encounter frequencies of EACH individual in EACH trap
    ## As simulated here it includes the "all 0" observations.  We want
    ## to delete those to mimic real data.
    if(discard0){
      totalcaps<-apply(Y,1,sum)
      Y<-Y[totalcaps>0,]
    }
    
    dimnames(Y)<-list(1:nrow(Y),paste("trap",1:ncol(Y),sep=""))
    
    if(array3d){
      
      ## Here we demonstrate how to simulate the full 3-dimensional
      ## encounter history
      # now generate the encounters of every individual in every trap
      
      Y<-array(NA,dim=c(N,ntraps, K))
      for(i in 1:nrow(Y)){
        for(j in 1:ntraps){
          Y[i,j,1:K]<-rbinom(K,1,probcap[i,j])
        }
      }
      
      if(discard0){
        ## create nind x ntraps array
        Y2d<- apply(Y,c(1,2),sum)
        ## which individuals were captured?
        ncaps<-apply(Y2d,1,sum)
        ## keep those ones that were captured
        Y<-Y[ncaps>0,,]
      }
      
    }
    
    
    list(Y=Y,traplocs=traplocs,xlim=c(Xl,Xu),ylim=c(Yl,Yu),N=N,alpha0=alpha0,alpha1=alpha1,sigma=sigma,K=K)
  }


### MLE
intlik1 <-
  function(parm,y=y,delta=.2,X=traplocs,ssbuffer=2){
    
    Xl<-min(X[,1])-ssbuffer
    Xu<-max(X[,1])+ ssbuffer
    Yu<-max(X[,2])+ ssbuffer
    Yl<-min(X[,2])-ssbuffer
    
    ## These commands set up the integration grid
    xg<-seq(Xl+delta/2,Xu-delta/2,by=delta)
    yg<-seq(Yl+delta/2,Yu-delta/2,by=delta)
    npix<-length(xg)
    G<-cbind(rep(xg,npix),sort(rep(yg,npix)))
    nG<-nrow(G)
    
    D<- e2dist(X,G)
    
    alpha0<-parm[1]
    ## alpha1 should probably be restricted to be positive here using exp()
    alpha1<-parm[2]
    probcap<- plogis(alpha0)*exp(-alpha1*D*D)
    
    
    Pm<-matrix(NA,nrow=nrow(probcap),ncol=ncol(probcap))
    # n0 here is the number of individuals with all-zero encounter histories
    # the next 3 lines replace ALL of the "all-0" encounter histories with a
    # single one, so that the calculation is only done once.
    n0<-sum(apply(y,1,sum)==0)
    ymat<-y[apply(y,1,sum)>0,]
    ymat<-rbind(ymat,rep(0,ncol(ymat)))
    lik.marg<-rep(NA,nrow(ymat))
    # for each encounter history, compute the likelihood for each value of s
    # then average over all possible values
    for(i in 1:nrow(ymat)){
      Pm[1:length(Pm)]<- (dbinom(rep(ymat[i,],nG),K,probcap[1:length(Pm)],log=TRUE))
      lik.cond<- exp(colSums(Pm))
      lik.marg[i]<- sum( lik.cond*(1/nG))
    }
    # nv is a vector of the frequency of each encounter history. Here it is
    #   a vector of 1's with n0 being the number of all-0 histories.
    nv<-c(rep(1,length(lik.marg)-1),n0)
    -1*( sum(nv*log(lik.marg)) )
    
  }

# make an array from encounter data
SCR23darray <-
  function(edf, tdf){
    ### Returns 3-d array "ind x trap x occasion"
    
    # Check for dups
    uq<- paste(edf[,2],edf[,3],edf[,4])
    uq<- unique(uq)
    if(any(table(uq)>1)) cat("Duplicate captures (same individual, trap, occasion) present in data set, these are not used",fill=TRUE)
    
    nind<-max(edf[,2])
    ntraps<-nrow(tdf)
    nperiods<-ncol(tdf)-3
    per.id<- as.numeric(dimnames(tdf)[[2]][4:ncol(tdf)])
    
    ind.id<- edf[,2]
    trap.id<- edf[,4]
    
    if( length(per.id) != length(min(per.id):max(per.id)) ){
      x<- 1:nperiods
      names(x)<-as.character(per.id)
      per.id <- x[as.character(edf[,3])]
    }
    else{
      per.id<-edf[,3]
    }
    
    y<-array(0,c(nind,ntraps, nperiods))
    
    tmp<-cbind(ind.id,trap.id,per.id)
    y[tmp]<-1
    y
  }

## 3rd version of the MLE 
intlik3 <-
  function (start = NULL, y = y, K = NULL, delta = 0.3, X = traplocs,
            ssbuffer = 2,model="B",predict=FALSE)
  {
    Xl <- min(X[, 1]) - ssbuffer
    Xu <- max(X[, 1]) + ssbuffer
    Yu <- max(X[, 2]) + ssbuffer
    Yl <- min(X[, 2]) - ssbuffer
    #SSarea <- (Xu - Xl) * (Yu - Yl)
    if (is.null(K))
      return("need sample size")
    xg <- seq(Xl + delta/2, Xu - delta/2, delta)
    yg <- seq(Yl + delta/2, Yu - delta/2, delta)
    npix.x <- length(xg)
    npix.y <- length(yg)
    area <- (Xu - Xl) * (Yu - Yl)/((npix.x) * (npix.y))
    G <- cbind(rep(xg, npix.y), sort(rep(yg, npix.x)))
    nG <- nrow(G)
    D <- e2dist(X, G)
    SSarea<- (delta*delta)*nrow(G)
    if (is.null(start))
      start <- c(0, 0, 0)
    alpha0 <- start[1]
    alpha1 <- exp(start[2])
    n0 <- exp(start[3])
    if(model=="B")
      probcap <- plogis(alpha0) * exp(-alpha1 * D * D)
    if(model=="P")
      probcap <- exp(alpha0) * exp(-alpha1 * D * D)
    Pm <- matrix(NA, nrow = nrow(probcap), ncol = ncol(probcap))
    ymat <- y
    ymat <- rbind(y, rep(0, ncol(y)))
    lik.marg <- rep(NA, nrow(ymat))
    for (i in 1:nrow(ymat)) {
      if(model=="B")
        Pm[1:length(Pm)] <- (dbinom(rep(ymat[i, ], nG), rep(K, nG), probcap[1:length(Pm)], log = TRUE))
      if(model=="P")
        Pm[1:length(Pm)] <- (dpois(rep(ymat[i, ], nG), rep(K, nG)*probcap[1:length(Pm)], log = TRUE))
      lik.cond <- exp(colSums(Pm))
      lik.marg[i] <- sum(lik.cond * (1/nG))
    }
    if(predict==FALSE){
      nv <- c(rep(1, length(lik.marg) - 1), n0)
      part1 <- lgamma(nrow(y) + n0 + 1) - lgamma(n0 + 1)
      part2 <- sum(nv * log(lik.marg))
      out <- -1 * (part1 + part2)
      attr(out, "SSarea") <- SSarea
      return(out)
    }
    if(predict==TRUE){
      
      posterior<-matrix(NA,nrow=nG,ncol=nrow(ymat))
      for(i in 1:nrow(ymat)){
        if(model=="B")
          Pm[1:length(Pm)] <- (dbinom(rep(ymat[i, ], nG), rep(K, nG), probcap[1:length(Pm)], log = TRUE))
        if(model=="P")
          Pm[1:length(Pm)] <- (dpois(rep(ymat[i, ], nG), rep(K, nG)*probcap[1:length(Pm)], log = TRUE))
        
        lik.cond <- exp(colSums(Pm))*(1/nG)
        posterior[,i]<- lik.cond/lik.marg[i]
      }
      
      return(cbind(G,posterior))
    }
    
    
  }

# Convert long format occasion data to 3-d 'ind x trap x occasion'
SCR23darray <-
  function(edf, tdf){
    ### Returns 3-d array "ind x trap x occasion"
    
    # Check for dups
    uq<- paste(edf[,2],edf[,3],edf[,4])
    uq<- unique(uq)
    if(any(table(uq)>1)) cat("Duplicate captures (same individual, trap, occasion) present in data set, these are not used",fill=TRUE)
    
    nind<-max(edf[,2])
    ntraps<-nrow(tdf)
    nperiods<-ncol(tdf)-3
    per.id<- as.numeric(dimnames(tdf)[[2]][4:ncol(tdf)])
    
    ind.id<- edf[,2]
    trap.id<- edf[,4]
    
    if( length(per.id) != length(min(per.id):max(per.id)) ){
      x<- 1:nperiods
      names(x)<-as.character(per.id)
      per.id <- x[as.character(edf[,3])]
    }
    else{
      per.id<-edf[,3]
    }
    
    y<-array(0,c(nind,ntraps, nperiods))
    
    tmp<-as.data.frame(cbind(ind.id,trap.id,per.id))
    twerp <- matrix(unlist(tmp), nrow = nrow(tmp), ncol = 3 )
    
    y[twerp]<-1
    y
  }

# Convert long format occasion data to 3-d 'ind x trap x occasion'
# where number of encounters are preserved
SCRENC23darray <-
  function(edf, tdf){
    ### Returns 3-d array "ind x trap x occasion"
    
    # Check for dups
    uq<- paste(edf[,2],edf[,3],edf[,4])
    uq<- unique(uq)
    if(any(table(uq)>1)) cat("Duplicate captures (same individual, trap, occasion) present in data set, these are not used",fill=TRUE)
    
    nind<-max(edf[,2])
    ntraps<-nrow(tdf)
    nperiods<-ncol(tdf)-3
    per.id<- as.numeric(dimnames(tdf)[[2]][4:ncol(tdf)])
    
    ind.id <- edf[,2]
    trap.id <- edf[,4]
    
    if( length(per.id) != length(min(per.id):max(per.id)) ){
      x<- 1:nperiods
      names(x)<-as.character(per.id)
      per.id <- x[as.character(edf[,3])]
    }
    else{
      per.id<-edf[,3]
    }
    
    y<-array( 0,c(nind, ntraps, nperiods) )
    
    tmp <- as.data.frame( cbind(ind.id,trap.id,per.id) )
    
    # Tally number of encounters
    tmpp <- as_tibble(tmp) %>% 
      group_by( ind.id, trap.id ) %>% 
      summarize_at( dplyr::vars(`per.id`), count ) %>% 
      mutate( sampp.id = per.id$x,
              enc = per.id$freq ) %>% 
      ungroup( ) 
    
    tmprar <-as.data.frame(cbind( tmpp$ind.id, tmpp$trap.id, 
                                  tmpp$sampp.id, tmpp$enc) )
    
    twerp <- matrix(unlist(tmprar), nrow = nrow(tmprar), ncol = 4 )
    
    y[ twerp[,c(1:3)] ] <- twerp[,4]
    y
  }
######################################################################
## spiderplot function
######################################################################
spiderplot<-function (y, traplocs) 
{
  dither<-FALSE
  dx<- max(traplocs[,1])-min(traplocs[,1])
  dy<- max(traplocs[,2])-min(traplocs[,2])
  dx<- .01*dx
  dy<- .01*dy
  
  
  
  
  
  if (length(dim(y)) == 3) {
    
    if(dim(y)[2]==nrow(traplocs)){
      # need to transpose here
      nind<-dim(y)[1]
      ntraps<-dim(y)[2]
      nocc<- dim(y)[3]
      newy<-array(NA,dim=c(nind,nocc,ntraps))
      for(i in 1:nind){
        newy[i,1:nocc,1:ntraps]<- t(y[i,,])
      }
      y<-newy
      
    }
    
    
    y3d <- y
    J <- dim(y3d)[3]
    T <- dim(y3d)[2]
    nind <- dim(y3d)[1]
    plot(traplocs, pch = 20, xlab = " ", ylab = " ",cex=1.5)
    avg.s <- matrix(NA, nrow = nind, ncol = 2)
    for (i in 1:nind) {
      tmp <- NULL
      for (t in 1:T) {
        aa <- y3d[i, t, ]
        if (sum(aa) > 0) {
          aa <- traplocs[aa > 0, ]
          tmp <- rbind(tmp, aa)
        }
      }
      avg.s[i, ] <- c(mean(tmp[, 1]), mean(tmp[, 2]))
      delta<- c(runif(1,-dx,dx),runif(1,-dy,dy))*ifelse(dither,1,0)
      points(avg.s[i, 1]+delta, 
             avg.s[i, 2]+delta, 
             pch = "S", cex = 1, 
             col = "red")
      for (m in 1:nrow(tmp)) {
        if (nrow(tmp) > 1) 
          lines(c(avg.s[i, 1], tmp[m, 1]), c(avg.s[i, 
                                                   2], tmp[m, 2]))
      }
    }
  }
  if (length(dim(y)) == 2) {
    y2d <- y
    J <- nrow(traplocs)
    T <- dim(y2d)[2]
    nind <- dim(y2d)[1]
    plot(traplocs, pch = 20, xlab = " ", ylab = " ",cex=1.5)
    avg.s <- matrix(NA, nrow = nind, ncol = 2)
    for (i in 1:nind) {
      tmp <- NULL
      for (t in 1:T) {
        aa <- y2d[i, t]
        if (aa <= J) {
          aa <- traplocs[aa, ]
          tmp <- rbind(tmp, aa)
        }
      }
      avg.s[i, ] <- c(mean(tmp[, 1]), mean(tmp[, 2]))
      points(avg.s[i, 1], avg.s[i, 2], pch = "S", cex = 1, 
             col = "red")
      for (m in 1:nrow(tmp)) {
        if (nrow(tmp) > 1) 
          lines(c(avg.s[i, 1], tmp[m, 1]), c(avg.s[i, 
                                                   2], tmp[m, 2]))
      }
    }
  }
  points(traplocs, pch = 20)
  Cx <- mean(traplocs[, 1])
  Cy <- mean(traplocs[, 2])
  xcent <- sqrt((avg.s[, 1] - Cx)^2 + (avg.s[, 2] - Cy)^2)
  list(xcent = xcent, avg.s = avg.s, center = c(Cx, Cy))
}