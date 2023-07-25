############################################################################
### Simulate spatial Jolly Seber
############################################################################
library(TMB)
library(tidyverse)
library(here)

source('SCR_Functions.R')

############################################
### Run a bunch of times to assess performance
# compile models outside of loop
# change working directory - TMB gets funny with here
setwd("CPP files")

# If you have problems compiling models below, it could be related to Rtools
# Might be helpful -> https://github.com/kaskr/adcomp/issues/215
# Try: system("gcc -v"), then: Sys.setenv(PATH="%PATH%;C:/Rtools/mingw_64/bin;c:/Rtools/bin")

compile("SCR_SearchEncounter_PoissonIntLik_multisession.cpp")
dyn.load(dynlib("SCR_SearchEncounter_PoissonIntLik_multisession"))

compile("SCR_SearchEncounter_PoissonIntLik_multisession_OPEN.cpp")
dyn.load(dynlib("SCR_SearchEncounter_PoissonIntLik_multisession_OPEN"))

compile("SCR_SearchEncounter_PoissonIntLik_multisession_OPEN_TagInt.cpp")
dyn.load(dynlib("SCR_SearchEncounter_PoissonIntLik_multisession_OPEN_TagInt"))

compile("SCR_SearchEncounter_PoissonIntLik_multisession_OPEN_NonSpatial.cpp")
dyn.load(dynlib("SCR_SearchEncounter_PoissonIntLik_multisession_OPEN_NonSpatial"))

compile("SCR_SearchEncounter_PoissonIntLik_multisession_OPEN_NoTE.cpp")
dyn.load(dynlib("SCR_SearchEncounter_PoissonIntLik_multisession_OPEN_NoTE"))

## Set up storage for each iteration
# to store simulated values
truePs <- matrix(NA, nrow=100, ncol=3)
truebs <- matrix(NA, nrow=100, ncol=19)
truePHIs <- matrix(NA, nrow=100, ncol=18)
trueIMMs <- matrix(NA, nrow=100, ncol=17)
trueNs <- matrix(NA, nrow=100, ncol=19)
trueSups <- matrix(NA, nrow=100, ncol=2)

K_sims <- NULL
eff_sims <- NULL
Y_sims <- NULL
CH_sims <- NULL
Tag_sims <- NULL
StateObs_sims <- NULL

# to store estimates from standard SCR version
estPs <- matrix(NA, nrow=100, ncol=3)
estNs <- matrix(NA, nrow=100, ncol=19)
# and sds
sdPs <- matrix(NA, nrow=100, ncol=3)
sdNs <- matrix(NA, nrow=100, ncol=19)
# convergence status
convg <- matrix(NA, nrow=100, ncol=1)

# to store estimates from open SCR version
estPsOP <- matrix(NA, nrow=100, ncol=3)
estbsOP <- matrix(NA, nrow=100, ncol=19)
estPHIsOP <- matrix(NA, nrow=100, ncol=18)
estIMMsOP <- matrix(NA, nrow=100, ncol=17)
estNsOP <- matrix(NA, nrow=100, ncol=19)
estSupsOP <- matrix(NA, nrow=100, ncol=2)
# and sds
sdPsOP <- matrix(NA, nrow=100, ncol=3)
sdbsOP <- matrix(NA, nrow=100, ncol=19)
sdPHIsOP <- matrix(NA, nrow=100, ncol=18)
sdIMMsOP <- matrix(NA, nrow=100, ncol=17)
sdNsOP <- matrix(NA, nrow=100, ncol=19)
sdSupsOP <- matrix(NA, nrow=100, ncol=2)
# convergence status
convgOP <- matrix(NA, nrow=100, ncol=1)

# from tag integrated version
estPsTI <- matrix(NA, nrow=100, ncol=3)
estbsTI <- matrix(NA, nrow=100, ncol=19)
estPHIsTI <- matrix(NA, nrow=100, ncol=18)
estIMMsTI <- matrix(NA, nrow=100, ncol=17)
estNsTI <- matrix(NA, nrow=100, ncol=19)
estSupsTI <- matrix(NA, nrow=100, ncol=2)
# and sds
sdPsTI <- matrix(NA, nrow=100, ncol=3)
sdbsTI <- matrix(NA, nrow=100, ncol=19)
sdPHIsTI <- matrix(NA, nrow=100, ncol=18)
sdIMMsTI <- matrix(NA, nrow=100, ncol=17)
sdNsTI <- matrix(NA, nrow=100, ncol=19)
sdSupsTI <- matrix(NA, nrow=100, ncol=2)
# convergence status
convgTI <- matrix(NA, nrow=100, ncol=1)

# from no temporary immigration version
estPsNO <- matrix(NA, nrow=100, ncol=3)
estbsNO <- matrix(NA, nrow=100, ncol=19)
estPHIsNO <- matrix(NA, nrow=100, ncol=18)
estNsNO <- matrix(NA, nrow=100, ncol=19)
estSupsNO <- matrix(NA, nrow=100, ncol=2)
# and sds
sdPsNO <- matrix(NA, nrow=100, ncol=3)
sdbsNO <- matrix(NA, nrow=100, ncol=19)
sdPHIsNO <- matrix(NA, nrow=100, ncol=18)
sdNsNO <- matrix(NA, nrow=100, ncol=19)
sdSupsNO <- matrix(NA, nrow=100, ncol=2)
# convergence status
convgNO <- matrix(NA, nrow=100, ncol=1)

# to store estimates from standard nonspatial version
estPsNOSp <- matrix(NA, nrow=100, ncol=2)
estbsNOSp <- matrix(NA, nrow=100, ncol=19)
estPHIsNOSp <- matrix(NA, nrow=100, ncol=18)
estIMMsNOSp <- matrix(NA, nrow=100, ncol=17)
estNsNOSp <- matrix(NA, nrow=100, ncol=19)
estSupsNOSp <- matrix(NA, nrow=100, ncol=2)
# and sds
sdPsNOSp <- matrix(NA, nrow=100, ncol=2)
sdbsNOSp <- matrix(NA, nrow=100, ncol=19)
sdPHIsNOSp <- matrix(NA, nrow=100, ncol=18)
sdIMMsNOSp <- matrix(NA, nrow=100, ncol=17)
sdNsNOSp <- matrix(NA, nrow=100, ncol=19)
sdSupsNOSp <- matrix(NA, nrow=100, ncol=2)
# convergence status
convgNOSp <- matrix(NA, nrow=100, ncol=1)

# start with a low number to see how it goes before bumping up
for (sims in 1:100){ 
  #### Simulate data based on JS model
  ## code for nonspatial JS from https://www.vogelwarte.ch/en/projects/publications/bpa/duplicate-of-complete-code-and-data-files-of-the-book
  # 10.4. Models with constant survival and time-dependent entry
  # Define parameter values
  n.occasions <- primoc                         # Number of capture occasions
  N <- 1000 #round( runif(1, min = 400, max = 4000) )      # Superpopulation size
  # Survival probabilities
  phi <- rep(0.6, n.occasions-1) #TIest[23:40] #runif( n.occasions-1, 0.1, 0.9 )
  
  # generate entry probs with a dirichlet
  #alphd <- runif( n.occasions, 10, 100)
  #b <- MCMCpack::rdirichlet( 1, alphd )    # Entry probabilities
  b <- c(0.2, rep( 0.8/17, n.occasions-1)) #3TIest[4:22]# from first model fit c(0.2, rep(0.08888889, n.occasions-1))
  imm <- rep(0.2, n.occasions-2) #TIest[41:57] #runif( n.occasions-2, 0.1, 0.5 ) #rep(runif(1, 0.05, 0.5), n.occasions-1) # Prob moves into surveyed region
  
  alpha0 <- -3.0 #TIest[1] #-2.5
  sigma <- 1.5 #TIest[2] #1.5
  alpha1 <- 1/(2*sigma*sigma)
  alpha2 <- 0.5 #TIest[3] #0.5 # runif(1, 0.5, 2) # coefficient for log effort
  
  PHI <- matrix(rep(phi, N), ncol = n.occasions-1, nrow = N, byrow = T)
  IMM <- matrix(rep(imm, N), ncol = n.occasions-2, nrow = N, byrow = T)
  
  # Function to simulate capture-recapture data under the JS model
  B <- rmultinom(1, N, b) # Generate no. of entering ind. per occasion
  # set up matrix for survival/staying and tagging 
  CH.sur <- matrix(0, ncol = n.occasions, nrow = N)
  
  # Define a vector with the occasion of entering the population
  ent.occ <- numeric()
  for (t in 1:n.occasions){
    ent.occ <- c(ent.occ, rep(t, B[t]))
  }
  
  # Simulating survival/staying in area and movement back in
  for (i in 1:N){
    CH.sur[i, ent.occ[i]] <- 1   # Write 1 when ind. enters the pop.
    # if only enters at last occasion, nothing else to do -> move to next individual
    if (ent.occ[i] == n.occasions) next
    # otherwise, use phi to determine whether stays in state space or leaves after entry occasion
    for (t in (ent.occ[i]+1):n.occasions){
      # after first occasion, stay or go
      if( t == ent.occ[i]+1 ){
        # Bernoulli trial: does individual survive/stay?
        CH.sur[i,t] <- rbinom(1, 1, PHI[i,t-1])
        # for subsequent time steps
      } else {
        # if there in previous time step, stay with prob phi
        if( CH.sur[i,t-1] == 1 ){
          CH.sur[i,t] <- rbinom(1, 1, PHI[i,t-1])
          # else come back with prob imm
        } else {
          CH.sur[i,t] <- rbinom(1, 1, IMM[i,t-2])
        }
      }
    } #t
  } #i
  
  # Simulating spatial capture histories over our trapping grid
  # our trapping grid equates to cells we spent time sampling in
  traplocs <- cbind(otraps$x, otraps$y)/10000
  
  ###########################################
  ## if you want to see if this is just related to your sampling design
  # rando <- sample(seq(1,213,1), 36, replace = F)
  # subtraps <- otraps %>% 
  #             dplyr::select( 'TrapID', 'x', 'y' ) %>% 
  #             filter( TrapID %in% c(rando) )
  #  
  # traplocs <- cbind(subtraps$x, subtraps$y)/10000
  #######################################
  
  # things related to grid cells
  G <- cbind(gridc$coords.x1, gridc$coords.x2)/10000 
  nG <- nrow(G)
  
  # calculate distance between each grid cell and each trap location
  Dmat <- e2dist(G,traplocs)
  ntraps <- nrow(traplocs)
  
  # # simulate activity centers over our state space
  # # do this by selecting grid cells IDs from a uniform distribution
  sgrid <- round( runif(N, min(hexsf$GridID), max(hexsf$GridID) ) )
  # get coordinates for each grid cell
  sgrid <- as_tibble(sgrid) %>% dplyr::rename( GridID = value )
  
  S <- left_join(sgrid,
                 hex_coords,
                 by = "GridID")
  # extract coords and scale for model fitting
  S <- S[,c(4:5)]/10000
  
  # how far is each individual's COA from each trap?
  D <- e2dist(S,traplocs)
  
  ################################################################
  #### Try just a grid with no mask
  # Here "delta" just adds
  # a fixed buffer to the outer extent of the traps.
  # buffer <- 30000
  # Xl <- min(traplocs[,1]*10000 - buffer)
  # Xu <- max(traplocs[,1]*10000 + buffer)
  # Yl <- min(traplocs[,2]*10000 - buffer)
  # Yu <- max(traplocs[,2]*10000 + buffer)
  # 
  # # simulate activity centers
  # sx <- runif(N,Xl,Xu)
  # sy <- runif(N,Yl,Yu)
  # S <- cbind(sx,sy)/10000
  # 
  # # how far is each individual's COA from each trap?
  # D <- e2dist(S,traplocs)
  
  ###########################################################################################
  ### If you want to simulate activity centers that move in each time step
  # simulate activity centers over our state space
  # have COAs for encountered individuals from step4
  # str(coas)
  # # generate for other individuals - do this by selecting grid cells IDs from a uniform distribution
  # nsgrid <- matrix( NA, nrow=N-nind, ncol=primoc)
  # nsgrid <- runif(nrow(nsgrid)*ncol(nsgrid), min(hexsf$GridID), max(hexsf$GridID) )
  # 
  # # 
  # # # merge with encountered
  #  sgrid <- rbind(coas,nsgrid)
  # # # 
  # # # # get coordinates for each grid cell and store in an array
  #  S <- array( NA, dim = c(N, 2, n.occasions) )
  #  tmp <- c(0,0)
  # # 
  # for (i in 1:N) {
  #   for (t in 1:n.occasions) {
  #     subgrid <- ifelse( !is.na(sgrid[i,t]), sgrid[i,t],
  #                        round(runif(1, min(hexsf$GridID), max(hexsf$GridID) )))
  #     tmp[1] <- as.numeric(hex_coords[subgrid,4])
  #     tmp[2] <- as.numeric(hex_coords[subgrid,5])
  # 
  #     # extract coords and scale
  #     S[i,,t] <- c(tmp[1],tmp[2])/10000
  #     }
  #   }
  # # 
  # # # How far is each individual from each trap?
  #  D <- array( NA, dim = c(N, ntraps, n.occasions) )
  # # 
  #  for (t in 1:n.occasions) {
  # #     # extract coords and scale
  #      D[,,t] <- e2dist(S[,,t],traplocs)
  #    }
  ###########################################################################################
  
  # use K from step3 
  K_sims[[sims]] <- K
  # if want to look at impact of K, use 25 trials the way they do in SCR book
  #K_sims[[sims]] <- matrix(10, nrow=19, ncol=11)
  
  #K_sims[[sims]] <- matrix( round(runif(180,0,12)), nrow=n.occasions, ncol=ntraps)
  
  # use effort data from step3
  eff_sims[[sims]] <- SEff
  #eff_sims[[sims]] <- matrix(runif(180,0,20), nrow=ntraps, ncol=n.occasions)
  
  # convert effort to log scale
  logeff <- log( eff_sims[[sims]] + 0.1 )
  #logeff[logeff == '-Inf'] <- 0
  logeff <- logeff-mean(as.vector(logeff)) #scale(logeff, center = T, scale = F)
  
  ## Chose your encounter model
  # cloglog since accounting for sampling effort
  logpois <- array( NA, dim = c(N, ntraps, n.occasions) )
  probcap <- array( NA, dim = c(N, ntraps, n.occasions) )
  for (i in 1:N){
    for (t in 1:n.occasions) {
      #probcap[i,,t] <- plogis(alpha0)*exp(-alpha1*D[i,]*D[i,])
      logpois[i,,t] <- alpha0 + alpha2*logeff[,t] - alpha1*D[i,]*D[i,]
      probcap[i,,t] <- 1-exp( -exp(logpois[i,,t]) ) 
    }
  }
  summary(probcap)
  
  # now generate the encounters of every individual in every trap in each occasion if present
  Y <- array( NA, dim = c(N, ntraps, n.occasions) )
  
  for (i in 1:N) {
    for (t in 1:n.occasions) {
      # if present, encountered with prob specified above
      if (CH.sur[i,t] == 1) {
        Y[i,,t] <- rbinom(ntraps, K_sims[[sims]][t,], probcap[i,,t])
        # if not, can't encounter so fill with zeros
      } else {
        Y[i,,t] <- rep(0, ntraps)
      }
    }
  }
  
  # Full capture-recapture array
  CH <- Y
  # Remove individuals never captured
  ## create nind x ntraps array
  cap.sum <- apply(Y,c(1,2),sum)
  ## which individuals were captured?
  ncaps <- apply(cap.sum,1,sum)
  ## keep those ones that were captured
  CH <- CH[ncaps>0,,]
  str(CH)
  str(Y)
  sum(CH)
  # Actual population size in each time step is those that were present, even if not seen
  Nt <- colSums(CH.sur)    
  
  # set up matrix for tagging for encountered individuals
  CH.tag <- matrix(0, ncol = n.occasions, nrow = nrow(CH))
  tagProb <- 0.12 # specify prob an individual is tagged given it was encountered - our rough calc was 12% based on #s tagged to untagged in each trip
  # For captured individuals, generate tagging histories
  for ( i in 1:nrow(CH) ){
    # If encountered, tagged with specified probability but only tagged once
    for (t in 1:n.occasions){
      # on first survey, tagged with specified probability if encountered
      if( t == 1 ){
        # Bernoulli trial: if encountered, was individual tagged?
        ifelse ( sum(CH[i,,t]) > 0, CH.tag[i,t] <- rbinom(1, 1, tagProb), 0 )
        # for subsequent time steps
      } else {
        # if tagged in previous time step, stay tagged
        if( CH.tag[i,t-1] == 1 ){
          CH.tag[i,t] <- 1
        } else {
          ifelse ( sum(CH[i,,t]) > 0, CH.tag[i,t] <- rbinom(1, 1, tagProb), 0 )
        }
      }
    } #t
  } #i
  
  ## Create matrix that specifies state for tagged individuals
  # first boil down CH.sur to only encountered individuals
  CH.sub <- CH.sur[ncaps>0,]
  str(CH.sub)
  str(CH)
  str(Y)
  
  # then specify which time steps we knew what the shark was up to based on tagging
  CH.known <- matrix( 2, ncol = n.occasions, nrow = nrow(CH.sub) )
  # loop through to determine if state known or not
  for( i in 1:nrow(CH.known) ){
    for( k in 1:ncol(CH.known) ){
      
      CH.known[i,k] <- ifelse( CH.tag[i,k] == 1, CH.sub[i,k], 2)
      
    }
  }
  
  
  ## Save values for comparisons later
  truePs[sims,] <- c(alpha0, sigma, alpha2)
  truebs[sims,] <- b
  truePHIs[sims,] <- phi
  trueIMMs[sims,] <- imm
  trueNs[sims,] <- Nt
  trueSups[sims,] <- c(N, N)
  
  Y_sims[[sims]] <- Y
  CH_sims[[sims]] <- CH
  Tag_sims[[sims]] <- CH.tag
  StateObs_sims[[sims]] <- CH.known
}

#########################################################################
### Fit open SCR model based on Poisson integrated likelihood
for (sims in 6:100){ 
  
  # generate grid for fitting
  # setup a grid to approximate the marginalization over activity centers
  # here we use grid the area specified above
  # we wouldn't know the true state space in real life
  # but makes things simpler for figuring out the impact of integrating tagging data
  #D <- e2dist(G, traplocs) 
  
  # number of individuals encountered in each primary period
  CH <- CH_sims[[sims]]
  eff <- eff_sims[[sims]]
  TH <- Tag_sims[[sims]]
  SH <- StateObs_sims[[sims]]
  
  # convert effort to log scale
  logeff <- log(eff + 0.1)
  #logeff[logeff == '-Inf'] <- 0
  
  # things related to grid cells
  G <- cbind(gridc$coords.x1, gridc$coords.x2)/10000 
  nG <- nrow(G)
  
  #distance from each trap location to every grid cell where an activity center could be
  D <- e2dist(G, traplocs) #flopped from way it was before
  
  # individuals encountered per primary period
  nenc_surv <- NULL 
  for (p in 1:n.occasions) {
    tpsub <- apply(CH[,,p],c(1),sum)
    nenc_surv[p] <- length( tpsub[tpsub>0] )
    
  }
  nenc_surv
  
  # matrix of Ks for each time step
  Kmat <- K_sims[[sims]]
  Kmax <- apply(Kmat,1,max)
  
  ## format encounter history
  Data <- list( "nenc" = nrow(CH), 
                "nsurv" = nenc_surv,
                "ntrap" = nrow(traplocs),  
                "ngrid" = nG,
                "ntime" = n.occasions,
                "K" = Kmat, 
                "Kmax" = Kmax,
                "log_effort" = logeff-mean(as.vector(logeff)), #scale(logeff, center=T, scale=F),
                "dist" = D,
                "y" = CH,
                "tagged" = TH,
                "obs_state" = SH, 
                "taghist" = TH, # this isn't tag hist proper but doesn't matter because we're not using
                "surveyarea" = as.numeric( sf::st_area(sa)/1000/1000 ),
                "garea_prop" = as.vector(hexsf$gridprop),
                "nstate" = 3,
                "year" = c(0,0,0,0,
                           1,1,1,1,1,
                           2,2,2,2,2,
                           3,3,3,3,3),
                "month" = c(2,3,4,5,
                            1,2,3,4,5,
                            1,2,3,4,5,
                            1,2,3,4,5)
  )
  
  Params <- list(
    "alpha0" = c(0),
    "log_alpha1" = c(0),
    "alpha2" = c(0),
    "log_lambda" = rep(0, n.occasions)
  )
  
  ###### Build and run standard SCR model
  Obj0 = MakeADFun(
    data = Data,
    parameters = Params,
    DLL = "SCR_SearchEncounter_PoissonIntLik_multisession"
  )
  
  # Prove that function and gradient calls work
  # Obj0$fn( Obj0$par )
  # Obj0$gr( Obj0$par )
  
  # Optimize
  Opt0 <- nlminb( start = Obj0$par,
                  objective = Obj0$fn,
                  gradient = Obj0$gr,
                  control = list( trace = 1, eval.max = 1e4, iter.max = 1e4),
                  verbose = TRUE
  )
  
  # Extract derived values specified in the reporting section of you .cpp file
  # rep0 <- Obj0$report()
  # rep0
  # # rep0$inidist
  # #
  # # # Convergence diagnostics - final gradient should be small
  # # Obj0$diagnostics = data.frame("name" = names(Obj0$par),
  # #                                "Est" = Opt0$par,
  # #                                "final_gradient" = as.vector(Obj0$gr(Opt0$par)))
  # # Obj0$diagnostics
  #
  # Calculate standard errors and extract parameter estimates
  sd_Obj0 <- sdreport( Obj0 ) #, bias.correct=T )
  
  pest0 <- sd_Obj0$value
  psd0 <- sd_Obj0$sd
  
  # Save estimates
  estPs[sims,] <- pest0[1:3]
  estNs[sims,] <- pest0[4:22]
  # and sds
  sdPs[sims,] <- psd0[1:3]
  sdNs[sims,] <-  psd0[4:22]
  # convergence status
  convg[sims,1] <- Opt0$convergence
  
  #########################################################################
  ### Build and run open SCR model
  Params <- list(
    "alpha0" = rep(0,1),
    "log_alpha1" = rep(0,1),
    "alpha2" = c(0),
    "log_lambda" = rep(0, primoc),
    "logit_gamma" = rep(0, primoc), 
    "logit_phi" = rep(0, primoc-1),
    "logit_imm" = rep(0, primoc-2)
  )
  
  # # Build and run standard SCR model
  Obj1 = MakeADFun(
    data = Data,
    parameters = Params,
    DLL = "SCR_SearchEncounter_PoissonIntLik_multisession_OPEN"
  )
  
  # Optimize
  Opt1 <- nlminb( start = Obj1$par,
                  objective = Obj1$fn,
                  gradient = Obj1$gr,
                  control = list( trace = 1, eval.max = 1e4, iter.max = 1e4),
                  verbose = TRUE
  )
  
  # Calculate standard errors and extract parameter estimates
  sd_Obj1 <- sdreport( Obj1 ) #, bias.correct=T )
  
  pest1 <- sd_Obj1$value
  psd1 <- sd_Obj1$sd
  
  # Save estimates
  estPsOP[sims,] <- pest1[1:3]
  estbsOP[sims,] <- pest1[4:22]
  estPHIsOP[sims,] <- pest1[23:40]
  estIMMsOP[sims,] <- pest1[41:57]
  estNsOP[sims,] <- pest1[58:76]
  estSupsOP[sims,] <- pest1[77:78]
  # and sds
  sdPsOP[sims,] <- psd1[1:3]
  sdbsOP[sims,] <- psd1[4:22]
  sdPHIsOP[sims,] <- psd1[23:40]
  sdIMMsOP[sims,] <- psd1[41:57]
  sdNsOP[sims,] <-  psd1[58:76]
  sdSupsOP[sims,] <- psd1[77:78]
  # convergence status
  convgOP[sims,1] <- Opt1$convergence
  
  
  #########################################################################
  ### Build and run tag integrated
  Params <- list(
    "alpha0" = c(0),
    "log_alpha1" = c(0),
    "alpha2" = c(0),
    "log_lambda" = rep(0, n.occasions),
    "logit_gamma" = rep(0, n.occasions),
    "logit_phi" = rep(0, n.occasions-1),
    "logit_imm" = rep(0, n.occasions-2)
  )
  
  ObjTag = MakeADFun(
    data = Data,
    parameters = Params,
    DLL = "SCR_SearchEncounter_PoissonIntLik_multisession_OPEN_TagInt"
  )
  
  # # # Prove that function and gradient calls work
  # # ObjTag$fn( ObjTag$par )
  # # ObjTag$gr( ObjTag$par )
  # #
  # # # Optimize
  OptTag = nlminb( start = ObjTag$par,
                   objective = ObjTag$fn,
                   gradient = ObjTag$gr,
                   control = list( trace = 1, eval.max = 1e4, iter.max = 1e4),
                   verbose = TRUE
  )
  # #
  # # # Calculate standard errors and extract parameter estimates
  sd_ObjTag <- sdreport( ObjTag ) #, bias.correct=T )
  
  pest2 <- sd_ObjTag$value
  psd2 <- sd_ObjTag$sd
  #
  # # Save estimates
  estPsTI[sims,] <- pest2[1:3]
  estbsTI[sims,] <- pest2[4:22]
  estPHIsTI[sims,] <- pest2[23:40]
  estIMMsTI[sims,] <- pest2[41:57]
  estNsTI[sims,] <- pest2[58:76]
  estSupsTI[sims,] <- pest2[77:78]
  # # and sds
  sdPsTI[sims,] <- psd2[1:3]
  sdbsTI[sims,] <- psd2[4:22]
  sdPHIsTI[sims,] <- psd2[23:40]
  sdIMMsTI[sims,] <- psd2[41:57]
  sdNsTI[sims,] <-  psd2[58:76]
  sdSupsTI[sims,] <- psd2[77:78]
  # # convergence status
  convgTI[sims,1] <- OptTag$convergence
  
  #########################################################################
  ### Build and run standard nonspatial
  Params <- list(
    "alpha0" = c(0),
    #  "log_alpha1" = c(0),
    #  "alpha2" = c(0),
    "log_lambda" = rep(0, n.occasions),
    "logit_gamma" = rep(0, n.occasions),
    "logit_phi" = rep(0, n.occasions-1),
    "logit_imm" = rep(0, n.occasions-2)
  )
  
  # # Build and run standard SCR model
  Obj3 = MakeADFun(
    data = Data,
    parameters = Params,
    DLL = "SCR_SearchEncounter_PoissonIntLik_multisession_OPEN_NonSpatial"
  )
  
  # Optimize
  Opt3 <- nlminb( start = Obj3$par,
                  objective = Obj3$fn,
                  gradient = Obj3$gr,
                  control = list( trace = 1, eval.max = 1e4, iter.max = 1e4),
                  verbose = TRUE
  )
  
  # Calculate standard errors and extract parameter estimates
  sd_Obj3 <- sdreport( Obj3 ) #, bias.correct=T )
  
  pest3 <- sd_Obj3$value
  psd3 <- sd_Obj3$sd
  
  # Save estimates
  estPsNOSp[sims,] <- pest3[1]
  estbsNOSp[sims,] <- pest3[2:20]
  estPHIsNOSp[sims,] <- pest3[21:38]
  estIMMsNOSp[sims,] <- pest3[39:55]
  estNsNOSp[sims,] <- pest3[56:74]
  estSupsNOSp[sims,] <- pest3[75:76]
  # and sds
  sdPsNOSp[sims,] <- psd3[1]
  sdbsNOSp[sims,] <- psd3[2:20]
  sdPHIsNOSp[sims,] <- psd3[21:38]
  sdIMMsNOSp[sims,] <- psd3[39:55]
  sdNsNOSp[sims,] <-  psd3[56:74]
  sdSupsNOSp[sims,] <- psd3[75:76]
  # convergence status
  convgNOSp[sims,1] <- Opt3$convergence
  
  ########################################################################
  ### Build and run model without temporary emigration
  Params <- list(
    "alpha0" = c(0),
    "log_alpha1" = c(0),
    "alpha2" = c(0),
    "log_lambda" = rep(0, n.occasions),
    "logit_gamma" = rep(0, n.occasions),
    "logit_phi" = rep(0, n.occasions-1)
  )
  
  # Build and run standard model
  Obj4 = MakeADFun(
    data = Data,
    parameters = Params,
    DLL = "SCR_SearchEncounter_PoissonIntLik_multisession_OPEN_NoTE"
  )
  
  # Optimize
  Opt4 = nlminb( start = Obj4$par,
                 objective = Obj4$fn,
                 gradient = Obj4$gr,
                 control = list( trace = 1, eval.max = 1e4, iter.max = 1e4),
                 verbose = TRUE
  )
  #
  # # Calculate standard errors and extract parameter estimates
  sd_Obj4 <- sdreport( Obj4 ) #, bias.correct=T )
  
  pest4 <- sd_Obj4$value
  psd4 <- sd_Obj4$sd
  
  # # Save estimates
  estPsNO[sims,] <- pest4[1:3]
  estbsNO[sims,] <- pest4[4:22]
  estPHIsNO[sims,] <- pest4[23:40]
  estNsNO[sims,] <- pest4[41:59]
  estSupsNO[sims,] <- pest4[60:61]
  # # and sds
  sdPsNO[sims,] <- psd4[1:3]
  sdbsNO[sims,] <- psd4[4:22]
  sdPHIsNO[sims,] <- psd4[23:40]
  sdNsNO[sims,] <-  psd4[41:59]
  sdSupsNO[sims,] <- psd4[60:61]
  # # convergence status
  convgNO[sims,1] <- Opt4$convergence
  # #
  
  # rinse and repeat
}

# save workspace to keep these handy
#save.image(file='SimRun_PoisIntLik_TagInt_25kmbuff_100km2grid_IamanIdiot.RData')
load('data/SimRun_PoisIntLik_TagInt_25kmbuff_100km2grid.RData')

# how did it do??
#mapeP <- ( (estPs - truePs)/truePs )*100

# create dataframe for boxplots
diffPs <- tibble(
  alpha0 = estPs[,1] - truePs[,1],
  sigma = estPs[,2] - truePs[,2],
  alpha2 = estPs[,3] - truePs[,3],
  Model = '5' #SCR'
)

diffPsNoSp <- tibble(
  alpha0 = estPsNOSp[,1] - truePs[,1],
  sigma = rep(NA, 100),
  alpha2 = rep(NA, 100),
  Model = '4' #NonSCR'
)

diffPNO <- tibble(
  alpha0 = estPsNO[,1] - truePs[,1],
  sigma = estPsNO[,2] - truePs[,2],
  alpha2 = estPsNO[,3] - truePs[,3],
  Model = '3' #OPSCR_noTE'
)

diffPOP <- tibble(
  alpha0 = estPsOP[,1] - truePs[,1],
  sigma = estPsOP[,2] - truePs[,2],
  alpha2 = estPsOP[,3] - truePs[,3],
  Model = '2' #OPSCR'
)

diffPTI <- tibble(
  alpha0 = estPsTI[,1] - truePs[,1],
  sigma = estPsTI[,2] - truePs[,2],
  alpha2 = estPsTI[,3] - truePs[,3],
  Model = '1' #TI_OPSCR'
)

################################### If using MPE instead
#mapeP <- ( (estPs - truePs)/truePs )*100
# 
# diffPs <- tibble(
#   alpha0 = ( (estPs[,1] - truePs[,1])/truePs[,1] )*100,
#   sigma = ( (estPs[,2] - truePs[,2])/truePs[,2] )*100,
#   alpha2 = ( (estPs[,3] - truePs[,3])/truePs[,3] )*100,
#   model = '5' #SCR'
# )
# 
# diffPsNoSp <- tibble(
#   alpha0 = ( (estPsNOSp[,1] - truePs[,1])/truePs[,1] )*100,
#   sigma = rep(NA, 100),
#   alpha2 = rep(NA, 100),
#   model = '4' #NonSCR'
# )
# 
# diffPNO <- tibble(
#   alpha0 = ( (estPsNO[,1] - truePs[,1])/truePs[,1] )*100,
#   sigma = ( (estPsNO[,2] - truePs[,2])/truePs[,2] )*100,
#   alpha2 = ( (estPsNO[,3] - truePs[,3])/truePs[,3] )*100,
#   model = '3' #OPSCR_noTE'
# )
# 
# diffPOP <- tibble(
#   alpha0 = ( (estPsOP[,1] - truePs[,1])/truePs[,1] )*100,
#   sigma = ( (estPsOP[,2] - truePs[,2])/truePs[,2] )*100,
#   alpha2 = ( (estPsOP[,3] - truePs[,3])/truePs[,3] )*100,
#   model = '2' #OPSCR'
# )
# 
# diffPTI <- tibble(
#   alpha0 = ( (estPsTI[,1] - truePs[,1])/truePs[,1] )*100,
#   sigma = ( (estPsTI[,2] - truePs[,2])/truePs[,2] )*100,
#   alpha2 = ( (estPsTI[,3] - truePs[,3])/truePs[,3] )*100,
#   model = '1' #'TI_OPSCR'
# )


Pdiffs <- rbind(diffPs, diffPsNoSp, diffPNO, diffPOP, diffPTI)

# Boxplots for each
# a0 <- ggplot(Pdiffs, aes(x=model, y=alpha0)) + 
#   geom_boxplot() +
#   ylim(-20, 20) +
#   geom_hline( yintercept=0 )
# sig <- ggplot(Pdiffs, aes(x=model, y=sigma)) + 
#   geom_boxplot() +
#  ylim(-20, 20) +
#   geom_hline( yintercept=0 )
# a2 <- ggplot(Pdiffs, aes(x=model, y=alpha2)) + 
#   geom_boxplot() +
#   ylim(-20, 20)  +
#   geom_hline( yintercept=0 )
# 
# ggpubr::ggarrange(a0, sig, a2, 
#           labels = c("A", "B", "C"),
#           ncol = 1, 
#           nrow = 3)

## To do a grouped boxplot, make params a field
LPdiffs <- Pdiffs %>% 
  pivot_longer(cols=c('alpha0', 'sigma','alpha2'),
               names_to = 'param',
               values_to = 'value') %>% 
  mutate( ParamN = ifelse(param == 'alpha0','beta[0]', 
                          ifelse(param == 'sigma', 'sigma', 'beta[1]')) )

LPdiffs$Param <- factor( LPdiffs$ParamN,
                         labels =  c(bquote( "(a)" ~ beta[0] ),
                                     bquote( "(b)" ~ sigma ),
                                     bquote( "(c)" ~ beta[1] )) )
LPdiffs

LPdiffs <- LPdiffs %>% 
  filter( value < 100)

# grouped boxplot
ggplot(LPdiffs, aes(x=Param, y=value, fill=Model)) + 
  geom_boxplot() +
  #  ylim(-100, 100) +
  theme( axis.text.x = element_blank(),
         axis.title.x = element_blank() ) +
  geom_hline( yintercept=0 ) + 
  facet_wrap(~Param, scales = 'free', labeller = label_parsed) +
  theme(strip.text = element_text(hjust = 0)) +
  ylab("Difference") 

############################# bs - won't have them for all
## first time step
d0bsNoSp <- estbsNOSp[,1] - truebs[,1]
d0bsNoSp <- as_tibble(d0bsNoSp)
names(d0bsNoSp) <- 1           
# convert into long format
d0bNoSp <- gather(d0bsNoSp)
d0bNoSp$Model <- '4'
d0bNoSp$param <- 'b0'

d0bsNo <- estbsNO[,1] - truebs[,1]
d0bsNo <- as_tibble(d0bsNo)
names(d0bsNo) <- 1          
# convert into long format
d0bNo <- gather(d0bsNo)
d0bNo$Model <- '3'
d0bNo$param <- 'b0'

d0bsOP <- estbsOP[,1] - truebs[,1]
d0bsOP <- as_tibble(d0bsOP)
names(d0bsOP) <- 1          
# convert into long format
d0bOP <- gather(d0bsOP)
d0bOP$Model <- '2'
d0bOP$param <- 'b0'

d0bsTI <- estbsTI[,1] - truebs[,1]
d0bsTI <- as_tibble(d0bsTI)
names(d0bsTI) <- 1           
# convert into long format
d0bTI <- gather(d0bsTI)
d0bTI$Model <- '1'
d0bTI$param <- 'b0'

bd0s <- rbind(d0bNoSp, d0bNo, d0bOP, d0bTI)
bd0s <- bd0s %>% 
  dplyr::select( 'Model', 'param', 'value') %>% 
  mutate( ParamN = 'alpha[0]' )

bd0s$Param <- factor( bd0s$ParamN,
                      labels =  c(bquote( "(d)" ~ alpha[0] ) ) )
bd0s

## To do a grouped boxplot, make params a field
Pcs <- rbind(LPdiffs, bd0s)

# grouped boxplot
ggplot(Pcs, aes(x=Param, y=value, fill=Model)) + 
  geom_boxplot() +
  # ylim(-100, 100) +
  theme( axis.text.x = element_blank(),
         axis.title.x = element_blank() ) +
  geom_hline( yintercept=0 ) + 
  facet_wrap(~Param, scales = 'free', 
             labeller = label_parsed ) +
  theme(strip.text = element_text(hjust = 0)) +
  ylab("Difference") 

### all others - same so lump
diffbsNoSp <- estbsNOSp[,2:19] - truebs[,2:19]
diffbsNoSp <- as_tibble(diffbsNoSp)
names(diffbsNoSp) <- seq(2,19,1)           
# convert into long format
dbNoSp <- gather(diffbsNoSp)
dbNoSp$Model <- '4'
dbNoSp$param <- 'b'

diffbsNo <- estbsNO[,2:19] - truebs[,2:19]
diffbsNo <- as_tibble(diffbsNo)
names(diffbsNo) <- seq(2,19,1)           
# convert into long format
dbNo <- gather(diffbsNo)
dbNo$Model <- '3'
dbNo$param <- 'b'

diffbsOP <- estbsOP[,2:19] - truebs[,2:19]
diffbsOP <- as_tibble(diffbsOP)
names(diffbsOP) <- seq(2,19,1)           
# convert into long format
dbOP <- gather(diffbsOP)
dbOP$Model <- '2'
dbOP$param <- 'b'

diffbsTI <- estbsTI[,2:19] - truebs[,2:19]
diffbsTI <- as_tibble(diffbsTI)
names(diffbsTI) <- seq(2,19,1)           
# convert into long format
dbTI <- gather(diffbsTI)
dbTI$Model <- '1'
dbTI$param <- 'b'

bdiffs <- rbind(dbNoSp, dbNo, dbOP, dbTI)
bdiffs <- bdiffs %>% 
  dplyr::select( 'Model', 'param', 'value') %>% 
  mutate( ParamN = 'alpha[t]' )

bdiffs$Param <- factor( bdiffs$ParamN,
                        labels =  c(bquote( "(e)" ~ alpha[t] ) ) )
bdiffs
## To do a grouped boxplot, make params a field
Pcs <- rbind(Pcs, bdiffs)

# grouped boxplot
ggplot(Pcs, aes(x=Param, y=value, fill=Model)) + 
  geom_boxplot() +
  # ylim(-100, 100) +
  theme( axis.text.x = element_blank(),
         axis.title.x = element_blank() ) +
  geom_hline( yintercept=0 ) + 
  facet_wrap(~Param, scales = 'free', labeller = label_parsed) +
  theme(strip.text = element_text(hjust = 0)) +
  ylab("Difference") 


# # Boxplots for each
# bNoSP <- ggplot(dbNoSp, aes(x=key, y=value) ) + 
#   geom_boxplot() +
#   ylim(-0.15, 0.15) +
#   geom_hline( yintercept=0 )
# bNo <- ggplot(dbNo, aes(x=key, y=value) ) + 
#   geom_boxplot() +
#   ylim(-0.15, 0.15) +
#   geom_hline( yintercept=0 )
# bOP <- ggplot(dbOP, aes(x=key, y=value) ) + 
#   geom_boxplot() +
#   ylim(-0.15, 0.15) +
#   geom_hline( yintercept=0 )
# bTI <- ggplot(dbTI, aes(x=key, y=value) ) + 
#   geom_boxplot() +
#   ylim(-0.15, 0.15) +
#   geom_hline( yintercept=0 )
# 
# ggpubr::ggarrange(bNoSP, bNo, bOP, bTI,
#                   labels = c("A", "B", "C", "D"),
#                   ncol = 1, 
#                   nrow = 4)
# 
# 
# ### aggregated
# bNoSP <- ggplot(dbNoSp, aes(y=value) ) + 
#   geom_boxplot() +
#   ylim(-0.15, 0.15) +
#   geom_hline( yintercept=0 )
# bNo <- ggplot(dbNo, aes(y=value) ) + 
#   geom_boxplot() +
#   ylim(-0.15, 0.15) +
#   geom_hline( yintercept=0 )
# bOP <- ggplot(dbOP, aes(y=value) ) + 
#   geom_boxplot() +
#   ylim(-0.15, 0.15) +
#   geom_hline( yintercept=0 )
# bTI <- ggplot(dbTI, aes(y=value) ) + 
#   geom_boxplot() +
#   ylim(-0.15, 0.15) +
#   geom_hline( yintercept=0 )
# 
# ggpubr::ggarrange(bNoSP, bNo, bOP, bTI,
#                   labels = c("A", "B", "C", "D"),
#                   ncol = 1, 
#                   nrow = 4)


############################# phis - won't have them for all
dPHIsNoSp <- estPHIsNOSp - truePHIs
dPHIsNoSp <- as_tibble(dPHIsNoSp)
names(dPHIsNoSp) <- seq(1,18,1)           
# convert into long format
dPHIsNoSp <- gather(dPHIsNoSp)
dPHIsNoSp$Model <- '4'
dPHIsNoSp$param <- 'phi'

dPHIsNo <- estPHIsNO - truePHIs
dPHIsNo <- as_tibble(dPHIsNo)
names(dPHIsNo) <- seq(1,18,1)           
# convert into long format
dPHIsNo <- gather(dPHIsNo)
dPHIsNo$Model <- '3'
dPHIsNo$param <- 'phi'

dPHIsOP <- estPHIsOP - truePHIs
dPHIsOP <- as_tibble(dPHIsOP)
names(dPHIsOP) <- seq(1,18,1)           
# convert into long format
dPHIsOP <- gather(dPHIsOP)
dPHIsOP$Model <- '2'
dPHIsOP$param <- 'phi'

dPHIsTI <- estPHIsTI - truePHIs
dPHIsTI <- as_tibble(dPHIsTI)
names(dPHIsTI) <- seq(1,18,1)           
# convert into long format
dPHIsTI <- gather(dPHIsTI)
dPHIsTI$Model <- '1'
dPHIsTI$param <- 'phi'

phidiffs <- rbind(dPHIsNoSp, dPHIsNo, dPHIsOP, dPHIsTI)
phidiffs <- phidiffs %>% 
  dplyr::select( 'Model', 'param', 'value') %>% 
  mutate( ParamN = 'phi[t]' )


phidiffs$Param <- factor( phidiffs$ParamN,
                          labels =  c(bquote( "(f)" ~ phi[t] ) ) )

## To do a grouped boxplot, make params a field
Pcs <- rbind(Pcs, phidiffs)

phidiffs
# grouped boxplot
ggplot(Pcs, aes(x=Param, y=value, fill=Model)) + 
  geom_boxplot() +
  # ylim(-100, 100) +
  theme( axis.text.x = element_blank(),
         axis.title.x = element_blank() ) +
  geom_hline( yintercept=0 ) + 
  facet_wrap(~Param, scales = 'free', labeller = label_parsed) +  
  theme(strip.text = element_text(hjust = 0)) +
  ylab("Difference") 
# 
# # Boxplots for each
# bNoSP <- ggplot(dPHIsNoSp, aes(x=key, y=value) ) + 
#   geom_boxplot() +
#  # ylim(-0.15, 0.15) +
#   geom_hline( yintercept=0 )
# bNo <- ggplot(dPHIsNo, aes(x=key, y=value) ) + 
#   geom_boxplot() +
# #  ylim(-0.15, 0.15) +
#   geom_hline( yintercept=0 )
# bOP <- ggplot(dPHIsOP, aes(x=key, y=value) ) + 
#   geom_boxplot() +
# #  ylim(-0.15, 0.15) +
#   geom_hline( yintercept=0 )
# bTI <- ggplot(dPHIsTI, aes(x=key, y=value) ) + 
#   geom_boxplot() +
# #  ylim(-0.15, 0.15) +
#   geom_hline( yintercept=0 )
# 
# ggpubr::ggarrange(bNoSP, bNo, bOP, bTI,
#                   labels = c("A", "B", "C", "D"),
#                   ncol = 1, 
#                   nrow = 4)

############################# imms - won't have them for all
dIMMsNoSp <- estIMMsNOSp - trueIMMs
dIMMsNoSp <- as_tibble(dIMMsNoSp)
names(dIMMsNoSp) <- seq(1,17,1)           
# convert into long format
dIMMsNoSp <- gather(dIMMsNoSp)
dIMMsNoSp$Model <- '4'
dIMMsNoSp$param <- 'imm'

dIMMsOP <- estIMMsOP - trueIMMs
dIMMsOP <- as_tibble(dIMMsOP)
names(dIMMsOP) <- seq(1,17,1)           
# convert into long format
dIMMsOP <- gather(dIMMsOP)
dIMMsOP$Model <- '2'
dIMMsOP$param <- 'imm'

dIMMsTI <- estIMMsTI - trueIMMs
dIMMsTI <- as_tibble(dIMMsTI)
names(dIMMsTI) <- seq(1,17,1)           
# convert into long format
dIMMsTI <- gather(dIMMsTI)
dIMMsTI$Model <- '1'
dIMMsTI$param <- 'imm'

immdiffs <- rbind(dIMMsNoSp, dIMMsOP, dIMMsTI)
immdiffs <- immdiffs %>% 
  dplyr::select( 'Model', 'param', 'value') %>%  
  mutate( ParamN = 'eta[t]' )

immdiffs$Param <- factor( immdiffs$ParamN,
                          labels =  c(bquote( "(g)" ~ eta[t] ) ) )


## To do a grouped boxplot, make params a field
Pcs <- rbind(Pcs, immdiffs)

# grouped boxplot
ggplot(Pcs, aes(x=Param, y=value, fill=Model)) + 
  geom_boxplot() +
  # ylim(-100, 100) +
  theme( axis.text.x = element_blank(),
         axis.title.x = element_blank() ) +
  geom_hline( yintercept=0 ) + 
  facet_wrap(~Param, scales = 'free', labeller = label_parsed) +
  ylab("Difference") 

# Boxplots for each
# dNoSP <- ggplot(dIMMsNoSp, aes(x=key, y=value) ) + 
#   geom_boxplot() +
#   ylim(-0.15, 0.15) +
#   geom_hline( yintercept=0 )
# dOP <- ggplot(dIMMsOP, aes(x=key, y=value) ) + 
#   geom_boxplot() +
#   ylim(-0.15, 0.15) +
#   geom_hline( yintercept=0 )
# dTI <- ggplot(dIMMsTI, aes(x=key, y=value) ) + 
#   geom_boxplot() +
#   ylim(-0.15, 0.15) +
#   geom_hline( yintercept=0 )
# 
# ggpubr::ggarrange(dNoSP, dOP, dTI,
#                   labels = c("A", "B", "C"),
#                   ncol = 1, 
#                   nrow = 3)

############################################# Ns
# create dataframe for boxplots
dNs <- estNs - trueNs
dNs <- as_tibble(dNs)
names(dNs) <- seq(1,19,1)           
# convert into long format
dNs <- gather(dNs)
dNs$Model <- '5'
dNs$param <- 'N'

dNsNoSp <- estNsNOSp - trueNs
dNsNoSp <- as_tibble(dNsNoSp)
names(dNsNoSp) <- seq(1,19,1)           
# convert into long format
dNsNoSp <- gather(dNsNoSp)
dNsNoSp$Model <- '4'
dNsNoSp$param <- 'N'

dNsNo <- estNsNO - trueNs
dNsNo <- as_tibble(dNsNo)
names(dNsNo) <- seq(1,19,1)           
# convert into long format
dNsNo <- gather(dNsNo)
dNsNo$Model <- '3'
dNsNo$param <- 'N'

dNsOP <- estNsOP - trueNs
dNsOP <- as_tibble(dNsOP)
names(dNsOP) <- seq(1,19,1)           
# convert into long format
dNsOP <- gather(dNsOP)
dNsOP$Model <- '2'
dNsOP$param <- 'N'

dNsTI <- estNsTI - trueNs
dNsTI <- as_tibble(dNsTI)
names(dNsTI) <- seq(1,19,1)           
# convert into long format
dNsTI <- gather(dNsTI)
dNsTI$Model <- '1'
dNsTI$param <- 'N'

Ndiffs <- rbind(dNs, dNsNoSp, dNsNo, dNsOP, dNsTI)
Ndiffs <- Ndiffs %>% 
  dplyr::select( 'Model', 'param', 'value') %>% 
  mutate( ParamN = 'N[t]' )

Ndiffs$Param <- factor( Ndiffs$ParamN,
                        labels =  c(bquote( "(h)" ~ N[t] ) ) )

## To do a grouped boxplot, make params a field
Pcs <- rbind(Pcs, Ndiffs)

# grouped boxplot
ggplot(Pcs, aes(x=Param, y=value, fill=Model)) + 
  geom_boxplot() +
  # ylim(-100, 100) +
  theme( axis.text.x = element_blank(),
         axis.title.x = element_blank() ) +
  geom_hline( yintercept=0 ) + 
  facet_wrap(~Param, scales = 'free', labeller = label_parsed) +
  ylab("Difference") 

# ### boxplots
# dN <- ggplot(dNs, aes(x=key,y=value) ) + 
#   geom_boxplot() +
#  # ylim(-400, 600) +
#   geom_hline( yintercept=0 )
# dNoSP <- ggplot(dNsNoSp, aes(x=key,y=value) ) + 
#   geom_boxplot() +
#  # ylim(-400, 600) +
#   geom_hline( yintercept=0 )
# dNo <- ggplot(dNsNo, aes(x=key,y=value) ) + 
#   geom_boxplot() +
# #  ylim(-400, 600) +
#   geom_hline( yintercept=0 )
# dOP <- ggplot(dNsOP, aes(x=key,y=value) ) + 
#   geom_boxplot() +
# #  ylim(-400, 600) +
#   geom_hline( yintercept=0 )
# dTI <- ggplot(dNsTI, aes(x=key,y=value) ) + 
#   geom_boxplot() +
#  # ylim(-400, 600) +
#   geom_hline( yintercept=0 )
# 
# ggpubr::ggarrange(dN, dNoSP, dNo, dOP, dTI,
#                   labels = c("A", "B", "C", "D", "E"),
#                   ncol = 1, 
#                   nrow = 5)

# create dataframe for superpop size
diffSupsNoSp <- tibble(
  # Nsup = estSupsNOSp[,1] - trueSups[,1],
  Nold = estSupsNOSp[,2] - trueSups[,2],
  Model = '4' #NonSCR'
)

diffSupsNo <- tibble(
  # Nsup = estSupsNO[,1] - trueSups[,1],
  Nold = estSupsNO[,2] - trueSups[,2],
  Model = '3' #OPSCR_noTE'
)

diffSupsOP <- tibble(
  # Nsup = estSupsOP[,1] - trueSups[,1],
  Nold = estSupsOP[,2] - trueSups[,2],
  Model = '2' #OPSCR'
)

diffSupsTI <- tibble(
  # Nsup = estSupsTI[,1] - trueSups[,1],
  Nold = estSupsTI[,2] - trueSups[,2],
  Model = '1' #TI_OPSCR'
)

NSUPdiffs <- rbind(diffSupsNoSp, diffSupsNo, diffSupsOP, diffSupsTI)
## To do a grouped boxplot, make params a field
NSUPdiffs <- NSUPdiffs %>% 
  pivot_longer(cols=c( 'Nold'),
               names_to = 'param',
               values_to = 'value') %>% 
  mutate( ParamN = 'N[super]' )


NSUPdiffs$Param <- factor( NSUPdiffs$ParamN,
                           labels =  c(bquote( "(i)" ~ N[super] ) ) )

## To do a grouped boxplot, make params a field
Pcs <- rbind(Pcs, NSUPdiffs)

# grouped boxplot
ggplot(Pcs, aes(x=Param, y=value, fill=Model)) + 
  geom_boxplot() +
  # ylim(-100, 100) +
  theme( axis.text.x = element_blank(),
         axis.title.x = element_blank() ) +
  geom_hline( yintercept=0 ) + 
  facet_wrap(~Param, scales = 'free', labeller = label_parsed) +
  ylab("Difference") 


## Save after all this, lordamighty!
ggsave("Figure_5.pdf", 
       device = 'pdf', width = 7, height = 7, units = "in")


#####################################################################
#### MAPE version
#####################################################################
################################### If using MPE instead
#mapeP <- ( (estPs - truePs)/truePs )*100
# 
diffPs <- tibble(
  alpha0 = ( (estPs[,1] - truePs[,1])/truePs[,1] )*100,
  sigma = ( (estPs[,2] - truePs[,2])/truePs[,2] )*100,
  alpha2 = ( (estPs[,3] - truePs[,3])/truePs[,3] )*100,
  Model = '5' #SCR'
)

diffPsNoSp <- tibble(
  alpha0 = ( (estPsNOSp[,1] - truePs[,1])/truePs[,1] )*100,
  sigma = rep(NA, 100),
  alpha2 = rep(NA, 100),
  Model = '4' #NonSCR'
)

diffPNO <- tibble(
  alpha0 = ( (estPsNO[,1] - truePs[,1])/truePs[,1] )*100,
  sigma = ( (estPsNO[,2] - truePs[,2])/truePs[,2] )*100,
  alpha2 = ( (estPsNO[,3] - truePs[,3])/truePs[,3] )*100,
  Model = '3' #OPSCR_noTE'
)

diffPOP <- tibble(
  alpha0 = ( (estPsOP[,1] - truePs[,1])/truePs[,1] )*100,
  sigma = ( (estPsOP[,2] - truePs[,2])/truePs[,2] )*100,
  alpha2 = ( (estPsOP[,3] - truePs[,3])/truePs[,3] )*100,
  Model = '2' #OPSCR'
)

diffPTI <- tibble(
  alpha0 = ( (estPsTI[,1] - truePs[,1])/truePs[,1] )*100,
  sigma = ( (estPsTI[,2] - truePs[,2])/truePs[,2] )*100,
  alpha2 = ( (estPsTI[,3] - truePs[,3])/truePs[,3] )*100,
  Model = '1' #'TI_OPSCR'
)

Pdiffs <- rbind(diffPs, diffPsNoSp, diffPNO, diffPOP, diffPTI)

## To do a grouped boxplot, make params a field
LPdiffs <- Pdiffs %>% 
  pivot_longer(cols=c('alpha0', 'sigma','alpha2'),
               names_to = 'param',
               values_to = 'value') %>% 
  mutate( Param = ifelse(param == 'alpha0','beta[0]', 
                         ifelse(param == 'sigma', 'sigma', 'beta[1]')) )

LPdiffs <- LPdiffs %>% 
  filter( value < 25)

# grouped boxplot
ggplot(LPdiffs, aes(x=Param, y=value, fill=Model)) + 
  geom_boxplot() +
  #  ylim(-100, 100) +
  theme( axis.text.x = element_blank(),
         axis.title.x = element_blank() ) +
  geom_hline( yintercept=0 ) + 
  facet_wrap(~Param, scales = 'free', labeller = label_parsed) +
  ylab("Difference") 

############################# bs - won't have them for all
( (estPsTI[,1] - truePs[,1])/truePs[,1] )*100

## first time step
d0bsNoSp <- ( (estbsNOSp[,1] - truebs[,1])/truebs[,1] )*100
d0bsNoSp <- as_tibble(d0bsNoSp)
names(d0bsNoSp) <- 1           
# convert into long format
d0bNoSp <- gather(d0bsNoSp)
d0bNoSp$Model <- '4'
d0bNoSp$param <- 'b0'

d0bsNo <- ( (estbsNO[,1] - truebs[,1])/truebs[,1] )*100
d0bsNo <- as_tibble(d0bsNo)
names(d0bsNo) <- 1          
# convert into long format
d0bNo <- gather(d0bsNo)
d0bNo$Model <- '3'
d0bNo$param <- 'b0'

d0bsOP <- ( (estbsOP[,1] - truebs[,1])/truebs[,1] )*100
d0bsOP <- as_tibble(d0bsOP)
names(d0bsOP) <- 1          
# convert into long format
d0bOP <- gather(d0bsOP)
d0bOP$Model <- '2'
d0bOP$param <- 'b0'

d0bsTI <- ( (estbsTI[,1] - truebs[,1])/truebs[,1] )*100
d0bsTI <- as_tibble(d0bsTI)
names(d0bsTI) <- 1           
# convert into long format
d0bTI <- gather(d0bsTI)
d0bTI$Model <- '1'
d0bTI$param <- 'b0'

bd0s <- rbind(d0bNoSp, d0bNo, d0bOP, d0bTI)
bd0s <- bd0s %>% 
  dplyr::select( 'Model', 'param', 'value') %>% 
  mutate( Param = 'alpha[0]' )

## To do a grouped boxplot, make params a field
Pcs <- rbind(LPdiffs, bd0s)

# grouped boxplot
ggplot(Pcs, aes(x=Param, y=value, fill=Model)) + 
  geom_boxplot() +
  # ylim(-100, 100) +
  theme( axis.text.x = element_blank(),
         axis.title.x = element_blank() ) +
  geom_hline( yintercept=0 ) + 
  facet_wrap(~Param, scales = 'free', 
             labeller = label_parsed ) +
  ylab("Difference") 

### all others - same so lump
diffbsNoSp <- ( (estbsNOSp[,2:19] - truebs[,2:19])/truebs[,2:19] )*100
diffbsNoSp <- as_tibble(diffbsNoSp)
names(diffbsNoSp) <- seq(2,19,1)           
# convert into long format
dbNoSp <- gather(diffbsNoSp)
dbNoSp$Model <- '4'
dbNoSp$param <- 'b'

diffbsNo <- ( (estbsNO[,2:19] - truebs[,2:19])/truebs[,2:19] )*100
diffbsNo <- as_tibble(diffbsNo)
names(diffbsNo) <- seq(2,19,1)           
# convert into long format
dbNo <- gather(diffbsNo)
dbNo$Model <- '3'
dbNo$param <- 'b'

diffbsOP <- ( (estbsOP[,2:19] - truebs[,2:19])/truebs[,2:19] )*100
diffbsOP <- as_tibble(diffbsOP)
names(diffbsOP) <- seq(2,19,1)           
# convert into long format
dbOP <- gather(diffbsOP)
dbOP$Model <- '2'
dbOP$param <- 'b'

diffbsTI <- ( (estbsTI[,2:19] - truebs[,2:19])/truebs[,2:19] )*100
diffbsTI <- as_tibble(diffbsTI)
names(diffbsTI) <- seq(2,19,1)           
# convert into long format
dbTI <- gather(diffbsTI)
dbTI$Model <- '1'
dbTI$param <- 'b'

bdiffs <- rbind(dbNoSp, dbNo, dbOP, dbTI)
bdiffs <- bdiffs %>% 
  dplyr::select( 'Model', 'param', 'value') %>% 
  mutate( Param = 'alpha[t]' )

bdiffs <- bdiffs %>% 
  filter(value < 400)

## To do a grouped boxplot, make params a field
Pcs <- rbind(Pcs, bdiffs)

# grouped boxplot
ggplot(Pcs, aes(x=Param, y=value, fill=Model)) + 
  geom_boxplot() +
  # ylim(-100, 100) +
  theme( axis.text.x = element_blank(),
         axis.title.x = element_blank() ) +
  geom_hline( yintercept=0 ) + 
  facet_wrap(~Param, scales = 'free', labeller = label_parsed) +
  ylab("Difference") 


############################# phis - won't have them for all
dPHIsNoSp <- ( (estPHIsNOSp - truePHIs)/truePHIs)*100
dPHIsNoSp <- as_tibble(dPHIsNoSp)
names(dPHIsNoSp) <- seq(1,18,1)           
# convert into long format
dPHIsNoSp <- gather(dPHIsNoSp)
dPHIsNoSp$Model <- '4'
dPHIsNoSp$param <- 'phi'

dPHIsNo <- ( (estPHIsNO - truePHIs)/truePHIs)*100
dPHIsNo <- as_tibble(dPHIsNo)
names(dPHIsNo) <- seq(1,18,1)           
# convert into long format
dPHIsNo <- gather(dPHIsNo)
dPHIsNo$Model <- '3'
dPHIsNo$param <- 'phi'

dPHIsOP <- ( (estPHIsOP - truePHIs)/truePHIs)*100
dPHIsOP <- as_tibble(dPHIsOP)
names(dPHIsOP) <- seq(1,18,1)           
# convert into long format
dPHIsOP <- gather(dPHIsOP)
dPHIsOP$Model <- '2'
dPHIsOP$param <- 'phi'

dPHIsTI <- ( (estPHIsTI - truePHIs)/truePHIs)*100
dPHIsTI <- as_tibble(dPHIsTI)
names(dPHIsTI) <- seq(1,18,1)           
# convert into long format
dPHIsTI <- gather(dPHIsTI)
dPHIsTI$Model <- '1'
dPHIsTI$param <- 'phi'

phidiffs <- rbind(dPHIsNoSp, dPHIsNo, dPHIsOP, dPHIsTI)
phidiffs <- phidiffs %>% 
  dplyr::select( 'Model', 'param', 'value') %>% 
  mutate( Param = 'phi[t]' )

## To do a grouped boxplot, make params a field
Pcs <- rbind(Pcs, phidiffs)

# grouped boxplot
ggplot(Pcs, aes(x=Param, y=value, fill=Model)) + 
  geom_boxplot() +
  # ylim(-100, 100) +
  theme( axis.text.x = element_blank(),
         axis.title.x = element_blank() ) +
  geom_hline( yintercept=0 ) + 
  facet_wrap(~Param, scales = 'free', labeller = label_parsed) +
  ylab("Difference") 
# 
# # Boxplots for each
# bNoSP <- ggplot(dPHIsNoSp, aes(x=key, y=value) ) + 
#   geom_boxplot() +
#  # ylim(-0.15, 0.15) +
#   geom_hline( yintercept=0 )
# bNo <- ggplot(dPHIsNo, aes(x=key, y=value) ) + 
#   geom_boxplot() +
# #  ylim(-0.15, 0.15) +
#   geom_hline( yintercept=0 )
# bOP <- ggplot(dPHIsOP, aes(x=key, y=value) ) + 
#   geom_boxplot() +
# #  ylim(-0.15, 0.15) +
#   geom_hline( yintercept=0 )
# bTI <- ggplot(dPHIsTI, aes(x=key, y=value) ) + 
#   geom_boxplot() +
# #  ylim(-0.15, 0.15) +
#   geom_hline( yintercept=0 )
# 
# ggpubr::ggarrange(bNoSP, bNo, bOP, bTI,
#                   labels = c("A", "B", "C", "D"),
#                   ncol = 1, 
#                   nrow = 4)

############################# imms - won't have them for all
dIMMsNoSp <- ( (estIMMsNOSp - trueIMMs)/trueIMMs )*100
dIMMsNoSp <- as_tibble(dIMMsNoSp)
names(dIMMsNoSp) <- seq(1,17,1)           
# convert into long format
dIMMsNoSp <- gather(dIMMsNoSp)
dIMMsNoSp$Model <- '4'
dIMMsNoSp$param <- 'imm'

dIMMsOP <- ( (estIMMsOP - trueIMMs)/trueIMMs )*100
dIMMsOP <- as_tibble(dIMMsOP)
names(dIMMsOP) <- seq(1,17,1)           
# convert into long format
dIMMsOP <- gather(dIMMsOP)
dIMMsOP$Model <- '2'
dIMMsOP$param <- 'imm'

dIMMsTI <- ( (estIMMsTI - trueIMMs)/trueIMMs )*100
dIMMsTI <- as_tibble(dIMMsTI)
names(dIMMsTI) <- seq(1,17,1)           
# convert into long format
dIMMsTI <- gather(dIMMsTI)
dIMMsTI$Model <- '1'
dIMMsTI$param <- 'imm'

immdiffs <- rbind(dIMMsNoSp, dIMMsOP, dIMMsTI)
immdiffs <- immdiffs %>% 
  dplyr::select( 'Model', 'param', 'value') %>%  
  mutate( Param = 'eta[t]' )

## To do a grouped boxplot, make params a field
Pcs <- rbind(Pcs, immdiffs)

# grouped boxplot
ggplot(Pcs, aes(x=Param, y=value, fill=Model)) + 
  geom_boxplot() +
  # ylim(-100, 100) +
  theme( axis.text.x = element_blank(),
         axis.title.x = element_blank() ) +
  geom_hline( yintercept=0 ) + 
  facet_wrap(~Param, scales = 'free', labeller = label_parsed) +
  ylab("Difference") 

# Boxplots for each
# dNoSP <- ggplot(dIMMsNoSp, aes(x=key, y=value) ) + 
#   geom_boxplot() +
#   ylim(-0.15, 0.15) +
#   geom_hline( yintercept=0 )
# dOP <- ggplot(dIMMsOP, aes(x=key, y=value) ) + 
#   geom_boxplot() +
#   ylim(-0.15, 0.15) +
#   geom_hline( yintercept=0 )
# dTI <- ggplot(dIMMsTI, aes(x=key, y=value) ) + 
#   geom_boxplot() +
#   ylim(-0.15, 0.15) +
#   geom_hline( yintercept=0 )
# 
# ggpubr::ggarrange(dNoSP, dOP, dTI,
#                   labels = c("A", "B", "C"),
#                   ncol = 1, 
#                   nrow = 3)

############################################# Ns
# create dataframe for boxplots
dNs <- ( (estNs - trueNs)/trueNs )*100
dNs <- as_tibble(dNs)
names(dNs) <- seq(1,19,1)           
# convert into long format
dNs <- gather(dNs)
dNs$Model <- '5'
dNs$param <- 'N'

dNsNoSp <- ( (estNsNOSp - trueNs)/trueNs )*100
dNsNoSp <- as_tibble(dNsNoSp)
names(dNsNoSp) <- seq(1,19,1)           
# convert into long format
dNsNoSp <- gather(dNsNoSp)
dNsNoSp$Model <- '4'
dNsNoSp$param <- 'N'

dNsNo <- ( (estNsNO - trueNs)/trueNs )*100
dNsNo <- as_tibble(dNsNo)
names(dNsNo) <- seq(1,19,1)           
# convert into long format
dNsNo <- gather(dNsNo)
dNsNo$Model <- '3'
dNsNo$param <- 'N'

dNsOP <- ( (estNsOP - trueNs)/trueNs )*100
dNsOP <- as_tibble(dNsOP)
names(dNsOP) <- seq(1,19,1)           
# convert into long format
dNsOP <- gather(dNsOP)
dNsOP$Model <- '2'
dNsOP$param <- 'N'

dNsTI <- ( (estNsTI - trueNs)/trueNs )*100
dNsTI <- as_tibble(dNsTI)
names(dNsTI) <- seq(1,19,1)           
# convert into long format
dNsTI <- gather(dNsTI)
dNsTI$Model <- '1'
dNsTI$param <- 'N'

Ndiffs <- rbind(dNs, dNsNoSp, dNsNo, dNsOP, dNsTI)
Ndiffs <- Ndiffs %>% 
  dplyr::select( 'Model', 'param', 'value') %>% 
  mutate( Param = 'N[t]' )

## To do a grouped boxplot, make params a field
Pcs <- rbind(Pcs, Ndiffs)

# grouped boxplot
ggplot(Pcs, aes(x=Param, y=value, fill=Model)) + 
  geom_boxplot() +
  # ylim(-100, 100) +
  theme( axis.text.x = element_blank(),
         axis.title.x = element_blank() ) +
  geom_hline( yintercept=0 ) + 
  facet_wrap(~Param, scales = 'free', labeller = label_parsed) +
  ylab("Difference") 

# ### boxplots
# dN <- ggplot(dNs, aes(x=key,y=value) ) + 
#   geom_boxplot() +
#  # ylim(-400, 600) +
#   geom_hline( yintercept=0 )
# dNoSP <- ggplot(dNsNoSp, aes(x=key,y=value) ) + 
#   geom_boxplot() +
#  # ylim(-400, 600) +
#   geom_hline( yintercept=0 )
# dNo <- ggplot(dNsNo, aes(x=key,y=value) ) + 
#   geom_boxplot() +
# #  ylim(-400, 600) +
#   geom_hline( yintercept=0 )
# dOP <- ggplot(dNsOP, aes(x=key,y=value) ) + 
#   geom_boxplot() +
# #  ylim(-400, 600) +
#   geom_hline( yintercept=0 )
# dTI <- ggplot(dNsTI, aes(x=key,y=value) ) + 
#   geom_boxplot() +
#  # ylim(-400, 600) +
#   geom_hline( yintercept=0 )
# 
# ggpubr::ggarrange(dN, dNoSP, dNo, dOP, dTI,
#                   labels = c("A", "B", "C", "D", "E"),
#                   ncol = 1, 
#                   nrow = 5)

# create dataframe for superpop size
diffSupsNoSp <- tibble(
  # Nsup = estSupsNOSp[,1] - trueSups[,1],
  Nold = ( (estSupsNOSp[,2] - trueSups[,2])/trueSups[,2] )*100,
  Model = '4' #NonSCR'
)

diffSupsNo <- tibble(
  # Nsup = estSupsNO[,1] - trueSups[,1],
  Nold = ( (estSupsNO[,2] - trueSups[,2])/trueSups[,2] )*100,
  Model = '3' #OPSCR_noTE'
)

diffSupsOP <- tibble(
  # Nsup = estSupsOP[,1] - trueSups[,1],
  Nold = ( (estSupsOP[,2] - trueSups[,2])/trueSups[,2] )*100,
  Model = '2' #OPSCR'
)

diffSupsTI <- tibble(
  # Nsup = estSupsTI[,1] - trueSups[,1],
  Nold = ( (estSupsTI[,2] - trueSups[,2])/trueSups[,2] )*100,
  Model = '1' #TI_OPSCR'
)

NSUPdiffs <- rbind(diffSupsNoSp, diffSupsNo, diffSupsOP, diffSupsTI)
## To do a grouped boxplot, make params a field
NSUPdiffs <- NSUPdiffs %>% 
  pivot_longer(cols=c( 'Nold'),
               names_to = 'param',
               values_to = 'value') %>% 
  mutate( Param = 'N[super]' )

## To do a grouped boxplot, make params a field
Pcs <- rbind(Pcs, NSUPdiffs)
# reorder levels for facet wrap
Pcs$Param <- factor(Pcs$Param, levels=c('beta[0]', 'sigma', 'beta[1]',
                                        'alpha[0]', 'alpha[t]',
                                        'phi[t]', 'eta[t]',
                                        'N[t]', 'N[super]' ) )

# grouped boxplot
ggplot(Pcs, aes(x=Param, y=value, fill=Model)) + 
  geom_boxplot() +
  # ylim(-100, 100) +
  theme( axis.text.x = element_blank(),
         axis.title.x = element_blank() ) +
  geom_hline( yintercept=0 ) + 
  facet_wrap(~Param, 
             scales = 'free', 
             labeller = label_parsed) +
  theme(strip.text = element_text(hjust = 0)) +
  ylab("Difference") 


## End ##