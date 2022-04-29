##  Code associated with paper on
## 
##      Surviving environmental change: when increasing
##        population size can increase extinction risk
## 
##       --- by Mark M Tanaka and Lindi M Wahl ---
##
##   Evolutionary rescue with standing variation given by Wright's
##   steady-state distribution assuming a fitness cost to the survival
##   allele before the change in environment, and assuming reversible
##   mutation at equal rates.
## 
##   - Plot Wright's distribution and probability an allele is present 
##   - Plot mean allele frequency at steady state
##   - Plot Wright's distribution vs Wright-Fisher simulation
##   - Simulate surviving env change with Wright steady state as init cond
##       - include sims where surv allele cannot be fixed before env change
##       - include Ricker version of the simulation
##   - Plot probabilities of survival as function of population size
##      - analytical and simulated
##      - components of P_survival
##   - Plot locally worst population size
## 
##====================================================================



## Set the default parameters for most of the analyses
setp <- function(sb = 0.05,      # rescue selective effect
                 delta = 0.01,   # population decline rate after change
                 sd = 0.001,     # fitness cost of rescue allele
                 N0 = 200,       # initial population size at change
                 mu = 10^(-5),   # mutation rate
                 resc0 = 5,      # initial number of rescue allele copies
                 simnum = 5000,  # how many simulations per point
                 p0.threshold = 1.1, # take init freq under this value
                                     # if >1 that means use all freqs
                 density = 12,   # num values of N to sample in sims
                 model = 1,      # 1=classic (Orr-Unckless), 2=Ricker
                 r = 0.05    # intrinsic growth parameter [for Ricker only]
                 ) {
   list(sb=sb, delta=delta, sd=sd, 
        N0=N0, mu=mu,
        resc0=resc0, simnum=simnum,
        p0.threshold=p0.threshold, density=density,
        model=model, r=r
        )
}
## Print with line ending
catn <- function(...){
   cat(...,"\n")
}


##==========================================================
## ----- Wright's (1931) stationary distribution -----
##==========================================================

## Wright: probability density function
wright.pdf <- function (x, N, sd, mu) { # Density function
   require(gsl)       # use GNU Scientific Library 
   alpha <- 2*N*mu    # scaled mutation rate
   beta <- 2*N*mu     # let mutation rates be equal for now
   gamma <- -2*N*sd   # make effect of sd negative 
   b <- beta(alpha,beta)  # beta function 
   f1 <- hyperg_1F1(beta,alpha+beta,gamma)  # confluent hypergeom fn from gsl
   phiC <- exp(gamma*(x)) * x^(alpha-1) * (1-x)^(beta-1)
   phiC/(b*f1)  # normalise to return pdf 
}
## Discretised Wright dist, integrating in binned freqs;
##   bins are (0,0.5/N), (0.5/N, 1.5/N), etc 
tabulate.wright.hist <- function(N, sd, mu){
   ##  N+1 bins centred on 1/N, ..., 1-1/N but the two outer bins
   ##  integrated over 0 to 1/(2N) and over 1-1/(2N) to 1
   bins <- N
   x <- seq(1/N, 1, length.out=bins) -1/(2*N)
   from <- c(0,x)
   to <- c(x,1) 
   wpmf <- c()
   for (i in 1:(bins+1)) {
      ## integrate between breakpoints
      wpmf[i] <- integrate(wright.pdf, from[i], to[i], N, sd, mu)$value
   }
   ## centres of the bins: 
   cen <- seq(0, 1, length.out=bins+1) 
   list(freq=cen, prob=wpmf)
}
## Kimura's (1963) approach to numerical integration of Wright dist
##  with equal mutation rates and for haploids
kimura.integrates.wright <- function(N, sd, mu){
   require(gsl)       # GNU Scientific Library 
   alpha <- 2*N*mu
   gamma <- -2*N*sd   
   b <- beta(alpha,alpha)
   f1 <- hyperg_1F1(alpha,2*alpha,gamma)
   Qx <- function(x, alpha, gamma){ # special fn to integrate
      px <- exp(x*gamma)*(x*(1-x))^(alpha-1)
      px - x^(alpha-1) - exp(gamma)*(1-x)^(alpha-1)
   }
   y <- seq(1/N, 1, length.out=N) -1/(2*N)  # bins 
   y <- c(y,1) # attach 1 at end 
   cdf <- c()  # vector of integrated densities 
   for (i in 1:(N+1)) {
      qint <- integrate(Qx, 0, y[i], alpha, gamma)$value
      num <- qint + y[i]^alpha/alpha +(1-(1-y[i])^alpha)*exp(gamma)/alpha
      cdf[i] <- num / (b*f1)
   }
   cdf <- c(0,cdf)  # add zero to front of the cdf vector
   cen <- seq(0, 1, length.out=N+1) # centres of bins
   prob <- cdf[2:(N+2)]-cdf[1:(N+1)] # N+1 elements
   list(freq=cen, prob=prob) 
}

## Draw random from Wright integrated into bins 
##  - try the kimura-integration version for small N
rwright <- function(n, N, sd, mu){
   ## An ad hoc criterion for when to switch to Kimura integration:
   ##   if (N*sd>5 && N*mu>0.0001) # trial and error thresholds
   ##      dist <- tabulate.wright.hist(N,sd,mu)
   ##   else # this one is better for smaller N*sd and N*mu: 
   ##      dist <- kimura.integrates.wright(N,sd,mu)

   ## Alternative approach: try kimura first
   dist <- kimura.integrates.wright(N,sd,mu)
   ## Switch to integration if Kimura method causes a problem: 
   if (min(dist$prob)<0)  # numerical problem so try direct integration
      dist <- tabulate.wright.hist(N,sd,mu)
   ## Randomise (though maybe not necessary): 
   sample(x=dist$freq, size=n, replace=T, prob=dist$prob)
}

## ----- Mutation-selection-drift balance -----

## (Deterministic) mean freq of allele with
##   fitness cost s and equal mutation rates
##   at equilibrium (for large N) 
mean.freq.determ <- function(s,mu){
   (s+2*mu - sqrt(s^2 + 4*mu^2))/(2*s)  # "exact"
   mu/((s+2*mu))   # approx 
}
## (Stochastic) mean freq of allele with
##   fitness cost s and equal mutation rates
##   under the Wright steady state distribution 
##  - Kimura et al '63 version with negative selection 
mean.freq <- function(N,s,mu){
   require(gsl)      # GNU Scientific Library 
   alpha <- 2*N*mu   # again symmetric mutation for now; can generalise later
   gamma <- -2*N*s   # 
   f1.top <- hyperg_1F1(alpha+1, 2*alpha+1,gamma) 
   f1.bot <- hyperg_1F1(alpha, 2*alpha,gamma)
   (1/2)*f1.top/f1.bot
}

## Plot density of Wright's steady state distribution
curve.wright.pdf <- function(from,to,N, sd, mu, ...){
   L <- 1000
   freq <- seq(from,to, length.out=L)
   wpdf <- c()
   for (i in 1:L) {
      wpdf[i] <- wright.pdf(freq[i], N, sd, mu)
   }
   lines(freq,wpdf , ...)
}
## Plot Pr(an allele is present in the population)
##  integrate Wright pdf from 1/N for a range of N values
curve.resc.present <- function(sd, mu, ...){
   L <- 200
   N_ <- 10^seq(2, 5.5, length=L)
   Pp <- c()
   for (i in 1:L)
      Pp[i] <- integrate(wright.pdf, 1/N_[i], 1, N_[i], sd, mu)$value
   lines(N_, Pp, ...)
}

## Make a 2x2 figure showing Wright distribution and features
plot.wright.aspects <- function(){
   require("sfsmisc")
   p <- setp()   
   #---
   pdf("fig-Wright-dist.pdf",width=6.5, height=6.5)
   par(mfrow=c(2,2), cex=0.99, mar=c(4.1, 4, 1.2, 1) )
   ## Set the layout
   layout(
      matrix(c(1,2,3,4), ncol=2, byrow=TRUE), 
      widths=c(1,1), heights=c(3,4) )

   ## ===== Top row =====
   sd <- p$sd
   mu <- p$mu

   ## Plot Wright PDF
   ## Top left: near left boundary for different N values
   from <- 0.0001; to <- 0.05
   plot(1,1, type="n", xlim=c(0,to), ylim=c(10^-6, 80) #, log="y"
      , xlab="(Low) frequency of allele", ylab="Probability density")
   title("A", adj=0)
   curve.wright.pdf(from,to,1000, sd, mu, col="darkorange",lty=5)
   curve.wright.pdf(from,to,10000, sd, mu, col="forestgreen")
   curve.wright.pdf(from,to,100000, sd, mu, col="purple", lty=6)
   ## Top right: near right boundary 
   from <- 0.95; to <- 0.99999
   plot(1,1, type="n", xlim=c(from,1), ylim=c(0, 80) #, log="y"
      , xlab="(High) frequency of allele", ylab="Probability density")

   curve.wright.pdf(from,to,1000, sd, mu, col="darkorange",lty=5)
   curve.wright.pdf(from,to,10000, sd, mu, col="forestgreen")
   curve.wright.pdf(from,to,100000, sd, mu, col="purple",lty=6)
   legend("topleft", c("N = 1000", "N = 10,000", "N = 100,000"),
          col=c("darkorange","forestgreen","purple"), lty=c(5,1,6)
        , bty="n", cex=0.8)

   ## ===== Bottom row =====
   ## Bottom left:  mean freq -- vary s and mu 
   s=c(0.001, 0.02)   # two values of selective effect s
   mu=c(10^-4, 10^-6) # two values of mutation rate mu 
   N.arr <- 10^seq(1,5.3,length=50)  # array of pop size N
   mf1 <- mean.freq(N.arr, s=s[1], mu=mu[1]) 
   mf3 <- mean.freq(N.arr, s=s[1], mu=mu[2]) 
   mf2 <- mean.freq(N.arr, s=s[2], mu=mu[1]) 
   mf4 <- mean.freq(N.arr, s=s[2], mu=mu[2])
   h1  <- mean.freq.determ(s[1], mu[1])
   h3  <- mean.freq.determ(s[1], mu[2])
   h2  <- mean.freq.determ(s[2], mu[1])
   h4  <- mean.freq.determ(s[2], mu[2])
   L <- c() # labels for legend 
   L[1] <- expression(paste(s[d]," = 0.001, ", mu, " = ", 10^-4))
   L[2] <- expression(paste(s[d]," = 0.001, ", mu, " = ", 10^-6))
   L[3] <- expression(paste(s[d]," = 0.02,  ", mu, " = ", 10^-4))
   L[4] <- expression(paste(s[d]," = 0.02,  ", mu, " = ", 10^-6))
   
   plot(N.arr, mf1, t="l", log="xy", ylim=c(0.00005,2)
      , xaxt="n", yaxt="n"
      , xlab="Population size N"
      , ylab="Expected frequency") # under Wright's distribution")
   title("B", adj=0)
   eaxis(1, at.small=F)
   axt1 <- c(0.1, 0.01, 0.001, 0.0001) 
   axt2 <- c(0.5, 0.1, 0.01, 0.001, 0.0001) 
   eaxis(2, at.small=F, at=axt2, labels=c(0.5, pretty10exp(axt1,drop.1=T)))
   lines(N.arr, mf3, col=1, lty=2)
   lines(N.arr, mf2, col=2)
   lines(N.arr, mf4, col=2, lty=2)
   legend("topright", L
        , col=c(1,1,2,2), lty=c(1,2,1,2), cex=0.85
        , bty="n")
   c2 <- rgb(red=255/255,green=0,blue=0, alpha=0.3) # transparent red
   c1 <- rgb(red=0,green=0,blue=0, alpha=0.3)       # transparent black
   abline(h=c(h1,h3,h2,h4), lty=c(1,2,1,2),col=c(c1,c1,c2,c2))

   ## Bottom right
   ## Now plot Prob(allele present) vs N
   sd <- p$sd  # redefine these 
   mu <- p$mu
   plot(1,1, type="n", ylim=c(0,1), xlim=c(10^2,10^5.5), log="x"
      , xaxt="n"
      , xlab="Population size", ylab="Probability allele is present")
   eaxis(1, at.small=F, n.axp=1) 
   curve.resc.present(sd=0.0001, mu, col="blue",lty=4)
   curve.resc.present(sd=0.001, mu, col="black",lty=1)
   curve.resc.present(sd=0.01, mu, col="red",lty=2)
   legend("topleft", c(expression(s[d]*" = 0.0001"),
                       expression(s[d]*" = 0.001"),
                       expression(s[d]*" = 0.01")),
          col=c(4,1,2), lty=c(4,1,2), bty="n", cex=0.85)
   title("C", adj=0)

   dev.off()
}

## (expected) absolute number of allele copies 
plot.resc.number <- function(sd, mu, ...) {
   s=c(1e-4,0.001, 0.01, 0.05)   # values of selective effect s

   mu=10^-5
   N.arr <- 10^seq(2,5.3,length=100)  # array of pop size N
   mf1 <- mean.freq(N.arr, s=s[1], mu=mu) 
   mf2 <- mean.freq(N.arr, s=s[2], mu=mu) 
   mf3 <- mean.freq(N.arr, s=s[3], mu=mu) 

   an1 <- mf1*N.arr # absolute num = mean freq x N 
   an2 <- mf2*N.arr # 
   an3 <- mf3*N.arr # 

   pdf("exp-num-allelles.pdf",width=3.3,height=3.3)
   par(mar=c(4,4,1,1), cex=0.9)
   plot(N.arr, an1, t="l", log="xy", ylim=c(0.5,10000)
      , xaxt="n", yaxt="n"
      , xlab="Population size N" , ylab="Expected number of alleles" )
   lines(N.arr, an2, col=2) 
   lines(N.arr, an3, col=3) 

   eaxis(1, at.small=F, n.axp=1)
   eaxis(2, at.small=F, n.axp=1)

   legend("topleft", c(expression(paste(s[d]," = 0.0001")),
                       expression(paste(s[d]," = 0.001")),
                       expression(paste(s[d]," = 0.01"))
                       )
        , col=c(1,2,3,4), lty = 1, cex=0.85, bty="n"
          )
   dev.off()
}


##----------------------------------------------------------
## Compare Wright's classic diffusion-based distribution
##  against simulation of Wright-Fisher process
##   for small populations
##----------------------------------------------------------

## Simulate two-type WF model with reversible mutation
##    ssize = sample size
##       N = population size
## fitness = VECTOR of fitnesses
##      mu = VECTOR of two mutation rates
##    skip = how many generations to skip before sampling again.
## Note: this samples the frequency of the FIRST allele p[1] only
## 
##  Note: Link and Eaton 2011 Methods Ecol Evol argue against thinning
##   to avoid autocorrelation; precision and efficiency more important. 
wright.fisher.sample <- function(ssize,N,
                                 fitness=c(1,1),mu=c(10^-4,10^-4),
                                 skip=0){
   p <- c(0.5, 0.5) #  hard-coding initial frequency
   ## thin: sample by skipping every <skip> steps
   ## burnin: start sampling after burnin steps 
   ##  - set this to avg substitution time for neutral mutant
   ##    (using the first mutation rate supplied)
   burnin <- round(1/mu[1] + 2*N)
   tmax <- burnin + (skip+1)*(ssize)  ## how long to run
   sample <- c() # initialise sample vector
   for (i in 1:tmax) {
      ns <- t(rmultinom(1,size=N, prob=p*fitness))
      muts <- rbinom(2, ns, c(mu[1],mu[2]))       # generate mutants
      ## muts[1] means num mutants from type 1 --> type 2
      ns[1] <- ns[1] - muts[1] + muts[2]
      ns[2] <- ns[2] + muts[1] - muts[2]
      p <- ns/N  # normalise
      if (i>burnin) { # then start sampling 
         if ((i-burnin) %% (skip+1) == 0) {
            sample <- c(sample, p[1]) 
         }
      } 
   }
   sample      # return the FIRST frequency
}

## Run WF sim for particular params, and collect data into a file
## Here, let mu be the single mutation rate used for both directions.
## seed : for random number generator
## pset : label for parameter set
run.wf <- function(seed=1,pset,ssize,N,sd,mu,skip){
   set.seed(seed)
   catn("Chain #",seed) # each chain gets a random seed 

   ## simulate Wright-Fisher with symmetric mutation
   ##  - fitness cost sd applies to first allele: 
   wfs <- wright.fisher.sample(ssize,N,c(1-sd,1), c(mu,mu),skip)
   b <- (seq(-1,N)+0.5)/N  # breaks to make freq table 
   hwfs <- hist(wfs, breaks=b, plot=F) # make freq table with sim data
   ## Store results, label file with seed number
   fn <- paste("dat-wf-ps",pset,"-",seed,sep="")
   sink(fn,split=T)
   catn("#  ssize N sd mu") 
   catn("# ", ssize, N, sd, mu)
   sink()
   write.table(cbind(hwfs$mids, hwfs$counts)
             , file=fn, quote=F, col.names=F, row.names=F
             , append=T)
}

## Run simulations for a single parameter set (not the default ones)
##     on multiple cores [set run.sims=T]
## Read in data files and plot 
##  chains = number of chains of simulations
##   ssize = sample size per simulated chain 
##    skip = thinning parameter
collect.and.plot <- function(run.sims=F, cores=8, chains=16, ssize=10^5,skip=0,
                             N=10,sd=0.001,mu=10^-4) {
   require("parallel")
   f <- function(s){run.wf(seed=s,pset=1,ssize,N,sd,mu)}
   if (run.sims)  # if simulations need to be (re-)generated
      mclapply(c(1:chains), f, mc.cores=cores) 

   ## Plot Wright's theoretical distribution (diffusion approx)
   twh <- tabulate.wright.hist(N, sd, mu)
   pdf("WrightFisher-diffusion-vs-sim.pdf", height=3,width=3)
   par(mar=c(4,4,1,1), cex=0.7)
   plot(twh$freq, twh$prob, log="y" 
      , t="l" , yaxt="n"
      , xlab="Frequency", ylab="Probability") 
   eaxis(2, f.smalltcl=0, n.axp=1)
   ## Read and combine sim data
   cumul <- rep(0,N+1)   # accumulate counts in data file 
   for (i in 1:chains) {
      name <- paste("dat-wf",i,sep="") 
      d <- read.table(name)
      cumul <- cumul + d$V2
   }
   ## Plot combined data 
   tot <- sum(cumul)
   points(d$V1, cumul/tot, pch=15, col=2)   
   legend("top",c("Diffusion approx", "Wright-Fisher sim")
        , col=c(1,2), pch=c(NA,15), lty=c(1,0))
   dev.off()
}

## For plotting Wright's distribution and Wright-Fisher sims,
##  instead of using a default parameter set define
##   many sets of parameters 
choose.pset <- function(ps){
   if      (ps==1) p <- list(N=100, sd=0.001,mu=10^-4 )
   else if (ps==2) p <- list(N=100, sd=0.01, mu=10^-4 )
   else if (ps==3) p <- list(N=100, sd=0.001,mu=10^-5 )
   else if (ps==4) p <- list(N=100, sd=0.01, mu=10^-5 )
   else if (ps==5) p <- list(N=200, sd=0.001,mu=10^-4 )
   else if (ps==6) p <- list(N=200, sd=0.01, mu=10^-4 )
   else if (ps==7) p <- list(N=200, sd=0.001,mu=10^-5 )
   else if (ps==8) p <- list(N=200, sd=0.01, mu=10^-5 )
   else {p <- NULL}
   p
}
## Run multiple sets and make multi-panel fig
## Run each of 8 sets separately
run.wf.multi <- function(run.ps=1, chains=32){
   ## run.ps is a label to specify which set to run.
   require("parallel")  # to use multiple cores 
   ## Some hard-coded simulation params:
   cores <- 8;    ssize <- 10^5;   skip=0
   p <- choose.pset(run.ps) 
   f <- function(s){ #redefine function with new param set pset
      run.wf(seed=s,pset=run.ps, ssize, N=p$N, sd=p$sd, mu=p$mu, skip=skip)}
   ## If out of bounds then don't run: 
   if (run.sims>0 && run.sims<9)  # if simulations need to be (re-)generated
      mclapply(c(1:chains), f, mc.cores=cores)
}

## Plot the WF sim data; assume all data files are present
plot.wf.multi <- function(chains=c(32,32,32,64,
                                   32,32,32,64)){
   require("sfsmisc") 
   ## Set up the plot 
   pdf("wf-diffusion-v-sim.pdf", height=4.3)
   par(mfrow=c(2,4), oma=c(3,3,0.2,0.7), mar=c(2,2,1,1), cex=0.55
     , xpd= TRUE)
   
   for (j in 1:8) { # each parameter set
      p <- choose.pset(j)   # Get parameter values for set j
      mt <- paste("sd=",p$sd,", mu=", p$mu, sep="")

      ## Plot Wright's theoretical distribution (diffusion approx)
      twh <- tabulate.wright.hist(p$N, p$sd, p$mu)
      plot(twh$freq, twh$prob,  log="y", t="n" # don't plot yet 
         , ylim=c(10^-5.4, 1)
         , xlab="", ylab="", yaxt="n")
      eaxis(2, f.smalltcl=0, n.axp=1)
      mtext(mt,3,cex=0.7, line=-1.2)
      ## Read and combine sim data
      catn("ps number: ",j)
      cumul <- rep(0,p$N+1)   # accumulate counts in data file 
      for (i in 1:chains[j]) {
         name <- paste("dat-wf-ps",j,"-",i,sep="")
         d <- read.table(name)
         cumul <- cumul + d$V2
      }
      catn(" cumul freqs", cumul) 
      ## Plot combined data 
      tot <- sum(cumul)
      print(tot)
      points(d$V1, cumul/tot, pch=15, col=2)
      lines(twh$freq, twh$prob)
   }
   ## Label axes in margins
   mtext("Frequency",1, outer=T, cex=0.9, line=1.8)
   mtext("Probability",2, outer=T, cex=0.9, line=1.8)
   mtext("N=100", 4, outer=T, cex=0.8, line=-0.5, adj=0.82)
   mtext("N=200", 4, outer=T, cex=0.8, line=-0.5, adj=0.28)
   dev.off()
}


##==========================================================
##   Theoretical probability of rescue:
##       Orr&Unckless + Wright
##   Below, these probabilities will be called "rescue"
##    even though in some cases the survival allele is already
##    fixed in the population before the environment changes.
##==========================================================

## Numerically evaluate exact fixation probability
##  -this is more accurate than 2*sr where sr=sb-delta
pfix.exact <- function(s.r){
   f <- function(x,s.r){
      1 - exp(-(1+s.r)*x) - x }
   sol <- uniroot(f, c(10^-6,0.9), s.r)
   sol$root
}

## The rescue probability averaged over phi(p)
## This one uses the confluent hypergeometric function
##  N = vector of N values 
## pfixmod = model for prob(fixation) 
prob.TW <- function(N,sb,delta,mu, sd, pfixmod=3){
   require("gsl")
   if (pfixmod==1) {         # 1st order approx [classic result]
      pfix <- 2*(sb-delta) 
   } else if (pfixmod==2) {  # 2nd order approx
      pfix <- (sqrt(24*(sb-delta)+9)-3)/2  
   } else {    # "exact" fixation probability [do this by default]
      pfix <- pfix.exact(sb-delta)
   }
   ptw <- c() # vector of probs to return
   for (i in 1:length(N)) {
         alpha <- 2*N[i]*mu  
         gamma <- 2*N[i]*sd 
         lambda <- N[i]*pfix*(1-mu/delta) + gamma
         f1.top <- hyperg_1F1(alpha,2*alpha,-lambda)
         f1.bot <- hyperg_1F1(alpha,2*alpha,-gamma)
         ptw[i] <- 1 - exp(-pfix*N[i]*mu/delta)*f1.top/f1.bot
   }
   ptw  ## return vector of probabilities
}
## This time compute components as marginal probabilities
## P_stand (regardless of de novo) [see Eq 9 in Orr&Unckless 08]
## P_denovo (regardless of standing var)
##  use exact fixation probability 
prob.new.old <- function(N,sb,delta,mu, sd){
   require("gsl")
   pfix <- pfix.exact(sb-delta)
   pst <- pdn <- c() # initialise prob vectors
   for (i in 1:length(N)) {
         alpha <- 2*N[i]*mu  
         gamma <- 2*N[i]*sd 
         lambda.st <- gamma + N[i]*pfix 
         f1.top.st <- hyperg_1F1(alpha,2*alpha,-lambda.st)
         f1.bot <- hyperg_1F1(alpha,2*alpha,-gamma)
         pst[i] <- 1 - f1.top.st/f1.bot
         pdn[i] <- 1 - exp(-pfix*N[i]*mu/delta) 
   }
   list(pst=pst, pdn=pdn)
}


##==========================================================
##    --- RESCUE/SURVIVAL SIMULATIONS ---
##==========================================================

## Run one sim using parameters in p
## Classic rescue model with geometric increase of mutant lineage
onesim.classic <- function(p){
   stopsize <- 2000  # at this size, deem escape successful
   with(p, {
      tcap  <- round(-log(N0)/log(1-delta)) # deterministic time till N=1
      win <- 0            # indicate escape from extinction 
      ## The next one collects num rescue type at each generation
      resc <- c(resc0)    # start array from [1]: resc0
      wt <- c(N0-resc0)   # initial number of wildtype
      for (i in 2:tcap) {  # go up to tcap generations

         if (wt[i-1]>0) { # if any wildtype left 
            wt.os <- rpois(1, wt[i-1]*(1-delta)) # wt offspring for next gen
            muts.f <- rbinom(1, wt.os, mu)       # rescue mutants (wt-->resc)
         } else {
            wt.os <- 0
            muts.f <- 0
         }
         if (resc[i-1]>0) { # if any rescue left
            ## carriers of rescue allele in next gen
            resc.os <- rpois(1, resc[i-1]*(1+sb-delta))
            muts.b <- rbinom(1, resc.os, mu)    # back mutants (resc-->wt)
         } else {
            resc.os <- 0
            muts.b <- 0
         }
         
         wt[i] <- wt.os     -muts.f +muts.b  # move mutants
         resc[i] <- resc.os -muts.b +muts.f  # move mutants
         
         if (resc[i]>stopsize) { # stop if many rescue 
            win <- 1
            break
         } else if (resc[i]==0 && wt[i]==0) { # stop if pop extinct
            win <- 0
            break
         }         
      }
      list(win=win,  wt=wt, resc=resc)
   })
}

## Run one sim using parameters in p
## Ricker model 
onesim.ricker <- function(p){
   with(p, {
      tcap <- 10*round(-log(N0)/log(1-delta)) # 10x deterministic time till N=1
      win <- 0            # indicate escape from extinction 
      ## The next one collects num rescue type at each generation
      resc <- c(resc0)    # start array from [1]: resc0
      kk <- N0  # carrying capacity of wildtype is initially initial pop size
      wt <- c(N0-resc0)   # initial number of wildtype

      for (i in 2:tcap) {  # go up to tcap generations
         kk <- kk*(1-delta)        # reduction in carrying capacity
         NN <- resc[i-1] + wt[i-1] # pop size last gen
         if (wt[i-1]>0) { # if any wildtype left 
            wt.os <- rpois(1, wt[i-1]*exp(r*(1-NN/kk))) # wt offspr next gen
            muts.f <- rbinom(1, wt.os, mu)       # rescue mutants (wt-->resc)
         } else {
            wt.os <- 0
            muts.f <- 0
         }
         if (resc[i-1]>0) { # if any rescue left
            ## carriers of rescue allele in next gen
            resc.os <- rpois(1, resc[i-1]*exp(r*(1-NN/N0)))
            muts.b <- rbinom(1, resc.os, mu)    # back mutants (resc-->wt)
         } else {
            resc.os <- 0
            muts.b <- 0
         }
         wt[i] <- wt.os     -muts.f +muts.b  # move mutants
         resc[i] <- resc.os -muts.b +muts.f  # move mutants

         if (resc[i] > N0 + sqrt(N0)) { # stop if many rescue 
            win <- 1
            break    
         } else if (resc[i]==0 && wt[i]==0) { # stop if pop extinct
            win <- 0
            break
         }
      }
      list(win=win,  wt=wt, resc=resc)
   })
}

onesim <- function(p){
   if (p$model==1)
      onesim.classic(p)
   else if (p$model==2)
      onesim.ricker(p)
   else
      catn("Unknown model: ",p$model)
}

## run sim once with arbitrary resc0, and plot result [not used]
plot.sim <- function(model=1){
   require("sfsmisc")
   ## To run denovo only, let resc0=0
   out <- onesim(setp(resc0=2, sb=0.02, N0=50000,model=model)) 
   pdf("sim-env-change.pdf", width=3.5,height=3)
   par(mar=c(4,4,1,1), cex=0.5)
   t <- c(1:length(out$wt))
   plot(t,out$resc, t="l", col=2, ylim=c(0.1,max(c(out$resc,out$wt))), 
        ylab="Number of wildtype and mutants", xlab="Time (generations)"
        , log="y", yaxt="n")
   lines(t, out$wt)
   eaxis(2, f.smalltcl=0, n.axp=1)
   dev.off()
   ## Report the result
   if (out$win==1) {
      catn("Evolutionary rescue!")
   } else {
      catn("Population extinction")
   }
}
## Simulate many (p$simnum) times the rescue model with wright standing var
##   dist as initial p, for one particular parameter set
##   - exclude initial frequencies above a threshold
sim.orr.wright <- function (p){
   wfs <- rwright(p$simnum, p$N0, p$sd, p$mu) # draw random Wright freqs
   wfs <- wfs[wfs<p$p0.threshold]  # select ones less than threshold
   len <- length(wfs) 
   wins <- 0  # for counting successful rescue 
   for (i in 1:len) {
      p$resc0 <- round(wfs[i]*p$N0) # initial number of rescue allele copies
      out <- onesim(p)
      if (out$win==1)
         wins <- wins+1
   }
   r <- wins/len                 # proportion rescued
   se <- sqrt(r*(1-r)/len)       # approx SE
   list(r=r, se=se)    
}

## Simulate for many N values and store in external file
##  whose name specifies which param has changed
##    density : number of values of N to sample 
sim.orr.wright.over.N <- function(p, ext="default"){
   N <-  round(10^seq(2,5,len=p$density))

   ## Do simulations and store results 
   fname <- paste("dat-orr-wright-sim-",ext,sep="")
   sink(fname, split=T)
   cat("# ", names(p), "\n")   # store names of params in comment
   cat("# ", as.vector(unlist(p)),"\n") # store values of params in comment
   set.seed(1)   # set seed
   presc <- se <-  c()
   catn("N presc se")  # write header 
   for (i in 1:length(N)) {
      p$N0 <- N[i]    # set pop size
      out <- sim.orr.wright(p)   # run simulation
      presc[i] <- out$r
      se[i] <- out$se
      catn(N[i], presc[i], se[i]) 
   }
   sink()
}

## Uncomment lines to run simulations
many.sims <- function(){

#   sim.orr.wright.over.N(setp(sb=0.4),ext="sb4")  # 
#   sim.orr.wright.over.N(setp(sb=0.05),ext="sb05")  # default 
#   sim.orr.wright.over.N(setp(sb=0.02),ext="sb02") #
#   sim.orr.wright.over.N(setp(sb=0.011),ext="sb011") #
   
#   sim.orr.wright.over.N(setp(mu=1e-4),ext="mu-1e-4") #
#   sim.orr.wright.over.N(setp(mu=1e-6),ext="mu-1e-6")  #  

#   sim.orr.wright.over.N(setp(sd=1e-1),ext="sd-1e-1")  #
#   sim.orr.wright.over.N(setp(sd=1e-2),ext="sd-1e-2")  #
#   sim.orr.wright.over.N(setp(sd=1e-4),ext="sd-1e-4")  #
#   sim.orr.wright.over.N(setp(sd=10^-3.5, density=12),ext="sd-1e-3.5")
#   sim.orr.wright.over.N(setp(sd=10^-3.7, density=14),ext="sd-1e-3.7")
#   sim.orr.wright.over.N(setp(sd=1e-5),ext="sd-1e-5")

#   sim.orr.wright.over.N(setp(delta=0.002),ext="del002")  #
#   sim.orr.wright.over.N(setp(delta=0.02),ext="del02")

   ## ---Survival vs rescue; rare alleles ---> 
#   sim.orr.wright.over.N(setp(p0.threshold=1.1, density=18) ) # default again (incl p0=1)
#   sim.orr.wright.over.N(setp(p0.threshold=1,density=16),ext="p1") # all but p0=1
#   sim.orr.wright.over.N(setp(p0.threshold=0.5, density=14),ext="p0.5")
#   sim.orr.wright.over.N(setp(p0.threshold=0.2, density=17),ext="p0.2")

   ## Another set with sb=0.011
#   sim.orr.wright.over.N(setp(
#      p0.threshold=1.1, sb=0.011, density=18),ext="sb011-p11" ) 
#   sim.orr.wright.over.N(setp(
#      p0.threshold=1,   sb=0.011, density=16),ext="sb011-p1") 
#   sim.orr.wright.over.N(setp(
#      p0.threshold=0.5, sb=0.011, density=15),ext="sb011-p0.5")
#   sim.orr.wright.over.N(setp(
#      p0.threshold=0.909, sb=0.011, density=17),ext="sb011-p0.909")

   ## --- Ricker --->  
#   sim.orr.wright.over.N(setp(model=2,delta=0.01),ext="m2-d01") 
#   sim.orr.wright.over.N(setp(model=2,delta=0.05),ext="m2-d05") 
#   sim.orr.wright.over.N(setp(model=2,delta=0.1),ext="m2-d1")
  
#   sim.orr.wright.over.N(setp(model=2,r=0.001),ext="m2-r001")
#   sim.orr.wright.over.N(setp(model=2,r=0.01),ext="m2-r01")
#   sim.orr.wright.over.N(setp(model=2,r=0.05),ext="m2-r05")
#   sim.orr.wright.over.N(setp(model=2,r=0.1),ext="m2-r1") 
}


##==========================================================
##           PLOT SIM and THEORY 
##==========================================================

## Plot rescue probabilities -- vary parameters
## --- Multi-panel figure --- 
plot.vary.params <- function(){
   setup.panel <- function(lo=2,hi=5,panel.lab){
      require("sfsmisc") 
      plot(1,1, pty="n", xaxt="n", log="x", xlim=c(10^lo,10^hi), ylim=c(0,1.12)
         , xlab="Population size N"   
         , ylab="", las=1, cex.lab=1.2) 
      eaxis(1, f.smalltcl=0, n.axp=1)     # nice x axis
      title(panel.lab, adj=0)
   }
   ## Add a series from theoretical calculation 
   one.series <- function(p,col=1,lo=1,hi=5){
      N_ <- 10^seq(lo,hi, length.out=100)
      presc.tw <- prob.TW(N_, p$sb,p$delta,p$mu,p$sd)
      lines(N_, presc.tw, col=col)
   }
   ## Add points from simulation files 
   add.pts <- function(fname,col=1,bars=0){ 
      d <- read.table(fname,header=T) 
      points(d$N,d$presc, pch=16, cex=0.8,col=col)
      if (bars==1) 
         arrows(d$N,d$presc-1.96*d$se,d$N,d$presc+1.96*d$se,
                length=0.03, angle=90, code=3, cex=0.8, col=col)
   }
   pdf("prob-surv-vary-params.pdf",height=2.6, width=7)
   par(mfrow=c(1,3), mar=c(4, 3, 1.2, 0.5), oma=c(.1, 1.2, .1,.1)
     , cex=0.75) 
   
   ## Vary sb in a panel 
   setup.panel(panel.lab="A")
   mtext("Probability of survival",2,outer=TRUE, cex=0.9, adj=0.8)
   one.series(setp(sb=0.4),4)  # add theoretical curve 
   one.series(setp(sb=0.05),1)
   one.series(setp(sb=0.02),3)
   one.series(setp(sb=0.011),2)
   add.pts("dat-orr-wright-sim-sb4",4) # add simulation points 
   add.pts("dat-orr-wright-sim-sb05",1)
   add.pts("dat-orr-wright-sim-sb02",3)
   add.pts("dat-orr-wright-sim-sb011",2)
   legend("topleft",c( # 
      expression("s"[b]*" = 0.4"),
      expression("s"[b]*" = 0.05"),
      expression("s"[b]*" = 0.02"),
      expression("s"[b]*" = 0.011"))
    , col=c(4,1,3,2),lty=1, cex=0.84, bty="n", seg.len=0.95, lwd=1.2)
   ## Vary mu in a panel 
   setup.panel(panel.lab="B")
   one.series(setp(mu=1e-4),4)   
   one.series(setp())
   one.series(setp(mu=1e-6),2)
   add.pts("dat-orr-wright-sim-mu-1e-4",4)
   add.pts("dat-orr-wright-sim-sb05",1)  # the default again 
   add.pts("dat-orr-wright-sim-mu-1e-6",2)
   legend("topleft",c(expression(mu*"=10"^-4),
                      expression(mu*"=10"^-5),
                      expression(mu*"=10"^-6))
        , col=c(4,1,2),lty=1, cex=0.84, bty="n", seg.len=0.95, lwd=1.2)
   ## Vary sd in a panel 
   setup.panel(panel.lab="C")
   one.series(setp(sd=1e-4),4)
   one.series(setp(sd=10^-3.5),6)
   one.series(setp())
   one.series(setp(sd=1e-2),3)
   one.series(setp(sd=0.1),2)
   add.pts("dat-orr-wright-sim-sd-1e-1",2)
   add.pts("dat-orr-wright-sim-sd-1e-2",3)
   add.pts("dat-orr-wright-sim-sb05",1)  # the default again 
   add.pts("dat-orr-wright-sim-sd-1e-3.5",6)   
   add.pts("dat-orr-wright-sim-sd-1e-4",4)  
#   add.pts("dat-orr-wright-sim-sd-1e-3.5",6)  # extra 
#   add.pts("dat-orr-wright-sim-sd-1e-3.7",7)  # extra
   legend("topleft",c(
      expression("s"[d]*" = 10"^-4),
      expression("s"[d]*" = 10"^-3.5),
      expression("s"[d]*" = 10"^-3),
      expression("s"[d]*" = 10"^-2),
      expression("s"[d]*" = 0.1")
      ),
      col=c(4,6,1,3,2),lty=1, cex=0.84, bty="n", seg.len=0.95, lwd=1.2)
   dev.off()
   
   ## And delta separately 
   pdf("prob-surv-vary-delta.pdf",height=2.6, width=2.6)
   par(mfrow=c(1,1), mar=c(4, 3, 1.2, 0.5)
     , oma=c(.1, 1.2, .1,.1)
     , cex=0.75) 
   
   setup.panel(panel.lab="")
   mtext("Probability of survival",2,outer=T, cex=0.9, adj=0.8)
   one.series(setp(delta=0.002),4)  # add theoretical curve 
   one.series(setp(delta=0.01),1)
   one.series(setp(delta=0.02),3)
   add.pts("dat-orr-wright-sim-del002",4) # add simulation points 
   add.pts("dat-orr-wright-sim-sb05",1)   # default set 
   add.pts("dat-orr-wright-sim-del02",3)
   legend("topleft",c( # 
      expression(delta*" = 0.002"),
      expression(delta*" = 0.01"),
      expression(delta*" = 0.02"))
    , col=c(4,1,3),lty=1, cex=0.84, bty="n", seg.len=0.95)
   dev.off()
}

plot.ricker.case <- function(){
   require("sfsmisc")
   cola=c("tomato","purple","navy","forestgreen") # colours 
   ## Vary delta 
   pdf("prob-surv-ricker-delta.pdf", height=2.6,width=2.6)
   d <- read.table("dat-orr-wright-sim-m2-d01",header=T)
   e <- read.table("dat-orr-wright-sim-m2-d05",header=T)
   f <- read.table("dat-orr-wright-sim-m2-d1",header=T)
   par(mar=c(4,4,1,1) , cex=0.7)
   plot(d$N,d$presc, cex=0.9, las=1, xaxt="n", log="x"
      , xlab="Population size N", ylab="Probability of survival"
      , pch =16, col=cola[1], ylim=c(0,1.12),
        main="Ricker model", cex.main=0.99
        )
   points(e$N, e$presc, pch=16, col=cola[2])
   points(f$N, f$presc, pch=16, col=cola[3])
   eaxis(1, f.smalltcl=0, n.axp=1)
   legend("topleft",c(
      expression(delta*" = 0.01"),
      expression(delta*" = 0.05"),
      expression(delta*" = 0.1"))
    , col=cola, pch=16, cex=0.9, bty="n") 
   dev.off()
   ## vary r
   pdf("prob-surv-ricker-r.pdf", height=2.6,width=2.6)
   d <- read.table("dat-orr-wright-sim-m2-r1",header=T)
   e <- read.table("dat-orr-wright-sim-m2-r05",header=T)
   f <- read.table("dat-orr-wright-sim-m2-r01",header=T)
   g <- read.table("dat-orr-wright-sim-m2-r001",header=T)
   par(mar=c(4,4,1,1)  , cex=0.75)
   plot(d$N,d$presc,  las=1, xaxt="n", log="x"
      , xlab="Population size N", ylab="Probability of survival"
      , cex.lab=1.2
      , pch =16, col=cola[1], ylim=c(0,1.12),
        main="Ricker model" , cex.main=0.9
        )
   lines(d$N, d$presc, pch=16, col=cola[1])
   lines(e$N, e$presc, pch=16, col=cola[2])
   lines(f$N, f$presc, pch=16, col=cola[3])
   lines(g$N, g$presc, pch=16, col=cola[4])
   points(e$N, e$presc, pch=16, col=cola[2])
   points(f$N, f$presc, pch=16, col=cola[3])
   points(g$N, g$presc, pch=16, col=cola[4])
   eaxis(1, f.smalltcl=0, n.axp=1)
   legend("topleft",c(
      expression(r*" = 0.1"),
      expression(r*" = 0.05"),
      expression(r*" = 0.01"),
      expression(r*" = 0.001")
   )
    , col=cola, pch=16, cex=0.9, bty="n") 
   dev.off()
}

## This time, plot runs where only rare initial frequencies are used
plot.init.rare.allele <- function(){
   require("sfsmisc") 
   pdf("prob-rescue-v-survival.pdf", height=2.7,width=5)
   co <- c("black", "purple", "dodgerblue","limegreen") 
   par(mfrow=c(1,2), mar=c(4,4,1.2,1) , cex=0.7)
   ## first with sb=0.05 
   fnm <- c("dat-orr-wright-sim-default",
           "dat-orr-wright-sim-p1",
           "dat-orr-wright-sim-p0.5",
           "dat-orr-wright-sim-p0.2")
   plot(1,1, t="n", log="x", cex=0.95, xlim=c(10^2,10^5), ylim=c(0,1.00)
      , las=1, xaxt="n", main=expression(s[b]*" = 0.05")
      , xlab="Population size N", ylab="Probability of rescue or survival"
      , cex.lab=1.2) 
      eaxis(1, f.smalltcl=0, n.axp=1)     # nice x axis
   for (i in 1:length(co)) {
      d <- read.table(fnm[i],header=T)
      points(d$N,d$presc, pch=16, col=co[i] , cex=0.8)      
      lines(d$N,d$presc, col=co[i], lty=1, lwd=0.7)
   }
   legend("topleft", c("all p", "p<1", "p<0.5", "p<0.2"),
          lty=1, lwd=0.7, pch=16, col=co, cex=0.9, bty="n")
   ## second with sb=0.011
   co <- c("black", "purple", "darkorange", "dodgerblue") 
   fnm <- c("dat-orr-wright-sim-sb011-p11",
           "dat-orr-wright-sim-sb011-p1",
           "dat-orr-wright-sim-sb011-p0.909",
           "dat-orr-wright-sim-sb011-p0.5")
   plot(1,1, t="n", log="x", cex=0.95, xlim=c(10^2,10^5), ylim=c(0,1.00)
      , las=1, xaxt="n", main=expression(s[b]*" = 0.011")
      , xlab="Population size N", ylab="Probability of rescue or survival"
      , cex.lab=1.2) 
      eaxis(1, f.smalltcl=0, n.axp=1) 
   for (i in 1:length(co)) {
      d <- read.table(fnm[i],header=T)
      points(d$N,d$presc, pch=16, col=co[i] , cex=0.8)      
      lines(d$N,d$presc, col=co[i], lty=1, lwd=0.7)
   }
   legend("topleft", c("all p", "p<1", "p<0.909", "p<0.5"),
          lty=1, lwd=0.7, pch=16, col=co, cex=0.9, bty="n")
   
   dev.off()
}

## Plot the components of the rescue probability 
plot.components <- function(lo=2, hi=6){
   require("sfsmisc")
   p <- setp()
   N_ <- 10^seq(lo,hi, length.out=100)
   pdf("p-rescue-components.pdf", height=2.6,width=3.4
     , title="Rescue components")
   par(mar=c(4,4,1,1),cex=0.65)
   plot(1,1, pty="n", xaxt="n", log="x", xlim=c(10^lo,10^hi), ylim=c(0,1)
      , xlab="Population size N", ylab="Probability of survival")
   eaxis(1, at.small=F, n.axp=1) 
   cola <- c("black","dodgerblue2","orange","magenta")

   ## A version with marginal probs which don't add to 1
   pno <- prob.new.old(N_, p$sb,p$delta,p$mu,p$sd)
   ptw <- prob.TW(N_, p$sb,p$delta,p$mu,p$sd)
   lines(N_, ptw, col=cola[1] )
   lines(N_, pno$pst, col=cola[2], lty=4)
   lines(N_, pno$pdn, col=cola[3], lty=2)
   abline(v=c(1/(1*p$mu), 1/p$sd,
              p$delta*log(2)/2/(p$sb-p$delta)/p$mu), lty=3, col=c(6,6,6)
            , lwd=c(1,1,1))
   legend("topleft", c("Total","standing variation","de novo")
        , lty=c(1,4,2), col=cola
        , bg="white", cex=0.9)
   ## Add special text to indicate regimes 
   text(12,0.2,"weak\n rescue", cex=0.8)  # <--- off the screen
   text(125,0.54,"strong\n drift", cex=0.8)
   text(3*10^3,0.08,"strong\n selection", cex=0.8)
   text(2*10^5,0.88,"strong\n mutation", cex=0.8)
   
   dev.off()
}


##==========================================================
##    (Locally) worst population size 
##==========================================================

## Find a min in vector starting from the end of the vector
min.from.right <- function(x,y){
   l <- length(y)
   argmin <- -999  # send to a negative num if no local min
   min <- y[1] 
   for (i in (l-1):1) {
      ## Strict inequality (do not count inflection point): 
      if (y[i+1]<y[i]) {  
         min <- y[i+1]
         argmin <- x[i+1]
         break
      } 
   }
   c(argmin,min)
}
## Find the first local minimum of rescue prob wrt N from the right
## In this version, do it in multiple steps from coarse to fine [faster]
get.min.N <- function(p, res=500, say=FALSE, interv=c(10,10^5)){
   ## Round 1
   ##N <- 10^seq(1,5, length.out=res) # logarithmic
   N <- seq(interv[1],interv[2],length.out=res) # linear
   ptw <- prob.TW(N, p$sb,p$delta,p$mu,p$sd) # compute probabilities
   mfr <- min.from.right(N,ptw) # get argmin (and min)
   if (mfr[1]>0) { # if local min found in round 1
      ## Round 2 
      N <- seq(mfr[1]/2, mfr[1]*2, length.out=res) # assume within these bounds
      ptw <- prob.TW(N, p$sb,p$delta,p$mu,p$sd)
      mfr <- min.from.right(N,ptw)
      ## Round 3 - let's go crazy! 
      N <- seq(mfr[1]*(0.7), mfr[1]*(1.3), length.out=res) 
      ptw <- prob.TW(N, p$sb,p$delta,p$mu,p$sd)
      mfr <- min.from.right(N,ptw)
   }
   if (say) { # if you need to know 
      catn(mfr[1], mfr[2])
   }
   mfr[1]
}

## Generate one series of worst-N for a sequence of parameter values
wN.series <- function(seq,sb=0.05,vary="sd",ext="sd-sb05"){
   p <- setp() # default params
   fname <- paste("dat-worst-N-",ext,sep="")
   sink(fname, split=T)
   catn("# ", names(p))   # store names of params in comment
   catn("# ", as.vector(unlist(p))) # store values of params in comment
   catn("x worstN")
   for (i in 1:length(seq) ) {
      if (vary=="sd") 
         p <- setp(sd=seq[i], sb=sb)   # change params
      else # vary mu
         p <- setp(mu=seq[i], sb=sb)
      mn <- get.min.N(p)
      if (mn>0) {  # keep if min N is in range
         catn(seq[i], mn)
      }
   }
   sink()
}
## Generate many series of worst-N
all.wN.series <- function(){
   ## sequences of mu and of sd
   mu.seq <- 10^seq(-6,-4.2, length.out=80) 
   sd.seq <- 10^seq(-4,-1, length.out=80)
   wN.series(mu.seq,sb=0.011,vary="mu",ext="mu-sb011")
   wN.series(mu.seq,sb=0.02,vary="mu",ext="mu-sb02")
   wN.series(mu.seq,sb=0.05,vary="mu",ext="mu-sb05")
   
   wN.series(sd.seq,sb=0.011,vary="sd",ext="sd-sb011")
   wN.series(sd.seq,sb=0.02,vary="sd",ext="sd-sb02")
   wN.series(sd.seq,sb=0.05,vary="sd",ext="sd-sb05")
}

## Plot worst N vs mu, and vs sd -- get data from files 
plot.worst.N <- function(){
   require("sfsmisc") 
   pdf("worst-N-vs-mu.pdf", width=5, height=2.5)
   par(mfrow=c(1,2), mar=c(4.5,4,1.2,1),cex=0.6)
   ## First panel -- worst N vs mu
   wN011 <- read.table("dat-worst-N-mu-sb011", header=T)
   wN02 <- read.table("dat-worst-N-mu-sb02", header=T)
   wN05 <- read.table("dat-worst-N-mu-sb05", header=T)
   plot(wN05$x, wN05$worstN, xlim=c(10^-6, 10^-4)
      , t="n" , log="x", xaxt="n"
      , ylim=c(0,max(wN05$worstN,wN02$worstN,wN011$worstN))
      , ylab="Locally worst population size" 
      , xlab=expression("Mutation rate, "*mu), cex.lab=1.2)
   eaxis(1, n.axp=1, at.small=F) 
   title("A", adj=0)
   lines(wN011$x, wN011$worstN, col=2)
   lines(wN02$x, wN02$worstN, col=4)
   lines(wN05$x, wN05$worstN, col=1)
   legend("topright", c(expression(s[b]*" = 0.011"),
                        expression(s[b]*" = 0.02"),
                        expression(s[b]*" = 0.05")), 
          col=c(2,4,1), lty=1, bty="n") 
   ## Second panel -- worst N vs sd
   wN011 <- read.table("dat-worst-N-sd-sb011", header=T)
   wN02 <- read.table("dat-worst-N-sd-sb02", header=T)
   wN05 <- read.table("dat-worst-N-sd-sb05", header=T)
   plot(wN05$x,wN05$worstN
      , t="n", log="x", xaxt="n"
      , ylim=c(0, max(wN05$worstN,wN02$worstN,wN011$worstN))
      , xlim=c(10^-4, 10^-1)  # setting manually 
      , ylab="Locally worst population size"
      , xlab=expression("Deleterious effect, "*s[d]), cex.lab=1.2)
   eaxis(1, n.axp=1, at.small=F) 
   title("B", adj=0)
   lines(wN011$x, wN011$worstN, col=2)
   lines(wN02$x, wN02$worstN, col=4)
   lines(wN05$x, wN05$worstN, col=1)
   sd.seq <- 10^seq(-4,-1, length.out=100)

   legend("topright", c(expression(s[b]*" = 0.011"),
                        expression(s[b]*" = 0.02"),
                        expression(s[b]*" = 0.05")
                        ), 
          col=c(2,4,1), lty=c(1,1,1), bty="n") 

   dev.off()
}

