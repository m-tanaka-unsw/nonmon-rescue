##  Code associated with paper on
##
##      Evolutionary rescue is non-monotonic in population size
##       --- by Mark Tanaka and Lindi Wahl ---
##
##   Evolutionary rescue with standing variation given by Wright's
##   steady-state distribution assuming a fitness cost to the rescue
##   allele before the change in environment, and assuming reversible
##   mutation at equal rates.
## 
##   - Plot Wright's PDF and probability an allele is present 
##   - Plot mean allele frequency at steady state
##   - Simulate evolutionary rescue with Wright steady state as init cond
##   - Plot probabilities of rescue as function of population size
##      - analytical and simulated
##   - Plot locally worst population size
## 
##====================================================================



## Set the default parameters 
setp <- function(sb = 0.05,      # rescue selective effect
                 delta = 0.01,   # population decline rate after change
                 sd = 0.001,     # fitness cost of rescue allele
                 N0 = 500,       # initial population size at change
                 mu = 10^(-5),   # mutation rate
                 resc0 = 5,      # initial number of rescue allele copies
                 simnum = 5000   # how many simulations per point
                 ) {
   ## Find time till pop size gets to less than 1 under the assumption of
   ## deterministic decline at rate delta regardless of whether mutants
   ## appear; use a multiple of this in the simulations. 
   tmax <- round(-log(N0)/log(1-delta))+1 
   list(
      tmax=tmax,
      sb=sb, delta=delta, sd=sd, 
      N0=N0, mu=mu,
      resc0=resc0, simnum=simnum
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
   alpha <- 2*N*mu
   beta <- 2*N*mu     # for now let mutation rates be equal
   gamma <- -2*N*sd   # make effect of sd negative 
   b <- beta(alpha,beta)  # beta function 
   f1 <- hyperg_1F1(beta,alpha+beta,gamma)  # confluent hypergeom fn from gsl
   phiC <- exp(gamma*(x)) * x^(alpha-1) * (1-x)^(beta-1)
   phiC/(b*f1)  # normalise to return pdf 
}
## Discretised Wright dist, integrating in binned freqs;
##   bins are (0,0.5/N), (0.5/N, 1.5/N), etc 
tabulate.wright.hist <- function(N, sd, mu){
   ##  N+1 bins centred on 1/N, ..., 1-1/N but the
   ##  two outer bins integrating over 0 to 1/(2N) and over 1-1/(2N) to 1
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
   ##   print(min(dist$prob))
   ## Switch to integration if Kimura method causes a problem: 
   if (min(dist$prob)<0)  # numerical problem so try straight integration
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
## Plot mean frequencies of deleterious allele under Wright model and
## under deterministic model
plot.mean.freq <- function(){
   require("sfsmisc") 
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
   
   pdf("mean-freq-Wright-dist.pdf",height=3, width=3
       , title="Mean frequency")
   par(mar=c(4,4,1,1), cex=0.7)
   plot(N.arr, mf1, t="l", log="xy", ylim=c(0.00005,2)
      , xaxt="n", yaxt="n"
      , xlab="Population size N"
      , ylab="Expected frequency under Wright's distribution")
   eaxis(1, at.small=F)
   axt1 <- c(0.1, 0.01, 0.001, 0.0001) 
   axt2 <- c(0.5, 0.1, 0.01, 0.001, 0.0001) 
   eaxis(2, at.small=F, at=axt2, labels=c(0.5, pretty10exp(axt1,drop.1=T)))
   lines(N.arr, mf3, col=1, lty=2)
   lines(N.arr, mf2, col=2)
   lines(N.arr, mf4, col=2, lty=2)
   legend("topright", L
        , col=c(1,1,2,2), lty=c(1,2,1,2), cex=0.8
        , bty="n")
   c2 <- rgb(red=255/255,green=0,blue=0, alpha=0.3) # transparent red
   c1 <- rgb(red=0,green=0,blue=0, alpha=0.3)       # transparent black
   abline(h=c(h1,h3,h2,h4), lty=c(1,2,1,2),col=c(c1,c1,c2,c2))
   dev.off()
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
plot.wright.and.rescue.present <- function(){
   require("sfsmisc") # for nice axis in plot
   p <- setp()
   sd <- p$sd
   mu <- p$mu
   pdf("wright-allele-present.pdf",width=6.7,height=2.5)
   par(mfrow=c(1,3), cex=0.7, mar=c(4,4,1.2,1) )
   ## Plot Wright PDF near left boundary for different N values
   from <- 0.0001; to <- 0.05
   plot(1,1, type="n", xlim=c(0,to), ylim=c(10^-6,10^2) #, log="y"
      , xlab="(Low) frequency of allele", ylab="Probability density")
   title("A", adj=0)
   curve.wright.pdf(from,to,1000, sd, mu, col="darkorange",lty=5)
   curve.wright.pdf(from,to,10000, sd, mu, col="forestgreen")
   curve.wright.pdf(from,to,100000, sd, mu, col="purple", lty=6)
   legend("topright", c("N = 1000", "N = 10,000", "N = 100,000"),
          col=c("darkorange","forestgreen","purple"), lty=c(5,1,6)
        , bty="n", cex=0.8)
   ## Now plot near right boundary 
   from <- 0.95; to <- 0.99999
   plot(1,1, type="n", xlim=c(from,1), ylim=c(0, 5) #, log="y"
      , xlab="(High) frequency of allele", ylab="Probability density")
   curve.wright.pdf(from,to,1000, sd, mu, col="darkorange",lty=5)
   curve.wright.pdf(from,to,10000, sd, mu, col="forestgreen")
   curve.wright.pdf(from,to,100000, sd, mu, col="purple",lty=6)
   title("B", adj=0)
   ## Now plot Prob(allele present) vs N
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
          col=c(4,1,2), lty=c(4,1,2), bty="n", cex=0.8)
   title("C", adj=0)

   dev.off()
}


##==========================================================
##   Theoretical probability of rescue:
##       Orr&Unckless + Wright
##==========================================================

## Numerically evaluated exact fixation probability
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
prob.TW <- function(N,sb,delta,mu, sd, pfixmod=3){
   ## pfixmod is the model for prob(fixation) 
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
         lambda.dn <- gamma - N[i]*pfix*mu/delta
         f1.top.st <- hyperg_1F1(alpha,2*alpha,-lambda.st)
         f1.top.dn <- hyperg_1F1(alpha,2*alpha,-lambda.dn)
         f1.bot <- hyperg_1F1(alpha,2*alpha,-gamma)
         pst[i] <- 1 - f1.top.st/f1.bot
         pdn[i] <- 1 - exp(-pfix*N[i]*mu/delta)*f1.top.dn/f1.bot
   }
   list(pst=pst, pdn=pdn)
}


##==========================================================
##    --- RESCUE SIMULATIONS ---
##==========================================================

## Run one sim using parameters in p 
onesim <- function(p){
   stopsize <- 2000  # at this size, deem escape successful
   with(p, {
      tcap <- tmax*10     # run for long time unless extinction or escape
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

## run sim once with arbitrary resc0, and plot result [not used]
plot.sim <- function(){
   ## To run denovo only, let resc0=0
   out <- onesim(setp(resc0=90, sb=0.05)) 
   pdf("sim-env-change.pdf", width=3.5,height=3)
   par(mar=c(4,4,1,1), cex=0.5)
   t <- c(1:length(out$wt))
   plot(t,out$resc, t="l", col=2, ylim=c(0,max(c(out$resc,out$wt))), 
        ylab="Number of wildtype and mutants", xlab="Time (generations)")
   lines(t, out$wt) 
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
sim.orr.wright <- function (p){
   wfs <- rwright(p$simnum, p$N0, p$sd, p$mu) # draw random Wright freqs
   wins <- 0  # for counting successful rescue 
   for (i in 1:p$simnum) {
      p$resc0 <- round(wfs[i]*p$N0) # initial number of rescue allele copies
      out <- onesim(p)
      if (out$win==1)
         wins <- wins+1
   }
   r <- wins/p$simnum            # proportion rescued
   se <- sqrt(r*(1-r)/p$simnum)  # approx SE
   list(r=r, se=se)    
}

## Simulate for many N values and store in external file
##  whose name specifies which param has changed
sim.orr.wright.over.N <- function(p, ext="mu1e-4"){
   N <-  10*2^(seq(1:14)-1)

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
##  - doing all takes a long time
many.sims <- function(){

   sim.orr.wright.over.N(setp(sb=0.4),ext="sb4")  # 
#   sim.orr.wright.over.N(setp(sb=0.05),ext="sb05")  # default 
   sim.orr.wright.over.N(setp(sb=0.02),ext="sb02") #
   sim.orr.wright.over.N(setp(sb=0.011),ext="sb011") #
   
   sim.orr.wright.over.N(setp(mu=1e-4),ext="mu-1e-4") #
   sim.orr.wright.over.N(setp(mu=1e-6),ext="mu-1e-6")  #  

   sim.orr.wright.over.N(setp(sd=1e-1),ext="sd-1e-1")  #
   sim.orr.wright.over.N(setp(sd=1e-2),ext="sd-1e-2")  #
   sim.orr.wright.over.N(setp(sd=1e-5),ext="sd-1e-5")
}


##==========================================================
##           PLOT SIM and THEORY 
##==========================================================

## Plot rescue probabilities -- vary parameters
## --- Multi-panel figure --- 
plot.vary.params <- function(){
   setup.panel <- function(lo=1,hi=5,panel.lab){
      require("sfsmisc") 
      plot(1,1, pty="n", xaxt="n", log="x", xlim=c(10^lo,10^hi), ylim=c(0,1)
         , xlab="Population size N", ylab="") 
      eaxis(1, f.smalltcl=0)     # nice x axis
      title(panel.lab, adj=0)
      title(ylab="Probability of rescue", line=2.1) 
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
   pdf("prob-rescue-vary-params.pdf",height=2.6, width=7)
   par(mfrow=c(1,3), mar=c(4,3.1, 1.2, 1), cex=0.75)

   ## Vary sb in a panel 
   setup.panel(panel.lab="A")
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
    , col=c(4,1,3,2),lty=1, cex=0.87, bty="n")
   ## Vary mu in a panel 
   setup.panel(panel.lab="B")
   one.series(setp(mu=1e-4),4)   
   one.series(setp())
   one.series(setp(mu=1e-6),2)
   add.pts("dat-orr-wright-sim-mu-1e-4",4)
   add.pts("dat-orr-wright-sim-sb05",1)  # the default again 
   add.pts("dat-orr-wright-sim-mu-1e-6",2)
   legend("topleft",c(expression(mu*" = 10"^-4),
                      expression(mu*" = 10"^-5),
                      expression(mu*" = 10"^-6))
         , col=c(4,1,2),lty=1, cex=0.87, bty="n")
   ## Vary sd in a panel 
   setup.panel(panel.lab="C")
   one.series(setp(sd=1e-5),4)
   one.series(setp())
   one.series(setp(sd=1e-2),3)
   one.series(setp(sd=0.1),2)
   add.pts("dat-orr-wright-sim-sd-1e-1",2)
   add.pts("dat-orr-wright-sim-sd-1e-2",3)
   add.pts("dat-orr-wright-sim-sb05",1)  # the default again 
   add.pts("dat-orr-wright-sim-sd-1e-5",4)
   legend("topleft",c(
      expression("s"[d]*" = 10"^-5),
      expression("s"[d]*" = 10"^-3),
      expression("s"[d]*" = 10"^-2),
      expression("s"[d]*" = 0.1")
   ),
      col=c(4,1,3,2),lty=1, cex=0.87, bty="n")
   dev.off()
}

## Plot the components of the rescue probability 
plot.components <- function(lo=1, hi=5){
   require("sfsmisc")
   p <- setp()
   N_ <- 10^seq(lo,hi, length.out=100)
   pdf("p-rescue-components.pdf", height=2.7,width=3.4
     , title="Rescue components")
   par(mar=c(4,4,1,1),cex=0.65)
   plot(1,1, pty="n", xaxt="n", log="x", xlim=c(10^lo,10^hi), ylim=c(0,1)
      , xlab="Population size N", ylab="Probability of rescue")
   eaxis(1, at.small=F) 
   cola <- c("black","dodgerblue2","orange","magenta")

   ## Newer version with marginal probs which don't add to 1
   pno <- prob.new.old(N_, p$sb,p$delta,p$mu,p$sd)
   ptw <- prob.TW(N_, p$sb,p$delta,p$mu,p$sd)
   lines(N_, ptw, col=cola[1] )
   lines(N_, pno$pst, col=cola[2], lty=4)
   lines(N_, pno$pdn, col=cola[3], lty=2)
   abline(v=c(1/(1*p$mu), 1/p$sd), lty=3, col=c(6,6))
   legend("topleft", c("Total","standing variation","de novo")
        , lty=c(1,4,2), col=cola
        , bg="white", cex=0.9)
   ## Add special text to indicate regimes 
   text(12,0.2,"weak\n rescue", cex=0.8)
   text(50,0.54,"strong\n drift", cex=0.8)
   text(3*10^3,0.08,"strong\n selection", cex=0.8)
   text(5*10^4,0.76,"strong\n mutation", cex=0.8)
   
   dev.off()
}



##==========================================================
##    (Locally) worst population size 
##==========================================================

## Find a min in vector starting from the end of the vector
min.from.right <- function(x,y){
   l <- length(y)
#   argmin <- x[1] # send to left if no local min
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
## resolution: try 10000 (it is slow)
get.min.N.slow <- function(p, res=10000, say=FALSE, interv=c(10,10^5)){
   N <- seq(interv[1],interv[2],length.out=res) # linear
   ##N <- 10^seq(1,5, length.out=res)  # exponent 
   ptw <- prob.TW(N, p$sb,p$delta,p$mu,p$sd)
   mfr <- min.from.right(N,ptw)
   if (say) { # if you need to know 
      catn(mfr[1], mfr[2])
   }
   mfr[1]
}
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
      , xlab=expression("Mutation rate, "*mu))
   eaxis(1, n.axp=1, at.small=F) 
   title("A", adj=0)
   lines(wN011$x, wN011$worstN, col=2)
   lines(wN02$x, wN02$worstN, col=4)
   lines(wN05$x, wN05$worstN, col=1)
   legend("topright", c(expression(s[b]*" = 0.011"),
                        expression(s[b]*" = 0.02"),
                        expression(s[b]*" = 0.05")), #pch=16,
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
      , xlab=expression("Deleterious effect, "*s[d]))
   eaxis(1, n.axp=1, at.small=F) 
   title("B", adj=0)
   lines(wN011$x, wN011$worstN, col=2)
   lines(wN02$x, wN02$worstN, col=4)
   lines(wN05$x, wN05$worstN, col=1)
   sd.seq <- 10^seq(-4,-1, length.out=100)
   lines(sd.seq, 1/sd.seq, lty=2, col=3) # plot 1/sd with dashed line
   legend("topright", c(expression(s[b]*" = 0.011"),
                        expression(s[b]*" = 0.02"),
                        expression(s[b]*" = 0.05")
                       ,expression(" 1/"*s[d])
                        ), 
          col=c(2,4,1,3), lty=c(1,1,1,2), bty="n") 

   dev.off()
}

