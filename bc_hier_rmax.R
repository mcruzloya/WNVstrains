setwd('/Users/cruzloya/git/WNVstrains')

set.seed(42) # Set seed for reproducibility.

library('R2jags')
library('mcmcplots')
library('MCMCvis')
library('scales')

quad_model_rmax <- function(T, Tmin, Tmax, rmax) {
  c <- 4 * rmax / (Tmax - Tmin)^2
  result <- (T > Tmin) * (T < Tmax) * c * (T - Tmin) * (Tmax - T)
  return(pmin(result, 1))
}

data <- read.csv('retrodata_max.csv')
data

historic <- subset(data, data$group == 0)
recent <- subset(data, data$group == 1)

historic

sink("bc_quad_hier_rmax.txt")
cat("
    model{
    
    ## Priors
    rmax_mu <- 0.35
    rmax_nprior <- 4.0
    rmax_m ~ dbeta(rmax_nprior * rmax_mu, rmax_nprior * (1 - rmax_mu))
    
    inv_sqrt_n_rmax ~ dnorm(0, 1^-2)T(0,)
    n_rmax <- 1.0 / inv_sqrt_n_rmax^2
    
    Tmin_m ~ dnorm(16.8, 1/2^2) # 95% likely in interval
    s_Tmin ~ dnorm(0, 2^-2)T(0,)
    
    Tmax_m ~ dnorm(38.9, 1/3^2 ) # 95% likely in interval 
    s_Tmax ~ dnorm(0, 2^-2)T(0,)
    
    c <- 4 * rmax / (Tmax - Tmin)^2
    
    for(st in 1:n.strains) {
      Tmin[st] ~ dnorm(Tmin_m, 1 / s_Tmin^2)
      Tmax[st] ~ dnorm(Tmax_m, 1 / s_Tmax^2)
      rmax[st] ~ dbeta(n_rmax * rmax_m, n_rmax * (1 - rmax_m))
      ## Derived Quantities and Predictions
      Topt[st] <- (Tmin[st] + Tmax[st]) / 2
    }
    
    ## Likelihood
    for(i in 1:N.obs){
      p[i] <- max(min(c[strain[i]] * (temp[i] - Tmin[strain[i]]) * (Tmax[strain[i]] - temp[i]) * (temp[i] < Tmax[strain[i]]) * (Tmin[strain[i]] < temp[i]),
      1), 10^-6)
      n[i] ~ dbin(p[i], N[i])
    }
    
    
    } # close model
    ",fill=TRUE)
sink()

##### Calculate initial values for MCMC.
## We are picking random values so that every chain will start at a different place. 
inits<-function(){list(
  Tmin_m = runif(1, 10, 15),
  Tmax_m = runif(1, 45, 55),
  rmax_m = 0.3,
  Tmin=runif(4, 12, 17),
  Tmax=runif(4, 35, 45),
  rmax = rep(0.3, 4))}

##### Parameters to Estimate
parameters <- c("Tmin_m", "Tmax_m", "rmax_m",
                "Tmin", "Tmax", "Topt", "rmax",
                "s_Tmin", "s_Tmax", "inv_sqrt_n_rmax")

##### MCMC Settings
# Number of posterior dist elements = [(ni - nb) / nt ] * nc = [ (25000 - 5000) / 8 ] * 3 = 7500
ni <- 200000 # number of iterations in each chain
nb <- 40000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 4 # number of chains


##### Organize Data for JAGS
N <- historic$total
n <- historic$disseminated
temp <- historic$temperature
strain <- historic$strain_id
N.obs <- length(N)
n.strains <- max(strain)

##### Bundle Data
jag.data<-list(N=N, n=n, temp = temp, N.obs=N.obs, strain=strain, 
               n.strains=n.strains)
print(jag.data)
##### Run JAGS
historic.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                 model.file="bc_quad_hier_rmax.txt", n.thin=nt, n.chains=nc, 
                 n.burnin=nb, n.iter=ni, DIC=TRUE, working.directory=getwd())
historic.out
mcmcplot(historic.out)

par(mfrow=c(2, 2))
temps <- seq(10, 55, 0.1)
stidx <- 1
col='black'
alpha=0.3
for(st in 1:4) {
  df <- subset(historic, (historic$strain_id == st))
  #print(df)
  p <- df$disseminated / df$total
  se <- sqrt(p * (1 - p) / df$total)
  stname <- df$strain[1]
  plot(df$temp, p, pch=20, xlim=c(10, 55), ylim=c(0,1),
       main=stname, xlab="temperature", ylab='bc (disseminated/exposed)', cex=1.5)
  arrows(df$temp, p, y1=p+se, col='black',
         angle=90, length=0.03, lwd=2)
  arrows(df$temp, p, y1=p-se, col='black',
         angle=90, length=0.03, lwd=2)
  
  chains <- MCMCchains(historic.out, params=c(paste("Tmin[", as.character(st), "]", sep=""), 
                                              paste("Tmax[", as.character(st), "]", sep=""),
                                              paste("rmax[", as.character(st), "]", sep="")),
                       ISB=FALSE)
  curves <- apply(chains, 1, 
                  function(x) quad_model_rmax(temps, x[1], x[2], x[3]))
  
  # Find mean curve and credible intervals.
  meancurve <- apply(curves, 1, mean)
  CI <- apply(curves, 1, quantile, c(0.025, 0.975))
  
  lines(temps, meancurve, col=col, linewidth=3)
  polygon(c(temps, rev(temps)), c(CI[1,], rev(CI[2,])), 
          col=alpha(col, alpha), lty=0)
  stidx <- stidx + 1
}

recent

##### Organize Data for JAGS
N <- recent$total
n <- recent$disseminated
temp <- recent$temperature
strain <- recent$strain_id
N.obs <- length(N)
n.strains <- max(strain)

##### Bundle Data
jag.data<-list(N=N, n=n, temp = temp, N.obs=N.obs, strain=strain, 
               n.strains=n.strains)
print(jag.data)
##### Run JAGS
recent.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                     model.file="bc_quad_hier_rmax.txt", n.thin=nt, n.chains=nc, 
                     n.burnin=nb, n.iter=ni, DIC=TRUE, working.directory=getwd())
recent.out
mcmcplot(recent.out)


par(mfrow=c(2, 2))
temps <- seq(10, 55, 0.1)
stidx <- 1
col='black'
alpha=0.3
for(st in 1:4) {
  df <- subset(recent, (recent$strain_id == st))
  #print(df)
  p <- df$disseminated / df$total
  se <- sqrt(p * (1 - p) / df$total)
  stname <- df$strain[1]
  plot(df$temp, p, pch=20, xlim=c(10, 55), ylim=c(0,1),
       main=stname, xlab="temperature", ylab='bc (disseminated/exposed)', cex=1.5)
  arrows(df$temp, p, y1=p+se, col='black',
         angle=90, length=0.03, lwd=2)
  arrows(df$temp, p, y1=p-se, col='black',
         angle=90, length=0.03, lwd=2)
  
  chains <- MCMCchains(recent.out, params=c(paste("Tmin[", as.character(st), "]", sep=""), 
                                              paste("Tmax[", as.character(st), "]", sep=""),
                                              paste("rmax[", as.character(st), "]", sep="")),
                       ISB=FALSE)
  curves <- apply(chains, 1, 
                  function(x) quad_model_rmax(temps, x[1], x[2], x[3]))
  
  # Find mean curve and credible intervals.
  meancurve <- apply(curves, 1, mean)
  CI <- apply(curves, 1, quantile, c(0.025, 0.975))
  
  lines(temps, meancurve, col=col, linewidth=3)
  polygon(c(temps, rev(temps)), c(CI[1,], rev(CI[2,])), 
          col=alpha(col, alpha), lty=0)

  stidx <- stidx + 1
}

recent.out
historic.out

saveRDS(historic.out, "bc_hier_historic.RDS")
saveRDS(recent.out, "bc_hier_recent.RDS")


