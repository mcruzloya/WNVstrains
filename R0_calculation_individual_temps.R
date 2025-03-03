setwd('/Users/cruzloya/git/WNVstrains')

# R0 calculation
library('R2jags')
library('MCMCvis')
library('coda')
library('scales')

# Relative R0 parametrization as in Eq. (A1) in Shocket et al. 
# (removing N and r).
R0 <- function(a, bc, lf, PDR, EFGC, EV, pLA, MDR) {
  eps <- 1e-10 # Small constant to prevent denominators being zero 
  # Calculate as log of R0^2 for numerical stability.
  logR02 <- (3 * log(a) + log(bc) - (1 / (PDR * lf + eps)) + log(EFGC)
             + log(EV) + log(pLA) + log(MDR) + 3 * log(lf))
  return(sqrt(exp(logR02)))
}

# Briere TPC model.
briere <- function(T, Tmin, Tmax, c) {
  result <- c * T * (T - Tmin) * sqrt((Tmax - T) * (T > Tmin) * (T < Tmax))
  return(result)
}

# Quadratic TPC model.
quad <- function(T, Tmin, Tmax, c) {
  result <- (T > Tmin) * (T < Tmax) * c * (T - Tmin) * (Tmax - T)
  return(result)
}

# Quadratic TPC model with upper limit at 1 for proportions.
quad_lim <- function(T, Tmin, Tmax, c) {
  result <- (T > Tmin) * (T < Tmax) * c * (T - Tmin) * (Tmax - T)
  return(pmin(result, 1))
}

quad_model_rmax <- function(T, Tmin, Tmax, rmax) {
  c <- 4 * rmax / (Tmax - Tmin)^2
  result <- (T > Tmin) * (T < Tmax) * c * (T - Tmin) * (Tmax - T)
  return(pmin(result, 1))
}

# Linear TPC model with upper limit at 15C.
linear_lim <- function(T, m, b) {
  result <- (-m * T + b) * (T < b/m)
  result[T < 15] <- (-m * 15 + b) * (T < b/m) # Set maximum value at 15 degrees
  return(result)
}

plot_curves <- function(temps, curves, xlim=c(0, 50), ylim=c(0, 1)) {
  meancurve <- apply(curves, 1, mean)
  CI <- apply(curves, 1, quantile, c(0.025, 0.975))
  plot(temps, meancurve, pch=20, col="black", xlim=xlim, ylim=ylim)
  points(temps, CI[1,], col="gray")
  points(temps, CI[2,], col="gray")
}

## 1. Load chains with posterior estimates of parameters and calculate trait
## values.

temps <- c(20, 24, 28) 

## Load data for biting rate
#load("./shocket_fits/jagsout_a_Cpip_inf.Rdata")
#a.Cpip.out.inf$BUGSoutput

# Only keep chains for TPC model parameters to save memory.
#a.chains <- MCMCchains(a.Cpip.out.inf, params=c("cf.T0", "cf.Tm", "cf.q"))
a.chains <- readRDS("./shocket_fits/a_Cpip_shocket.RDS")

# Remove data to save memory
#rm(a.Cpip.out.inf)

# Calculate trait values for each sample at all temperatures.
a.curves <-  apply(a.chains, 1, function(x) briere(temps, x[1], x[2], x[3]))
plot_curves(temps, a.curves, ylim=c(0, 0.5))

## Vector competence 

bc.hist.chains <- readRDS("bc_hier_historic.RDS")
bc.hist.chains
#View(bc.hist.chains$BUGSoutput$summary)
bc.hist.curves <- list()
for(st in 1:4) {
  chains <- MCMCchains(bc.hist.chains, params=c(paste("Tmin[", as.character(st), "]", sep=""), 
                                                paste("Tmax[", as.character(st), "]", sep=""),
                                                paste("rmax[", as.character(st), "]", sep="")),
                       ISB=FALSE)
  bc.hist.curves[[st]] <- apply(chains, 1, 
                                function(x) quad_model_rmax(temps, x[1], x[2], x[3]))
}
str(bc.hist.curves)
plot_curves(temps, bc.hist.curves[[2]], ylim=c(0,1))
bc.rec.chains <- readRDS("bc_hier_recent.RDS")
#View(bc.rec.chains$BUGSoutput$summary)

bc.rec.chains
bc.rec.curves <- list()
for(st in 1:4) {
  chains <- MCMCchains(bc.rec.chains, params=c(paste("Tmin[", as.character(st), "]", sep=""), 
                                               paste("Tmax[", as.character(st), "]", sep=""),
                                               paste("rmax[", as.character(st), "]", sep="")),
                       ISB=FALSE)
  bc.rec.curves[[st]] <- apply(chains, 1, 
                               function(x) quad_model_rmax(temps, x[1], x[2], x[3]))
}
plot_curves(temps, bc.rec.curves[[1]], ylim=c(0,1))

## Lifespan
#load("./shocket_fits/jagsout_lf_Cpip_inf.Rdata")
#lf.Cpip.out.inf

#lf.chains <- MCMCchains(lf.Cpip.out.inf, params=c("cf.m", "cf.b"))
#rm(lf.Cpip.out.inf)
lf.chains <- readRDS("./shocket_fits/lf_Cpip_shocket.RDS")

lf.curves <- apply(lf.chains, 1, function(x) linear_lim(temps, x[1], x[2]))

plot_curves(temps, lf.curves, ylim=c(0, 100))

## Pathogen development rate
#load("./shocket_fits/jagsout_PDR_CpipWNV_inf.Rdata")
#PDR.CpipWNV.out.inf$BUGSoutput$summary

#PDR.chains <- MCMCchains(PDR.CpipWNV.out.inf, params=c("cf.T0", "cf.Tm", "cf.q"))
#rm(PDR.CpipWNV.out.inf)
PDR.chains <- readRDS("./shocket_fits/PDR_CpipWNV_shocket.RDS")
PDR.curves <- apply(PDR.chains, 1, function(x) briere(temps, x[1], x[2], x[3]))

plot_curves(temps, PDR.curves, ylim=c(0, 0.3))

## Fecundity
#load("./shocket_fits/jagsout_EFOC_Cpip_inf.Rdata")
#EFOC.Cpip.out.inf
#EFOC.chains <- MCMCchains(EFOC.Cpip.out.inf, params=c("cf.T0", "cf.Tm", "cf.q"))
EFOC.chains <- readRDS("./shocket_fits/EFOC_Cpip_shocket.RDS")
EFOC.curves <- apply(EFOC.chains, 1, function(x) quad(temps, x[1], x[2], x[3]))
plot_curves(temps, EFOC.curves, ylim=c(0, 400))


## Egg viability 
#load("./shocket_fits/jagsout_EV_Cpip_inf.Rdata")
#EV.chains <- MCMCchains(EV.Cpip.out.inf, params=c("cf.T0", "cf.Tm", "cf.q"))
#EV.Cpip.out.inf$BUGSoutput$summary
#rm(EV.Cpip.out.inf)
EV.chains <- readRDS("./shocket_fits/EV_Cpip_shocket.RDS")
EV.curves <- apply(EV.chains, 1, function(x) quad_lim(temps, x[1], x[2], x[3]))
plot_curves(temps, EV.curves, ylim=c(0, 1))


## Larval survival
#load("./shocket_fits/jagsout_pLA_Cpip_inf.Rdata")
#pLA.Cpip.out.inf$BUGSoutput$summary
#pLA.chains <- MCMCchains(pLA.Cpip.out.inf, params=c("cf.T0", "cf.Tm", "cf.q"))
#rm(pLA.Cpip.out.inf)
pLA.chains <- readRDS("./shocket_fits/pLA_Cpip_shocket.RDS")

pLA.curves <- apply(pLA.chains, 1, function(x) quad_lim(temps, x[1], x[2], x[3]))
plot_curves(temps, pLA.curves, ylim=c(0, 1))

## Mosquito development rate
#load("./shocket_fits/jagsout_MDR_Cpip_inf.Rdata")
#MDR.Cpip.out.inf$BUGSoutput$summary
#MDR.chains <- MCMCchains(MDR.Cpip.out.inf, params=c("cf.T0", "cf.Tm", "cf.q"))

MDR.chains <- readRDS("./shocket_fits/MDR_Cpip_shocket.RDS")
MDR.curves <- apply(MDR.chains, 1, function(x) briere(temps, x[1], x[2], x[3]))
plot_curves(temps, MDR.curves, ylim=c(0, 0.3))

## R0
find_Topt <- function(temps, curves) {
  idx <- apply(curves, 2, which.max)
  return(temps[idx])
}
find_Tmin <- function(temps, curves) {
  idx <- apply(curves, 2, function(x) which(x > 0)[1] - 1)
  return(temps[idx])
}

find_Tmax <- function(temps, curves) {
  idx <- apply(curves, 2, function(x) tail(which(x > 0), 1) + 1)
  return(temps[idx])
}

R0.summary <- data.frame(group = c(rep('historic', 4), rep('recent', 4)),
                         strain= rep(1:4, 2),
                         R0_20=rep(0, 8), R0_20l=rep(0,8), R0_20h=rep(0,8),
                         R0_24=rep(0, 8), R0_24l=rep(0,8), R0_24h=rep(0,8),
                         R0_28=rep(0, 8), R0_28l=rep(0,8), R0_28h=rep(0,8))
R0.summary

R0.hist.curves <- list()

str(a.curves)
str(t(bc.curves.hist[[1]]))

str(PDR.curves)

# Resample indices from PDR to simulate another strain (with same curve).
#(a, bc, lf, PDR, EFGC, EV, pLA, MDR)
for(st in 1:4) {
  cond <- (R0.summary$group == "historic") & (R0.summary$strain == st) 
  R0.hist.curves[[st]] <- R0(a.curves, t(bc.curves.hist[[st]]), lf.curves, 
                             PDR.curves[, sample(1:80000, 80000)], EFOC.curves,
                             EV.curves, pLA.curves, MDR.curves)
  
  R0.summary$R0_20[cond] <- mean(R0.hist.curves[[st]][1, ])
  R0.summary$R0_20l[cond] <- quantile(R0.hist.curves[[st]][1, ], 0.025)
  R0.summary$R0_20h[cond] <- quantile(R0.hist.curves[[st]][1, ], 0.975)
  
  R0.summary$R0_24[cond] <- mean(R0.hist.curves[[st]][2, ])
  R0.summary$R0_24l[cond] <- quantile(R0.hist.curves[[st]][2, ], 0.025)
  R0.summary$R0_24h[cond] <- quantile(R0.hist.curves[[st]][2, ], 0.975)
  
  R0.summary$R0_28[cond] <- mean(R0.hist.curves[[st]][3, ])
  R0.summary$R0_28l[cond] <- quantile(R0.hist.curves[[st]][3, ], 0.025)
  R0.summary$R0_28h[cond] <- quantile(R0.hist.curves[[st]][3, ], 0.975)
  
}

R0.summary

R0.hist.curves[[1]][1,1:100]

as.matrix(t(bc.curves.recent[[1]][1:100, ]))

R0.rec.curves <- list()
for(st in 1:4) {
  cond <- (R0.summary$group == "recent") & (R0.summary$strain == st) 
  R0.rec.curves[[st]] <- R0(a.curves, t(bc.curves.recent[[st]]), lf.curves, 
                            PDR.curves[, sample(1:80000, 80000)], EFOC.curves, EV.curves, 
                             pLA.curves, MDR.curves)
  R0.summary$R0_20[cond] <- mean(R0.rec.curves[[st]][1, ])
  R0.summary$R0_20l[cond] <- quantile(R0.rec.curves[[st]][1, ], 0.025)
  R0.summary$R0_20h[cond] <- quantile(R0.rec.curves[[st]][1, ], 0.975)
  
  R0.summary$R0_24[cond] <- mean(R0.rec.curves[[st]][2, ])
  R0.summary$R0_24l[cond] <- quantile(R0.rec.curves[[st]][2, ], 0.025)
  R0.summary$R0_24h[cond] <- quantile(R0.rec.curves[[st]][2, ], 0.975)
  
  R0.summary$R0_28[cond] <- mean(R0.rec.curves[[st]][3, ])
  R0.summary$R0_28l[cond] <- quantile(R0.rec.curves[[st]][3, ], 0.025)
  R0.summary$R0_28h[cond] <- quantile(R0.rec.curves[[st]][3, ], 0.975)
}

R0.hist.curves[[1]]
View(R0.summary)

R0.summary
R0.summary.quot

R0.mean.rec <- (R0.rec.curves[[1]] + R0.rec.curves[[2]] + R0.rec.curves[[3]] + R0.rec.curves[[4]]) / 4
R0.mean.hist <- (R0.hist.curves[[1]] + R0.hist.curves[[2]] + R0.hist.curves[[3]] + R0.hist.curves[[4]]) / 4

R0.mean.diff <- R0.mean.rec - R0.mean.hist
R0.mean.summary <-data.frame(temperature=c(20, 24, 28), R0_diff=apply(R0.mean.diff, 1, mean),
                             R0_diff_CI_lo=apply(R0.mean.diff, 1, quantile, 0.025),
                             R0_diff_CI_hi=apply(R0.mean.diff, 1, quantile, 0.975))


R0.mean.quot <- R0.mean.rec / (R0.mean.hist + 1e-10)

R0.mean.summary <-data.frame(temperature=c(20, 24, 28), R0_quot=apply(R0.mean.quot, 1, mean),
                             R0_quot_CI_lo=apply(R0.mean.quot, 1, quantile, 0.025),
                             R0_quot_CI_hi=apply(R0.mean.quot, 1, quantile, 0.975),
                             prob_rec_R0_higher=apply(R0.mean.quot > 1, 1, mean))


R0.mean.summary.hist <-data.frame(temperature=c(20, 24, 28), R0=apply(R0.mean.hist, 1, mean),
                             R0_CI_lo=apply(R0.mean.hist, 1, quantile, 0.025),
                             R0_CI_hi=apply(R0.mean.hist, 1, quantile, 0.975))
R0.mean.summary.rec <-data.frame(temperature=c(20, 24, 28), R0=apply(R0.mean.rec, 1, mean),
                                  R0_CI_lo=apply(R0.mean.rec, 1, quantile, 0.025),
                                  R0_CI_hi=apply(R0.mean.rec, 1, quantile, 0.975))

R0.mean.summary.hist

strain_names

# Comparison plot.
strain_names <- c("2003.1", "2003.2", "2004.1", "2004.2",
                  "2017.1", "2017.2", "2018.1", "2018.2")

## Old labels:
## c("DS03", "US03", "DS04", "US04", 
## "DS17", "US17", "DS18", "US18")
strain_plot_order_hist <- c(1,2,3,4)
strain_plot_order_rec <- c(2,1,4,3)

R0.summary
#par(mfrow=c(1, 3), mar = c(5.1, 4.1, 4.1, 2.1), mgp=c(3,1,0))
par(mfrow=c(1, 3), mar = c(7.1, 4.1, 4.1, 2.1), mgp=c(5, 1, 0))
for(temp in c(20, 24, 28)) {
  plot("", ylim=c(0, 130), xlab="strain", xlim=c(0.5, 8.5),
           ylab="relative R0", main=paste(temp, '°C', sep=''),
           xaxt='n', yaxt='n', col="steelblue")
  axis(side=1, at=1:8, labels=strain_names, las=2)
  axis(side=2, at=c(0, 20, 40, 60 ,80, 100, 120), las=2)
  
  for(st in 1:4) {
    hist.strain <- strain_plot_order_hist[st]
    rec.strain <- strain_plot_order_rec[st]
    points(st, R0.summary[(R0.summary$group == 'historic') & 
                        (R0.summary$strain == hist.strain), 
                        paste('R0', temp, sep='_')], pch=20, col="steelblue")
    arrows(st, R0.summary[(R0.summary$group == 'historic') & 
                            (R0.summary$strain == hist.strain), 
                          paste('R0', temp, sep='_')],
            y1=R0.summary[(R0.summary$group == 'historic') & 
                            (R0.summary$strain == hist.strain), 
                          paste(paste('R0', temp, sep='_'), 'l', sep='')],
           angle=90, length=0.1, lwd=1, col="steelblue")
    arrows(st, R0.summary[(R0.summary$group == 'historic') & 
                            (R0.summary$strain == hist.strain), 
                          paste('R0', temp, sep='_')],
           y1=R0.summary[(R0.summary$group == 'historic') & 
                           (R0.summary$strain == hist.strain), 
                         paste(paste('R0', temp, sep='_'), 'h', sep='')],
           angle=90, length=0.1, lwd=1, col="steelblue")
    
    points(st + 4, R0.summary[(R0.summary$group == 'recent') & 
                              (R0.summary$strain == rec.strain), 
                              paste('R0', temp, sep='_')], pch=20, col="firebrick")
    arrows(st + 4, R0.summary[(R0.summary$group == 'recent') & 
                            (R0.summary$strain == rec.strain), 
                          paste('R0', temp, sep='_')],
           y1=R0.summary[(R0.summary$group == 'recent') & 
                           (R0.summary$strain == rec.strain), 
                         paste(paste('R0', temp, sep='_'), 'l', sep='')],
           angle=90, length=0.1, lwd=1, col="firebrick")
    arrows(st + 4, R0.summary[(R0.summary$group == 'recent') & 
                            (R0.summary$strain == rec.strain), 
                          paste('R0', temp, sep='_')],
           y1=R0.summary[(R0.summary$group == 'recent') & 
                           (R0.summary$strain == rec.strain), 
                         paste(paste('R0', temp, sep='_'), 'h', sep='')],
           angle=90, length=0.1, lwd=1, col="firebrick")
    #arrows(st + 4, R0.mean.summary.rec[, 'R0'], y1=R0.mean.summary.rec[, 'R0_CI_hi'], 
    #       angle=90, length=0.1, lwd=1, col="firebrick")
    #arrows(st + 4, R0.mean.summary.rec[, 'R0'], y1=R0.mean.summary.rec[, 'R0_CI_lo'], 
    #       angle=90, length=0.1, lwd=1, col="firebrick")
  }
  
}

par(mfrow=c(1, 3), mar = c(7.1, 4.1, 4.1, 2.1), mgp=c(5, 1, 0))
for(temp in c(20, 24, 28)) {
  plot("", ylim=c(20, 130), xlab="strain", xlim=c(0.5, 8.5),
       ylab="relative R0", main=paste(temp, '°C', sep=''),
       xaxt='n', yaxt='n', col="steelblue")
  axis(side=1, at=1:8, labels=strain_names, las=2)
  for(st in 1:4) {
    hist.strain <- strain_plot_order_hist[st]
    rec.strain <- strain_plot_order_rec[st]
    points(st, R0.summary[(R0.summary$group == 'historic') & 
                            (R0.summary$strain == hist.strain), 
                          paste('R0', temp, sep='_')], pch=20, col="steelblue")
    arrows(st, R0.summary[(R0.summary$group == 'historic') & 
                            (R0.summary$strain == hist.strain), 
                          paste('R0', temp, sep='_')],
           y1=R0.summary[(R0.summary$group == 'historic') & 
                           (R0.summary$strain == hist.strain), 
                         paste(paste('R0', temp, sep='_'), 'l', sep='')],
           angle=90, length=0.1, lwd=1, col="steelblue")
    arrows(st, R0.summary[(R0.summary$group == 'historic') & 
                            (R0.summary$strain == hist.strain), 
                          paste('R0', temp, sep='_')],
           y1=R0.summary[(R0.summary$group == 'historic') & 
                           (R0.summary$strain == hist.strain), 
                         paste(paste('R0', temp, sep='_'), 'h', sep='')],
           angle=90, length=0.1, lwd=1, col="steelblue")
    
    points(st + 4, R0.summary[(R0.summary$group == 'recent') & 
                                (R0.summary$strain == rec.strain), 
                              paste('R0', temp, sep='_')], pch=20, col="firebrick")
    arrows(st + 4, R0.summary[(R0.summary$group == 'recent') & 
                                (R0.summary$strain == rec.strain), 
                              paste('R0', temp, sep='_')],
           y1=R0.summary[(R0.summary$group == 'recent') & 
                           (R0.summary$strain == rec.strain), 
                         paste(paste('R0', temp, sep='_'), 'l', sep='')],
           angle=90, length=0.1, lwd=1, col="firebrick")
    arrows(st + 4, R0.summary[(R0.summary$group == 'recent') & 
                                (R0.summary$strain == rec.strain), 
                              paste('R0', temp, sep='_')],
           y1=R0.summary[(R0.summary$group == 'recent') & 
                           (R0.summary$strain == rec.strain), 
                         paste(paste('R0', temp, sep='_'), 'h', sep='')],
           angle=90, length=0.1, lwd=1, col="firebrick")
    #arrows(st + 4, R0.mean.summary.rec[, 'R0'], y1=R0.mean.summary.rec[, 'R0_CI_hi'], 
    #       angle=90, length=0.1, lwd=1, col="firebrick")
    #arrows(st + 4, R0.mean.summary.rec[, 'R0'], y1=R0.mean.summary.rec[, 'R0_CI_lo'], 
    #       angle=90, length=0.1, lwd=1, col="firebrick")
  }
  
}

R0.summary.quot <- data.frame(group = c(rep('historic', 4), rep('recent', 4)),
                         strain= rep(1:4, 2),
                         R0_20=rep(0, 8), R0_20l=rep(0,8), R0_20h=rep(0,8),
                         R0_24=rep(0, 8), R0_24l=rep(0,8), R0_24h=rep(0,8),
                         R0_28=rep(0, 8), R0_28l=rep(0,8), R0_28h=rep(0,8))


#ctrl.R0.curve <- R0.hist.curves[[2]]
ctrl.R0.curve <- R0.mean.hist

for(st in 1:4) {
  cond <- (R0.summary$group == "historic") & (R0.summary$strain == st) 
  curr.curves <- R0.hist.curves[[st]] / ctrl.R0.curve
  
  R0.summary.quot$R0_20[cond] <- mean(curr.curves[1, ])
  R0.summary.quot$R0_20l[cond] <- quantile(curr.curves[1, ], 0.025)
  R0.summary.quot$R0_20h[cond] <- quantile(curr.curves[1, ], 0.975)
  
  R0.summary.quot$R0_24[cond] <- mean(curr.curves[2, ])
  R0.summary.quot$R0_24l[cond] <- quantile(curr.curves[2, ], 0.025)
  R0.summary.quot$R0_24h[cond] <- quantile(curr.curves[2, ], 0.975)
  
  R0.summary.quot$R0_28[cond] <- mean(curr.curves[3, ])
  R0.summary.quot$R0_28l[cond] <- quantile(curr.curves[3, ], 0.025)
  R0.summary.quot$R0_28h[cond] <- quantile(curr.curves[3, ], 0.975)
}



for(st in 1:4) {
  cond <- (R0.summary$group == "recent") & (R0.summary$strain == st) 
  curr.curves <- R0.rec.curves[[st]] / ctrl.R0.curve
  R0.summary.quot$R0_20[cond] <- mean(curr.curves[1, ])
  R0.summary.quot$R0_20l[cond] <- quantile(curr.curves[1, ], 0.025)
  R0.summary.quot$R0_20h[cond] <- quantile(curr.curves[1, ], 0.975)
  
  R0.summary.quot$R0_24[cond] <- mean(curr.curves[2, ])
  R0.summary.quot$R0_24l[cond] <- quantile(curr.curves[2, ], 0.025)
  R0.summary.quot$R0_24h[cond] <- quantile(curr.curves[2, ], 0.975)
  
  R0.summary.quot$R0_28[cond] <- mean(curr.curves[3, ])
  R0.summary.quot$R0_28l[cond] <- quantile(curr.curves[3, ], 0.025)
  R0.summary.quot$R0_28h[cond] <- quantile(curr.curves[3, ], 0.975)
}

#View(R0.summary.quot)

#mgp=c(3,1,0)
par(mfrow=c(1, 3), mar = c(7.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0))
for(temp in c(20, 24, 28)) {
  plot("", ylim=c(0.2, 1.8), xlab="", xlim=c(0.5, 8.5),
       ylab="rel R0 / mean rel R0 historic", main=paste(temp, '°C', sep=''),
       xaxt='n', col="steelblue")
  axis(side=1, at=1:8, labels=strain_names, las=2)
  abline(h=1.0, lty=2)
  # Historic
  for(st in 1:4) {
    hist.strain <- strain_plot_order_hist[st]
    rec.strain <- strain_plot_order_rec[st]
    points(st, R0.summary.quot[(R0.summary.quot$group == 'historic') & 
                            (R0.summary.quot$strain == hist.strain), 
                          paste('R0', temp, sep='_')], pch=20, col="steelblue")
    arrows(st, R0.summary.quot[(R0.summary$group == 'historic') & 
                            (R0.summary$strain == hist.strain), 
                          paste('R0', temp, sep='_')],
           y1=R0.summary.quot[(R0.summary.quot$group == 'historic') & 
                           (R0.summary.quot$strain == hist.strain), 
                         paste(paste('R0', temp, sep='_'), 'l', sep='')],
           angle=90, length=0.1, lwd=1, col="steelblue")
    arrows(st, R0.summary.quot[(R0.summary.quot$group == 'historic') & 
                            (R0.summary.quot$strain == hist.strain), 
                          paste('R0', temp, sep='_')],
           y1=R0.summary.quot[(R0.summary.quot$group == 'historic') & 
                           (R0.summary.quot$strain == hist.strain), 
                         paste(paste('R0', temp, sep='_'), 'h', sep='')],
           angle=90, length=0.1, lwd=1, col="steelblue")
    
    points(st + 4, R0.summary.quot[(R0.summary.quot$group == 'recent') & 
                                (R0.summary.quot$strain == rec.strain), 
                              paste('R0', temp, sep='_')], pch=20, col="firebrick")
    arrows(st + 4, R0.summary.quot[(R0.summary.quot$group == 'recent') & 
                                (R0.summary.quot$strain == rec.strain), 
                              paste('R0', temp, sep='_')],
           y1=R0.summary.quot[(R0.summary.quot$group == 'recent') & 
                           (R0.summary.quot$strain == rec.strain), 
                         paste(paste('R0', temp, sep='_'), 'l', sep='')],
           angle=90, length=0.1, lwd=1, col="firebrick")
    arrows(st + 4, R0.summary.quot[(R0.summary.quot$group == 'recent') & 
                                (R0.summary.quot$strain == rec.strain), 
                              paste('R0', temp, sep='_')],
           y1=R0.summary.quot[(R0.summary.quot$group == 'recent') & 
                           (R0.summary.quot$strain == rec.strain), 
                         paste(paste('R0', temp, sep='_'), 'h', sep='')],
           angle=90, length=0.1, lwd=1, col="firebrick")
    #arrows(st + 4, R0.mean.summary.rec[, 'R0'], y1=R0.mean.summary.rec[, 'R0_CI_hi'], 
    #       angle=90, length=0.1, lwd=1, col="firebrick")
    #arrows(st + 4, R0.mean.summary.rec[, 'R0'], y1=R0.mean.summary.rec[, 'R0_CI_lo'], 
    #       angle=90, length=0.1, lwd=1, col="firebrick")
    
  }
  
}
R0.summary
R0.summary.quot

par(mfrow=c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1), mgp=c(3,1,0))

plot(R0.mean.summary.hist[, 'temperature'] - 0.5, R0.mean.summary.hist[, 'R0'], pch=20, ylim=c(20, 120),
     xlim=c(18, 30), xlab="Temperature [°C]", ylab="relative R0",
     xaxt='n', yaxt='n', col="steelblue")
axis(side=1, at=c(20, 24, 28))
arrows(R0.mean.summary.hist[, 'temperature'] - 0.5, R0.mean.summary.hist[, 'R0'], y1=R0.mean.summary.hist[, 'R0_CI_hi'], 
       angle=90, length=0.1, lwd=1, col="steelblue")
arrows(R0.mean.summary.hist[, 'temperature'] - 0.5, R0.mean.summary.hist[, 'R0'], y1=R0.mean.summary.hist[, 'R0_CI_lo'], 
       angle=90, length=0.1, lwd=1, col="steelblue")

points(R0.mean.summary.rec[, 'temperature'] + 0.5, R0.mean.summary.rec[, 'R0'], pch=20, col="firebrick")
arrows(R0.mean.summary.rec[, 'temperature'] + 0.5, R0.mean.summary.rec[, 'R0'], y1=R0.mean.summary.rec[, 'R0_CI_hi'], 
       angle=90, length=0.1, lwd=1, col="firebrick")
arrows(R0.mean.summary.rec[, 'temperature'] + 0.5, R0.mean.summary.rec[, 'R0'], y1=R0.mean.summary.rec[, 'R0_CI_lo'], 
       angle=90, length=0.1, lwd=1, col="firebrick")



plot(R0.mean.summary[, 'temperature'], R0.mean.summary[, 'R0_quot'], pch=20, ylim=c(0.5, 1.5),
     xlim=c(18, 30), xlab="Temperature [°C]", ylab="contemporary / historic", main="Relative R0 ratio",
     xaxt='n')
axis(side=1, at=c(20, 24, 28))
arrows(R0.mean.summary[, 'temperature'], R0.mean.summary[, 'R0_quot'], y1=R0.mean.summary[, 'R0_quot_CI_hi'], 
       angle=90, length=0.1, lwd=1)
arrows(R0.mean.summary[, 'temperature'], R0.mean.summary[, 'R0_quot'], y1=R0.mean.summary[, 'R0_quot_CI_lo'], 
       angle=90, length=0.1, lwd=1)
abline(h=1, lty=2)
