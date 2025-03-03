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
  plot(temps, meancurve, type="l", col="black", xlim=xlim, ylim=ylim)
  lines(temps, CI[1,], col="gray")
  lines(temps, CI[2,], col="gray")
}

## 1. Load chains with posterior estimates of parameters and calculate trait
## values.

temps <- seq(0, 50, 0.05) # Test with small vector first to debug.

## Load data for biting rate
load("./shocket_fits/jagsout_a_Cpip_inf.Rdata")
a.Cpip.out.inf

# Only keep chains for TPC model parameters to save memory.
a.chains <- MCMCchains(a.Cpip.out.inf, params=c("cf.T0", "cf.Tm", "cf.q"))

# Remove data to save memory
rm(a.Cpip.out.inf)

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
load("./shocket_fits/jagsout_lf_Cpip_inf.Rdata")
lf.Cpip.out.inf

lf.chains <- MCMCchains(lf.Cpip.out.inf, params=c("cf.m", "cf.b"))
rm(lf.Cpip.out.inf)

lf.curves <- apply(lf.chains, 1, function(x) linear_lim(temps, x[1], x[2]))

plot_curves(temps, lf.curves, ylim=c(0, 100))

## Pathogen development rate
load("./shocket_fits/jagsout_PDR_CpipWNV_inf.Rdata")
PDR.CpipWNV.out.inf$BUGSoutput$summary

PDR.chains <- MCMCchains(PDR.CpipWNV.out.inf, params=c("cf.T0", "cf.Tm", "cf.q"))
rm(PDR.CpipWNV.out.inf)

PDR.curves <- apply(PDR.chains, 1, function(x) briere(temps, x[1], x[2], x[3]))

plot_curves(temps, PDR.curves, ylim=c(0, 0.3))

## Fecundity
load("./shocket_fits/jagsout_EFOC_Cpip_inf.Rdata")
EFOC.Cpip.out.inf
EFOC.chains <- MCMCchains(EFOC.Cpip.out.inf, params=c("cf.T0", "cf.Tm", "cf.q"))

EFOC.curves <- apply(EFOC.chains, 1, function(x) quad(temps, x[1], x[2], x[3]))

plot_curves(temps, EFOC.curves, ylim=c(0, 400))
#plot_curves(temps, ER.alb.curves)


## Egg viability 
load("./shocket_fits/jagsout_EV_Cpip_inf.Rdata")
EV.chains <- MCMCchains(EV.Cpip.out.inf, params=c("cf.T0", "cf.Tm", "cf.q"))
EV.Cpip.out.inf$BUGSoutput$summary
rm(EV.Cpip.out.inf)

EV.curves <- apply(EV.chains, 1, function(x) quad_lim(temps, x[1], x[2], x[3]))

## Larval survival
load("./shocket_fits/jagsout_pLA_Cpip_inf.Rdata")
pLA.Cpip.out.inf$BUGSoutput$summary
pLA.chains <- MCMCchains(pLA.Cpip.out.inf, params=c("cf.T0", "cf.Tm", "cf.q"))
rm(pLA.Cpip.out.inf)

pLA.curves <- apply(pLA.chains, 1, function(x) quad_lim(temps, x[1], x[2], x[3]))

plot_curves(temps, pLA.curves)

## Mosquito development rate
load("./shocket_fits/jagsout_MDR_Cpip_inf.Rdata")
MDR.Cpip.out.inf$BUGSoutput$summary

MDR.chains <- MCMCchains(MDR.Cpip.out.inf, params=c("cf.T0", "cf.Tm", "cf.q"))

MDR.chains
MDR.curves <- apply(MDR.chains, 1, function(x) briere(temps, x[1], x[2], x[3]))

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
                         Tmin=rep(0, 8), Tminl=rep(0,8), Tminh=rep(0,8),
                         Topt=rep(0, 8), Toptl=rep(0,8), Topth=rep(0,8),
                         Tmax=rep(0, 8), Tmaxl=rep(0,8), Tmaxh=rep(0,8))
R0.summary

R0.hist.curves <- list()
#(a, bc, lf, PDR, EFGC, EV, pLA, MDR)
for(st in 1:4) {
  R0.hist.curves[[st]] <- R0(a.curves, bc.hist.curves[[st]], lf.curves, 
                             PDR.curves[, sample(1:80000, 80000)], EFOC.curves, EV.curves, 
                             pLA.curves, MDR.curves)
  #print(st)
  Topt.chains <- find_Topt(temps, R0.hist.curves[[st]])
  R0.summary$Topt[(R0.summary$group == 'historic') & 
                    (R0.summary$strain == st)] <- mean(Topt.chains)
  Topt.CI <- quantile(Topt.chains, c(0.025, 0.975))
  R0.summary$Toptl[(R0.summary$group == 'historic') & 
                     (R0.summary$strain == st)] <- Topt.CI[1]
  R0.summary$Topth[(R0.summary$group == 'historic') & 
                     (R0.summary$strain == st)] <- Topt.CI[2]
  rm(Topt.chains)
  
  Tmax.chains <- find_Tmax(temps, R0.hist.curves[[st]])
  R0.summary$Tmax[(R0.summary$group == 'historic') & 
                    (R0.summary$strain == st)] <- mean(Tmax.chains)
  Tmax.CI <- quantile(Tmax.chains, c(0.025, 0.975))
  R0.summary$Tmaxl[(R0.summary$group == 'historic') & 
                     (R0.summary$strain == st)] <- Tmax.CI[1]
  R0.summary$Tmaxh[(R0.summary$group == 'historic') & 
                     (R0.summary$strain == st)] <- Tmax.CI[2]
  rm(Tmax.chains)
  
  
  Tmin.chains <- find_Tmin(temps, R0.hist.curves[[st]])
  R0.summary$Tmin[(R0.summary$group == 'historic') & 
                    (R0.summary$strain == st)] <- mean(Tmin.chains)
  Tmin.CI <- quantile(Tmin.chains, c(0.025, 0.975))
  R0.summary$Tminl[(R0.summary$group == 'historic') & 
                     (R0.summary$strain == st)] <- Tmin.CI[1]
  R0.summary$Tminh[(R0.summary$group == 'historic') & 
                     (R0.summary$strain == st)] <- Tmin.CI[2]
  rm(Tmin.chains)
}
R0.hist.mean <- (R0.hist.curves[[1]] + R0.hist.curves[[2]] + R0.hist.curves[[3]]
                 + R0.hist.curves[[4]]) / 4

R0.rec.curves <- list()

for(st in 1:4) {
  R0.rec.curves[[st]] <- R0(a.curves, bc.rec.curves[[st]], lf.curves, 
                            PDR.curves[, sample(1:80000, 80000)], EFOC.curves, EV.curves, 
                             pLA.curves, MDR.curves)
  #print(st)
  Topt.chains <- find_Topt(temps, R0.rec.curves[[st]])
  R0.summary$Topt[(R0.summary$group == 'recent') & 
                    (R0.summary$strain == st)] <- mean(Topt.chains)
  Topt.CI <- quantile(Topt.chains, c(0.025, 0.975))
  R0.summary$Toptl[(R0.summary$group == 'recent') & 
                     (R0.summary$strain == st)] <- Topt.CI[1]
  R0.summary$Topth[(R0.summary$group == 'recent') & 
                     (R0.summary$strain == st)] <- Topt.CI[2]
  rm(Topt.chains)
  
  Tmax.chains <- find_Tmax(temps, R0.rec.curves[[st]])
  R0.summary$Tmax[(R0.summary$group == 'recent') & 
                    (R0.summary$strain == st)] <- mean(Tmax.chains)
  Tmax.CI <- quantile(Tmax.chains, c(0.025, 0.975))
  R0.summary$Tmaxl[(R0.summary$group == 'recent') & 
                     (R0.summary$strain == st)] <- Tmax.CI[1]
  R0.summary$Tmaxh[(R0.summary$group == 'recent') & 
                     (R0.summary$strain == st)] <- Tmax.CI[2]
  rm(Tmax.chains)
  
  
  Tmin.chains <- find_Tmin(temps, R0.rec.curves[[st]])
  R0.summary$Tmin[(R0.summary$group == 'recent') & 
                    (R0.summary$strain == st)] <- mean(Tmin.chains)
  Tmin.CI <- quantile(Tmin.chains, c(0.025, 0.975))
  R0.summary$Tminl[(R0.summary$group == 'recent') & 
                     (R0.summary$strain == st)] <- Tmin.CI[1]
  R0.summary$Tminh[(R0.summary$group == 'recent') & 
                     (R0.summary$strain == st)] <- Tmin.CI[2]
  rm(Tmin.chains)
}


rm(a.chains)
rm(EFOC.chains)
rm(EV.chains)
rm(lf.chains)
rm(MDR.chains)
rm(PDR.chains)
rm(pLA.chains)

rm(a.curves) 
rm(EFOC.curves)
rm(EV.curves)
rm(lf.curves)
rm(MDR.curves)
rm(PDR.curves)
rm(pLA.curves)


R0.rec.mean <- (R0.rec.curves[[1]] + R0.rec.curves[[2]] + R0.rec.curves[[3]]
                 + R0.rec.curves[[4]]) / 4


R0.ratio <- R0.rec.mean / R0.hist.mean
str(R0.ratio)

idx <- (temps >= 17.5) & (temps <= 30.5) 

R0.ratio[idx, ]
mean.R0.ratio <- apply(R0.ratio[idx, ], 1, mean)
median.R0.ratio <- apply(R0.ratio[idx, ], 1, median, na.rm=TRUE)

CI.R0.ratio <- apply(R0.ratio[idx, ], 1, quantile, c(0.025, 0.975), na.rm=TRUE)
plot(temps[idx], mean.R0.ratio, type='l', ylim=c(0.5, 1.5), xlim=c(18.4, 29.6),
     xlab='Temperature [°C]', ylab='contemporary / historic', main='Relative R0 ratio')
polygon(c(temps[idx], rev(temps[idx])), c(CI.R0.ratio[1,], rev(CI.R0.ratio[2,])), 
        col=alpha("black", 0.2), lty=0)
abline(h=1.0, lty=2)
#View(R0.summary)

# Comparison plot.
strain_names <- c("2003.1", "2003.2", "2004.1", "2004.2",
                  "2017.1", "2017.2", "2018.1", "2018.2")

## Old labels:
## c("DS03", "US03", "DS04", "US04", 
## "DS17", "US17", "DS18", "US18")
strain_plot_order_hist <- c(1,2,3,4)
strain_plot_order_rec <- c(2,1,4,3)

y <- 1:8
plot("", xlim=c(10, 40), ylim=c(-0.5, 8), xlab="Temperature [°C]", ylab="",
     main="Relative R0", xaxt="n", yaxt = "n")
axis(1, at = seq(10, 40, by = 10))
rug(x = seq(15, 40, 5), ticksize = -0.03, side = 1)
rug(x = 10:40, ticksize = -0.01, side = 1)

axis(2, at=1:8, labels=strain_names, las=2)
for(st in strain_plot_order_hist) {
  df <- subset(R0.summary, (R0.summary$group == "historic") & 
                 (R0.summary$strain == st))
  points(df$Topt, st, pch=20, col="steelblue")
  lines(c(df$Toptl, df$Topth), rep(st, 2), col="steelblue", lwd=2)
  
  points(df$Tmax, st, pch=20, col="steelblue")
  lines(c(df$Tmaxl, df$Tmaxh), rep(st, 2), col="steelblue", lwd=2)
  
  points(df$Tmin, st, pch=20, col="steelblue")
  lines(c(df$Tminl, df$Tminh), rep(st, 2), col="steelblue", lwd=2)
}

for(st in strain_plot_order_rec) {
  df <- subset(R0.summary, (R0.summary$group == "recent") & 
                 (R0.summary$strain == st))
  points(df$Topt, st+4, pch=20, col="firebrick")
  lines(c(df$Toptl, df$Topth), rep(st+4, 2), col="firebrick", lwd=2)
  
  points(df$Tmax, st+4, pch=20, col="firebrick")
  lines(c(df$Tmaxl, df$Tmaxh), rep(st+4, 2), col="firebrick", lwd=2)
  
  points(df$Tmin, st+4, pch=20, col="firebrick")
  lines(c(df$Tminl, df$Tminh), rep(st+4, 2), col="firebrick", lwd=2)
}

text(15, 0, "Tmin")
text(25, 0, "Topt")
text(35, 0, "Tmax")

## bc
bc.hist.chains$BUGSoutput$summary["Tmax[1]", "mean"]
paste("Topt[", st, "]")
bc.rec.chains$BUGSoutput$summary

plot("", xlim=c(5, 55), ylim=c(-0.5, 8), xlab="Temperature [°C]", ylab="",
     main="bc (disseminated/exposed)", xaxt="n", yaxt = "n")
axis(1, at = seq(10, 50, by = 10))
rug(x = seq(5, 55, 5), ticksize = -0.03, side = 1)
rug(x = 4:58, ticksize = -0.01, side = 1)

axis(2, at=1:8, labels=strain_names, las=2)

for(st in strain_plot_order_hist) {
  vname <- paste("Topt[", st, "]", sep="")
  points(bc.hist.chains$BUGSoutput$summary[vname, "mean"], st, pch=20, col="steelblue")
  lines(c(bc.hist.chains$BUGSoutput$summary[vname, "2.5%"], 
          bc.hist.chains$BUGSoutput$summary[vname, "97.5%"]),
        rep(st, 2), col="steelblue", lwd=2)
  
  vname <- paste("Tmax[", st, "]", sep="")
  points(bc.hist.chains$BUGSoutput$summary[vname, "mean"], st, pch=20, col="steelblue")
  lines(c(bc.hist.chains$BUGSoutput$summary[vname, "2.5%"], 
          bc.hist.chains$BUGSoutput$summary[vname, "97.5%"]),
        rep(st, 2), col="steelblue", lwd=2)
  
  vname <- paste("Tmin[", st, "]", sep="")
  points(bc.hist.chains$BUGSoutput$summary[vname, "mean"], st, pch=20, col="steelblue")
  lines(c(bc.hist.chains$BUGSoutput$summary[vname, "2.5%"], 
          bc.hist.chains$BUGSoutput$summary[vname, "97.5%"]),
        rep(st, 2), col="steelblue", lwd=2)
}

for(st in strain_plot_order_rec) {
  vname <- paste("Topt[", st, "]", sep="")
  points(bc.rec.chains$BUGSoutput$summary[vname, "mean"], st+4, pch=20, col="firebrick")
  lines(c(bc.rec.chains$BUGSoutput$summary[vname, "2.5%"], 
          bc.rec.chains$BUGSoutput$summary[vname, "97.5%"]),
        rep(st+4, 2), col="firebrick", lwd=2)
  
  vname <- paste("Tmax[", st, "]", sep="")
  points(bc.rec.chains$BUGSoutput$summary[vname, "mean"], st+4, pch=20, col="firebrick")
  lines(c(bc.rec.chains$BUGSoutput$summary[vname, "2.5%"], 
          bc.rec.chains$BUGSoutput$summary[vname, "97.5%"]),
        rep(st+4, 2), col="firebrick", lwd=2)
  
  vname <- paste("Tmin[", st, "]", sep="")
  points(bc.rec.chains$BUGSoutput$summary[vname, "mean"], st+4, pch=20, col="firebrick")
  lines(c(bc.rec.chains$BUGSoutput$summary[vname, "2.5%"], 
          bc.rec.chains$BUGSoutput$summary[vname, "97.5%"]),
        rep(st+4, 2), col="firebrick", lwd=2)
}

text(15, 0, "Tmin")
text(30, 0, "Topt")
text(45, 0, "Tmax")


plot("", xlim=c(0.25, 0.85), ylim=c(-0.5, 8), xlab="bc(Topt)", ylab="",
     main="max bc (disseminated/exposed)", yaxt = "n")
axis(2, at=1:8, labels=strain_names, las=2)
for(st in strain_plot_order_hist) {
  vname <- paste("rmax[", st, "]", sep="")
  points(bc.hist.chains$BUGSoutput$summary[vname, "mean"], st, pch=20, col="steelblue")
  lines(c(bc.hist.chains$BUGSoutput$summary[vname, "2.5%"], 
          bc.hist.chains$BUGSoutput$summary[vname, "97.5%"]),
        rep(st, 2), col="steelblue", lwd=2)
}

for(st in strain_plot_order_rec) {
  vname <- paste("rmax[", st, "]", sep="")
  points(bc.rec.chains$BUGSoutput$summary[vname, "mean"], st+4, pch=20, col="firebrick")
  lines(c(bc.rec.chains$BUGSoutput$summary[vname, "2.5%"], 
          bc.rec.chains$BUGSoutput$summary[vname, "97.5%"]),
        rep(st+4, 2), col="firebrick", lwd=2)
}

bc.hist.chains
bc.rec.chains


# Find mean curves and credible intervals.

alpha <- 0.1
plot("", main="Relative R0", xlab="Temperature [°C]", ylab="Relative R0",
     xlim=c(10, 40), ylim=c(0, 120.0), xaxt="n", yaxt="n")
axis(1, at=seq(10, 40, 10))
rug(x = seq(15, 35, 5), ticksize = -0.03, side = 1)
rug(x = 10:40, ticksize = -0.01, side = 1)
for(st in 1:4) {
  meancurve.hist <- apply(R0.hist.curves[[st]], 1, mean)
  CI.hist <- apply(R0.hist.curves[[st]], 1, quantile, c(0.025, 0.975))
  
  lines(temps, meancurve.hist, col="steelblue", lwd=1)
  #polygon(c(temps, rev(temps)), c(CI.hist[1,], rev(CI.hist[2,])), 
  #        col=alpha("steelblue", alpha), lty=0)
  
  meancurve.rec <- apply(R0.rec.curves[[st]], 1, mean)
  CI.rec <- apply(R0.rec.curves[[st]], 1, quantile, c(0.025, 0.975))
  
  lines(temps, meancurve.rec, col="firebrick", lwd=1)
  #polygon(c(temps, rev(temps)), c(CI.rec[1,], rev(CI.rec[2,])), 
  #        col=alpha("firebrick", alpha), lty=0)  
  
}


bc.rec.chains$BUGSoutput$summary



alpha.val <- 0.1
par(mfrow=c(2,2))

strain_plot_order_hist <- c(1,2,3,4)
strain_plot_order_rec <- c(2,1,4,3)
strain_names <- c("2003.1", "2003.2", "2004.1", "2004.2",
                  "2017.1", "2017.2", "2018.1", "2018.2")
R0.summary
for(st in 1:4) {
  hist.strain <- strain_plot_order_hist[st]
  plot("", main=strain_names[[st]], xlab="Temperature [°C]", ylab="Relative R0",
       xlim=c(10, 40), ylim=c(0, 140.0), xaxt="n", yaxt="n")
  axis(1, at=seq(10, 40, 10))
  rug(x = seq(15, 35, 5), ticksize = -0.03, side = 1)
  rug(x = 10:40, ticksize = -0.01, side = 1)
  meancurve.hist <- apply(R0.hist.curves[[hist.strain]], 1, mean)
  CI.hist <- apply(R0.hist.curves[[hist.strain]], 1, quantile, c(0.025, 0.975))
  
  lines(temps, meancurve.hist, col="steelblue", lwd=1)
  polygon(c(temps, rev(temps)), c(CI.hist[1,], rev(CI.hist[2,])), 
          col=alpha("steelblue", alpha.val), lty=0)
  
  for(temp in c(20, 24, 28)) {
    points(temp, R0.summary[(R0.summary$group == 'historic') & 
                            (R0.summary$strain == hist.strain), 
                          paste('R0', temp, sep='_')], pch=20, col="steelblue")
    arrows(temp, R0.summary[(R0.summary$group == 'historic') & 
                            (R0.summary$strain == hist.strain), 
                          paste('R0', temp, sep='_')],
           y1=R0.summary[(R0.summary$group == 'historic') & 
                           (R0.summary$strain == hist.strain), 
                         paste(paste('R0', temp, sep='_'), 'l', sep='')],
           angle=90, length=0.1, lwd=1, col="steelblue")
    arrows(temp, R0.summary[(R0.summary$group == 'historic') & 
                            (R0.summary$strain == hist.strain), 
                          paste('R0', temp, sep='_')],
           y1=R0.summary[(R0.summary$group == 'historic') & 
                           (R0.summary$strain == hist.strain), 
                         paste(paste('R0', temp, sep='_'), 'h', sep='')],
           angle=90, length=0.1, lwd=1, col="steelblue")
    }
  
}


alpha.val <- 0.1
par(mfrow=c(2,2))

for(st in 1:4) {
  rec.strain <- strain_plot_order_rec[st]
  plot("", main=strain_names[[st + 4]], xlab="Temperature [°C]", ylab="Relative R0",
       xlim=c(10, 40), ylim=c(0, 140.0), xaxt="n", yaxt="n")
  axis(1, at=seq(10, 40, 10))
  rug(x = seq(15, 35, 5), ticksize = -0.03, side = 1)
  rug(x = 10:40, ticksize = -0.01, side = 1)
  meancurve.hist <- apply(R0.hist.curves[[rec.strain]], 1, mean)
  CI.hist <- apply(R0.hist.curves[[rec.strain]], 1, quantile, c(0.025, 0.975))
  
  lines(temps, meancurve.rec, col="firebrick", lwd=1)
  polygon(c(temps, rev(temps)), c(CI.rec[1,], rev(CI.rec[2,])), 
          col=alpha("firebrick", alpha.val), lty=0) 
  
  for(temp in c(20, 24, 28)) {
  points(temp, R0.summary[(R0.summary$group == 'recent') & 
                              (R0.summary$strain == rec.strain), 
                            paste('R0', temp, sep='_')], pch=20, col="firebrick")
  arrows(temp, R0.summary[(R0.summary$group == 'recent') & 
                              (R0.summary$strain == rec.strain), 
                            paste('R0', temp, sep='_')],
         y1=R0.summary[(R0.summary$group == 'recent') & 
                         (R0.summary$strain == rec.strain), 
                       paste(paste('R0', temp, sep='_'), 'l', sep='')],
         angle=90, length=0.1, lwd=1, col="firebrick")
  arrows(temp, R0.summary[(R0.summary$group == 'recent') & 
                              (R0.summary$strain == rec.strain), 
                            paste('R0', temp, sep='_')],
         y1=R0.summary[(R0.summary$group == 'recent') & 
                         (R0.summary$strain == rec.strain), 
                       paste(paste('R0', temp, sep='_'), 'h', sep='')],
         angle=90, length=0.1, lwd=1, col="firebrick")
  }
  
}

