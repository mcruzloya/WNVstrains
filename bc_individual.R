setwd('/Users/cruzloya/git/WNVstrains')

set.seed(42) # Set seed for reproducibility.

library('R2jags')
library('mcmcplots')
library('MCMCvis')
library('scales')

data <- read.csv('retrodata_max.csv')
data
historic
historic <- subset(data, data$group == 0)
recent <- subset(data, data$group == 1)

df
sink("bc_indiv_temp.txt")
cat("
    model{
    ## Likelihood
    for(Tidx in 1:N.temp){
      p[Tidx] ~ dunif(0, 1) 
      n[Tidx] ~ dbin(p[Tidx], N[Tidx])
    }
    
    } # close model
    ",fill=TRUE)
sink()

##### Calculate initial values for MCMC.
## We are picking random values so that every chain will start at a different place. 


##### Parameters to Estimate
parameters <- c("p")

##### MCMC Settings
# Number of posterior dist elements = [(ni - nb) / nt ] * nc = [ (25000 - 5000) / 8 ] * 3 = 7500
ni <- 200000 # number of iterations in each chain
nb <- 40000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 4 # number of chains

print('historic')
hist.model_output <- list()
for(st in 1:4) {
  df <- subset(historic, (historic$strain_id == st))
  
  ##### Organize Data for JAGS
  N <- df$total
  n <- df$disseminated
  temp <- df$temperature
  N.obs <- length(N)
  
  inits<-function(){list(
    p = runif(N.obs, 0, 1))}
  
  ##### Bundle Data
  jag.data<-list(N=N, n=n, N.temp=N.obs)
  print(jag.data)
  ##### Run JAGS
  jags.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                   model.file="bc_indiv_temp.txt", n.thin=nt, n.chains=nc, 
                   n.burnin=nb, n.iter=ni, DIC=TRUE, working.directory=getwd())
  hist.model_output[[st]] <- jags.out
  print(st)
  print(jags.out)
}

print('recent')
recent.model_output <- list()
for(st in 1:4) {
  df <- subset(recent, (recent$strain_id == st))
  
  ##### Organize Data for JAGS
  N <- df$total
  n <- df$disseminated
  temp <- df$temperature
  N.obs <- length(N)
  
  inits<-function(){list(
    p = runif(N.obs, 0, 1))}
  
  ##### Bundle Data
  jag.data<-list(N=N, n=n, N.temp=N.obs)
  print(jag.data)
  ##### Run JAGS
  jags.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                   model.file="bc_indiv_temp.txt", n.thin=nt, n.chains=nc, 
                   n.burnin=nb, n.iter=ni, DIC=TRUE, working.directory=getwd())
  recent.model_output[[st]] <- jags.out
  print(st)
  print(jags.out)
}
historic

hist.model_output[[2]]



saveRDS(hist.model_output, "bc_individual_historic.RDS")
saveRDS(recent.model_output, "bc_individual_recent.RDS")


