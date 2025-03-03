library('R2jags')
library('mcmcplots')
library('MCMCvis')
library('scales')

historic.out <- readRDS("bc_hier_historic.RDS")
recent.out <- readRDS("bc_hier_recent.RDS")

historic.chains <- MCMCchains(historic.out)
recent.chains <- MCMCchains(recent.out)

# Calculate posterior distribution of mean Topt of four historic and contemporary
# strains.
hist.Topt.chains <- (historic.chains[, "Topt[1]"] + historic.chains[, "Topt[2]"] + 
                    + historic.chains[, "Topt[3]"] + historic.chains[, "Topt[4]"]) / 4

rec.Topt.chains <- (recent.chains[, "Topt[1]"] + recent.chains[, "Topt[2]"] + 
                       + recent.chains[, "Topt[3]"] + recent.chains[, "Topt[4]"]) / 4

# Mean and 95% CI.
mean(rec.Topt.chains)
quantile(rec.Topt.chains, c(0.025, 0.975))

# Mean and 95% CI.
mean(hist.Topt.chains)
quantile(hist.Topt.chains, c(0.025, 0.975))


# Mean and 95% CI.
mean(rec.Topt.chains - hist.Topt.chains)
quantile(rec.Topt.chains - hist.Topt.chains, c(0.025, 0.975))

# Posterior probability that contemporary Topt mean is higher than historic
# Topt mean.
mean((rec.Topt.chains - hist.Topt.chains) > 0)


mean(recent.chains[, "rmax_m"] - historic.chains[, "rmax_m"])
quantile(recent.chains[, "rmax_m"] - historic.chains[, "rmax_m"], c(0.025, 0.975))

mean((recent.chains[, "rmax_m"] - historic.chains[, "rmax_m"]) > 0)
