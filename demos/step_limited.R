#step_limited.R
library(BranchingTime)
library(socialR)
tags <- c("adaptivedynamics", "simulations")

log <- gitlog()
rep <- 16
cpu <- 16
K <- 10


#vary_mu <- lapply(1:K, function(i){
#	mu <- seq(1e-4, 1e-3, length=K)
#	ensemble_sim(rep=rep, sigma_mu = 0.05, mu =mu[i], sigma_c2 = 0.5, sigma_k2 = 1, ko = 500, xo = 0.1, threshold = 30, cpu=cpu, maxtime=5e6) 
#	})

#save(list=ls(), file="step_limited.Rdat")
load("step_limited.Rdat")
source("phases_plots.R")
mu <- seq(1e-4, 1e-3, length=K)
plot_phases(vary_mu, mu, xlab="mu")


