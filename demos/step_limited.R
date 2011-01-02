#step_limited.R
library(BranchingTime)
library(socialR)
tags <- c("adaptivedynamics simulations")

log <- gitlog()
rep <- 16*5
cpu <- 16


a <- ensemble_sim(rep=rep, sigma_mu = 0.05, mu =.0001, sigma_c2 = 0.3, sigma_k2 = 1, ko = 500, xo = 0.5, threshold = 30, cpu=cpu, maxtime=1e7) 
b <-  ensemble_sim(rep=rep, sigma_mu = 0.02, mu =.005, sigma_c2 = 0.9, sigma_k2 = 1, ko = 500, xo = 0.5, threshold = 30, cpu=cpu, maxtime=1e7) 

save(list=ls(), file="step_limited.Rdat")

social_plot(plottries(a), file="step_limited.png", tags=tags)
social_plot(plottries(b), file="step_limited.png", tags=tags)
