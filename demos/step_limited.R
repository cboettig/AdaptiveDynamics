#step_limited.R
library(BranchingTime)
library(socialR)
tags <- c("adaptivedynamics simulations")
source("../R/branching_time.R")

log <- gitlog()
rep <- 16*5
cpu <- 16

a <- ensemble_sim(rep=rep, sigma_mu = 0.05, mu =.001, sigma_c2 = 0.1, sigma_k2 = 1, ko = 500, xo = 0.5, threshold = 30, cpu=cpu, maxtime=1e7) 
save(list=ls(), file="step_limited.Rdat")
social_plot(plot_waitingtimes(a), file="step_limited.png", tags=tags)
social_plot(plot_failures(a), file="step_limited.png", tags=tags)

b <- ensemble_sim(rep=rep, sigma_mu = 0.05, mu =.001, sigma_c2 = 0.3, sigma_k2 = 1, ko = 500, xo = 0.5, threshold = 30, cpu=cpu, maxtime=1e7) 
save(list=ls(), file="step_limited.Rdat")
social_plot(plot_waitingtimes(b), file="step_limited.png", tags=tags)
social_plot(plot_failures(b), file="step_limited.png", tags=tags)

c <-  ensemble_sim(rep=rep, sigma_mu = 0.02, mu =.01, sigma_c2 = 0.8, sigma_k2 = 1, ko = 500, xo = 0.5, threshold = 30, cpu=cpu, maxtime=1e7) 
save(list=ls(), file="step_limited.Rdat")
social_plot(plot_waitingtimes(c), file="step_limited.png", tags=tags)
social_plot(plot_failures(c), file="step_limited.png", tags=tags)

#load("step_limited.Rdat")

