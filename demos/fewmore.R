#fewmore.R, based on step_limited.R
library(BranchingTime)
library(socialR)
tags <- c("adaptivedynamics simulations")
source("../R/branching_time.R")

log <- gitlog()
rep <- 16*5
cpu <- 16

f <- ensemble_sim(rep=rep, sigma_mu = 0.05, mu =.0005, sigma_c2 = 0.3, sigma_k2 = 1, ko = 500, xo = 0.5, threshold = 30, cpu=cpu, maxtime=1e7) 
save(list=ls(), file="fewmore.Rdat")
social_plot(plot_waitingtimes(f), file="step_limited.png", tags=tags, comment="f) sigmac2=.3, mu=.0005")
social_plot(plot_waitingtimes(f, HISTOGRAM=T), file="step_limited.png", tags=tags, comment="f) sigmac2=.3, mu=.0005")
social_plot(plot_failures(f), file="step_limited.png", tags=tags, comment="f) sigmac2=.3, mu=.0005")

