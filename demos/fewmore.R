#fewmore.R, based on step_limited.R
library(BranchingTime)
library(socialR)
tags <- c("adaptivedynamics simulations")
source("../R/branching_time.R")

log <- gitlog()
rep <- 16*5
cpu <- 16

d <- ensemble_sim(rep=rep, sigma_mu = 0.05, mu =.0001, sigma_c2 = 0.3, sigma_k2 = 1, ko = 500, xo = 0.5, threshold = 30, cpu=cpu, maxtime=1e7) 
save(list=ls(), file="fewmore.Rdat")
social_plot(plot_waitingtimes(d), file="step_limited.png", tags=tags, comment="d) sigmac2=.3, mu=.0001")
social_plot(plot_waitingtimes(d, HISTOGRAM=T), file="step_limited.png", tags=tags, comment="d) sigmac2=.3, mu=.0001")
social_plot(plot_failures(d), file="step_limited.png", tags=tags, comment="d) sigmac2=.3, mu=.0001")

e <- ensemble_sim(rep=rep, sigma_mu = 0.05, mu =.001, sigma_c2 = 0.5, sigma_k2 = 1, ko = 500, xo = 0.5, threshold = 30, cpu=cpu, maxtime=1e7) 
save(list=ls(), file="fewmore.Rdat")
social_plot(plot_waitingtimes(e), file="step_limited.png", tags=tags, comment="e) sigmac2=.5, mu=.001")
social_plot(plot_failures(e), file="step_limited.png", tags=tags, comment="e) sigmac2=.5, mu=.001")

