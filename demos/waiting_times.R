library(BranchingTime)
library(socialR)
tags <- c("adaptivedynamics", "simulations")

log <- gitlog()
# phasetime[1] = First time two coexisting branches are established
#			2  = Last (most recent time) two coexisting branches were established 
#			3  = Time to finish phase 2: (satisfies invade_pair() test: trimorphic, third type can coexist (positive invasion), is above threshold, two coexisting types already stored in pair
#			4  = Time to Finishline: traits seperated by more than critical threshold 
#			5  = Number of times the dimporphism is lost in phase 2 (in a single run, until reaching theshhold)
#			6  = Number of times the dimoprhism is lost in phase 1 (in a single run, until reaching theshhold)

rep <- 16
cpu <- 16
K <- 10

vary_sigma_c2 <- lapply(1:K, function(i){
	sigma_c2 <- seq(.1, .8, length=K)
	ensemble_sim(rep=rep, sigma_mu = 0.05, mu = 5e-4, sigma_c2 = sigma_c2[i], sigma_k2 = 1, ko = 500, xo = 0.1, threshold = 30, cpu=cpu) 
	})
	sigma_c2 <- seq(.1, .8, length=K)

plot_phases(vary_sigma_c2, sigma_c2)

vary_sigma_mu <- lapply(1:K, function(i){
	sigma_mu  <- seq(.02, .1, length=K)
	ensemble_sim(rep=rep, sigma_mu = sigma_mu[i], mu = 5e-4, sigma_c2 = 0.3, sigma_k2 = 1, ko = 500, xo = 0.1, threshold = 30, cpu=cpu) 
	})
	sigma_mu  <- seq(.02, .1, length=K)
plot_phases(vary_sigma_mu, sigma_mu)



vary_mu <- lapply(1:K, function(i){
	mu <- seq(1e-4, 1e-2, length=K)
	ensemble_sim(rep=rep, sigma_mu = 0.05, mu =mu[i], sigma_c2 = 0.3, sigma_k2 = 1, ko = 500, xo = 0.1, threshold = 30, cpu=cpu) 
	})
	mu <- seq(1e-4, 1e-2, length=K)
plot_phases(vary_mu, mu)



save(list=ls(), file="waiting_time2.Rdat")

system("git add waiting_time2.Rdat")
gitcommit()
log <- gitlog()
system("git push")
tweet("Finished and saved data:", tags=tags, commit=log$commitID)

