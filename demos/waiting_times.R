library(BranchingTime)

# phasetime[1] = First time two coexisting branches are established
#			2  = Last (most recent time) two coexisting branches were established 
#			3  = Time to finish phase 2: (satisfies invade_pair() test: trimorphic, third type can coexist (positive invasion), is above threshold, two coexisting types already stored in pair
#			4  = Time to Finishline: traits seperated by more than critical threshold 
#			5  = Number of times the dimporphism is lost in phase 2 (in a single run, until reaching theshhold)
#			6  = Number of times the dimoprhism is lost in phase 1 (in a single run, until reaching theshhold)

rep <- 10
cpu <- 2

all <- vector(mode="list", length=4)

all[[1]] <- ensemble_sim(rep=rep, sigma_mu = 0.05, mu = 1e-3, sigma_c2 = .8, sigma_k2 = 1, ko = 500, xo = 0.1, threshold = 30, cpu=cpu)
all[[2]] <- ensemble_sim(rep=rep, sigma_mu = 0.05, mu = 1e-3, sigma_c2 = .1, sigma_k2 = 1, ko = 500, xo = 0.1, threshold = 30, cpu=cpu)
all[[3]] <- ensemble_sim(rep=rep, sigma_mu = 0.05, mu = 5e-4, sigma_c2 = .8, sigma_k2 = 1, ko = 500, xo = 0.1, threshold = 30, cpu=cpu)
all[[4]] <- ensemble_sim(rep=rep, sigma_mu = 0.05, mu = 1e-4, sigma_c2 = .3, sigma_k2 = 1, ko = 500, xo = 0.1, threshold = 30, cpu=cpu)

save(file = "run1.Rdat")

for(i in 1:4){
	png(file=paste("waitingtimes_ens_", i, ".png", sep=""))
	plot_waitingtimes(all[[i]])
	dev.off()

	png(file=paste("butterfly_ens_", i, ".png", sep=""))
	plot_butterfly(all[[i]])
	dev.off()
	
	png(file=paste("failures_ens_", i, ".png", sep=""))
	plot_failures(all[[i]])
	dev.off()
}
