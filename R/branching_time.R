# File: branching_time.R
# Author: Carl Boettiger <cboettig@gmail.com>
# License: GPL v3.0


branch_simulation <- function(sigma_mu = 0.05, mu = 1e-3, sigma_c2 = .1, sigma_k2 = 1, ko = 500, xo = 0.1, seed = NULL){
	phasetime <- double(3);
	if(is.null(seed)) { 
		seed = runif(1)
	}
	out  <- .C( "_Z17branch_simulationPdS_S_S_S_S_S_Pi",
				as.double(sigma_mu), 
				as.double(mu), 
				as.double(sigma_c2), 
				as.double(sigma_k2), 
				as.double(ko), 
				as.double(xo), 
				double(3),
				as.integer(seed)
			)
	out[[7]]
}

branching_time <- function(reps = 10, sigma_mu = 0.05, mu = 1e-3, sigma_c2 = .1, sigma_k2 = 1, ko = 500, xo = 0.1, cpus = 2){
	require(snowfall)
	if (cpus > 1){ 
		sfInit(parallel=TRUE, cpus=cpus) 
	} else { 
		sfInit() 
	}
	sfLibrary(BranchingTime)
## sfLibrary taxes care of this since the package NAMESPACE calls useDynLib
#	loc <- system.file(package="BranchingTime")
#	lib <- paste(loc, "/libs/BranchingTime.so", sep="")
#	sfExport("lib")
#	sfClusterEval(dyn.load(lib) )
	seeds <- 1e9*runif(reps)
	out <- sfSapply(1:reps, function(i){ branch_simulation(sigma_mu, mu, sigma_c2, sigma_k2, ko, xo, seed=seeds[i]) })
	sfStop()
	out
}


analytic_distribution <- function(maximum = 1e4, minimum = 0, sigma_mu = 0.05, mu = 1e-3, sigma_c2 = .1, sigma_k2 = 1, ko = 500, xo = 0.1, n_pts = 100){
	times <- seq(minimum, maximum, length= n_pts);
	out <- .C("_Z9analyticsPdS_S_S_S_S_S_S_Pi",
		as.double(sigma_mu), 
		as.double(mu), 
		as.double(sigma_c2), 
		as.double(sigma_k2), 
		as.double(ko), 
		as.double(xo),
		as.double(times),
		double(n_pts),
		as.integer(n_pts)
	)
	list(y = out[[8]], x = times);
}
		

