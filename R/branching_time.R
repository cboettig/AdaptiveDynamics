# File: branching_time.R
# Author: Carl Boettiger <cboettig@gmail.com>
# License: GPL v3.0


branch_simulation <- function(sigma_mu = 0.05, mu = 1e-3, sigma_c2 = .1, sigma_k2 = 1, ko = 1000, xo = 0.5, seed = NULL){
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

branching_time <- function(reps = 10, sigma_mu = 0.05, mu = 1e-3, sigma_c2 = .1, sigma_k2 = 1, ko = 1000, xo = 0.5, cpus = 2){
	require(snowfall)
	if (cpus > 1){ 
		sfInit(parallel=TRUE, cpus=cpus) 
	} else { 
		sfInit() 
	}
	sfLibrary(BranchingTime)
	loc <- system.file(package="BranchingTime")
	lib <- paste(loc, "/libs/BranchingTime.so", sep="")
	sfExport("lib")
	seeds <- 1e9*runif(reps)
	sfClusterEval(dyn.load(lib) )
	out <- sfSapply(1:reps, function(i){ branch_simulation(seed=seeds[i]) })
	sfStop()
	out
}




