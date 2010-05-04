# File: branching_time.R
# Author: Carl Boettiger <cboettig@gmail.com>
# License: GPL v3.0


branch_simulation <- function(sigma_mu = 0.05, mu = 1e-3, sigma_c2 = .1, sigma_k2 = 1, ko = 500, xo = 0.1, threshold = 30, seed = NULL){
	phasetime <- double(6);
	if(is.null(seed)) { 
		seed = runif(1)
	}
	out  <- .C( "_Z17branch_simulationPdS_S_S_S_S_S_PiS0_",
				as.double(sigma_mu), 
				as.double(mu), 
				as.double(sigma_c2), 
				as.double(sigma_k2), 
				as.double(ko), 
				as.double(xo), 
				as.double(phasetime),
				as.integer(seed),
				as.integer(threshold)
			)
	out[[7]]
}

branching_time <- function(reps = 100, sigma_mu = 0.02, mu = 0.005, sigma_c2 = .1, sigma_k2 = 1, ko = 500, xo = 0.5, threshold = 30, cpus = 2){
	require(snowfall)
	if (cpus > 1){ 
		sfInit(parallel=TRUE, cpus=cpus) 
	} else { 
		sfInit() 
	}
	sfLibrary(BranchingTime)
	seeds <- 1e9*runif(reps)
	out <- sfSapply(1:reps, function(i){ branch_simulation(sigma_mu, mu, sigma_c2, sigma_k2, ko, xo, seed=seeds[i], threshold = threshold) })
	sfStop()
	pars <- list(sigma_mu = sigma_mu, mu=mu, sigma_c2 = sigma_c2, sigma_k2 = sigma_k2, ko = ko, xo = xo, threshold = threshold)
	list(data = out, pars = pars, reps = reps)
}


