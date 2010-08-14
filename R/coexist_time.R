# File coexist_time.R
coexist_simulation <- function(sigma_mu = 0.03, mu = 1e-2, sigma_c2 = .1, sigma_k2 = 1, ko = 100, xo = 0.4, seed = NULL, threshold = 30){
	GRID = 20
	if(is.null(seed)) { 
		seed = runif(1)
	}
	out  <- .C( 
				"coexist_simulation",
				as.double(sigma_mu), 
				as.double(mu), 
				as.double(sigma_c2), 
				as.double(sigma_k2), 
				as.double(ko), 
				as.double(xo), 
				double(GRID^2),
				as.integer(seed),
				as.integer(threshold),
				as.double(GRID^2),
				as.double(GRID^2)
			)
	names(out) <- c("sigma_mu", "mu", "sigma_c2", "sigma_k2", "ko", "xo", "coexist_time", "seed", "threshold", "xval", "yval")
	class(out) <- "ad_coexist_times"
	out
	
}




