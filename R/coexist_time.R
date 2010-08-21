# File coexist_time.R
coexist_simulation <- function(sigma_mu = 0.03, mu = 1e-2, sigma_c2 = .5, sigma_k2 = 1, ko = 100, xo = 0.2, seed = NULL, threshold = 30, maxtime=5e4, samples=1e5){
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
				double(GRID^2),
				double(GRID^2),
				as.double(maxtime),
				as.integer(samples)
			)
	names(out) <- c("sigma_mu", "mu", "sigma_c2", "sigma_k2", "ko", "xo", "coexist_time", "seed", "threshold", "xval", "yval", "maxtime", "samples")
	class(out) <- "ad_coexist_times"
	out
	
}

ensemble_coexistence <-function(reps = 100, cpu=2, sigma_mu = 0.03, mu = 1e-2, sigma_c2 = .5, sigma_k2 = 1, ko = 100, xo = 0.2, seed = NULL, threshold = 30, maxtime=5e4, samples=1e5){
	require(snowfall)
	if(cpu > 1){
		sfInit(parallel=TRUE, cpu=cpu)
	} else { 
		sfInit() 
	}
	sfLibrary(BranchingTime)
	seeds <- 1e9*runif(reps)
	times <- sfSapply(	1:reps, 
						function(i){
							sim <- coexist_simulation(	sigma_mu=sigma_mu, mu=mu, sigma_c2 = sigma_c2, 
														sigma_k2 = sigma_k2, ko = ko, xo = xo, seed = seeds[i],
														threshold =threshold, maxtime=maxtime, samples=samples)
							sim$coexist_time
						})
	out <- list(times = times, reps=reps, sigma_mu=sigma_mu, mu=mu, sigma_c2 = sigma_c2, 
				sigma_k2 = sigma_k2, ko = ko, xo = xo, seed = seed, threshold =threshold, maxtime=maxtime, samples=samples)
}


ensemble_coexist_stats <- function( object, log=FALSE, nlevels=20 ){
	mean_times <- rowMeans(object$times)
	gridsize <- sqrt(dim(object$times)[1] )
	
	grid <- 20
	x <- seq(-object$xo, object$xo, length=grid)
	z <- matrix(mean_times, nrow = gridsize)
	if(log) z <- log(z)
	contour(x,x,z, lty=3, lwd=2, labcex=1.5, nlevels=nlevels)
	lines(x, bdry(x, object$sigma_k2, object$sigma_c2), lwd=3, col="darkblue" )
	lines(x, mirr(x, object$sigma_k2, object$sigma_c2), lwd=3, lty=1, col="darkblue" )
}

plot_contours <- function( object, log=FALSE ){
	gridsize <- sqrt(length(object$xval) )
	x <- seq(min(object$xval), max(object$xval), length = gridsize )
	z <- matrix(object$coexist_time, nrow = gridsize)
	if(log) z <- log(z)
	contour(x,x,z, lty=2, lwd=3)
	lines(x, bdry(x, object$sigma_k2, object$sigma_c2), lwd=3 )
	lines(x, mirr(x, object$sigma_k2, object$sigma_c2), lwd=3, lty=2 )
}


coexist_analytics <- function(sigma_mu = 0.03, mu = 1e-2, sigma_c2 = .1, sigma_k2 = 1, ko = 100, xo = 0.2){
	GRID = 20 # this is not optional, as it is not yet passed as a value but hardwired into the C code!!
	if(is.null(seed)) { 
		seed = runif(1)
	}
	out  <- .C( 
				"analytic_contours_wrapper",
				as.double(sigma_mu), 
				as.double(mu), 
				as.double(sigma_c2), 
				as.double(sigma_k2), 
				as.double(ko), 
				as.double(xo), 
				double(GRID^2),
				double(GRID^2),
				double(GRID^2)
			)
	names(out) <- c("sigma_mu", "mu", "sigma_c2", "sigma_k2", "ko", "xo", "coexist_time", "xval", "yval")
	class(out) <- "ad_coexist_times"
	out
	
}



