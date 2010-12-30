# File: branching_time.R
# Author: Carl Boettiger <cboettig@gmail.com>
# License: GPL v3.0


branch_simulation <- function(sigma_mu = 0.03, mu = 1e-2, sigma_c2 = .1, sigma_k2 = 1, ko = 100, xo = 0.2, seed = NULL, threshold = 30, maxtime=1e6, samples=1e5){
	phasetime <- double(6);
	xpair <- double(4);
	ypair <- double(4);

	if(is.null(seed)) { 
		seed = runif(1)
	}
	out  <- .C( 
				"branch_simulation",
				as.double(sigma_mu), 
				as.double(mu), 
				as.double(sigma_c2), 
				as.double(sigma_k2), 
				as.double(ko), 
				as.double(xo), 
				double(6),
				as.integer(seed),
				as.integer(threshold),
				double(4),
				double(4),
				as.double(maxtime),
				as.integer(samples)
			)
	names(out) <- c("sigma_mu", "mu", "sigma_c2", "sigma_k2", "ko", "xo", "phasetime", "seed", "threshold", "xpair", "ypair", "maxtime", "samples")
	class(out) <- "ad_simulation"
	out
	#phasetime <- out[[7]]  # returns phasetime, a vector of length 6.  
	#xpair <- out[[10]]
	#ypair <- out[[11]]

# phasetime[1] = First time two coexisting branches are established
#			2  = Last (most recent time) two coexisting branches were established 
#			3  = Time to finish phase 2: (satisfies invade_pair() test: trimorphic, third type can coexist (positive invasion), is above threshold, two coexisting types already stored in pair
#			4  = Time to Finishline: traits seperated by more than critical threshold 
#			5  = Number of times the dimporphism is lost in phase 2 (in a single run, until reaching theshhold)
#			6  = Number of times the dimoprhism is lost in phase 1 (in a single run, until reaching theshhold)


	
}

ensemble_sim <- function(reps = 100, sigma_mu = 0.03, mu = 1e-2, sigma_c2 = .1, sigma_k2 = 1, ko = 100, xo = 0.4, threshold = 30, cpus = 2){
	require(snowfall)
	if (cpus > 1){ 
		sfInit(parallel=TRUE, cpus=cpus) 
	} else { 
		sfInit() 
	}
	sfLibrary(BranchingTime)
	seeds <- 1e9*runif(reps)
	out <- sfLapply(1:reps, 
					function(i){
						sim<-branch_simulation(sigma_mu, mu, sigma_c2, sigma_k2, ko, xo, seed=seeds[i], threshold = threshold) 
						list(	pos_time=data.frame(times=sim$phasetime[1:4], x=sim$xpair[1:4], y=sim$ypair[1:4]), 
								dimorph_fails = sim$phasetime[6], # "sigma_limited"?  Fails before coexistence establishes successful invader
								trimorph_fails = sim$phasetime[5]	  # "mu_limited" ? Fails after coexistence extablishes successful invader
							)
					}
			)
	sfStop()
	pars <- list(sigma_mu = sigma_mu, mu=mu, sigma_c2 = sigma_c2, sigma_k2 = sigma_k2, ko = ko, xo = xo, threshold = threshold)
	output <- list(data = out, pars = pars, reps = reps)
	class(output) <- "ad_ensemble"
	output
}


plot_waitingtimes <- function(ensemble, differences=FALSE, max_time=NULL, HISTOGRAM=FALSE){
#reorganize the data as matrix of nreps by 3, where rows are: time, xval, yval
# repeat for each stage in branching
	pts <- vector(mode="list", length=4)
	pts[[1]] <- sapply(1:ensemble$reps, function(i) ensemble$data[[i]]$pos_time[1,])
	pts[[2]] <- sapply(1:ensemble$reps, function(i) ensemble$data[[i]]$pos_time[2,])
	pts[[3]] <- sapply(1:ensemble$reps, function(i) ensemble$data[[i]]$pos_time[3,])
	pts[[4]] <- sapply(1:ensemble$reps, function(i) ensemble$data[[i]]$pos_time[4,])

	# use difference in times, not absolute time of event
	if(differences) for( i in 2:4) 	pts[[i]][1,] <- as.list( as.numeric(pts[[i]][1,]) - as.numeric(pts[[i-1]][1,]) )

	if(is.null(max_time)) max_time <- max( unlist( pts[[1]][1,]),  unlist( pts[[2]][1,]), unlist( pts[[3]][1,]),  unlist( pts[[4]][1,]) )

	if(!HISTOGRAM){
		plot(density(unlist( pts[[1]][1,] )), xlab='Time', lwd=3, lty=1, col='black', xlim=c(0,max_time), main="Waiting Time Distributions")
		lines(density(unlist( pts[[2]][1,] )), lwd=3, lty=2, col='purple')
		lines(density(unlist( pts[[3]][1,] )), lwd=3, lty=3, col='blue')
		lines(density(unlist( pts[[4]][1,] )), lwd=3, lty=4, col='green')
		legend(	'topright', 
				c("First dimorphism established", "Most recent dimorphism", "Dimorphism invaded", "Finishline"), 
				col=c("black", "purple", "blue", "green"), lty = c(1,2,3,4) )

	}

	if(HISTOGRAM){
		par(mfrow=c(2,2) )
		hist(unlist( pts[[1]][1,] ), xlab='Time', border='gray', col='gray', xlim=c(0,max_time), main="Frist dimorphism established")
		hist(unlist( pts[[2]][1,] ), add=F, border='violet', col='violet', xlim=c(0,max_time), main="Most recent dimorphism", xlab='Time')
		hist(unlist( pts[[3]][1,] ),add=F, border='lightblue', col='lightblue', xlim=c(0,max_time), main="Dimophism invaded", xlab='Time')
		hist(unlist( pts[[4]][1,] ),add=F, border='lightgreen', col='lightgreen', xlim=c(0,max_time), main="Finishline", xlab='Time' )
	}
}


# Histogram of failure rates
plot_failures <- function(ensemble, stage=0 ){
	fails <- sapply(1:ensemble$reps, function(i) c(ensemble$data[[i]]$dimorph_fails, ensemble$data[[i]]$trimorph_fails) )

	if(stage==2){ 
		hist(fails[1,], xlab="Number of attempts", main="Fails while waiting for mutation", col="darkred", border="white", cex.main=.8  ) 
	}else if(stage==3){ 
		hist(fails[2,], xlab="Number of attempts", main="Fails after mutation establishes", col="darkred", border="white", cex.main=.8  ) 
	} else {
		par(mfrow=c(1,2) )
		hist(fails[1,], xlab="Number of attempts", main="Fails while waiting for mutation", col="darkred", border="white", cex.main=.8  ) 
		hist(fails[2,], xlab="Number of attempts", main="Fails after mutation establishes", col="darkred", border="white", cex.main=.8  ) 
	}
}



bdry <- function(x, sigma_k2, sigma_c2){ -x* (1+ sigma_k2/sigma_c2)/(1-sigma_k2/sigma_c2) }
mirr <- function(x, sigma_k2, sigma_c2){ -x* (1- sigma_k2/sigma_c2)/(1+sigma_k2/sigma_c2) }


# Butterfly plot
plot_butterfly <- function(ensemble, k=2, differences=FALSE, lim=NULL, all=TRUE){
#reorganize the data as matrix of nreps by 3, where rows are: time, xval, yval
# repeat for each stage in branching
	pts <- vector(mode="list", length=4)
	pts[[1]] <- sapply(1:ensemble$reps, function(i) ensemble$data[[i]]$pos_time[1,])
	pts[[2]] <- sapply(1:ensemble$reps, function(i) ensemble$data[[i]]$pos_time[2,])
	pts[[3]] <- sapply(1:ensemble$reps, function(i) ensemble$data[[i]]$pos_time[3,])
	pts[[4]] <- sapply(1:ensemble$reps, function(i) ensemble$data[[i]]$pos_time[4,])

	# use difference in times, not absolute time of event
	if(differences) for( i in 2:4) 	pts[[i]][1,] <- as.list( as.numeric(pts[[i]][1,]) - as.numeric(pts[[i-1]][1,]) )



	if(all){ par(mfrow=c(2,2) ); plots=1:4 } else { plots = k }
	for(k in plots){
		col <- gray(as.numeric(pts[[k]][1,])/max(as.numeric(pts[[k]][1,])))

		if(is.null(lim) ){
			m <- max(c(abs(ensemble$pars$xo), unlist(pts[[4]][2,]), unlist(pts[[4]][3,]) ))
			lim <-c( -m, m )
		}
		plot(pts[[k]][2,], pts[[k]][3,], pch=16, col = col, cex=1.5, xlim = lim, ylim = lim, 
		xlab="resident, x", ylab="invader, y", main=paste("Position and time at collapse from ", k) )
		points(pts[[k]][2,], pts[[k]][3,], cex=1.5)
		xrange <- seq(lim[1], lim[2], length=100) 
		lines(xrange, bdry(xrange, ensemble$pars$sigma_k2, ensemble$pars$sigma_c2), lwd=3 )
		lines(xrange, mirr(xrange, ensemble$pars$sigma_k2, ensemble$pars$sigma_c2), lwd=3, lty=2 )
	}
}




