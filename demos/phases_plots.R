colVars <- function(X) sapply(1:dim(X)[2], function(i) var(X[,i]))
library(Hmisc)

##
# plot the output of waiting_times.R
# all is one of vary_mu, vary_sigma_mu, etc

plot_phases <- function(all, parameter, xlab){


	sigma_c2 <- parameter
	K<- length(all)
	reps <- all[[1]]$reps

	phase1_attempts <- sapply(1:K, function(i)
		sapply(1:reps, function(j) all[[i]]$data[[j]]$dimorph_fails)
	)

	phase2_attempts <- sapply(1:K, function(i)
		sapply(1:reps, function(j) all[[i]]$data[[j]]$trimorph_fails)
	)

	phase1_times <-  sapply(1:K, function(i)
		sapply(1:reps, function(j) all[[i]]$data[[j]]$pos_time[1,1])
	)
	phase2_times <-  sapply(1:K, function(i)
		sapply(1:reps, function(j) all[[i]]$data[[j]]$pos_time[2,1])
	)
	phase3_times <-  sapply(1:K, function(i)
		sapply(1:reps, function(j) all[[i]]$data[[j]]$pos_time[3,1])
	)


	plot_err <- function(phase, xlab, add=F, col="black", ylab=ylab){
		cm1 <- colMeans(phase)
		cv1 <- sqrt(colVars(phase))
			errbar(parameter, cm1, cm1+cv1, cm1-cv1, xlab=xlab, ylab=ylab, add=add, col=col)
			lines(parameter, cm1)
	}

	plotme <- function(){
		plot_err(phase1_attempts, xlab, ylab="Num. of attempts")
		plot_err(phase2_attempts, add=T, col="blue")
	}

	
	plotme <- function(){
	plot_err(phase1_times, xlab)
	plot_err(phase2_times, xlab)
	plot_err(phase3_times, xlab)
	}

social_plot(plotme(),
		file="attempts.png", tags="adaptivedynamics", comment="Number of attempts from phase 1")
	
}
#social_plot(errbar(sigma_c2, cm1/cm2, cm1/cm2+(cv1+cv2), cm1/cm2-(cv1+cv2)),
#	file="sigma_c2.png", tags="adaptivedynamics", comment="ratio")





#var_phase1_attempts


