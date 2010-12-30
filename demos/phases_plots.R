colVars <- function(X) sapply(1:dim(X)[2], function(i) var(X[,i]))

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


	plot_err <- function(phase, xlab){
		cm1 <- colMeans(phase)
		cv1 <- sqrt(colVars(phase))
		social_plot(
		errbar(parameter, cm1, cm1+cv1, cm1-cv1, xlab=xlab, ylab="time/frequency"),
		file="phase.png", tags="adaptivedynamics", comment="Number of attempts from phase 1")
	}
	plot_err(phase1_attempts, xlab)
	plot_err(phase2_attempts, xlab)
	plot_err(phase1_times, xlab)
	plot_err(phase2_times, xlab)
	plot_err(phase3_times, xlab)

}
#social_plot(errbar(sigma_c2, cm1/cm2, cm1/cm2+(cv1+cv2), cm1/cm2-(cv1+cv2)),
#	file="sigma_c2.png", tags="adaptivedynamics", comment="ratio")





#var_phase1_attempts


