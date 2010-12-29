	
##
# plot the output of waiting_times.R
# all is one of vary_mu, vary_sigma_mu, etc

plot_phases <- function(all, parameter){


sigma_c2 <- parameter
K<- length(all)
reps <- all[[1]]$reps

phase1_attempts <- sapply(1:K, function(i)
	sapply(1:reps, function(j) all[[i]]$data[[j]]$dimorph_fails)
)

phase2_attempts <- sapply(1:K, function(i)
	sapply(1:reps, function(j) all[[i]]$data[[j]]$trimorph_fails)
)

colVars <- function(X) sapply(1:dim(X)[2], function(i) var(X[,i]))


cm1 <- colMeans(phase1_attempts)
cv1 <- sqrt(colVars(phase1_attempts))
social_plot(
	errbar(sigma_c2, cm1, cm1+cv1, cm1-cv1),
	file="sigma_c2.png", tags="adaptivedynamics", comment="Number of attempts from phase 1")

cm2 <- colMeans(phase2_attempts)
cv2 <- sqrt(colVars(phase2_attempts))
social_plot(
	errbar(sigma_c2, cm2, cm2+cv2, cm2-cv2),
	file="sigma_c2.png", tags="adaptivedynamics", comment="Number of attempts from phase 2")
}
#social_plot(errbar(sigma_c2, cm1/cm2, cm1/cm2+(cv1+cv2), cm1/cm2-(cv1+cv2)),
#	file="sigma_c2.png", tags="adaptivedynamics", comment="ratio")





#var_phase1_attempts


