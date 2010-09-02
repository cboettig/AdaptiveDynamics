# coexist_demo.R
require(BranchingTime)


pairplot <- function(ens, analytic, i=1, png=FALSE){
	if(png) png(paste("coexist_sim",i, sep=""), 1200, 600)
	par(mfrow=c(1,2))
	ensemble_coexist_stats(ens, cex.axis=1.8, labcex=1.4, cex=1.8, lwd=2, main=paste("Sim:", ens$ko, "ko, ", ens$sigma_c2, "sigma_c2, ", ens$sigma_k2, "sigma_k2") )
	plot_contours(analytic, cex.axis=1.8, labcex=1.4, cex=1.8, lwd=2, main=paste("Analytic:", analytic$ko, "ko, ", analytic$sigma_c2, "sigma_c2, ", analytic$sigma_k2, "sigma_k2") )
	if(png) dev.off()

}

ens1 <- ensemble_coexistence(reps=1000, cpu=16, sigma_c2 = .5, sigma_k2 = 1, ko = 500, xo = 0.2)
analytic1 <- coexist_analytics(sigma_c2 = .5, sigma_k2 = 1, ko = 500, xo = 0.2)
pairplot(ens1, analytic1, i=1, png=T)

ens2 <- ensemble_coexistence(reps=1000, cpu=16, sigma_c2 = .8, sigma_k2 = 1, ko = 500, xo = 0.2)
analytic2 <- coexist_analytics(sigma_c2 = .8, sigma_k2 = 1, ko = 500, xo = 0.2)
pairplot(ens2, analytic2, i=2, png=T)

ens3 <- ensemble_coexistence(reps=1000, cpu=16, sigma_c2 = .3, sigma_k2 = 1, ko = 500, xo = 0.2)
analytic3 <- coexist_analytics(sigma_c2 = .3, sigma_k2 = 1, ko = 500, xo = 0.2)
pairplot(ens3, analytic3, i=3, png=T)


ens4 <- ensemble_coexistence(reps=1000, cpu=16,  sigma_c2 = .8, sigma_k2 = 1, ko = 1000, xo = 0.2)
analytic4 <- coexist_analytics(sigma_c2 = .8, sigma_k2 = 1, ko = 1000, xo = 0.2)
pairplot(ens4, analytic4, i=4, png=T)


ens5 <- ensemble_coexistence(reps=1000, cpu=16,  sigma_c2 = .3, sigma_k2 = 1, ko = 1000, xo = 0.2)
analytic5 <- coexist_analytics(sigma_c2 = .3, sigma_k2 = 1, ko = 1000, xo = 0.2)
pairplot(ens5, analytic5, i=5, png=T)

save(list=ls(), file="coexist.Rdat")

