# coexist_demo.R
require(BranchingTime)


ens1 <- ensemble_coexistence(reps=1000, cpu=16, sigma_c2 = .5, sigma_k2 = 1, ko = 500, xo = 0.2)
analytic1 <- coexist_analytics(sigma_c2 = .5, sigma_k2 = 1, ko = 500, xo = 0.2)
png("coexist_sim1.png")
par(mfrow=c(1,2))
ensemble_coexist_stats(ens1)
plot_contours(analytic1)
dev.off()

ens2 <- ensemble_coexistence(reps=1000, cpu=16, sigma_c2 = .8, sigma_k2 = 1, ko = 500, xo = 0.2)
analytic2 <- coexist_analytics(sigma_c2 = .8, sigma_k2 = 1, ko = 500, xo = 0.2)
png("coexist_sim2.png")
par(mfrow=c(1,2))
ensemble_coexist_stats(ens2)
plot_contours(analytic2)
dev.off()

ens3 <- ensemble_coexistence(reps=1000, cpu=16, sigma_c2 = .3, sigma_k2 = 1, ko = 500, xo = 0.2)
analytic3 <- coexist_analytics(sigma_c2 = .3, sigma_k2 = 1, ko = 500, xo = 0.2)
png("coexist_sim3.png")
par(mfrow=c(1,2))
ensemble_coexist_stats(ens3)
plot_contours(analytic3)
dev.off()


ens4 <- ensemble_coexistence(reps=1000, cpu=16,  sigma_c2 = .8, sigma_k2 = 1, ko = 1000, xo = 0.2)
png("coexist_sim4.png")
par(mfrow=c(1,2))
analytic4 <- coexist_analytics(sigma_c2 = .8, sigma_k2 = 1, ko = 1000, xo = 0.2)
ensemble_coexist_stats(ens4)
plot_contours(analytic4)
dev.off()


ens5 <- ensemble_coexistence(reps=1000, cpu=16,  sigma_c2 = .3, sigma_k2 = 1, ko = 1000, xo = 0.2)
png("coexist_sim5.png")
par(mfrow=c(1,2))
analytic5 <- coexist_analytics(sigma_c2 = .3, sigma_k2 = 1, ko = 1000, xo = 0.2)
ensemble_coexist_stats(ens5)
plot_contours(analytic5)
dev.off()

save(list=ls(), file="coexist.Rdat")

