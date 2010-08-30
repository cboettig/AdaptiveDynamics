# coexist_demo.R
require(BranchingTime)
ens1 <- ensemble_coexistence(reps=1000, cpu=16, sigma_mu = 0.05, mu = 5e-3, sigma_c2 = .8, sigma_k2 = 1, ko = 500, xo = 0.1)
png("coexist_sim1.png")
ensemble_coexist_stats(ens1)
dev.off()

ens2 <- ensemble_coexistence(reps=1000, cpu=16, sigma_mu = 0.05, mu = 1e-3, sigma_c2 = .1, sigma_k2 = 1, ko = 500, xo = 0.1)
png("coexist_sim2.png")
ensemble_coexist_stats(ens2)
dev.off()

ens3 <- ensemble_coexistence(reps=1000, cpu=16, sigma_mu = 0.03, mu = 5e-4, sigma_c2 = .3, sigma_k2 = 1, ko = 500, xo = 0.1)
png("coexist_sim3.png")
ensemble_coexist_stats(ens3)
dev.off()

ens4 <- ensemble_coexistence(reps=1000, cpu=16, sigma_mu = 0.05, mu = 1e-4, sigma_c2 = .3, sigma_k2 = 1, ko = 500, xo = 0.1)
png("coexist_sim4.png")
ensemble_coexist_stats(ens4)
dev.off()

save(list=ls(), file="coexist.Rdat")

