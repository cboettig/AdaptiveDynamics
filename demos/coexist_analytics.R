# analytics 
require(BranchingTime)
source("../R/branching_time.R")
source("../R/analytics.R")
source("../R/coexist_time.R")

analytic1 <- coexist_analytics(sigma_c2 = .8, sigma_k2 = 1, ko = 500, xo = 0.2)
save(list=ls(), file="analytics.Rdat")
png("analytics1.png")
plot_contours(analytic1)
dev.off()


analytic2 <- coexist_analytics(sigma_c2 = .5, sigma_k2 = 1, ko = 500, xo = 0.2)
save(list=ls(), file="analytics.Rdat")
png("analytics2.png")
plot_contours(analytic2)
dev.off()

analytic3 <- coexist_analytics(sigma_c2 = .3, sigma_k2 = 1, ko = 500, xo = 0.2)
save(list=ls(), file="analytics.Rdat")
png("analytics3.png")
plot_contours(analytic3)
dev.off()

analytic4 <- coexist_analytics(sigma_c2 = .8, sigma_k2 = 1, ko = 5000, xo = 0.2)
save(list=ls(), file="analytics.Rdat")
png("analytics4.png")
plot_contours(analytic4)
dev.off()

analytic5 <- coexist_analytics(sigma_c2 = .3, sigma_k2 = 1, ko = 5000, xo = 0.2)
save(list=ls(), file="analytics.Rdat")
png("analytics5.png")
plot_contours(analytic5)
dev.off()


