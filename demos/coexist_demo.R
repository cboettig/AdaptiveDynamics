# coexist_demo.R
require(BranchingTime)
ens <- ensemble_coexistence(reps=100, cpu=16)
save(list=ls(), file="coexist.Rdat")
png("contour.png")
ensemble_coexist_stats(ens)
dev.off()
