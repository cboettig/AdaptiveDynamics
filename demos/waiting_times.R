library(BranchingTime)

# phasetime[1] = First time two coexisting branches are established
#			2  = Last (most recent time) two coexisting branches were established 
#			3  = Time to finish phase 2: (satisfies invade_pair() test: trimorphic, third type can coexist (positive invasion), is above threshold, two coexisting types already stored in pair
#			4  = Time to Finishline: traits seperated by more than critical threshold 
#			5  = Number of times the dimporphism is lost in phase 2 (increments over ensemble runs)
#			6  = Number of times the dimoprhism is lost in phase 1  (increments over ensemble runs)



# Distribution of waiting times by phase.  Ensemble distribution of times observed for 1:4 from above.  
sim <- branching_time(rep=1000, cpu=8)
save(file = "run1.Rdat")
pdf("distribution.pdf")
m <- max(sim$data)
plot(density(sim$data[1,]), xlim = c(0,m), col="black" )
lines(density(sim$data[2,] ), col = "blue")
lines(density(sim$data[3,] ), col = "green")
lines(density(sim$data[4,] ), col = "red")
legend('topright', c("1", "2", "3", "4"), col=c("black", "blue", "green", "red"), pch = 15)
dev.off()


