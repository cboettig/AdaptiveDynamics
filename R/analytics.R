# analytics.R
MINF <- -10 
INF <- 10
MAXRES <- 1000

pars <- list(sigma_mu = 0.02, sigma_c2 = .4, sigma_k2 = 1, ko = 100, mu = 0.001, xo = 0.5)


## Functions as defined in writeup
K <- function(x, pars){ pars$ko*exp(-x^2/(2*pars$sigma_k2) ) }
C <- function(y,x, pars){ exp( -(y-x)^2/(2*pars$sigma_c2) ) }
S <- function(y,x, pars){ 1 - C(y,x, pars)*K(x, pars)/K(y, pars) }
M <- function(y,x, pars){exp(- (y-x)^2/(2*pars$sigma_mu^2) ) / (pars$sigma_mu*sqrt(2*pi) ) }
phi <- function(x, pars){ x*(pars$sigma_k2/pars$sigma_c2+1)/(pars$sigma_k2/pars$sigma_c2-1) }
psi <- function(x, pars){ x*(pars$sigma_k2/pars$sigma_c2-1)/(pars$sigma_k2/pars$sigma_c2+1) }

# We only want to calculate these once
reals <- seq(MINF, INF, length.out = MAXRES)
S_x_M  <- function(x){ sapply(reals, function(y){ S(y,x, pars)*M(y,x,pars) } ) }
SxM_matrix <- sapply(reals, S_x_M)

P <- function(x, pars=pars){
	step <- reals[2]-reals[1]
	x_loc <- (reals > x & reals <= x + step)
	pars$mu*K(x,pars)*( sum(SxM_matrix[x_loc, reals < phi(x,pars) ] ) + sum(SxM_matrix[x_loc, reals > psi(x,pars) ] ) )
}

# canonical path
## names makes this easier to read, but may slow down execution.
times <- seq(0, 1e4, length.out=MAXRES);
f <- function(t, y, p)
  {
	list(-y * 0.5 * p$mu * p$sigma_mu^2 * K(y,p) / p$sigma_k2)
  }
require(odesolve)
meanpath <- lsoda(.5,times,f, pars, rtol=1e-4, atol= 1e-6)
#plot(meanpath)
tstep = times[2] - times[1]

P_path <- P(meanpath[,2], pars) 
integral <- cumsum(P_path*tstep)

myPi <- function(Time){
	Time_loc <- which(times > Time & times < Time + tstep)
	P_path[Time_loc] * exp(-integral[Time_loc])
}
ti <- seq(1,1e3, len=100)
dist <- sapply(ti, myPi)
plot(ti, dist)
sum(dist*ti)/sum(dist)











Pi <- function(Time, pars, meanpath){
	times <- meanpath[,1]
	step = times[2] - times[1]
	Time_loc <- (times > Time & times < Time + step)
	integrand <- function(x_time){
		sapply(meanpath[times <= x_time,2], function(x){P(x, pars) } )
	}
	P(meanpath[Time_loc,2], pars) * exp(-integrate( integrand, 0, Time ) ) 
}
Pi(100, pars, meanpath)

prob <- sapply(times, function(t){ Pi(t, pars, meanpath) } ) 
plot(times, prob)


