# analytics.R
MINF <- -1 
INF <- 1
MAXRES <- 1000

pars <- list(sigma_mu = 0.01, sigma_c2 = .1, sigma_k2 = 1, ko = 100, mu = 0.01, xo = 0.5)


## Functions as defined in writeup
K <- function(x, pars){ pars$ko*exp(-x^2/(2*pars$sigma_k2) ) }
C <- function(y,x, pars){ exp( -(y-x)^2/(2*pars$sigma_c2) ) }
S <- function(y,x, pars){ 1 - C(y,x, pars)*K(x, pars)/K(y, pars) }
M <- function(y,x, pars){exp(- (y-x)^2/(2*pars$sigma_mu^2) ) / (pars$sigma_mu*sqrt(2*pi) ) }
phi <- function(x, pars){ x*(pars$sigma_k2/pars$sigma_c2+1)/(pars$sigma_k2/pars$sigma_c2-1) }
psi <- function(x, pars){ x*(pars$sigma_k2/pars$sigma_c2-1)/(pars$sigma_k2/pars$sigma_c2+1) }

# We only want to calculate these once
coexist_with_x <- function(x, pars, reals){   reals[ (reals < phi(x,pars) | reals > psi(x,pars) ) ] }

P <- function(x, pars, reals){
	y <- coexist_with_x(x, pars, reals)
	deltaX <- y[2] - y[1]
	K(x,pars)*pars$mu*sum( sapply(y, function(y){ S(y,x,pars)*M(y,x,pars)*deltaX } ) )
}
reals <- seq(MINF, INF, length.out = MAXRES)
P(.04, pars, reals)


# canonical path
## names makes this easier to read, but may slow down execution.
times <- seq(0, 5e4, length.out=MAXRES);
f <- function(t, y, p)
  {
	list(-y * 0.5 * p$mu * p$sigma_mu^2 * K(y,p) / p$sigma_k2)
  }
require(odesolve)
meanpath <- lsoda(pars$xo,times,f, pars, rtol=1e-4, atol= 1e-6)
#plot(meanpath)
tstep = times[2] - times[1]

P_path <- sapply(meanpath[,2], function(x) P(x,pars, reals) )
#plot(P_path)


integral <- cumsum(P_path*tstep)

myPi <- function(Time){
	Time_loc <- which(times > Time & times < Time + tstep)
	P_path[Time_loc] * exp(-integral[Time_loc])
}

myPi(500)

ti <- seq(1,5e3, len=MAXRES/2)
dist <- sapply(ti, myPi)
#plot(ti, dist)
print( sum(dist*ti)/sum(dist) )


