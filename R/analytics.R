# analytics.R
MINF <- -.25 
INF <- .25
MAXRES <- 50

pars <- list(sigma_mu = 0.02, sigma_c2 = .1, sigma_k2 = 1, ko = 100, mu = 0.001, xo = 0.2)


## Functions as defined in writeup
K <- function(x, pars){ pars$ko*exp(-x^2/(2*pars$sigma_k2) ) }
C <- function(y,x, pars){ exp( -(y-x)^2/(2*pars$sigma_c2) ) }
S <- function(y,x, pars){ 1 - C(y,x, pars)*K(x, pars)/K(y, pars) }
M <- function(y,x, pars){exp(- (y-x)^2/(2*pars$sigma_mu^2) ) / (pars$sigma_mu*sqrt(2*pi) ) }
phi <- function(x, pars){ x*(pars$sigma_k2/pars$sigma_c2+1)/(pars$sigma_k2/pars$sigma_c2-1) }
psi <- function(x, pars){ x*(pars$sigma_k2/pars$sigma_c2-1)/(pars$sigma_k2/pars$sigma_c2+1) }

# We only want to calculate these once
reals <- seq(MINF, INF, length.out = MAXRES)
coexist_with_x <- function(x, pars){   reals[ (reals < phi(x,pars) | reals > psi(x,pars) ) ] }

P <- function(x, pars){
	y <- coexist_with_x(x, pars)
	deltaX <- y[2] - y[1]
	K(x,pars)*pars$mu*sum( sapply(y, function(y){ S(y,x,pars)*M(y,x,pars)*deltaX } ) )
}

P(.04, pars)


P_matrix 
S_x_M  <- function(x){ sapply(reals, function(y){ S(y,x, pars)*M(y,x,pars) } ) }
SxM_matrix <- sapply(reals, S_x_M)

P <- function(x, pars=pars){
	step <- reals[2]-reals[1]
	x_loc <- (reals > x & reals <= x + step)
	pars$mu*K(x,pars)*( sum(SxM_matrix[reals < phi(x,pars),x_loc ] ) + sum(SxM_matrix[reals > psi(x,pars), x_loc ] ) )
}

# canonical path
## names makes this easier to read, but may slow down execution.
times <- seq(0, 1e4, length.out=MAXRES);
f <- function(t, y, p)
  {
	list(-y * 0.5 * p$mu * p$sigma_mu^2 * K(y,p) / p$sigma_k2)
  }
require(odesolve)
meanpath <- lsoda(2,times,f, pars, rtol=1e-4, atol= 1e-6)
#plot(meanpath)
tstep = times[2] - times[1]

P_path <- P(meanpath[,2], pars) 
integral <- cumsum(P_path*tstep)

myPi <- function(Time){
	Time_loc <- which(times > Time & times < Time + tstep)
	P_path[Time_loc] * exp(-integral[Time_loc])
}
ti <- seq(1,1e3, len=50)
dist <- sapply(ti, myPi)
plot(ti, dist)
sum(dist*ti)/sum(dist)


