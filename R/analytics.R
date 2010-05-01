# analytics.R


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

# canonical path
f <- function(t, y, p)
  {
	list(-y * 0.5 * p$mu * p$sigma_mu^2 * K(y,p) / p$sigma_k2)
  }

waiting_time <- function(	reals = seq(-2, 2, length.out = 1000),
							times = seq(0, 5e4, length.out = 1000), 
							pars = list(sigma_mu = 0.03, mu = 0.01, sigma_c2 = .1, sigma_k2 = 1, ko = 100, xo = 0.5)
						){
	meanpath <- lsoda(pars$xo,times,f, pars, rtol=1e-4, atol= 1e-6)
	P_path <- sapply(meanpath[,2], function(x) P(x,pars, reals) )
	tstep = times[2] - times[1]
	integral <- cumsum(P_path*tstep)

	Pi <- function(Time){
		Time_loc <- which(times > Time & times < Time + tstep)
		P_path[Time_loc] * exp(-integral[Time_loc])
	}

	ti <- seq(1,max(times)/2, len=length(times)/2)
	dist <- sapply(ti, Pi)
	print( sum(dist*ti)/sum(dist) )
	dist
}

