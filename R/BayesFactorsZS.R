
library(BAS)
library(gsl)
library(pracma)

trapint <- function(xgrid, fgrid) {
	ng <- length(xgrid)
	xvec <- xgrid[2:ng] - xgrid[1:(ng - 1)]
	fvec <- fgrid[1:(ng - 1)] + fgrid[2:ng]
	integ <- sum(xvec * fvec) / 2
	return(integ)
}


####################################################################################################

# Function for calculating the confluent hypergeometric function M on the log-scale

log.M.fun <- function(a,b,z)
{
	vv = c()
	vs = c()
	vv[1] = 0
	vs[1] = 1
	for (s in 1:1000) {
		vv[s+1] = log(a + s - 1) + log(z) - log(abs(b + s - 1)) - log(s)
		vs[s+1] = vs[s]*sign(b+s-1)
	}
	M = max(vv[vs==1])
	val = M + log(sum(vs*exp(vv - M)))

	return(list(vv=vv,vs=vs,val=val))
}

####################################################################################################
 
#log.U.fun.quad1 <- function(a,b,z,quadrule) 
#{
#	log.vals = (b - a - 1)*log(1 + quadrule$x/z)
#	M = max(log.vals)	
#	val = -a*log(z) - lgamma(a) + M +  log(sum(quadrule$w*exp(quadrule$logGammaAp1 + log.vals - M)))
#	return(val)
#}

####################################################################################################

log.U.fun.quad <- function(a,b,z,quadrule) 
{
	log.vals = (b - a - 1)*log(1 + quadrule$x/z) + (a - 1)*log(quadrule$x)
	M = max(log.vals)	
	val = -a*log(z) - lgamma(a) + M +  log(sum(quadrule$w*exp(quadrule$logGammaAp1 + log.vals - M)))
	return(val) 
}

####################################################################################################

#log.U.fun.quad3 <- function(a,b,z,M) {
#
#	A = z
#	B = -(b - z - 2)
#	C = -(a - 1)
#	disc = sqrt( B*B - 4*A*C )
#	val1 = (-B + disc)/(2*A)
#	val2 = (-B - disc)/(2*A)
#	x = max(c(val1,val2))
#	
#	h = -(a - 1)/(x^2) - (b - a - 1)/((1 + x)^2)
#	s = sqrt(-1/h)
#	
#	xg = seq(max(c(x-10*s,0)),x+30*s,,M)
#	fg = (xg^(a - 1))*((1 + xg)^(b - a - 1))*exp(-z*xg)
#	
	#plot(xg,fg,type="l")
#	
#	ng <- length(xg)
#	dvec <- xg[2:ng] - xg[1:(ng - 1)]
#	fvec <- fg[1:(ng - 1)] + fg[2:ng]
# 
#	integ <- sum(dvec * fvec) / 2
#	
#	vals =   log( integ  ) - lgamma(a)
#	
#	return(vals)
#} 

####################################################################################################

log.U.fun.Laplace <- function(a,b,z) 
{
	A = z
	B = -(b - z - 2)
	C = -(a - 1)
	disc = sqrt( B*B - 4*A*C )
	val1 = (-B + disc)/(2*A)
	val2 = (-B - disc)/(2*A)
	x = max(c(val1,val2))
	
	#EPS = 1.0E-8
	#if (x<EPS) {
	#	x = EPS
	#}
	
	#cat("x=",x,"\n")
	
	h = -(a - 1)/(x^2) - (b - a - 1)/((1 + x)^2)
	s2 = (-1/h)
	
	#cat("s2=",s2,"\n")
	#cat("a=",a,"\n")
	#cat("b - a - 1=",b - a - 1,"\n")
	
	log.f = (a - 1)*log(x) + (b - a - 1)*log(1 + x) - z*x
	
	if (FALSE) {
		s = sqrt(s2)
		xg = seq(max(c(x-6*s,0)),x+6*s,,1000)
		log.fg = (a - 1)*log(xg) + (b - a - 1)*log(1 + xg) - z*xg
	 
		plot(xg,exp(log.fg),type="l")
		points(x,exp(log.f))
		
		#print(x)
		#print(xg[which.max(log.fg)])
	}

	val = log.f + 0.5*log(2*pi*s2)  - lgamma(a) 
	
	#print(val)
	#ans <- readline()
	
	return(val)
}

####################################################################################################

log.U.fun <- function(a,b,z,method=c("Exact","Laplace","Quad"),quadrule) 
{
	if (a==1) {
		#val = log(hyperg_U(a, b, z))
		val = log.U.fun.quad(a,b,z, quadrule)
	} else {

		if (method=="Exact") {
			val = log(hyperg_U(a, b, z))
		}
		#if (method=="Asymptotic") {
		#	val = -a*log(z) + log(hyperg_2F0(a, a - b + 1, x=-1/z, give=FALSE, strict=TRUE))
		#}
		if (method=="Laplace") {
			val = log.U.fun.Laplace(a,b,z)
		}
		#if (method=="Quad1") {
		#	quadrule = my.gaussLaguerre(M, a = a - 1)
		#	val = log.U.fun.quad1(a,b,z, quadrule)
		#}
		if (method=="Quad") {
			val = log.U.fun.quad(a,b,z, quadrule)
		}
		#if (method=="Quad3") {
		#	val = log.U.fun.quad3(a,b,z, M)
		#}
	}
	#cat(method,"\n")
	#print(val)
	
	return(val)
}

####################################################################################################

# Returns the real roots of the equation
#  x^3 + a*x^2 + b*x + c = 0

solve.cubic <- function(a,b,c) {

	Q <- (a^2 - 3*b)/9
	R <- (2*(a^3) - 9*a*b + 27*c)/54
	
	if ((R^2)<=(Q^3)) {
		theta <- acos(R/sqrt(Q^3))
		x1 <- -2*sqrt(Q)*cos(theta/3) - a/3
		x2 <- -2*sqrt(Q)*cos((theta+2*pi)/3) - a/3
		x3 <- -2*sqrt(Q)*cos((theta-2*pi)/3) - a/3
		return(c(x1,x2,x3))
	} else {
		A <- - sign(R)*(abs(R) + sqrt(R^2 - Q^3))^(1/3)
		if (A!=0) {
			B <- Q/A
		} else {
			B <- 0
		}
		return( A + B - a/3 )
	}        
}

log.BF.ZS.Laplace <- function(vR2,n,vp)
{
	M <- length(vR2) 
	
	vals <- c()
	for (i in 1:M) {
		
		a = 0.5*(vp[i] - 1)
		b = 0.5*(n - vp[i] - 1)
		c = 0.5*(n - 1)
		d = 1/(1- vR2[i])
		e = -0.5*n
		
		ed = e*d
		
		A = (a*d + b*d - c*d + e + e*d)/ed
		B = (a + a*d + b - c*d + e)/ed
		C = a/ed
		
		x = solve.cubic(A,B,C)
		
		#print(res)
		#x = seq(0,0.01,,1000)
		
		log.f = a*log(x) + b*log(1+x) - c*log(1+d*x) + e*x
		
		#g = a/x + b/(1 + x) - c*d/(1 + d*x) + e	
		h = -a/(x*x) - b/((1+x)^2) + c*d*d/((1 + d*x)^2)
		
		s2 = -1/h
		s  = sqrt(s2)
		
	 
		val = log.f + 0.5*log(2*pi*s2) - c*log(1 - vR2[i]) + 0.5*log(n/2/pi)  
		vals[i] <- val
	}
		
	return(vals)
}

####################################################################################################

log.BF.ZellnerSiow.MC  = function(vR2, n, vp, N = 1.0E6)
{
	M = length(vR2)
	
	# Type checking
	if (length(vp)==1) {
		vp = rep(vp,M)
	} else {
		if (length(vp)!=M) {
			stop("The length of vR2 and vp should be equal or the length of vp should be 1")
		}
	}

	vx = 1/rgamma(N,0.5,0.5*n)
	vsigma2 = 1 - vR2
 
	#######################################
	
	log.f = (0.5*(n - vp - 1))%*%t(log(1 + vx)) 
	log.f = log.f - 0.5*(n-1)*log(1 + vsigma2%*%t(vx))
	
	vmaxVal = apply(log.f,1,max)
	
	vvals = vmaxVal  
	vvals = vvals + log(apply(exp(log.f - matrix(vmaxVal,M,N)),1,mean))

	return(vvals)
}	

####################################################################################################

my.gaussLaguerre = function (n, a = 0) 
{
	stopifnot(is.numeric(n), length(n) == 1, n >= 2)
	stopifnot(is.numeric(a), length(a) == 1, a >= 0)
	i <- 1:n
	d <- (2 * i - 1) + a
	b <- sqrt(i[1:(n - 1)] * ((1:(n - 1)) + a))
	B <- Diag(d) + Diag(b, 1) + Diag(b, -1)
	E <- eigen(B)
	L <- E$values
	V <- E$vectors
	inds <- order(L)
	x <- L[inds]
	V <- t(V[, inds])
	w <- V[, 1]^2
	logGammaAp1 = lgamma(a + 1)
	return(list(x = x, w = w, logGammaAp1=logGammaAp1))
}

####################################################################################################

log.BF.ZellnerSiow.trap  = function(R2, n, p, N = 1.0E5, R=8000)
{
	sigma2 = 1 - R2
	vx = 1/rgamma(N,0.5,0.5*n)
	vx = seq(1.0E-8,R,,N)
	
	log.f = 0.5*(n - p - 1)*log(1 + vx) - 0.5*(n - 1)*log(1 + vx*sigma2)
	log.f = log.f + 0.5*log(n/2/pi) - 1.5*log(vx) - 0.5*n/vx
	
	plot(vx,exp(log.f),type="l")
	
	M = max(log.f)
	xg = vx
    fg = exp(log.f - M)
	
	ng <- length(xg)
	dvec <- xg[2:ng] - xg[1:(ng - 1)]
	fvec <- fg[1:(ng - 1)] + fg[2:ng]
	
	integ <- sum(dvec * fvec) / 2
	
	val =  M +  log( integ  ) 
	
	return(val)
}	

####################################################################################################

log.BF.ZellnerSiow.quad  = function(vR2, n, vp, quadrule=gaussLaguerre(1000, a = 0))
{
	q = quadrule
	N = length(q$x)
	vx = q$x
	vw = q$w
	
	vsigma2 = 1 - vR2
	M = length(vR2)
	mW = matrix(vw,M,N,byrow=TRUE)
	vz = 0.5*n*vsigma2
	
	#######################################
	
	# Type checking
	if (length(vp)==1) {
		vp = rep(vp,M)
	} else {
		if (length(vp)!=M) {
			stop("The length of vR2 and vp should be equal or the length of vp should be 1")
		}
	}
	
	#######################################
	
	log.f = (0.5*(vp-1))%*%t(log(2*vx/n)) 
	log.f = log.f + (0.5*(n - vp - 1))%*%t(log(1 + 2*vx/n)) 
	log.f = log.f - 0.5*(n-1)*log(1 + (1/vz)%*%t(vx))
	log.f = log.f - (0.5*(n - 1)*log(vsigma2))%*%matrix(1,1,N)
	
	vmaxVal = apply(log.f,1,max)
	
	vvals = 0.5*log(n/2/pi) - log(n/2) 
	vvals = vvals + vmaxVal
	vvals = vvals + log(apply(mW*exp(log.f - matrix(vmaxVal,M,N)),1,sum))

	return(vvals)
}

####################################################################################################

log.BF.ZellnerSiow.integer  = function(R2, n, p, method=c("Exact","Laplace","Quad"),M=1000) 
{
	N = round((n - p - 1)/2)
	sigma2 = 1 - R2
	log.vals = rep(0,N+1)
	
	quadrule = my.gaussLaguerre(M, a = 0)
	
	for (k in 0:N) 
	{
		a = k + 0.5*(p + 1)
		b = k - N + 1.5
		z = 0.5*n*sigma2
		
		log.val = (b - 1)*log(sigma2)
		#cat(k,1,log.val,"\n")
		log.val = log.val + lchoose(N,k)
		#cat(k,2,log.val,"\n")		
		log.val = log.val + lgamma(a) 
		#cat(k,3,log.val,"\n")
		#cat(a,b,z,method,"\n")
		log.val = log.val + log.U.fun(a,b,z,method,quadrule) 
		#cat(k,4,log.val,"\n")
		
		log.vals[k+1] = log.val
		
		#cat(k,5,val,log.val,valU,"\n")
	}
	
	#print(log.vals)
	
	M = max(log.vals)
	val = M + log(sum(exp(log.vals - M))) - 0.5*log(pi)  +  0.5*log(n/2)
	 
	return(val)
}

####################################################################################################

