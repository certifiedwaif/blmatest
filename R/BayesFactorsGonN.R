
library(BAS)
library(gsl)
library(appell)
library(pracma)

####################################################################################################

# Naive approach using direct call of appellf1

log.BF.g.on.n.naive <- function(vR2, n, vp, flag) 
{
  a <- 3
  vlogBF <- rep(0, length(vR2))
  errors = c()
  for (i in 1:length(vR2)) {
  
  	if (vp[i]==0) {
  		vlogBF[i] <- 0
  		errors[i] <- FALSE
  	} else {
	    vlogBF[i] <- log(a - 2) - log(n) - log(vp[i] + a - 2)
	    res = try({ 
	    	appellf1(1.0, 0.5*a, 0.5*(n-1), 0.5*(vp[i]+a), (1-1/n), vR2[i], debug=FALSE, userflag = flag, hyp2f1 = "michel.stoitsov")
	    },silent=TRUE)
	  	#print(res)
	  	
	  	if (class(res)=="try-error") {
	  		vlogBF[i] <- NA
	  		errors[i] <- TRUE
	  	} else {
	  		val = Re(res$val)
	  		errors[i] = (val<0)|(is.nan(val))
	  		vlogBF[i] <- vlogBF[i] + log(val)
	  	}
	  }
  }
  vlogBF[vR2 < 1.0E-12] <- 0
  return(list(vals=vlogBF,errors=errors))
}

####################################################################################################

log.BF.g.on.n.integrand <- function(vu, R2, n, p, a) {

	vals = 0
	vals = vals + log(a - 2)
	vals = vals - log(2*n)
	vals = vals + 0.5*(p + a - 4)*log(1 - vu)
	vals = vals - 0.5*a*log(1 - vu*(1 - 1/n))
	vals = vals - 0.5*(n - 1)*log(1 - vu*R2)
	
	return(vals)
}

####################################################################################################

log.BF.g.on.n.quad <- function(R2, n, p, a, quadrule) 
{
	log.f = log.BF.g.on.n.integrand(quadrule$x, R2, n, p, a)
	log.f.star = max(log.f)
	log.f.til = log.f - log.f.star
	val = log.f.star + log(sum(quadrule$w*exp(log.f.til)))	
	return(val)
}


log.BF.g.on.n.quad.vec <- function(vR2, n, vp, a, quadrule=gaussLegendre(n=1000,a=0,b=1)) 
{
   	vu <- quadrule$x
	vw <- quadrule$w

	vals = 0
	vals = vals + log(a - 2)
	vals = vals - log(2*n)
	vals = vals + 0.5*(vp + a - 4)%*%t(log(1 - vu))
	vals = vals - matrix(0.5*a,length(vR2),1)%*%t(log(1 - vu*(1 - 1/n)))
	vals = vals - 0.5*(n - 1)*log(1 - vR2%*%t(vu))
	
	log.f.star = max(vals)
	log.f.til = vals - log.f.star
	logBF = log.f.star + log( exp(log.f.til)%*%vw )	
	
	return(logBF)
}

####################################################################################################



log.f.fn = function(x,a,n,p,R2) 
{
	s2 = 1 - R2
	log.f = -0.5*a*log(1 + x/n)
	log.f = log.f + 0.5*(n - p - 1)*log(1 + x)
	log.f = log.f - 0.5*(n - 1)*log(1 + x*s2)
	return(log.f)
}

g.fn = function(x,a,n,p,R2) 
{
	s2 = 1 - R2
	g = -0.5*a/(n + x)
	g = g + 0.5*(n - p - 1)/(1 + x)
	g = g - 0.5*(n - 1)*s2/(1 + x*s2)
	return(g)
}

h.fn = function(x,a,n,p,R2)
{
	s2 = 1 - R2
	h = 0.5*a/((n + x)^2)
	h = h - 0.5*(n - p - 1)/((1 + x)^2)
	h = h + 0.5*(n - 1)*s2*s2/((1 + x*s2)^2)
	return(h)
}



solve.cubic <- function(a,b,c) 
{
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

####################################################################################################

log.BF.g.on.n.Laplace <- function(vR2, n, vp) 
{
	M = length(vR2)
	vals = c()
	
	a = 3
	
	for (m in 1:M) 
	{
		p = vp[m]
		R2 = vR2[m]
		
		
		A = a
		B = (n - p - 1)
		C = (n - 1)
		s = 1 - R2
		
		c2 = (A - B + C)*s
		c2 = a + p  
		c1 = A - B + A*s  - B*n*s + C*s*n + C*s
		c0 = A - B*n + C*s*n
		
		solp = (-c1 + sqrt(c1*c1 - 4*c2*c0))/(2*c2)
		#solm = (-c1 - sqrt(c1*c1 - 4*c2*c0))/(2*c2)
		
		x.star = max(c(0,solp))
		
		log.f.star = log.f.fn(x.star,a,n,p,R2) 
		sigma2 = -1/h.fn(x.star,a,n,p,R2) 
		
		if (sigma2<=0) {
			val = NA
		} else {
			val = log(a - 2) - log(2*n) + log.f.star + 0.5*log(2*pi*sigma2)
		}
		
		vals[m] = val
		
		
		#ans <- readline()
	}
	
	return(vals)
}


