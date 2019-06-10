
library(BAS)
library(gsl)
library(appell)
library(pracma)

####################################################################################################

CalculateProbabilities <- function(marg, mGamma) {
  mod.prob <- exp(marg - max(marg)) / sum(exp(marg - max(marg)))
  var.prob <- t(mGamma) %*% mod.prob
  return(list(mod.prob = mod.prob, var.prob = var.prob))
}

####################################################################################################

log.BF.to.modelprobs <- function(vlog.VB,p,mprior=c("uniform","beta-binomial"),mpriorvec) 
{
	log.marg <- vlog.VB
	if (mprior=="beta-binomial") {
		a <- mpriorvec[1]
		b <- mpriorvec[2]
		log.marg <- log.marg + lbeta(a + vp, b + p - vp)
	}
	modelprobs <- exp(log.marg - max(log.marg)) / sum(exp(log.marg - max(log.marg)))
	return(modelprobs)
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

####################################################################################################

log.BF.cake <- function(vR2, n, vp) {
  vBIC <- n * log(1 - vR2) + vp * log(n)
  vlogBF <- -0.5 * vBIC
  vlogBF[vR2 < 1.0E-12] <- 0
  return(vlogBF)
}

####################################################################################################

log.BF.ZE <- function(vR2, n, vp) {
  a <- -0.75
  vb <- 0.5 * (n - vp - 5) - a
  vd <- 0.5 * vp + a
  inds <- which(((vd + 1) > 0) & ((vb + 1) > 0))
  vlogBF <- -(vb[inds] + 1) * log(1 - vR2[inds]) + lbeta(vd[inds] + 1, vb[inds] + 1) - lbeta(a + 1, vb[inds] + 1)
  vlogBF[-inds] <- -Inf
  vlogBF[vR2 < 1.0E-12] <- 0
  return(vlogBF)
}

####################################################################################################

log.hyperg_2F1 <- function(b, c, x) {
  val <- 0
  val <- val + log(c - 1)
  val <- val + (1 - c) * log(x)
  val <- val + (c - b - 1) * log(1 - x)
  val <- val + lbeta(c - 1, b - c + 1)
  val <- val + pbeta(x, shape1 = (c - 1), shape2 = (b - c + 1), log = TRUE)
  val[x == 0] <- 0
  return(val)
}

####################################################################################################

log.BF.g.naive <- function(vR2, n, vp) {
  library(gsl)
  a <- 3
  vlogBF <- log(a - 2) - log(vp + a - 2) + log(hyperg_2F1(0.5 * (n - 1), 1, 0.5 * (vp + a), vR2, give = FALSE, strict = TRUE))
  vlogBF[vR2 < 1.0E-12] <- 0
  return(vlogBF)
}

####################################################################################################

log.BF.g.safe1 <- function(vR2, n, vp) {
  a <- 3
  vlogBF <- log(a - 2) - log(vp + a - 2) + log.hyperg_2F1(0.5 * (n - 1), 0.5 * (vp + a), vR2)
  vlogBF[vR2 < 1.0E-12] <- 0
  return(vlogBF)
}

####################################################################################################

log.BF.g.safe2 <- function(vR2, n, vp) {
  a <- 3
  shape1 <- 0.5 * (vp + a - 2)
  shape2 <- 0.5 * (n - vp - a + 1)
  val <- -log(vR2)
  val <- val - log(1 - vR2)
  val <- val + pbeta(vR2, shape1 = shape1, shape2 = shape2, log = TRUE)
  val <- val - dbeta(vR2, shape1 = shape1, shape2 = shape2, log = TRUE)
  vlogBF <- log(a - 2) - log(2) + val
  vlogBF[vR2 < 1.0E-12] <- 0
  return(vlogBF)
}

####################################################################################################

# Equation 26 of Bayarri et al (2012) with gsl implementation of 2F1
log.BF.robust.naive1 <- function(vR2, n, vp) {
  vsigma2 <- 1 - vR2
  vlogBF <- -0.5 * vp * log(n + 1)
  vlogBF <- vlogBF + 0.5 * vp * log(vp + 1)
  vlogBF <- vlogBF - 0.5 * (n - 1) * log(vsigma2)
  vlogBF <- vlogBF - log(vp + 1)
  vlogBF <- vlogBF + log(hyperg_2F1(0.5 * (n - 1), 0.5 * (vp + 1), 0.5 * (vp + 3), (1 - 1 / vsigma2) * (vp + 1) / (n + 1), give = FALSE, strict = TRUE))
  vlogBF[vR2 < 1.0E-12] <- 0
  return(vlogBF)
}

####################################################################################################

# Equation 26 of Bayarri et al (2012) with BAS implementation of 2F1
log.BF.robust.naive2 <- function(vR2, n, vp) {
  vsigma2 <- 1 - vR2
  vlogBF <- -0.5 * vp * log(n + 1)
  vlogBF <- vlogBF + 0.5 * vp * log(vp + 1)
  vlogBF <- vlogBF - 0.5 * (n - 1) * log(vsigma2)
  vlogBF <- vlogBF - log(vp + 1)
  for (i in 1:length(vR2)) {
    vlogBF[i] <- vlogBF[i] + hypergeometric2F1(0.5 * (n - 1), 0.5 * (vp[i] + 1), 0.5 * (vp[i] + 3), (1 - 1 / vsigma2[i]) * (vp[i] + 1) / (n + 1), method = "Cephes", log = TRUE)
  }
  vlogBF[vR2 < 1.0E-12] <- 0
  return(vlogBF)
}

####################################################################################################

# Equation 26 of Bayarri et al (2012)
# with BAS implementation of 2F1 using Laplace method
log.BF.robust.naive3 <- function(vR2, n, vp) {
  vsigma2 <- 1 - vR2
  vlogBF <- -0.5 * vp * log(n + 1)
  vlogBF <- vlogBF + 0.5 * vp * log(vp + 1)
  vlogBF <- vlogBF - 0.5 * (n - 1) * log(vsigma2)
  vlogBF <- vlogBF - log(vp + 1)
  for (i in 1:length(vR2)) {
    vlogBF[i] <- vlogBF[i] + hypergeometric2F1(0.5 * (n - 1), 0.5 * (vp[i] + 1), 0.5 * (vp[i] + 3), (1 - 1 / vsigma2[i]) * (vp[i] + 1) / (n + 1), method = "Laplace", log = TRUE)
  }
  vlogBF[vR2 < 1.0E-12] <- 0
  return(vlogBF)
}

####################################################################################################

trapint <- function(xgrid, fgrid) {
  ng <- length(xgrid)
  xvec <- xgrid[2:ng] - xgrid[1:(ng - 1)]
  fvec <- fgrid[1:(ng - 1)] + fgrid[2:ng]
  integ <- sum(xvec * fvec) / 2
  return(integ)
}

####################################################################################################

# Using Equation 17 and quadrature
log.BF.robust.quad1 <- function(vR2, n, vp) {
  vr <- (1 + n) / (1 + vp)
  vL <- vr - 1
  vlogBF <- rep(0, length(vR2))
  for (i in 1:length(vR2))
  {
    vg <- seq(vL[i], 10000, , 10000)
    log.j <- -log(2) + 0.5 * log(vr[i]) + 0.5 * (n - vp[i] - 4) * log(1 + vg) - 0.5 * (n - 1) * log(1 + vg * (1 - vR2[i]))
    vlogBF[i] <- log(trapint(vg, exp(log.j)))
  }
  vlogBF[vR2 < 1.0E-12] <- 0
  return(vlogBF)
}

####################################################################################################

# Using Equation 17 after the transform x = g - L
log.BF.robust.quad2 <- function(vR2, n, vp) {
  vr <- (1 + n) / (1 + vp)
  vL <- vr - 1
  vx <- seq(0, 10000, , 10000)
  vsigma2 <- 1 - vR2
  vbeta <- (1 + vsigma2 * vL) / vsigma2
  vlogBF <- rep(0, length(vR2))
  for (i in 1:length(vR2))
  {
    log.j <- -log(2) + 0.5 * log(vr[i]) - 0.5 * (n - 1) * log(vsigma2[i]) + 0.5 * (n - vp[i] - 4) * log(vr[i] + vx) - 0.5 * (n - 1) * log(vbeta[i] + vx)
    vlogBF[i] <- log(trapint(vx, exp(log.j)))
  }
  vlogBF[vR2 < 1.0E-12] <- 0

  return(vlogBF)
}

####################################################################################################

# Equation 19 using gsl implementation of 2F1
log.BF.robust.naive4 <- function(vR2, n, vp) {
  vsigma2 <- 1 - vR2
  vr <- (1 + n) / (1 + vp)
  vL <- vr - 1
  vz <- (1 - vsigma2) / (1 + vsigma2 * vL)
  vlogBF <- 0.5 * (n - vp - 1) * log(vr)
  vlogBF <- vlogBF - log(1 + vp)
  vlogBF <- vlogBF - 0.5 * (n - 1) * log(1 + vL * vsigma2)
  vlogBF <- vlogBF + log(hyperg_2F1(0.5 * (n - 1), 1, 0.5 * (vp + 3), vz, give = FALSE, strict = TRUE))
  vlogBF[vR2 < 1.0E-12] <- 0
  return(vlogBF)
}

####################################################################################################

# Equation 19 using safe implementation 1 of 2F1 special case
log.BF.robust.safe1 <- function(vR2, n, vp) {
  vsigma2 <- 1 - vR2
  vr <- (1 + n) / (1 + vp)
  vL <- vr - 1
  vz <- (1 - vsigma2) / (1 + vsigma2 * vL)
  vlogBF <- 0.5 * (n - vp - 1) * log(vr)
  vlogBF <- vlogBF - log(1 + vp)
  vlogBF <- vlogBF - 0.5 * (n - 1) * log(1 + vL * vsigma2)
  vlogBF <- vlogBF + log.hyperg_2F1(0.5 * (n - 1), 0.5 * (vp + 3), vz)
  vlogBF[vR2 < 1.0E-12] <- 0
  return(vlogBF)
}

####################################################################################################

# Equation 19 using safe implementation 2 of 2F1 special case
log.BF.robust.safe2 <- function(vR2, n, vp) {
  vsigma2 <- 1 - vR2
  vr <- (1 + n) / (1 + vp)
  vL <- vr - 1
  vz <- (1 - vsigma2) / (1 + vsigma2 * vL)
  vlogBF <- rep(0, length(vR2))
  vlogBF <- 0.5 * (n - vp - 1) * log(vr)
  vlogBF <- vlogBF - log(2)
  vlogBF <- vlogBF - 0.5 * (n - 1) * log(1 + vL * vsigma2)
  shape1 <- 0.5 * (vp + 1)
  shape2 <- 0.5 * (n - vp - 2)
  vlogBF <- vlogBF - log(vz)
  vlogBF <- vlogBF - log(1 - vz)
  vlogBF <- vlogBF + pbeta(vz, shape1 = shape1, shape2 = shape2, log = TRUE)
  vlogBF <- vlogBF - dbeta(vz, shape1 = shape1, shape2 = shape2, log = TRUE)
  vlogBF[vR2 < 1.0E-12] <- 0
  return(vlogBF)
}

####################################################################################################
 