#
#
#LogBF_ZS_Laplace <- function(R2, n, p) {
#	M = length(R2)
#	res <- R_LogBF_ZS_Laplace(
#		r2curr=R2, n=as.integer(n), dim=as.integer(p), nmodels=as.integer(M), logmarg=as.double(rep(0,M)))
#	return(res)
#}
#
#
#LogBF_ZS_QUAD <- function(R2, n, p) {
#	M = length(R2)
#	res <- R_LogBF_ZS_QUAD(
#		r2curr=R2, n=as.integer(n), dim=as.integer(p), nmodels=as.integer(M), logmarg=as.double(rep(0,M)))
#	return(res)
#}
#
#LogBF_GonN_Laplace <- function(R2, n, p, alpha) {
#	M = length(R2)
#	res <- R_LogBF_Hg_null_vect(
#		r2curr=R2, n=as.integer(n), dim=as.integer(p), nmodels=as.integer(M), 
#		logmarg=as.double(rep(0,M)), alpha=as.double(alpha), gpower=0)
#	return(res)
#}
#
#n = 100
#p = 5
#R2 = 0.75
#alpha = 3
#
#val1 = LogBF_ZS_Laplace(R2, n, p)
#val2 = LogBF_ZS_QUAD(R2, n, p)
#val3 = LogBF_GonN_Laplace(R2, n, p, alpha)
#
