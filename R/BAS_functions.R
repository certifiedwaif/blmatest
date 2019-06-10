

LogBF_ZS_Laplace <- function(R2, n, p) {
      M = length(R2)
      res <- R_LogBF_ZS_Laplace(
              r2curr=R2, n=as.integer(n), dim=as.integer(p), nmodels=as.integer(M), logmarg=as.double(rep(0,M)))
      return(res)
}


LogBF_ZS_QUAD <- function(R2, n, p) {
      M = length(R2)
      res <- R_LogBF_ZS_QUAD(
              r2curr=R2, n=as.integer(n), dim=as.integer(p), nmodels=as.integer(M), logmarg=as.double(rep(0,M)))
      return(res)
}

LogBF_GonN_Laplace <- function(R2, n, p, alpha) {
      M = length(R2)
      res <- R_LogBF_Hg_null_vect(
              r2curr=R2, n=as.integer(n), dim=as.integer(p), nmodels=as.integer(M), 
              logmarg=as.double(rep(0,M)), alpha=as.double(alpha), gpower=0)
      return(res)
}
