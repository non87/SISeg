#### Bunch of functions that are used in the calculation of the delta method

#For HT, calculate the var-covar matrix from a multinomial
mult_vcv <- function(p, n){
  v <- diag(p)
  cv <- -p %*% t(p)
  return(n*(v + cv))
}

#For HT, calculate the empirical var-covar matrix
mult_vcv_est <- function(v, N){
  n <- sum(v)
  v <- v/sum(n)
  return(mult_vcv(v,N))
}


##Fast version of Matrix :: .bdiag() -- for the case of *many*  (k x k) matrices:
##lmat: list(<mat1>, <mat2>, ....., <mat_N>)  where each mat_j is a  k x k 'matrix'
##return a sparse (N*k x N*k) matrix of class  dgCMatrix.
bdiag_m <- function(lmat) {
  ## Copyright (C) 2016 Martin Maechler, ETH Zurich
  if(!length(lmat)){
    m <- Matrix:: sparseMatrix(c(), c(), dims=c(0,0))
    return(methods::as(m, "dgCMatrix"))
  }
  stopifnot(is.list(lmat), is.matrix(lmat[[1]]),
            (k <- (d <- dim(lmat[[1]]))[1]) == d[2], # k x k
            all(vapply(lmat, dim, integer(2)) == k)) # all of them
  N <- length(lmat)
  if(N * k > .Machine$integer.max)
    stop("resulting matrix too large; would be  M x M, with M=", N*k)
  M <- as.integer(N * k)
  ## result: an   M x M  matrix
  ## 'i :' maybe there's a faster way (w/o matrix indexing), but elegant?
  i <- as.vector(matrix(0L:(M-1L), nrow=k)[, rep(seq_len(N), each=k)]) + 1
  j <- unlist(lapply(seq(M), function(x) rep(x, k)))
  x <- as.double(unlist(lmat, recursive=FALSE, use.names=FALSE))
  Matrix::sparseMatrix(i = i, j = j, x = x, dims = c(M,M), index1 = TRUE)
}


#From a data matrix to a Sigma matrix of se
#For the dm sampling scheme
sigma_dm <- function(env){
  N <- sum(env)
  n_m <- sum(env[1,])
  n_f <- sum(env[2,])
  N <- sum(n_m, n_f)
  Ms_s <- list(mult_vcv_est(env[1,], 1) * N/n_m,  mult_vcv_est(env[2,], 1) * N/n_f)
  return(bdiag_m(Ms_s))
}

#For the fm sampling scheme
sigma_fm <- function(env){
  N <- sum(env)
  env <- c(env[1,], env[2,])
  return(mult_vcv(env/N, 1))
}

#For the mult.-binomial sampling scheme
sigma_mb <- function(env){
  Ms <- lapply(seq(ncol(env)), function(x) mult_vcv_est(c(env[1,x],env[2,x]), 1) )
  return(bdiag_m(Ms))
}
