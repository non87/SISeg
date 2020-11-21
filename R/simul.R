# Contains the sampling functions. One per framework


#' Simulate Independent Group Framework
#'
#' \code{r_independent_group} samples from two independent multinomial
#' distribution or, following a simple conjugate model, two Dirichlet-
#' multinomial. Returns the sample in the form of a matrix. At least, one
#' among p_1 and dir_1 and one among p_2 and dir_2 must be specified.
#' The length of p_1 (dir_1) and p_2 (dir_2) must be identical.
#'
#' Simulate one sample from two independent distributions, one per group. The
#' user can sample from two independent multinomial, two independent
#' Dirichlet-multinomial, a multinomial for group 1 and Dirichlet-multinomial
#' for group 2, or vice-versa. The multinomial distributions will take as
#' probability vectors p_1 and p_2. The Dirichlet-multinomial distributions will
#' take as alpha vectors dir_1 and dir_2. In case p_1 (p_2) is specified,
#' dir_1 (dir_2) will be ignored and group 1 (2) will be sampled from a multinomial
#' distribution. p_1 and p_2 do not have to be normalized.
#'
#' At least, one among p_1 and dir_1 and one among p_2 and dir_2 must be
#' specified. The length of p_1 (dir_1) and p_2 (dir_2) must be identical.
#'
#' @export
#'
#' @param n_1 \code{numeric} representing the number of individuals from group 1.
#' @param n_1 \code{numeric} representing the number of individuals from group 2.
#' @param p_1 \code{vector} of probability for group 1.
#' @param p_2 \code{vector} of probability for group 2.
#' @param dir_1 \code{vector} the shape parameter of the Dirichlet-multinomial
#'    for group 1.
#' @param dir_2 \code{vector} the shape parameter of the Dirichlet-multinomial
#'    for group 2.
#'
#' @return \code{matrix} containing the generated sample
#' @family frameworks
#' @examples
#' env_ <- matrix(c(1,2,3,4,5,6,7,8), nrow = 2, byrow = T)
#' n_1 <- sum(env[1,])
#' n_2 <- sum(env[2,])
#' p_1 <- env[1,]
#' p_2 <- env[2,]
#' samp <- r_independent_group(n_1 = n_1, n_2 = n_2, p_1 = p_1, p_2 = p_2)
#' samp2 <- r_independent_group(n_1 = n_1, n_2 = n_2, dir_1 = p_1, dir_2 = p_2)
#' all.equal(rowSums(env), rowSums(samp), rowSums(samp2))
r_independent_group <- function(n_1, n_2, p_1 = F, p_2 = F, dir_1 = F, dir_2 = F){
  # There is nothing to sample if we only have one column
  if (length(p_1) == 1){
    if (length(dir_1) > 1){
      p_1 <- gtools::rdirichlet(1, dir_1)
    } else {
      stop("If you do not specify p_1, you need to specify dir_1")
    }}
  if (length(p_2) == 1){
    if (length(dir_2) > 1){
      p_2 <- gtools::rdirichlet(1, dir_2)
    } else {
      stop("If you do not specify p_2, you need to specify dir_2")
    }}
  observed_1 <- stats::rmultinom(1, n_1, p_1)
  observed_2 <- stats::rmultinom(1, n_2, p_2)
  return(matrix(c(observed_1, observed_2), byrow = T, nrow = 2))
}

#' Simulate Full Multinomial Framework
#'
#' \code{r_full_multinomial} samples from one
#' distribution or, following a simple conjugate model, a Dirichlet-
#' multinomial. Returns the sample in the form of a matrix. At least, one
#' among p and dir must be specified.
#'
#' Simulate one sample from from a multinomial distribution or
#' Dirichlet-multinomial distribution. The multinomial distribution will take as
#' probability vectors p. The Dirichlet-multinomial distribution will
#' take as alpha vectors dir_alpha. In case p is specified,
#' dir_alpha will be ignored and the sample will be simulated from a multinomial
#' distribution. p does not have to be normalized.
#'
#'
#' At least, one among p and dir_alpha must be specified.
#'
#' @export
#'
#' @param n \code{numeric} representing the number of individuals in the sample.
#' @param p \code{vector} of probability for the multinomial distribution.
#' @param dir_alpha \code{vector} the shape parameter of the Dirichlet-multinomial
#'    distribution.
#'
#' @return \code{matrix} containing the generated sample
#' @family frameworks
#' @examples
#' env_ <- matrix(c(1,2,3,4,5,6,7,8), nrow = 2, byrow = T)
#' n <- sum(env)
#' p <- c(t(env))
#' samp <- r_full_multinomial(n = n, p = p)
#' samp2 <- r_full_multinomial(n = n, dir_alpha = p)
#' all.equal(sum(env), sum(samp), sum(samp2))
r_full_multinomial <- function(n, p = F, dir_alpha = F){
  # We use all.equal to enforce some tolerance
  if (length(p) == 1){
    if (length(dir_alpha) > 1){
      p <-gtools::rdirichlet(1, dir_alpha)
    } else {
      stop("If you do not specify p, you need to specify dir")
    }}
  observed <- stats::rmultinom(1, n, p)
  return(matrix(observed, byrow = T, nrow = 2))
}

#This is used in the simulation to create random
#multiplicative-binomial samples. It implies nrow = 2
#weights: the # of trials for the various binomials
#ps: the probabilities for the binomials
#alfas: the alfas for the beta (in case you want beta-binomial)
#betas: the betas for the beta (in case you want beta-binomial)
#If you specify ps, alfas and betas are ignored
#ps and weights must be equal in length.
#alfas and betas must be of equal length.
#They can be either of length = 1 (all  probabilities from the same beta)
#or length = length(weights)

#' Simulate Independent-unit Framework
#'
#' \code{r_independent_unit} samples from k independent multinomial
#' distributions or, following a simple Bayesian  conjugate model, k beta-
#' binomial distributions -- where k is the number of units in an environment.
#' Returns the sample in the form of a 2xk matrix. At least, one
#' among ps and both alfas and betas must be specified.
#'
#' Simulate a sample from k independent distributions, one per unit. The
#' user can sample from independent binomial, or independent
#' beta-binomial. The kth binomial distributions will take as
#' probability vectors the kth element of ps. The kth beta-bininomial
#' distributions will take as alfa (beta) the kth element in alfas (betas) if
#' alfas is a \code{vector}. Otherwise all beta-binomial distributions will
#' have alfas (betas) as the alfa (beta) parameter, when alfas (betas) is
#' \code{numeric}.
#' In case ps is specified, alfas and betas will be ignored and
#' the units' composition will be sampled from binomial distributions.
#'
#' At least, one among ps and both alphas and betas must be specified.
#' The length of ns and ps (or alphas and betas) must be identical.
#'
#' @export
#'
#' @param ns A k-length \code{vector} containing the trial parameter for each
#'     binomial distribution.
#' @param ps A k-length \code{vector} containing the probability
#'    parameter for each binomial distribution.
#' @param alfas A k-length \code{vector} containing the alfa parameters for each
#'     beta-binomial distribution. Alternatively a \code{numeric} showing the
#'     alfa parameter for all beta-binomial distribution.
#' @param betas A k-length \code{vector} containing the beta parameters for each
#'     beta-binomial distribution. Alternatively a \code{numeric} showing the
#'     beta parameter for all beta-binomial distribution.
#'
#' @return A 2xk \code{matrix} containing the generated sample.
#' @family frameworks
#' @examples
#' env_ <- matrix(c(1,2,3,4,5,6,7,8), nrow = 2, byrow = T)
#' ns <- colSums(env)
#' ps <- env[1,]/ns
#' samp <- r_independent_unit(ns = ns, ps = ps)
#' alfas <- env[1,]
#' betas <- env[2,]
#' samp2 <- r_independent_unit(ns = ns, alfas = alfas, betas = betas)
#' alfas <- 2
#' betas <- 3
#' samp3 <- r_independent_unit(ns = ns, alfas = alfas, betas = betas)
#' all.equal(colSums(env), colSums(samp), colSums(samp2), colSums(samp3))
r_independent_unit <- function(ns, ps = F, alfas = F, betas = F){
  if (length(ps) == 1){
    # Beta and alfa cannot be 0 in any case. So this is a valid test
    if ((betas) & (alfas)){
      alfas <- rep(alfas, length(ns) - length(alfas) + 1)
      betas <- rep(alfas, length(ns) - length(alfas) + 1)
      ps <- unlist(lapply(seq(length(ns)), function(x) stats::rbeta(1, alfas[x]
                                                                  , betas[x])))
    } else {
      stop("If you do not specify ps, you need to specify alfas and betas for
           the beta-binomial model")
    }}
  g1_observed <- mapply(function(x, N) stats::rbinom(1, N, x), x = ps, N = ns)
  g2_observed <- ns - g1_observed
  return(matrix(c(g1_observed, g2_observed), byrow = T, ncol = length(ns)))
}
