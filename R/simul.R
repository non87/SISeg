# Contains the sampling functions. One per framework


r_ig <- function(n_1, n_2, p_1 = F, p_2 = F, dir_1 = F, dir_2 = F){
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
  return(matrix(c(observed_1, observed_2), byrow = TRUE, nrow = 2))
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
#' env <- matrix(c(1,2,3,4,5,6,7,8), nrow = 2, byrow = TRUE)
#' n <- sum(env)
#' p <- c(t(env))
#' samp <- r_full_multinomial(n = n, p = p)
#' samp2 <- r_full_multinomial(n = n, dir_alpha = p)
#' all.equal(sum(env), sum(samp), sum(samp2))
r_fm <- function(n, p = F, dir_alpha = F){
  # We use all.equal to enforce some tolerance
  if (length(p) == 1){
    if (length(dir_alpha) > 1){
      p <-gtools::rdirichlet(1, dir_alpha)
    } else {
      stop("If you do not specify p, you need to specify dir")
    }}
  observed <- stats::rmultinom(1, n, p)
  return(matrix(observed, byrow = TRUE, nrow = 2))
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
#' env <- matrix(c(1,2,3,4,5,6,7,8), nrow = 2, byrow = TRUE)
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
r_iu <- function(ns, ps = F, alfas = F, betas = F){
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
  return(matrix(c(g1_observed, g2_observed), byrow = TRUE, ncol = length(ns)))
}


#Resampling schemes for bootstrap
#If bay is true, they will use bayesian resampling.
#IN this case the output is a list with 1. a matrix of probabilities (one sample per row)
#and 2. a matrix of drawn tables based on the probabilities in 1.
#If bay is false, then the output is simply a matrix of drawn tables

#' Simulate Independent Group Framework
#'
#' \code{r_ig} samples from one independent multinomial distribution per
#' group or, following a simple conjugate model, one Dirichlet-
#' multinomial per group. Returns \code{n_b} simulated samples in the form of a
#' list. Will also include simulated probability vector when Bayesian option is
#' activated.
#'
#' Simulate one sample from one independent distribution per group. The
#' user can sample from independent multinomial or independent
#' Dirichlet-multinomial. The multinomial distributions will take as
#' probability vectors the MLE based on the input environment -- that is, the
#' simple proportions. By default, the Dirichlet-multinomial distributions will
#' use a Jeffrey's prior. But the user can input a different prior as a list
#' of k-long vectors, where k is \code{ncol(env)}.
#'
#' @export
#'
#' @param env nxk \code{matrix} of \code{numeric}.
#' @param n_b \code{numeric}. How many simulations should be returned?
#' @param bay \code{logical} Should the Bayesian model be used?
#' @param prior \code{list} of probability vectors for the Dirichlet
#'      distribution. Will be ignoder if \code{bay} is FALSE.
#'
#' @return A \code{matrix} or a \code{list}. If bay is FALSE, the return is
#'     a matrix. Each row represents a sampled environment (ordered by row).
#'     If bay is TRUE, the function returns
#'     a list of matrices where the first element is the matrix of generated
#'     samples and the second element is a matrix of probability vectors used
#'     in the dirichlet-multinomial sampling and generated from
#'     the posterior Dirichlet distribution.
#'
#' @family frameworks
#' @examples
#' env <- matrix(c(1,2,3,4,5,6,7,8), nrow = 2, byrow = TRUE)
#' samp <- r_ig(env, n_b = 1000)
#' samp2 <- r_ig(env, n_b = 1000, bay = TRUE)
r_ig <- function(env, n_b, bay = FALSE, prior = NULL){
  # Much of the convultion in this function is to generalize to more than 2 rows
  ns <- rowSums(env)
  if (!(bay)){
    resample <- lapply(seq(nrow(env)), function(x) stats::rmultinom(n_b, ns[x],
                                                             env[x,]/ns[x]))
    #THis is another horrible one-liner with a sapply inside a mapply
    resample <- mapply(function(x, y) resample[[x]][,y],
                rep(seq(nrow(env)), n_b), c(sapply(seq(n_b),
                function(x) rep(x, nrow(env)))))
  } else {
    ns_p <- ns/sum(ns)
    # Uninformative jeffrey's prior
    if (is.null(prior)){ prior <- lapply(seq(nrow(env)), function(x)
      rep(1/ncol(env), ncol(env))) }
    post <- lapply(seq(nrow(env)), function(x) env[x,] + prior[[x]])
    dirs <- lapply(post, function(x) gtools::rdirichlet(n_b, x))
    #Here we have to use the usual trick of sapply inside mapply
    resample <- mapply(function(x,y) stats::rmultinom(1, ns[x], dirs[[x]][y,]),
                            rep(seq(nrow(env)), n_b),
                       c(sapply(seq(n_b), function(x) rep(x, nrow(env)))))
    dirs <- mapply( function(x,y) dirs[[x]][y,], rep(seq(nrow(env)), n_b),
                    c(sapply(seq(n_b), function(x) rep(x, nrow(env)))))
    dirs <- t(t(dirs)*ns_p)
    dirs <- c(dirs)
    dirs <- matrix(dirs, byrow = T, nrow = n_b)
  }
  resample <- c(resample)
  resample <- matrix(resample, byrow = T, nrow = n_b)
  # We pass the sampled probability vectors alongside the
  # resampled environments
  if (bay){ resample <- list(resample = resample, ps = dirs) }
  return(resample)
}



#For full-multi bootstrap
boot_fm <- function(mat, n_b, bay = FALSE, prior = FALSE){
  mat <- unlist(lapply(seq(nrow(mat)), function(x) mat[x,]))
  N <- sum(mat)
  if (!(bay)){
    resample <- rmultinom(n_b, N, mat/N)
    resample <- t(resample)
  } else {
    if (!(prior)){ prior <- rep(1/length(mat), length(mat)) }
    post <- prior + mat
    dirs <- rdirichlet(n_b, post)
    resample <- matrix(unlist(lapply(seq(nrow(dirs)), function(x)  rmultinom(1, N, dirs[x,]) ) ),
                       byrow = T, nrow = n_b)
    resample <- list(resample = resample, ps = dirs)
  }
  return(resample)
}

#For the multiple-binomial sampling. ws is a vector of weights for each column
#The bayesian estimation can opt to a beta-mixture model instead of a saturated model
#TO do that, the prior must be given as a list, which must have a 'beta_mix' element alongside
#a 'alfas' and 'betas' element. The length of the three components must be the same
#That is, a prior for each beta distribution of the mix and a prior for the mixture!
#If no prior is specified, the function will default to a saturated model
#Only works with 2 rows
boot_mb <- function(mat, n_b, ws = FALSE, bay = FALSE, prior = FALSE){
  if (!(sum(ws))){ws <- colSums(mat)}
  if (!(bay)){
    #Here we sample, the rest of the lines are to put the sample in a nicer format
    resample <- lapply(seq(ncol(mat)) , function(x) rbinom(n_b, ws[x], mat[1,x]/ws[x]) )
    resample <- mapply(function(x, y) resample[[x]][y], x = rep(seq(ncol(mat)), n_b), y = c(sapply(seq(n_b), function(x) rep(x, ncol(mat)))) )
    resample <- matrix(resample, nrow = n_b , byrow = T)
  } else{
    #This case estimates a beta binomial model with a Gibbs sampling
    if ('beta_mix' %in% names(prior)) {
      mix = beta_mix_prior[['beta_mix']]
      alfas <- beta_mix_prior[['alfas']]
      betas <- beta_mix_prior[['betas']]
      dir_prior <- rep(1/length(mix), length(mix))
      rbetas <- estimate_bin_bay(n_b, mat, alfas, betas, burn_out, mix = mix, dir_prior = dir_prior)
      #We are only really interested into the probabilities coming out of the estimation
      #We reshape this to ease the commands later
      rbetas <- c(rbetas$Pas)
      #This case goes with a simpler saturated binomial model
    } else {
      #Uninformative Jeffrey prior on all alfas and betas
      if (!(prior)) {
        prior <- list()
        prior$alfas <- rep(1/nrow(mat), ncol(mat))
        prior$betas <- rep(1/nrow(mat), ncol(mat))
      }
      post <- list()
      post$alfas <- sapply(seq(ncol(mat)), function(x) prior$alfas[x] + mat[1,x])
      post$betas <- sapply(seq(ncol(mat)), function(x) prior$betas[x] + mat[2,x])
      #After the posterior for the betas, we get the actual probability draw
      rbetas <- unlist(lapply(seq(ncol(mat)), function(x) rbeta(n_b, post$alfas[x], post$betas[x])))
    }
    #Now we draw the various binomial
    resample <- mapply( function(x,y)  rbinom(1, x, y), c(sapply(ws, function(x) rep(x, n_b) )) , y = rbetas )
    resample <- matrix(resample, nrow = n_b)
    #This is for reworking the probability into a data-like matrix
    ws_p <- ws/sum(ws)
    rbetas <- matrix(rbetas, nrow = n_b, byrow = F)
    rbetas <- unlist(lapply(seq(n_b), function(x) rbetas[x,]*ws_p))
    #Finally, we pack the beta draws in a way similar to a data-matrix
    rbetas <- matrix(rbetas, nrow = n_b, byrow = T)
    be_diff <- lapply(seq(ncol(rbetas)), function(x) ws_p[x] - rbetas[,x])
    be_diff <- mapply(function(x, y) be_diff[[x]][y], x = rep(seq(ncol(mat)), n_b), y = c(sapply(seq(n_b), function(x) rep(x, ncol(mat)))) )
    be_diff <- matrix(be_diff, nrow = n_b , byrow = T)
    rbetas <- cbind(rbetas, be_diff)
  }
  #We only have the first line of the matrix. This provides also the second
  diff <- lapply(seq(ncol(resample)), function(x) ws[x] - resample[,x])
  diff <- mapply(function(x, y) diff[[x]][y], x = rep(seq(ncol(mat)), n_b), y = c(sapply(seq(n_b), function(x) rep(x, ncol(mat)))) )
  diff <- matrix(diff, nrow = n_b , byrow = T)
  resample <- cbind(resample, diff)
  if (bay){ resample <- list(resample = resample, ps = rbetas) }
  return(resample)
}

