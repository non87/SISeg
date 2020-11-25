##### Contains the sampling functions. One per framework


#' Simulate Independent Group Framework
#'
#' \code{r_ig} samples from one independent multinomial distribution per
#' group or, following a simple conjugate model, one Dirichlet-
#' multinomial per group. Returns \code{n_b} simulated samples in the form of a
#' list or a matrix. Will also include simulated probability vector when
#' Bayesian option is activated.
#'
#' \code{r_ig} simulates samples from k independent distributions, one per group.
#' The user can sample from independent multinomial distributions, or
#' independent Dirichlet-multinomial distributions.
#' If \code{bay} is \code{TRUE}, the function will sample from
#' Dirichlet-multinomial distributions following a simple conjugate model.
#' Whatever the case, It will sample \code{n_b} times and return the results
#' in a matrix or list format, depending on the \code{list_out} argument.
#' The user can specify a \code{prior} for the model, otherwise, the function
#' will default to a Jeffrey's prior.
#'
#' @export
#'
#' @param env nxk \code{matrix} of \code{numeric}.
#' @param n_b \code{numeric}. How many simulations should be returned?
#' @param bay \code{logical} Should the Bayesian model be used?
#' @param prior n-long \code{list} containing n probability vectors for the
#'      Dirichlet distributions. Will be ignored if \code{bay} is FALSE.
#' @param list_out \code{logical} should the output be in a list form?
#'
#' @return The return depends on the value of \code{list_out} and \code{bay}. By
#'     default, both are \code{FALSE}. In that case, the function returns a
#'     \code{matrix} where every row contains a simulated environment. Within
#'     each row, the values are ordered by row. In case of \code{list_out} set
#'     to \code{TRUE} and \code{bay} set to \code{FALSE}, the output will
#'     be a \code{n_b}-long list of \code{matrix} showing the simulated
#'     environments.
#'     When \code{bay} is \code{TRUE}, the output will always be a \code{list}.
#'     The first element of the list contains the simulated environments and
#'     the second element contains the Dirichlet samples used to create the
#'     simulated environment, scaled so that the total sum for each environment
#'     is one. The way these pieces of information are shaped depends on the
#'     value of \code{list_out}. If \code{list_out} is \code{FALSE},
#'     the information will be stored in \code{matrix}: one \code{matrix} for
#'     the simulated
#'     environment, one \code{matrix} for the Dirichlet samples. In the two
#'     matrices, every row represent a simulated environment, ordered by
#'     row. The two matrices are ordered so that the kth row in one corresponds
#'     to the kth row in the other. If \code{list_out} is \code{TRUE} (and
#'     \code{bay} is \code{TRUE}), the
#'     output will be shaped in a \code{list} of 2 lists. The first \code{list}
#'     contains the simulated environment. The second \code{list} contains the
#'     Dirichlet sample used to create the simulated environment. The lists
#'     are ordered so that the kth environment was created from the kth
#'     Dirichlet sample.
#'
#' @family frameworks
#' @examples
#' env <- matrix(seq(12), nrow = 3, byrow = TRUE)
#' samp <- r_ig(env, n_b = 1000, list_out = TRUE)
#' samp2 <- r_ig(env, n_b = 1000, bay = TRUE, list_out = TRUE)
#' # Very strong priors
#' priors <- list(c(0,0,30, 0))
#' priors <- rep(priors, nrow(env))
#' samp3 <- r_ig(env, n_b = 1000, list_out = TRUE, bay = TRUE, prior =
#' priors)
#' all.equal(rowSums(env), rowSums(samp[[1]]), rowSums(samp2[[1]][[1]]),
#'  rowSums(samp3[[1]][[1]]))
#' all.equal(rowSums(samp2[[2]][[1]]), rowSums(samp3[[2]][[1]]))
#' # How are the list and matrix output related?
#' set.seed(1234)
#' samp_list <- r_ig(env, n_b = 2, list_out = TRUE)
#' set.seed(1234)
#' samp_mat <- r_ig(env, n_b = 2, list_out = FALSE)
#' all.equal(samp_list[[2]], matrix(samp_mat[2,], byrow = TRUE,
#'   nrow = nrow(env)))
#' # For the Dirchlet samples is really the same
#' set.seed(1234)
#' samp_list <- r_ig(env, n_b = 2, list_out = TRUE, bay = TRUE)
#' set.seed(1234)
#' samp_mat <- r_ig(env, n_b = 2, list_out = FALSE, bay = TRUE)
#' all.equal(samp_list[[2]][[2]], matrix(samp_mat[[2]][2,], byrow = TRUE,
#'   nrow = nrow(env)))
r_ig <- function(env, n_b = 5000, bay = FALSE, prior = NULL, list_out = FALSE){
  # Much of the convultion in this function is to generalize to more than 2 rows
  ns <- rowSums(env)
  k <- ncol(env)
  if (!(bay)){
    resample <- lapply(seq(nrow(env)), function(x) stats::rmultinom(n_b, ns[x],
                                                             env[x,]/ns[x]))
    #THis is another horrible one-liner with a sapply inside a mapply
    resample <- mapply(function(x, y) resample[[x]][,y],
                x = rep(seq(nrow(env)), n_b),
                y = rep(seq(n_b), each = nrow(env)))
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
                       rep(seq(n_b), each= nrow(env)))
    dirs <- mapply( function(x,y) dirs[[x]][y,], rep(seq(nrow(env)), n_b),
                    rep(seq(n_b), each = nrow(env)))
    dirs <- t(t(dirs)*ns_p)
    dirs <- c(dirs)
    dirs <- matrix(dirs, byrow = T, nrow = n_b)
  }
  resample <- c(resample)
  resample <- matrix(resample, byrow = T, nrow = n_b)
  if (list_out){
    resample <- lapply(seq(nrow(resample)), function(x) matrix(resample[x, ],
                              ncol = k, byrow = T))
  }
  # We pass the sampled probability vectors alongside the
  # resampled environments
  if (bay){
    resample <- list(resample = resample, ps = dirs)
    if (list_out){
      resample[[2]] <- lapply(seq(nrow(resample[[2]])), function(x)
                              matrix(resample[[2]][x, ], ncol = k, byrow = T))
    }
  }
  return(resample)
}



#' Simulate Full Multinomial Framework
#'
#' \code{r_fm} samples an environment from one multinomial distribution or,
#' following a simple conjugate model, one Dirichlet-multinomial distribution.
#' Returns \code{n_b} simulated samples in the form of a
#' list or a matrix. If the Bayesian model is used, \code{r_fm} also returns
#' simulated the probability vectors from the Dirichlet distribution.
#'
#' \code{r_fm} simulates samples from one distribution for the entire
#' environment. The user can sample from a multinomial, or a Dirichlet-multinomial
#' distribution depending on the value of  \code{bay}.
#' If \code{bay} is \code{TRUE}, the function will sample from
#' Dirichlet-multinomial distributions following a simple conjugate model.
#' Whatever the case, \code{r_fm} will sample \code{n_b} times and return the
#' results in a matrix or list format, depending on the \code{list_out}
#' argument. The user can specify a \code{prior} for the model, otherwise, the
#' function will default to a Jeffrey's prior.
#'
#' @export
#'
#' @param env nxk \code{matrix} of \code{numeric}.
#' @param n_b \code{numeric}. How many simulations should be returned?
#' @param bay \code{logical} Should the Bayesian model be used?
#' @param prior n-long \code{list} containing n probability vectors for the
#'      Dirichlet distributions. Will be ignored if \code{bay} is FALSE.
#' @param list_out \code{logical} should the output be in a list form?
#'
#' @return The return depends on the value of \code{list_out} and \code{bay}. By
#'     default, both are \code{FALSE}. In that case, the function returns a
#'     \code{matrix} where every row contains a simulated environment. Within
#'     each row, the values are ordered by row. In case of \code{list_out} set
#'     to \code{TRUE} and \code{bay} set to \code{FALSE}, the output will
#'     be a \code{n_b}-long list of \code{matrix} showing the simulated
#'     environments.
#'     When \code{bay} is \code{TRUE}, the output will always be a \code{list}.
#'     The first element of the list contains the simulated environments and
#'     the second element contains the Dirichlet samples used to create the
#'     simulated environment, scaled so that the total sum for each environment
#'     is one. The way these pieces of information are shaped depends on the
#'     value of \code{list_out}. If \code{list_out} is \code{FALSE},
#'     the information will be stored in \code{matrix}: one \code{matrix} for
#'     the simulated
#'     environment, one \code{matrix} for the Dirichlet samples. In the two
#'     matrices, every row represent a simulated environment, ordered by
#'     row. The two matrices are ordered so that the kth row in one corresponds
#'     to the kth row in the other. If \code{list_out} is \code{TRUE} (and
#'     \code{bay} is \code{TRUE}), the
#'     output will be shaped in a \code{list} of 2 lists. The first \code{list}
#'     contains the simulated environment. The second \code{list} contains the
#'     Dirichlet sample used to create the simulated environment. The lists
#'     are ordered so that the kth environment was created from the kth
#'     Dirichlet sample.
#'
#' @family frameworks
#' @examples
#' env <- matrix(seq(12), nrow = 3, byrow = TRUE)
#' samp <- r_fm(env, n_b = 1000, list_out = TRUE)
#' samp2 <- r_fm(env, n_b = 1000, list_out = TRUE, bay = TRUE)
#' # Very strong priors
#' priors <- rep(0, nrow(env)*ncol(env))
#' priors[10] <- 50
#' samp3 <- r_fm(env, n_b = 1000, list_out = TRUE, bay = TRUE, prior =
#'  priors)
#' all.equal(sum(env), sum(samp[[1]]), sum(samp2[[1]][[1]]),
#'   sum(samp3[[1]][[1]]))
#' all.equal(sum(samp2[[2]][[1]]), sum(samp3[[2]][[1]]))
#' # How are the list and matrix output related?
#' set.seed(1234)
#' samp_list <- r_fm(env, n_b = 2, list_out = TRUE)
#' set.seed(1234)
#' samp_mat <- r_fm(env, n_b = 2, list_out = FALSE)
#' all.equal(samp_list[[2]], matrix(samp_mat[2,], byrow = TRUE,
#'   nrow = nrow(env)))
#' # For the Dirchlet samples is really the same
#' set.seed(1234)
#' samp_list <- r_fm(env, n_b = 2, list_out = TRUE, bay = TRUE)
#' set.seed(1234)
#' samp_mat <- r_fm(env, n_b = 2, list_out = FALSE, bay = TRUE)
#' all.equal(samp_list[[2]][[2]], matrix(samp_mat[[2]][2,], byrow = TRUE,
#'   nrow = nrow(env)))
r_fm <- function(env, n_b, bay = FALSE, prior = NULL, list_out = FALSE){
  k <- ncol(env)
  env <- unlist(lapply(seq(nrow(env)), function(x) env[x,]))
  N <- sum(env)
  if (!(bay)){
    resample <- stats::rmultinom(n_b, N, env/N)
    resample <- t(resample)
  } else {
    if (is.null(prior)){ prior <- rep(1/length(env), length(env)) }
    post <- prior + env
    dirs <- gtools::rdirichlet(n_b, post)
    resample <- matrix(unlist(lapply(seq(nrow(dirs)), function(x)
      stats::rmultinom(1, N, dirs[x,]) ) ), byrow = T, nrow = n_b)
  }
  if (list_out){
    resample <- lapply(seq(nrow(resample)), function(x) matrix(resample[x, ],
                                                          ncol = k, byrow = T))
  }
  if (bay){
    resample <- list(resample = resample, ps = dirs)
    if (list_out){
      resample[[2]] <- lapply(seq(nrow(resample[[2]])), function(x)
        matrix(resample[[2]][x, ], ncol = k, byrow = T))
    }
  }
  return(resample)
}




# For future
#The bayesian estimation can opt to a beta-mixture model instead of a saturated model
#TO do that, the prior must be given as a list, which must have a 'beta_mix' element alongside
#a 'alfas' and 'betas' element. The length of the three components must be the same
#That is, a prior for each beta distribution of the mix and a prior for the mixture!
#If no prior is specified, the function will default to a saturated model
#Only works with 2 rows

#' Simulate Independent-unit Framework
#'
#' \code{r_iu} samples n_b times from k independent multinomial
#' distributions or, following a simple Bayesian  conjugate model, k Dirichlet-
#' multinomial distributions -- where k is the number of units in an environment.
#' Returns the simulations in the form of a matrix or list.
#' If the Bayesian model is used, \code{r_iu} also returns the
#' simulated probability vectors from the Dirichlet distributions.
#'
#' \code{r_iu} simulates samples from k independent distributions, one per unit.
#' The user can sample from independent multinomial, or independent
#' Dirichlet-multinomial, depending on the value of  \code{bay}.
#' If \code{bay} is \code{TRUE}, the function will sample from
#' Dirichlet-multinomial distributions following a simple conjugate model.
#' Whatever the case,  \code{r_iu} will sample \code{n_b} times and return the
#' results in a matrix or list format, depending on the \code{list_out}
#' argument.
#' The user can specify a \code{prior} for the model, otherwise, the function
#' will default to a Jeffrey's prior. The parameter \code{ws} regulates the
#' number of individuals in each unit, regardless of the input environment.
#' If \code{ws} is not specified, the simulations will have the same number of
#' individuals per unit as the input environment.
#'
#'
#' @export
#'
#' @param env nxk \code{matrix} of \code{numeric}.
#' @param n_b \code{numeric}. How many simulations should be returned?
#' @param bay \code{logical} Should the Bayesian model be used?
#' @param prior k-long \code{list} containing n probability vectors for the
#'      Dirichlet distributions. Will be ignored if \code{bay} is FALSE.
#' @param ws A k-long \code{vector} containing the number of individual in each
#'     unit.
#' @param list_out \code{logical} should the output be in a list form?
#'      distributions -- one per column. Will be ignored if \code{bay} is FALSE.
#'
#' @return The return depends on the value of \code{list_out} and \code{bay}. By
#'     default, both are \code{FALSE}. In that case, the function returns a
#'     \code{matrix} where every row contains a simulated environment. Within
#'     each row, the values are ordered by row. In case of \code{list_out} set
#'     to \code{TRUE} and \code{bay} set to \code{FALSE}, the output will
#'     be a \code{n_b}-long list of \code{matrix} showing the simulated
#'     environments.
#'     When \code{bay} is \code{TRUE}, the output will always be a \code{list}.
#'     The first element of the list contains the simulated environments and
#'     the second element contains the Dirichlet samples used to create the
#'     simulated environment, scaled so that the total sum for each environment
#'     is one. The way these pieces of information are shaped depends on the
#'     value of \code{list_out}. If \code{list_out} is \code{FALSE},
#'     the information will be stored in \code{matrix}: one \code{matrix} for
#'     the simulated
#'     environment, one \code{matrix} for the Dirichlet samples. In the two
#'     matrices, every row represent a simulated environment, ordered by
#'     row. The two matrices are ordered so that the kth row in one corresponds
#'     to the kth row in the other. If \code{list_out} is \code{TRUE} (and
#'     \code{bay} is \code{TRUE}), the
#'     output will be shaped in a \code{list} of 2 lists. The first \code{list}
#'     contains the simulated environment. The second \code{list} contains the
#'     Dirichlet sample used to create the simulated environment. The lists
#'     are ordered so that the kth environment was created from the kth
#'     Dirichlet sample.
#'
#' @family frameworks
#' @examples
#' env <- matrix(seq(12), nrow = 3, byrow = TRUE)
#' samp <- r_iu(env, n_b = 1000, list_out = TRUE)
#' samp2 <- r_iu(env, n_b = 1000, list_out = TRUE, bay = TRUE)
#' # Very strong priors
#' priors <- list(c(0,0,30))
#' priors <- rep(priors, ncol(env))
#' samp3 <- r_iu(env, n_b = 1000, list_out = TRUE, bay = TRUE, prior =
#'  priors)
#' all.equal(colSums(env), colSums(samp[[1]]), colSums(samp2[[1]][[1]]),
#' colSums(samp3[[1]][[1]]))
#' all.equal(colSums(samp2[[2]][[1]]), colSums(samp3[[2]][[1]]))
#' # How are the list and matrix output related?
#' set.seed(1234)
#' samp_list <- r_iu(env, n_b = 2, list_out = TRUE)
#' set.seed(1234)
#' samp_mat <- r_iu(env, n_b = 2, list_out = FALSE)
#' all.equal(samp_list[[2]], matrix(samp_mat[2,], byrow = TRUE,
#'   nrow = nrow(env)))
#' # For the Dirchlet samples is really the same
#' set.seed(1234)
#' samp_list <- r_iu(env, n_b = 2, list_out = TRUE, bay = TRUE)
#' set.seed(1234)
#' samp_mat <- r_iu(env, n_b = 2, list_out = FALSE, bay = TRUE)
#' all.equal(samp_list[[2]][[2]], matrix(samp_mat[[2]][2,], byrow = TRUE,
#'   nrow = nrow(env)))
r_iu <- function(env, n_b = 5000, ws = NULL, bay = FALSE, prior = NULL,
                    list_out = FALSE){
  k <- ncol(env)
  n_r <- nrow(env)
  if (is.null(ws)){ws <- colSums(env)}
  if (!(bay)){
    #Here we sample, the rest of the lines are to put the sample in a nicer format
    resample <- lapply(seq(k), function(x) stats::rmultinom(n_b, ws[x],
                                                          env[,x]))
    resample <- mapply(function(x, y) resample[[y]][,x],
                       x = rep(seq(n_b), each = k), y = rep(seq(k), n_b))
    resample <- lapply(seq(0, n_b-1), function(x)  c(t( resample[,((x*k)+1):
                                                                  (k*(x+1))] )))
    resample <- do.call("rbind", resample)
  } else {
    # This case estimates a beta binomial model with a Gibbs sampling
    # Leave here for lter
    # if ('beta_mix' %in% names(prior)) {
    #   mix = beta_mix_prior[['beta_mix']]
    #   alfas <- beta_mix_prior[['alfas']]
    #   betas <- beta_mix_prior[['betas']]
    #   dir_prior <- rep(1/length(mix), length(mix))
    #   rbetas <- estimate_bin_bay(n_b, env, alfas, betas, burn_out, mix = mix, dir_prior = dir_prior)
    #   #We are only really interested into the probabilities coming out of the estimation
    #   #We reshape this to ease the commands later
    #   rbetas <- c(rbetas$Pas)
    #   #This case goes with a simpler saturated binomial model
    # } else {
    #Uninformative Jeffrey priors
    if (is.null(prior)) {
      prior <- list()
      prior <- lapply(seq(k), function(x) rep(1/(n_r), n_r))
    }
    post <- lapply(seq(k), function(x) prior[[x]] + env[,x])
    # The normalized weight for each column
    ws_p <- ws/sum(ws)
    #After the posterior for the dirichlet, we get the actual probability draw
    rdiri <- lapply(seq(ncol(env)), function(x) gtools::rdirichlet(n_b,
                                post[[x]] ))
    #Now we draw the various multinomial
    resample <- mapply( function(x,y)  stats::rmultinom (1, ws[x], rdiri[[x]][y,]),
                x = rep(seq(k), n_b),
                y = rep(seq(n_b), each = k)  )
    # resample <- matrix(c(resample), nrow = n_b, byrow = T)
    resample <- lapply(seq(0, n_b-1), function(x)
                       c(t( resample[,((x*k)+1):(k*(x+1))] )))
    resample <- do.call("rbind", resample)
    # Reshape the dirichlet samples as well
    # Normalize the total weight to 1
    rdiri <- lapply(seq(length(rdiri)), function(x) rdiri[[x]]*ws_p[x])
    rdiri <- do.call("cbind", rdiri)
    rdiri <- lapply(seq(nrow(rdiri)) , function(x) c(t(matrix(rdiri[x,],
                                                              ncol = k))))
    rdiri <- do.call("rbind", rdiri)
    #This is for reworking the probability into a data-like matrix
  }
  if (list_out){
    resample <- lapply(seq(nrow(resample)), function(x) matrix(resample[x, ],
                                                          ncol = k, byrow = T))
  }
  #We only have the first line of the matrix. This provides also the second
  if (bay){
    resample <- list(resample = resample, ps = rdiri)
    if (list_out){
      resample[[2]] <- lapply(seq(nrow(resample[[2]])), function(x)
        matrix(resample[[2]][x, ], ncol = k, byrow = T))
    }
  }
  return(resample)
}
