#For HT we need to have all the derivatives of the various
#indexes under the various sampling schemes. This amounts to number-indices
# times 3 (frameworks) different derivatives
##########DERIVATIVE FOR INDEPENDENT-GROUPS FRAMEWORK##########

#' Derivatives of Segregation Indices
#'
#' Functions that takes an environment as an input and output the derivative of
#' a segregation index for the given environment. Derivatives are necessary in
#' some part of hypothesis testing. The exact derivatives depend on both the
#' index and the framework used. There is a derivative available for all
#' indices-framework coupling.
#'
#'
#' Some part of hypothesis testing relies on the delta method, which requires
#' the use of derivatives with respect to the framework's parameters. This is
#' where these functions are used. All functions take as input a 2xN matrix
#' representing the analyzed environment. As for all other functions in SISeg,
#' the derivatives assume a 2-groups environment.
#' These functions output the derivative of a certain index for the given
#' environment under a certain statistical framework. The derivatives are
#' returned in the form of a \code{vector}. The vectors are
#' ordered in a specific way that depends on the framework used:
#' \describe{
#' \item{Independent-groups & Full-multinomial}{In the case of the
#'     independent-groups and full-multinomial frameworks,
#'     the partial derivatives are ordered by row. So, the first element
#'     of the returned vector is the derivative for the probability parameter
#'     of the first-group-first-unit cell; the second element
#'     of the returned vector is the derivative for the probability parameter
#'     of the first-group-second-unit cell; etc.}
#' \item{Independent-units}{In the case of the independent-units framework,
#'     the partial derivatives are ordered by column So, the first element
#'     of the returned vector is the derivative for the probability parameter
#'     of the first-group-first-unit cell; the second element
#'     of the returned vector can be interpreted as  the derivative for the
#'     probability parameter of the second-group-second-unit cell; the
#'     derivative for the probability parameter of the first-group-second-unit
#'     cell; etc.}
#' }
#'
#' The independent-units framework is particularly attuned to the binomial form
#' of the indices. As is the case for the indices in binomial form, the
#' environment can be passed in a binomial form when using these functions.
#' This means that the environment is represented as a 2xk matrix where the
#' first row show the probability of a group 1 member in each unit and the
#' second row shows the number of individuals in the unit. If the environment
#' is passed in binomial form, set the \code{in_form} to \code{TRUE} argument to
#'  signal it.
#'
#' Finally, notice that the independent-units framework only has
#' k free parameters, since one probability per unit is sufficient to fully specify
#' it -- the numbers of individuals in each unit are not considered free in the
#' framework. For this reason, half of the derivatives under the
#' independent-units framework will result to be 0. The 0s are substantially
#' place-holders.
#'
#' @section Naming convention:
#' All derivative functions' names start with a "d_". Then, they
#' have the complete name of the index function they take a derivative of, for
#' example "d_ind". Finally they have an abbreviated name of the statistical
#' framework, specifically:
#' \itemize{
#' \item{_ig for the independent-groups framework}
#' \item{_fm for the full-multinomial framework}
#' \item{_iu for the independent-units framework}
#' }
#' For example, the function ``d_d_ind_fm'' calculates the derivative of the
#' D index under the multinomial framework. Notice, everything is lowercase.
#'
#' @name seg_deriv
#'
#' @param env 2xk matrix of numeric.
#' @param in_form Is the \code{env} matrix in a binomial form?
#'
#' @return A 2 times k  \code{vector} showing the partial derivatives' value for
#'     each of the framework's parameters.
#' @examples
#' env <- matrix(c(1,2,3,4,5,6,7,8), nrow = 2, byrow = TRUE)
#' d_d_ind_ig(env)
#' d_gini_ind_ig(env)
#' env_bin_form <- env
#' env_bin_form[2,] <- colSums(env_bin_form)
#' env_bin_form[1,] <- env_bin_form[1,]/env_bin_form[2,]
#' env_bin_form[2,] <- env_bin_form[2,]/sum(env_bin_form[2,])
#' der_d <- d_d_ind_iu(env)
#' der_d_bin <-  d_d_ind_iu(env_bin_form, in_form = TRUE)
#' all.equal(der_d, der_d_bin)
NULL


#' @describeIn seg_deriv D Index, independent-groups framework.
#' @export
d_d_ind_ig <- function(env){
  env <- env/rowSums(env)
  der <- c(sign(env[1,] - env[2,]), -sign(env[1,] - env[2,]))
  return(1/2*der)
}

#' @describeIn seg_deriv Gini Index, independent-groups framework. As usual for
#' the Gini index, this is more complex than other derivatives under this
#' framework.
#' @export
d_gini_ind_ig <- function(env){
  env <- env/rowSums(env)
  abs_matrix <- t(t(rep(1, ncol(env)))) %*% t(env[1,]) * env[2,] %*% t(rep(1, ncol(env)))
  abs_matrix <- t(abs_matrix) - abs_matrix
  sign_matrix <- sign(abs_matrix)
  # abs_matrix <- abs(abs_matrix)
  # A <- sum(abs_matrix)
  val_mat_1 <- matrix(rep(env[2,],ncol(env)), byrow = TRUE,  nrow = ncol(env))
  val_mat_1 <- val_mat_1 * sign_matrix
  val_mat_2 <- matrix(rep(env[1,],ncol(env)), byrow = TRUE,  nrow = ncol(env))
  #Notice the minus
  val_mat_2 <- -val_mat_2 * sign_matrix
  val_mats <- list(val_mat_1, val_mat_2)
  ders <- rep(0, ncol(env)*2)
  for (i in seq(2*ncol(env))){
    sel = 1
    der_n <- i
    if (i > ncol(env)){
      i <- i - ncol(env)
      sel <- 2}
    ders[der_n] <- sum(val_mats[[sel]][i,])
  }
  return(ders)
}

#' @describeIn seg_deriv Mutual Information Index, independent-groups framework.
#' @export
d_mutual_info_ind_ig <- function(env){
  k <- ncol(env)
  p <- rowSums(env)/sum(env)
  env <- env/rowSums(env)
  mat_ext <- unlist(lapply(seq(nrow(env)), function(x) env[x, ]))
  #To avoid problems we add a little probability to 0 cells
  MN <- min(mat_ext[mat_ext!=0])
  mat_ext <- ifelse(mat_ext !=0, mat_ext, MN/length(mat_ext) )
  der <- unlist(lapply(seq(k*nrow(env)),function(x) p[ceiling(x/k)] * log(mat_ext[x]/(p %*% env[, x-k*(ceiling(x/k)-1)])[1]) ))
  return(der)
}

#' @describeIn seg_deriv Theil Index, independent-groups framework.
#' @export
d_theil_ig <- function(env){
  p <- rowSums(env)
  ent <- entropy(p)
  der <- 1/ent*d_mutual_info_ind_ig(env)
  return(der)
}


#' @describeIn seg_deriv Atkinson Index, independent-groups framework. Supposes
#'     that all weights are equal.
#' @export
d_atkinson_ind_ig<- function(env){
  env <- env/rowSums(env)
  k <- ncol(env)
  mat_ext <- unlist(lapply(seq(nrow(env)), function(x) env[x, ]))
  #To avoid problems we add a little probability to 0 cells
  MN <- min(mat_ext[mat_ext!=0])
  mat_ext <- ifelse(mat_ext !=0, mat_ext, MN/length(mat_ext) )
  der <- -0.5*unlist(lapply(seq(k*nrow(env)), function(x) sqrt(prod(env[, x-k*(ceiling(x/k)-1)])/(mat_ext[x]^2)) ))
  return(der)
}

#' @describeIn seg_deriv Isolation Index, independent-groups framework.
#' @export
d_isolation_ind_ig <- function(env){
  p <- rowSums(env)
  c_sm <- colSums(env)
  env <- env/rowSums(env)
  # derivative for the first group
  der <- p[1]*env[1,]/c_sm^2 * (env[1,]*p[1] + 2*env[2,]*p[2])
  # derivative for the second group
  der <- c(der, -p[1]*p[2]*(env[1,]/c_sm)^2  )
  return(der)
}

#' @describeIn seg_deriv Derivative for the Isolation Index (independent-groups
#' framework) applied to the second row of the environment matrix.
#' @export
d_isolationInv_ind_ig <- function(env){
  env <-matrix(c(env[2,], env[1,]), nrow = 2, byrow = TRUE)
  return(d_isolation_ind_ig(env))
}

#' @describeIn seg_deriv V Index, independent-groups framework.
#' @export
d_v_ind_ig <- function(env){
  p <- rowSums(env)[1]/sum(env)
  # It is very convenient to use d_xpx_dm
  der <-  d_isolation_ind_ig(env)/(1-p)
  return(der)
}


##########DERIVATIVE FOR FULL MULTINOMIAL FRAMEWORK##########

#' @describeIn seg_deriv D Index, full multinomial framework.
#' @export
d_d_ind_fm <- function(env){
  k <- ncol(env)
  sums <- rowSums(env)/(sum(env))
  env <- env/rowSums(env)
  signs <- matrix(c(sign(env[1,] - env[2,]), -sign(env[1,] - env[2,])), nrow = 2, byrow = TRUE)
  addendo <- signs %*% t(env)
  der <- 1/2 * unlist(lapply(seq(k*nrow(env)), function(x) 1/sums[ceiling(x/k)] *
                               (signs[ceiling(x/k), x-k*((ceiling(x/k))-1)]  -
                                  addendo[ceiling(x/k), ceiling(x/k)]) ))
  return(der)
}

#' @describeIn seg_deriv Gini Index, full multinomial framework. As usual for
#' the Gini index, this is more complex than other derivatives under this
#' framework.
d_gini_ind_fm <- function(env){
  p <- rowSums(env)[1]/sum(env)
  norm <- 1/2*1/(p*(1-p))
  env <- env/sum(env)
  abs_matrix <- t(t(rep(1, ncol(env)))) %*% t(env[1,]) * env[2,] %*% t(rep(1, ncol(env)))
  abs_matrix <- t(abs_matrix) - abs_matrix
  sign_matrix <- sign(abs_matrix)
  abs_matrix <- abs(abs_matrix)
  A <- sum(abs_matrix)
  val_mat_1 <- matrix(rep(env[2,],ncol(env)), byrow = TRUE,  nrow = ncol(env))
  val_mat_1 <- val_mat_1 * sign_matrix
  val_mat_2 <- matrix(rep(env[1,],ncol(env)), byrow = TRUE,  nrow = ncol(env))
  #Notice the minus
  val_mat_2 <- -val_mat_2 * sign_matrix
  val_mats <- list(val_mat_1, val_mat_2)
  ders <- rep(0, ncol(env)*2)
  for (i in seq(2*ncol(env))){
    sel = 1
    der_n <- i
    if (i > ncol(env)){
      i <- i - ncol(env)
      sel <- 2}
    ders[der_n] <- sum(val_mats[[sel]][i,])
  }
  addend <- c(rep(-A/p, ncol(env)), rep(-A/(1-p), ncol(env)))
  norm * (addend+ 2*ders)
}

#' @describeIn seg_deriv Mutual Information Index, full multinomial framework.
#' @export
d_mutual_info_ind_fm <- function(env){
  k <- ncol(env)
  env <- env/sum(env)
  mat_ext <- unlist(lapply(seq(nrow(env)), function(x) env[x, ]))
  #To avoid problems we add a little probability to 0 cells
  MN <- min(mat_ext[mat_ext!=0])
  mat_ext <- ifelse(mat_ext !=0, mat_ext, MN/length(mat_ext) )
  # mat_ext <- unlist(lapply(seq(k*nrow(env)), function(x) env[ceiling(x/k), x-k*(ceiling(x/k)-1)]))
  p <- rowSums(env)
  m1 <- colSums(env)
  # Assumes 2 rows
  d_ent_m1 <- rep(-(log(m1) + 1), nrow(env))
  der <- d_ent_m1 - unlist(lapply(seq(k*nrow(env)),function(x) log(p[ceiling(x/k)]) -log(mat_ext[x])  ))
  return(der)
}

#' @describeIn seg_deriv Theil Index, full multinomial framework.
#' @export
d_theil_ind_fm <- function(env){
  p <- rowSums(env)
  ent <- entropy(p)
  mi <- mutual_info_ind(env)
  der <- d_mutual_info_ind_fm(env)
  p <- unlist(lapply(seq(nrow(env)), function(x) rep(p[x], ncol(env))))
  der <- 1/ent * ( der + (1 + log(p))/ent * mi )
}
#Works

#' @describeIn seg_deriv Atkinson Index, full multinomial framework. All weights
#'     are supposed to be equal.
#' @export
d_atkinson_ind_fm<- function(env){
  k <- ncol(env)
  p <- rowSums(env)/sum(env)
  mat_ext <- unlist(lapply(seq(nrow(env)), function(x) env[x, ]))
  #To avoid problems we add a little probability to 0 cells
  MN <- min(mat_ext[mat_ext!=0])
  mat_ext <- ifelse(mat_ext !=0, mat_ext, MN/length(mat_ext) )
  norm <- 1/sqrt(prod(p))
  c_prod <- sum(unlist(lapply(seq(k), function(x) sqrt(prod(env[,x]/sum(env))))))
  # c_prod <- c_prod*norm
  # der <- unlist(lapply(seq(k*nrow(env)),function(x) 1/p[ceiling(x/k)] ))
  der <- -0.5*norm*unlist(lapply(seq(k*nrow(env)), function(x) sqrt(prod(env[,
              x-k*(ceiling(x/k)-1)])/(mat_ext[x]^2)) - c_prod/p[ceiling(x/k)] ))
  return(der)
}

#' @describeIn seg_deriv Isolation Index, full multinomial framework.
#' @export
d_isolation_ind_fm <- function(env){
  env <- env/sum(env)
  p <- rowSums(env)[1]
  # derivative for the first group
  addendo <- -1/p^2*sum(env[1, ]^2/colSums(env))
  der <- 1/p*(1/colSums(env))*env[1,]*(
    2 - env[1,]/ colSums(env)
  )
  der <- addendo + der
  # Add the derivative for the second group (very different)
  der <- c(der, -1/p*(env[1,]/colSums(env))^2)
  return(der)
}

#' @describeIn seg_deriv Derivative for the Isolation Index (full multinomial
#' framework) applied to the second row of the environment matrix.
#' @export
d_isolationInv_ind_fm <- function(env){
  env <- matrix(c(env[2,], env[1,]), nrow = 2, byrow = TRUE)
  return(d_isolation_ind_fm(env))
}

#' @describeIn seg_deriv V Index, full multinomial framework.
#' @export
d_v_ind_fm <- function(env){
  env <- env/sum(env)
  p <- rowSums(env)[1]
  # It is useful to use xpx
  der <- d_isolation_ind_fm(env)
  # Second group already done
  der <- 1/(1-p)*der
  # First group case
  der[1:ncol(env)] <- v_ind(env)/(1-p) + (der[1:ncol(env)] - 1/(1-p))
  return(der)
}


##########DERIVATIVE FOR INDEPENDENT-UNITS##########

#Each function has a 'in_form' parameter.
#This parameter check if the matrix is already in the same form as
#real_matrix is.

#' @describeIn seg_deriv D Index, independent-units framework.
#' @export
d_d_ind_iu <- function(env, in_form = FALSE){
  if (!(in_form)){
    env[2,] <- colSums(env)
    env[1,] <- env[1,]/env[2,]
    env[2,] <- env[2,]/sum(env[2,])
  }
  p <- (env[1,] %*% env[2,])[1]
  norm <- 1/2 * 1/(p*(1-p))
  signs <- sign(env[1,] - p)
  #There is a fixed part that does not need to be calculated for each p
  addendo <- -((1-2*p)/(p*(1-p)))*(env[2,] %*% (abs(env[1,] - p)))[1]
  addendo <- addendo - (env[2,] %*%  signs)[1]
  der <- env[2,]
  der <- norm * der * unlist(lapply( seq(ncol(env)), function(x) signs[x] + addendo ))
  der <- unlist(lapply(seq(length(der)), function(x) c(der[x],0) ))
  return(der)
}

#' @describeIn seg_deriv Gini Index, independent-units framework. As usual for
#' the Gini index, this is more complex than other derivatives under this
#' framework.
d_gini_ind_iu <- function(env, in_form = FALSE){
  if (!(in_form)){
    env[2,] <- colSums(env)
    env[1,] <- env[1,]/env[2,]
    env[2,] <- env[2,]/sum(env[2,])
  }
  p <- c(env[1,]%*% env[2,])
  norm <- (1/(p*(1-p)))
  abs_matrix <- t(t(rep(1, ncol(env)))) %*% t(env[1,])
  abs_matrix <- t(abs_matrix) - abs_matrix
  sign_matrix <- sign(abs_matrix)
  w_matrix <- matrix(rep(env[2,],ncol(env)), byrow = TRUE,  nrow = ncol(env))
  w_matrix <- sign_matrix*w_matrix
  ders <- rowSums(w_matrix)
  ders <- env[2,]*norm*(ders - norm*1/2*(1-2*p)*(t(env[2,]) %*% abs(abs_matrix) %*% t(t(env[2,])))[1])
  #We put the derivatives (= 0) for the probabilities of the second group
  ders <- unlist(lapply(seq(length(ders)), function(x) c(ders[x],0) ))
  return(ders)
}

#' @describeIn seg_deriv Mutual Information Index, independent-units framework.
#' @export
d_mutual_info_ind_iu <- function(env, in_form = FALSE) {
  if (!(in_form)){
    env[2,] <- colSums(env)
    env[1,] <- env[1,]/env[2,]
    env[2,] <- env[2,]/sum(env[2,])
  }
  p <- (env[1,]%*% env[2,])[1]
  mat_ext <- env[1,]
  #To avoid problems we add a little probability to 0 cells
  MN <- min(mat_ext[mat_ext!=0])
  mat_ext <- ifelse(mat_ext !=0, mat_ext, MN/length(mat_ext) )
  mat_ext <- ifelse(mat_ext !=1, mat_ext, 1- MN/length(mat_ext) )
  env[1,] <- mat_ext
  der <- env[2,]*(log(env[1,]/(1-env[1,])) -  log(p/(1-p)))
  der <- unlist(lapply(seq(length(der)), function(x) c(der[x],0) ))
  return(der)
}

#' @describeIn seg_deriv D Index, independent-units framework.
#' @export
d_theil_ind_iu <- function(env, in_form = FALSE) {
  if (!(in_form)){
    mi <- mutual_info_ind(env)
    env[2,] <- colSums(env)
    env[1,] <- env[1,]/env[2,]
    env[2,] <- env[2,]/sum(env[2,])
  }
  p <- (env[1,]%*% env[2,])[1]
  p <- c(p, 1-p)
  ent <- entropy(p)
  #We need to calculate mi if the matrix was in_form from the beginning
  if (in_form){
    mi <- mutual_info_ind(matrix(c(env[1,]*env[2,], (1-env[1,])*env[2,]),
                                 nrow =2, byrow = TRUE))
  }
  der <- d_mutual_info_ind_iu(env, in_form = TRUE)
  mat_2 <- unlist(lapply(seq(ncol(env)), function(x) c(env[2,x],0) ))
  der <- 1/ent * (der + mat_2/ent* log(p[1]/p[2])*mi)
  return(der)
}


#' @describeIn seg_deriv Atkinson Index, independent-units framework. Supposes
#'     that all weights are equal.
#' @export
d_atkinson_ind_iu<- function(env, in_form = FALSE){
  if (!(in_form)){
    env[2,] <- colSums(env)
    env[1,] <- env[1,]/env[2,]
    env[2,] <- env[2,]/sum(env[2,])
  }
  k <- ncol(env)
  p <- c(env[1,] %*% env[2,])
  norm <- 1/sqrt(p*(1-p))
  mat_ext <- env[1,]
  #To avoid problems we add a little probability to 0 cells
  MN <- min(mat_ext[mat_ext!=0])
  mat_ext <- ifelse(mat_ext !=0, mat_ext, MN/length(mat_ext) )
  mat_ext <- ifelse(mat_ext !=1, mat_ext, 1- MN/length(mat_ext) )
  env[1,] <- mat_ext
  #This is fixed for every term
  addendo <-norm^2*(1-2*p)*c(sqrt(env[1,]*(1- env[1,])) %*% env[2,])
  der <- 1/2*norm*env[2,]
  der <- der * (addendo -(1-2*env[1,])/sqrt(env[1,]*(1-env[1,])))
  #We add 0s for the derivative of the 1-p terms
  der <- unlist(lapply(seq(length(der)), function(x) c(der[x],0) ))
  return(der)
}

#' @describeIn seg_deriv Isolation Index, independent-units framework.
#' @export
d_isolation_ind_iu <- function(env, in_form = FALSE){
  if (!(in_form)){
    mi <- mutual_info_ind(env)
    env[2,] <- colSums(env)
    env[1,] <- env[1,]/env[2,]
    env[2,] <- env[2,]/sum(env[2,])
  }
  p <- (env[1,]%*% env[2,])[1]
  der <- env[2,]/p*(2*env[1,] - isolation_ind_bin(env, in_form = TRUE))
  #We add 0s for the derivative of the 1-p terms
  der <- unlist(lapply(seq(length(der)), function(x) c(der[x],0) ))
  return(der)
}

#' @describeIn seg_deriv Derivative for the Isolation Index (independent-units
#' framework) applied to the second row of the environment matrix.
#' @export
d_isolationInv_ind_iu <- function(env, in_form = FALSE){
  if (!(in_form)){
    mi <- mutual_info_ind(env)
    env[2,] <- colSums(env)
    env[1,] <- env[1,]/env[2,]
    env[2,] <- env[2,]/sum(env[2,])
  }
  env[1,] <- 1 - env[1,]
  return(d_isolation_ind_iu(env, in_form=TRUE))
}



#' @describeIn seg_deriv V Index, independent-units framework.
#' @export
d_v_ind_iu <- function(env, in_form = FALSE){
  if (!(in_form)){
    mi <- mutual_info_ind(env)
    env[2,] <- colSums(env)
    env[1,] <- env[1,]/env[2,]
    env[2,] <- env[2,]/sum(env[2,])
  }
  p <- (env[1,]%*% env[2,])[1]
  # convenient for the actual calculation
  p <- 1 - p
  addendo <- env[2,]/p * v_ind_bin(env, in_form = TRUE)
  # Useful to use d_isolation_ind_iu, but only non-zero values
  # Here we are making the assumption of only two rows`
  non_zero_d_isolation <- d_isolation_ind_iu(env, in_form = TRUE)[seq(1,ncol(env)*2, 2)]
  der <- addendo + 1/(p) * (non_zero_d_isolation - env[2,])
  #We add 0s for the derivative of the 1-p terms
  der <- unlist(lapply(seq(length(der)), function(x) c(der[x],0) ))
  return(der)
}
