# This file contains all definitions of indices.
# All functions work the same: they take a matrix in input, output a real number
# Row should represent groups; columns represent units.
# Naming convention: INDEXNAME_ind[_bin?]. bin indicates that the function
# uses the binomial version of the index as defined in the paper. However,
# the user should stick to the normal versions, which are more efficient.


#' Segregation Indices
#'
#' All these functions calculate an index for the input environment, which is
#' represented as a matrix. They output a positive \code{numeric} representing
#' the segregation in the environment,
#'
#' All functions have names end in "_ind". They take in input a 2xN matrix.
#' In the matix, each cell shows how many individuals from a given group are in
#' a certain unit. The rows of the matrix are considered the groups in the
#' environment. The columns are considered the different units. The functions
#'  returns the value of the desired index for the environment as a
#'  \code{numeric}.
#'
#' At this point, all indices suppose that there are only two groups in the
#' environment. This will change later.
#'
#' @references
#' For the index formulas and properties, see discussion in (sources
#'  are in order of relevance)
#'
#' \describe{
#'   \item{Atkinson Index}{
#'   James, David R., and Karl E. Taeuber. 1985. ``Measures of Segregation.''
#'   Sociological Methodology 15: 1–32.
#'
#'   Frankel, David M., and Oscar Volij. 2011. “Measuring School Segregation.”
#'   Journal of Economic Theory 146 (1): 1–38.
#'   }
#'   \item{D Index}{Duncan, Otis Dudley, and Beverly Duncan. 1955. “A Methodological Analysis of
#'  Segregation Indexes.” American Sociological Review 20 (2): 210–17.}
#'   \item{Gini Index}{Hutchens, Robert. 1991. “Segregation Curves, Lorenz
#'   Curves, and Inequality in the Distribution of People across Occupations.”
#'   Mathematical Social Sciences 21 (1): 31 – 51.}
#'   \item{Mutual Information Index and Theil Index}{
#'   Mora, Ricardo, and Javier Ruiz-Castillo. 2011. “Entropy-Based Segregation
#'   Indices.” Sociological Methodology 41 (1): 159–194.
#'
#'   Reardon, Sean F., and Glenn Firebaugh. 2002. “Measures of Multigroup
#'   Segregation.” Sociological Methodology 32 (1): 33–67.
#'
#'   Frankel, David M., and Oscar Volij. 2011. “Measuring School Segregation.”
#'   Journal of Economic Theory 146 (1): 1–38.
#'   }
#'   \item{V and Isolation Index}{
#'   Bell, Wendell. 1954. “A Probability Model for the Measurement of Ecological
#'    Segregation.” Social Forces 32 (4): 357–364.
#'
#'   Massey, Douglas S., and Nancy A. Denton. 1988. ``The Dimensions of
#'   Residential Segregation.'' Social Forces 67 (2): 281–315.
#'   }
#' }
#'
#' @name seg_ind
#'
#' @param env 2xN matrix of numeric.
#'
#' @return \code{numeric} representing the index value for the data.
#' @family indices
#' @examples
#' env <- matrix(c(1,2,3,4,5,6,7,8), nrow = 2, byrow = TRUE)
#' d_ind(env)
#' gini_ind(env)
#' mutual_info_ind(env)
#' theil_ind(env)
#' atkinson_ind(env)
#' v_ind(env)
#' isolation_ind(env)
#' isolationInv_ind(env)
#' isolation_ind(env[c(2,1),]) == isolationInv_ind(env)
NULL

#' @describeIn seg_ind D Index
#' @export
d_ind <- function(env){
  Mal <- sum(env[1,])
  Fem <- sum(env[2,])
  return(sum(abs(env[1,]/Mal - env[2,]/Fem))*.5)
}

#' @describeIn seg_ind Gini index. It shoud be noticed that the Gini index requires more complex
#'     calculations with respect to the other indices in the SISeg package. For this
#'     reason, it takes more time to compute and the difference may be noticeable
#'     in the construction of confidence intervals.
#' @export
gini_ind <- function(env){
  p <- rowSums(env)[1]/sum(env)
  env <- env/rowSums(env)
  abs_matrix <- t(t(rep(1, ncol(env)))) %*% t(env[1,]) * env[2,] %*%
    t(rep(1, ncol(env)))
  abs_matrix <- abs(abs_matrix- t(abs_matrix))
  norm <- 1/2
  norm * sum(abs_matrix)
}

# For MI below
# Not exported. v is a vector or matrix. returns entropy
entropy <- function(v){
  v <- v/sum(v)
  v <- v[v != 0]
  ent <- -sum(v*log(v))
  return(ent)
}

#' @describeIn seg_ind Mutual Information Index. It uses the exponential as a
#'     base for logarithms.
#' @export
mutual_info_ind <- function(env){
  p <-rowSums(env)/sum(env)
  ents <- unlist(lapply(1:nrow(env), function(x) entropy(env[x,])))
  base_ent <- entropy(colSums(env)/sum(env))
  ents <- base_ent-ents
  return(sum(p*ents))
}

#' @describeIn seg_ind Theil Index (One of the possible normalization of the
#'     Mutual Information Index). It uses the exponential as a
#'     base for logarithms.
#' @export
theil_ind <- function(env){
  p <-rowSums(env)/sum(env)
  ent <- entropy(p)
  return(1/ent * mutual_info_ind(env))
}

#' @describeIn seg_ind Atkinson Index as defined by Frankel and Volij (2011).
#'     The weights are assumed to be uniform over the units.
#' @export
atkinson_ind <- function(env){
  w <- nrow(env)^-1
  env <- env/rowSums(env)
  prods <- unlist(lapply(seq(ncol(env)), function(x) prod(env[,x])))
  prods <- sum(prods^w)
  return(1-prods)
}

#' @describeIn seg_ind V Index, also known as Eta-squared index.
#' @export
v_ind <- function(env){
  env <- env/sum(env)
  ps <- rowSums(env)
  Os <- colSums(env)
  env <- env/ps
  ps[1]*(sum(env[1,]^2/Os) - 1)/(ps[2])
}

#' @describeIn seg_ind Isolation Index. Notice, this is asymmetric.
#' @export
isolation_ind <- function(env){
  env <- env/sum(env)
  p <- rowSums(env)[1]
  m1 <- env/rowSums(env)
  return(p*sum(m1[1,]^2/colSums(env)))
}

#' @describeIn seg_ind A convenience function to calculate the Isolation Index for
#'    the second group in an environment.
#' @export
isolationInv_ind <- function(env){
  env <- env[c(2,1),]
  return(isolation_ind(env))
}

#We also define some indexes for the labor-market perspective
#That is the binomial product
#Here the matrix must be a vector of probabilities (first row)
#And a vector of weights (second row)


#' Binomial Form of Segregation Indices
#'
#' All these functions calculate the value of a segregation index in an
#' environment using the binomial form. They output a positive \code{numeric}
#' representing the segregation in the environment. These functions use the
#' binomial form of the indices. Results are identical to the use of standard
#' forms implemented in  \link[=seg_ind]{Segregation Indices}.
#'
#' In general, it is possible to write segregation indices in two different
#' forms. The binomial form uses the proportion of a group in each unit, as
#' opposed to the standard form, which uses the proportion of each unit in the
#' group.It is usually faster and convenient to use the standard form, which is
#' indeed implemented in the basic functions to calculate segregation indices.
#'
#' For a given index, the function calculating the index value through the
#' binomial form has a name that simply adds "_bin" to the basic function for
#' the index. For example, \code{d_ind()} and \code{d_ind_bin()}.
#' As for the base versions of the indices, these function take a 2xN matrix
#' as input. The rows of the matrix are considered the groups in the environment.
#' The columns are considered the different units. They returns the value of the
#' desired index for the environment as a \code{numeric}.
#'
#' The input matrix can be in the so-called ``binomial form''. In this form,
#' the first row represents the proportion of individuals from group one in
#' a unit. The second row represents the proportions of individual from a
#' unit in the environment. If the input matrix is in this form, then the
#' argument \code{in_form} should be set to \code{TRUE}
#'
#' The use of a binomial form is convenient, theoretically, to formulate
#' multi-group indices. See
#'
#' Reardon, Sean F., and Glenn Firebaugh. 2002. “Measures of Multigroup
#' Segregation.” Sociological Methodology 32 (1): 33–67.
#'
#' @name seg_ind_bin
#'
#' @param env 2xN \code{matrix} of \code{numeric}.
#' @param in_form \code{logical}, is \code{env} in binomial form?
#'
#' @return \code{numeric} representing the index value for the data. This result
#'     is identical to the basic segregation indices functions.
#' @family indices
#' @examples
#' env <- matrix(c(1,2,3,4,5,6,7,8), nrow = 2, byrow = TRUE)
#' basic <- d_ind(env)
#' bin <-  d_ind_bin(env)
#' all.equal(basic, bin)
#' env_bin_form <- env
#' env_bin_form[2,] <- colSums(env_bin_form)
#' env_bin_form[1,] <- env_bin_form[1,]/env_bin_form[2,]
#' env_bin_form[2,] <- env_bin_form[2,]/sum(env_bin_form[2,])
#' bin_bin_form <- d_ind_bin(env_bin_form, in_form = TRUE)
#' all.equal(bin_bin_form, bin)
NULL




#' @describeIn seg_ind_bin D Index
#' @export
d_ind_bin <- function(env, in_form = FALSE){
  if (!(in_form)){
    env[2,] <- colSums(env)
    env[1,] <- env[1,]/env[2,]
    env[2,] <- env[2,]/sum(env[2,])
  }
  #Calculate the row margin
  p_f <- sum(unlist(lapply(seq(ncol(env)), function(x) prod(env[,x]))))
  #Calculate the normalizing factor
  p_f <- c(env[1,]%*% env[2,])
  norm <- (2*p_f*(1-p_f))^(-1)
  #Proceed to sum
  su <- sum(env[2,]*abs(p_f - env[1,]))
  return(su*norm)
}

#' @describeIn seg_ind_bin Gini Index. As for its basic version, the Gini Index
#'     is quite cumbersome to calculate.
#' @export
gini_ind_bin <- function(env, in_form = FALSE){
  if (!(in_form)){
    env[2,] <- colSums(env)
    env[1,] <- env[1,]/env[2,]
    env[2,] <- env[2,]/sum(env[2,])
  }
  p <- c(env[1,]%*% env[2,])
  norm <- 1/2*(1/(p*(1-p)))
  abs_matrix <- t(t(rep(1, ncol(env)))) %*% t(env[1,])
  abs_matrix <- abs(abs_matrix - t(abs_matrix))
  c(norm * t(env[2,]) %*% abs_matrix %*% t(t(env[2,])))
}

#Mutual information in binomial version
#This is the classic form by Theil

#' @describeIn seg_ind_bin Mutual Information Index.
#' @export
mutual_info_ind_bin <- function(env, in_form = FALSE ){
  if (!(in_form)){
    env[2,] <- colSums(env)
    env[1,] <- env[1,]/env[2,]
    env[2,] <- env[2,]/sum(env[2,])
  }
  p <- (env[1,] %*% env[2,])[1]
  ents <- unlist(lapply(1:ncol(env), function(x) env[2,x]*entropy(c(env[1,x], 1-env[1,x]))))
  base_ent <- entropy(c(p, 1-p))
  return(base_ent-sum(ents))
}

#' @describeIn seg_ind_bin Theil index. (One of the possible normalization of the
#'     Mutual Information Index)
#' @export
theil_ind_bin <- function(env, in_form = FALSE ){
  if (!(in_form)){
    env[2,] <- colSums(env)
    env[1,] <- env[1,]/env[2,]
    env[2,] <- env[2,]/sum(env[2,])
  }
  p <- (env[1,] %*% env[2,])[1]
  ent <- entropy(c(p, 1-p))
  return(1/ent * mutual_info_ind_bin(env, in_form = TRUE))
}


#' @describeIn seg_ind_bin Atkinson Index as defined by Frankel and Volij (2011).
#'     The weights are assumed to be uniform over the units.
#' @export
atkinson_ind_bin <- function(env, in_form = FALSE){
  if (!(in_form)){
    env[2,] <- colSums(env)
    env[1,] <- env[1,]/env[2,]
    env[2,] <- env[2,]/sum(env[2,])
  }
  p <- (env[1,] %*% env[2,])[1]
  norm <- 1/(sqrt(p*(1-p)))
  env[1,] <- sqrt(env[1,]*(1-env[1,]))
  seg <- 1 - norm*(env[1,] %*% env[2,])[1]
  return(seg)
}


#' @describeIn seg_ind_bin V Index, also known as Eta-squared index.
#' @export
v_ind_bin <- function(env, in_form = FALSE){
  if (!(in_form)){
    env[2,] <- colSums(env)
    env[1,] <- env[1,]/env[2,]
    env[2,] <- env[2,]/sum(env[2,])
  }
  #Calculate the row margin
  p_f <- sum(unlist(lapply(seq(ncol(env)), function(x) prod(env[,x]))))
  tot <- (env[2,]%*%(env[1,]^2))[1]
  return(tot/(p_f*(1-p_f)) - p_f/(1-p_f))
}

#' @describeIn seg_ind_bin Isolation Index. Notice, this is asymmetric.
#' @export
isolation_ind_bin <- function(env, in_form = FALSE){
  if (!in_form){
    env[2,] <- colSums(env)
    env[1,] <- env[1,]/env[2,]
    env[2,] <- env[2,]/sum(env[2,])
  }
  p <- (env[1,]%*%env[2,])[1]
  return(1/p*sum(env[2,]*env[1,]^2))
}

#' @describeIn seg_ind_bin A convenience function to calculate the Isolation Index for
#'    the second group in an environment.
#' @export
isolationInv_ind_bin <- function(env, in_form = FALSE){
  if (!in_form){
    env[2,] <- colSums(env)
    env[1,] <- env[1,]/env[2,]
    env[2,] <- env[2,]/sum(env[2,])
  }
  env[1,]  <- 1-env[1,]
  return(isolation_ind_bin(env, in_form = TRUE))
}
