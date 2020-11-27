#### Bunch of functions that are used in the calculation of the delta method



#For HT, calculate the var-covar matrix from a multinomial
mult_vcv <- function(p, n){
  v <- diag(p)
  cv <- -p %*% t(p)
  return(n*(v + cv))
}

#For HT, calculate the empirical var-covar matrix
mult_vcv_est <- function(env){
  n <- sum(env)
  env <- env/sum(n)
  return(mult_vcv(env,n))
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


#' Create the Var-Covar Matrix of a Framework from Data
#'
#' These function produces the variance-covariance matrix for the statistical
#' frameworks starting from the data. The independence structures of the
#' frameworks are very different, resulting in more or less sparse
#' variance-covariance matrices.
#'
#' @name sigmas
#'
#' @param env 2xk \code{matrix} of numeric
#'
#' @return A kxk Variance-Covariance \code{matrix}
#'
#' @examples
#' env <- matrix(seq(8), nrow = 2, byrow = TRUE)
#' est1 <- sigma_ig(env)
#' est2 <- sigma_fm(env)
#' est3 <- sigma_iu(env)
NULL

#' @describeIn sigmas Variance-covariance matrix for the independent-groups
#'     framework
#' @export
sigma_ig <- function(env){
  n_m <- sum(env[1,])
  n_f <- sum(env[2,])
  N <- sum(n_m, n_f)
  # This is a way to avoid doing a multiplication and a division
  Ms_s <- list(mult_vcv_est(env[1,]/sum(env[1,])) * N/n_m,
               mult_vcv_est(env[2,]/sum(env[2,])) * N/n_f)
  return(bdiag_m(Ms_s))
}

#' @describeIn sigmas Variance-covariance matrix for the full multinomial
#'     framework
#' @export
sigma_fm <- function(env){
  N <- sum(env)
  env <- c(env[1,], env[2,])
  return(mult_vcv(env/N, 1))
}

#' @describeIn sigmas Variance-covariance matrix for the independent-units
#'     framework
#' @export
sigma_iu <- function(env){
  Ms <- lapply(seq(ncol(env)), function(x) mult_vcv_est(c(env[1,x],env[2,x])/
                                                          sum(env[,x])) )
  return(bdiag_m(Ms))
}

#####Here a delicate part. We create an environment to stock the indices/
#####derivatives/frameworks and their connection. It is also possible to add
#####an index anew.
#### The code is written thinking about the future addition of new indices and
#### frameworks
#https://www.r-bloggers.com/2013/04/package-wide-variablescache-in-r-packages/
.index_env <- new.env(parent=emptyenv())
# Will contain user defined indices
assign('other_indices', NULL, envir = .index_env)
#other_indices_initialize()



# Basic creation of vector of indices that are already supported
index <- function(){
  # This is the fundamental that is used to create names, retrieve all
  # functions, create lists
  supported_indices <- c("D","Gini", "Mutual_info", "Atkinson",
                          "Theil", "V", "isolation", "isolationinv")
  assign('index', c(supported_indices,
                    .index_env$other_indices), pos = .index_env)
}

# create a list of indices function based on the naming convention
create_index_list <- function(){
  assign('indexes', list(), pos = .index_env)
  for (ind in .index_env$index){
    # We lower all captals for the name
    ind <- tolower(ind)
    name <- paste0(ind, "_ind")
    .index_env$indexes[[ind]] <- match.fun(name)
  }
}

# Name the frameworks and their abbreviation
frameworks <- function(){
  assign('frmews', c("Independent_groups", "Full_Multinomial",
                         "Independent_Units"), pos = .index_env)
  assign('frmews_abbr', c("ig", "fm",
                         "iu"), pos = .index_env)
  names(.index_env$frmews_abbr) <- .index_env$frmews
  x <- vector("list", length=length(.index_env$frmews))
  names(x) <- .index_env$frmews_abbr
  # THis is more readable than *apply and it is only for a short list
  for (frm in  .index_env$frmews_abbr){
    fun_name <- paste('r', frm , sep = "_")
    x[[frm]] <- match.fun(fun_name)
  }
  assign('boots', x, pos = .index_env)
}

# Create a list of functions that create var-covar matrices in the
# different environments. Must run after frameworks()
sigmas <- function(){
  x <-vector("list", length=length(.index_env$frmews))
  names(x) <- .index_env$frmews_abbr
  # THis is more readable than *apply and it is only for a short list
  for (frm in  .index_env$frmews_abbr){
    fun_name <- paste('sigma', frm , sep = "_")
    x[[frm]] <- match.fun(fun_name)
  }
  assign('sigmas', x, pos = .index_env)
}

# Create a list of derivative for a given framework in .ind_environment
# Suppose framworks() and index() were run already
# Returns the name of the created list
create_derivative_list_framework <- function(framework = NULL){
  if (is.null(framework == 0)){
    stop("Must specify an existing framework.")
  }
  if (!(framework %in% .index_env$frmews)){
    stop("framework must be within existing frameworks")
  }
  abbr <- .index_env$frmews_abbr[framework]
  list_name <- paste('ders', abbr, sep = "_")
  # Create the derivative list outside and insert that last
  x <- vector("list", length=length(.index_env$index))
  names(x) <- tolower(.index_env$index)
  # Assign functions to list elements
  for (ind in .index_env$index){
    ind <- tolower(ind)
    fun_name <- paste("d", ind, "ind", abbr, sep = "_")
    x[[ind]] <- match.fun(fun_name)
  }
  assign(list_name, x, pos = .index_env)
  return(list_name)
}

# Create a list of lists containing the derivative for each framework/index
# In the list of lists the first level is the framework, the second the index
# e.g. ders[["fm"]][["D"]]
create_derivative_lists <- function(){
  list_names <- c()
  for (frmew in .index_env$frmews){
    list_names <- append(list_names, create_derivative_list_framework(
      framework = frmew))
    names(list_names)[length(list_names)] <- .index_env$frmews_abbr[frmew]
  }
  x <- vector("list", length=length(list_names))
  names(x) <- .index_env$frmews_abbr
  assign('ders', x, pos = .index_env)
  for (frmew in .index_env$frmews_abbr){
    .index_env$ders[[frmew]] <- get(list_names[frmew], pos = .index_env)
  }
}

# Create list of valid method names for inference
create_method_vector <- function(){
  x <- c("Boot", "Asymp", "TBoot", "CBoot", "BayModel", "CBayBoot")
  assign('methods', x, pos = .index_env)
}

# Wrapper to be called on load
create_env_lists <- function(){
  index()
  create_index_list()
  frameworks()
  sigmas()
  create_derivative_lists()
  create_method_vector()
}


# Function to create new index for inference
# Indicate the name. Derivatives must be created in advance

#' Inference on New Indices
#'
#' The user can define new segregation indices and SISeg will perform
#' inference on those automatically. However, the user must define the
#' all the needed functions for inference following the naming conventions.
#' In particular, the user must define an index and three (one per framework)
#' derivative functions (see details). Then the function \code{add_index} will
#' be used to signal to SISeg that it can perform inference on a new index.
#'
#' The inference techniques in SISeg are general and can be applied to any
#' (derivable) index. The user can specify a new index and then use the
#' inference functions in SISeg to get inference for the newly specified index.
#' However, the user must specify four functions returning to do this:
#' \itemize{
#' \item{An index function to calculate the index in an environment. This
#'    function must take exactly one positional argument -- \code{env} --
#'    for the environment whose segregation needs to be quantified.}
#' \item{A derivative function of the index under the independent-groups
#'     framework. This function must take exactly one positional argument
#'      -- \code{env} --  for the environment whose segregation needs to be
#'      quantified.}
#' \item{A derivative function of the index under the full multinomial framework.
#'     framework. This function must take exactly one positional argument
#'      -- \code{env} --  for the environment whose segregation needs to be
#'      quantified.}
#' \item{A derivative function of the index under the independent-units
#'     framework. This function must take exactly one positional argument
#'      -- \code{env} --  for the environment whose segregation needs to be
#'      quantified. It must have the optional argument \code{in_form} to
#'      handle environments in the binomial form.}
#' }
#' The four functions must follow the naming conventions of SISeg to add the
#' index available in SISeg. In itself, the index can have any name, but avoid
#' names of indices already present in the package.
#'
#' Derivatives are necessary to create bootstrap studentized confidence intervals
#' and asymptotic confindence interval. In case the user is interested in other
#' kinds of confidence intervals, the derivative functions can be simple
#' placeholders returning a numeric. However, it is still necessary to have the
#' derivative functions defined. In the future, it will be possible to avoid
#' the derivative in the studentized confidence intervals (at the expense of
#' computation time).
#'
#' @section Naming convention:
#' For the use of SISeg inferential techniques on new indices it is necessary
#' that the users follows the naming conventions of SISeg in their functions.
#' The naming conventions for the index function can be found in
#'  \code{\link{d_ind}}. The naming conventions for the derivative functions can
#' be found in \code{\link{d_d_ind_fm}}.
#'
#' @param index_name A \code{string} with the name of the index to be added.
#'
#' @return \code{NULL}
#'
#' @examples
#' # Bogus index
#' bogus_ind <- function(env) return(0.3)
#' # Its derivative
#' d_bogus_ind_ig <- function(env) return(0)
#' d_bogus_ind_fm <- function(env) return(0)
#' d_bogus_ind_iu <- function(env) return(0)
#' # We are good to go.
#' add_index("Bogus")
#'
#' @export
add_index <- function(index_name){
  # Add the new index name
  .index_env$other_indices <- append(.index_env$other_indices, index_name)
  # Now we re-create the lists with the new index. Supposes that
  # the derivative were there already
  index()
  create_index_list()
  # frameworks()
  create_derivative_lists()
}




