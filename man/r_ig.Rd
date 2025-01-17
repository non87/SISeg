% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simul.R
\name{r_ig}
\alias{r_ig}
\title{Simulate Independent Group Framework}
\usage{
r_ig(env, n_b = 5000, bay = FALSE, prior = NULL, list_out = FALSE)
}
\arguments{
\item{env}{nxk \code{matrix} of \code{numeric}.}

\item{n_b}{\code{numeric}. How many simulations should be returned?}

\item{bay}{\code{logical} Should the Bayesian model be used?}

\item{prior}{n-long \code{list} containing n probability vectors for the
Dirichlet distributions. Will be ignored if \code{bay} is FALSE.}

\item{list_out}{\code{logical} should the output be in a list form?}
}
\value{
The return depends on the value of \code{list_out} and \code{bay}. By
default, both are \code{FALSE}. In that case, the function returns a
\code{matrix} where every row contains a simulated environment. Within
each row, the values are ordered by row. In case of \code{list_out} set
to \code{TRUE} and \code{bay} set to \code{FALSE}, the output will
be a \code{n_b}-long list of \code{matrix} showing the simulated
environments.
When \code{bay} is \code{TRUE}, the output will always be a \code{list}.
The first element of the list contains the simulated environments and
the second element contains the Dirichlet samples used to create the
simulated environment, scaled so that the total sum for each environment
is one. The way these pieces of information are shaped depends on the
value of \code{list_out}. If \code{list_out} is \code{FALSE},
the information will be stored in \code{matrix}: one \code{matrix} for
the simulated
environment, one \code{matrix} for the Dirichlet samples. In the two
matrices, every row represent a simulated environment, ordered by
row. The two matrices are ordered so that the kth row in one corresponds
to the kth row in the other. If \code{list_out} is \code{TRUE} (and
\code{bay} is \code{TRUE}), the
output will be shaped in a \code{list} of 2 lists. The first \code{list}
contains the simulated environment. The second \code{list} contains the
Dirichlet sample used to create the simulated environment. The lists
are ordered so that the kth environment was created from the kth
Dirichlet sample.
}
\description{
\code{r_ig} samples from one independent multinomial distribution per
group or, following a simple conjugate model, one Dirichlet-
multinomial per group. Returns \code{n_b} simulated samples in the form of a
list or a matrix. Will also include simulated probability vector when
Bayesian option is activated.
}
\details{
\code{r_ig} simulates samples from k independent distributions, one per group.
The user can sample from independent multinomial distributions, or
independent Dirichlet-multinomial distributions.
If \code{bay} is \code{TRUE}, the function will sample from
Dirichlet-multinomial distributions following a simple conjugate model.
Whatever the case, It will sample \code{n_b} times and return the results
in a matrix or list format, depending on the \code{list_out} argument.
The user can specify a \code{prior} for the model, otherwise, the function
will default to a Jeffrey's prior.
}
\examples{
env <- matrix(seq(12), nrow = 3, byrow = TRUE)
samp <- r_ig(env, n_b = 1000, list_out = TRUE)
samp2 <- r_ig(env, n_b = 1000, bay = TRUE, list_out = TRUE)
# Very strong priors
priors <- list(c(0,0,30, 0))
priors <- rep(priors, nrow(env))
samp3 <- r_ig(env, n_b = 1000, list_out = TRUE, bay = TRUE, prior =
priors)
all.equal(rowSums(env), rowSums(samp[[1]]), rowSums(samp2[[1]][[1]]),
 rowSums(samp3[[1]][[1]]))
all.equal(rowSums(samp2[[2]][[1]]), rowSums(samp3[[2]][[1]]))
# How are the list and matrix output related?
set.seed(1234)
samp_list <- r_ig(env, n_b = 2, list_out = TRUE)
set.seed(1234)
samp_mat <- r_ig(env, n_b = 2, list_out = FALSE)
all.equal(samp_list[[2]], matrix(samp_mat[2,], byrow = TRUE,
  nrow = nrow(env)))
# For the Dirchlet samples is really the same
set.seed(1234)
samp_list <- r_ig(env, n_b = 2, list_out = TRUE, bay = TRUE)
set.seed(1234)
samp_mat <- r_ig(env, n_b = 2, list_out = FALSE, bay = TRUE)
all.equal(samp_list[[2]][[2]], matrix(samp_mat[[2]][2,], byrow = TRUE,
  nrow = nrow(env)))
}
\seealso{
Other frameworks: 
\code{\link{r_fm}()},
\code{\link{r_iu}()}
}
\concept{frameworks}
