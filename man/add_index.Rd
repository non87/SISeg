% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{add_index}
\alias{add_index}
\title{Inference on New Indices}
\usage{
add_index(index_name)
}
\arguments{
\item{index_name}{A \code{string} with the name of the index to be added.}
}
\value{
\code{NULL}
}
\description{
The user can define new segregation indices and SISeg will perform
inference on those automatically. However, the user must define the
all the needed functions for inference following the naming conventions.
In particular, the user must define an index and three (one per framework)
derivative functions (see details). Then the function \code{add_index} will
be used to signal to SISeg that it can perform inference on a new index.
}
\details{
The inference techniques in SISeg are general and can be applied to any
(derivable) index. The user can specify a new index and then use the
inference functions in SISeg to get inference for the newly specified index.
However, the user must specify four functions returning to do this:
\itemize{
\item{An index function to calculate the index in an environment. This
function must take exactly one positional argument -- \code{env} --
for the environment whose segregation needs to be quantified.}
\item{A derivative function of the index under the independent-groups
framework. This function must take exactly one positional argument
-- \code{env} --  for the environment whose segregation needs to be
quantified.}
\item{A derivative function of the index under the full multinomial framework.
framework. This function must take exactly one positional argument
-- \code{env} --  for the environment whose segregation needs to be
quantified.}
\item{A derivative function of the index under the independent-units
framework. This function must take exactly one positional argument
-- \code{env} --  for the environment whose segregation needs to be
quantified. It must have the optional argument \code{in_form} to
handle environments in the binomial form.}
}
The four functions must follow the naming conventions of SISeg to add the
index available in SISeg. In itself, the index can have any name, but avoid
names of indices already present in the package.

Derivatives are necessary to create bootstrap studentized confidence intervals
and asymptotic confindence interval. In case the user is interested in other
kinds of confidence intervals, the derivative functions can be simple
placeholders returning a numeric. However, it is still necessary to have the
derivative functions defined. In the future, it will be possible to avoid
the derivative in the studentized confidence intervals (at the expense of
computation time).
}
\section{Naming convention}{

For the use of SISeg inferential techniques on new indices it is necessary
that the users follows the naming conventions of SISeg in their functions.
The naming conventions for the index function can be found in
\code{\link{d_ind}}. The naming conventions for the derivative functions can
be found in \code{\link{d_d_ind_fm}}.
}

\examples{
# Bogus index
bogus_ind <- function(env) return(0.3)
# Its derivative
d_bogus_ind_ig <- function(env) return(0)
d_bogus_ind_fm <- function(env) return(0)
d_bogus_ind_iu <- function(env) return(0)
# We are good to go.
add_index("Bogus")

}
