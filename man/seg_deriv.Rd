% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/derivatives.R
\name{seg_deriv}
\alias{seg_deriv}
\alias{d_d_ind_ig}
\alias{d_gini_ind_ig}
\alias{d_mutual_info_ind_ig}
\alias{d_theil_ind_ig}
\alias{d_atkinson_ind_ig}
\alias{d_isolation_ind_ig}
\alias{d_isolationinv_ind_ig}
\alias{d_v_ind_ig}
\alias{d_d_ind_fm}
\alias{d_gini_ind_fm}
\alias{d_mutual_info_ind_fm}
\alias{d_theil_ind_fm}
\alias{d_atkinson_ind_fm}
\alias{d_isolation_ind_fm}
\alias{d_isolationinv_ind_fm}
\alias{d_v_ind_fm}
\alias{d_d_ind_iu}
\alias{d_gini_ind_iu}
\alias{d_mutual_info_ind_iu}
\alias{d_theil_ind_iu}
\alias{d_atkinson_ind_iu}
\alias{d_isolation_ind_iu}
\alias{d_isolationinv_ind_iu}
\alias{d_v_ind_iu}
\title{Derivatives of Segregation Indices}
\usage{
d_d_ind_ig(env)

d_gini_ind_ig(env)

d_mutual_info_ind_ig(env)

d_theil_ind_ig(env)

d_atkinson_ind_ig(env)

d_isolation_ind_ig(env)

d_isolationinv_ind_ig(env)

d_v_ind_ig(env)

d_d_ind_fm(env)

d_gini_ind_fm(env)

d_mutual_info_ind_fm(env)

d_theil_ind_fm(env)

d_atkinson_ind_fm(env)

d_isolation_ind_fm(env)

d_isolationinv_ind_fm(env)

d_v_ind_fm(env)

d_d_ind_iu(env, in_form = FALSE)

d_gini_ind_iu(env, in_form = FALSE)

d_mutual_info_ind_iu(env, in_form = FALSE)

d_theil_ind_iu(env, in_form = FALSE)

d_atkinson_ind_iu(env, in_form = FALSE)

d_isolation_ind_iu(env, in_form = FALSE)

d_isolationinv_ind_iu(env, in_form = FALSE)

d_v_ind_iu(env, in_form = FALSE)
}
\arguments{
\item{env}{2xk matrix of numeric.}

\item{in_form}{Is the \code{env} matrix in a binomial form?}
}
\value{
A 2 times k  \code{vector} showing the partial derivatives' value for
each of the framework's parameters.
}
\description{
Functions that takes an environment as an input and output the derivative of
a segregation index for the given environment. Derivatives are necessary in
some part of hypothesis testing. The exact derivatives depend on both the
index and the framework used. There is a derivative available for all
indices-framework coupling.
}
\details{
Some part of hypothesis testing relies on the delta method, which requires
the use of derivatives with respect to the framework's parameters. This is
where these functions are used. All functions take as input a 2xN matrix
representing the analyzed environment. As for all other functions in SISeg,
the derivatives assume a 2-groups environment.
These functions output the derivative of a certain index for the given
environment under a certain statistical framework. The derivatives are
returned in the form of a \code{vector}. The vectors are
ordered in a specific way that depends on the framework used:
\describe{
\item{Independent-groups & Full-multinomial}{In the case of the
independent-groups and full-multinomial frameworks,
the partial derivatives are ordered by row. So, the first element
of the returned vector is the derivative for the probability parameter
of the first-group-first-unit cell; the second element
of the returned vector is the derivative for the probability parameter
of the first-group-second-unit cell; etc.}
\item{Independent-units}{In the case of the independent-units framework,
the partial derivatives are ordered by column So, the first element
of the returned vector is the derivative for the probability parameter
of the first-group-first-unit cell; the second element
of the returned vector can be interpreted as  the derivative for the
probability parameter of the second-group-second-unit cell; the
derivative for the probability parameter of the first-group-second-unit
cell; etc.}
}

The independent-units framework is particularly attuned to the binomial form
of the indices. As is the case for the indices in binomial form, the
environment can be passed in a binomial form when using these functions.
This means that the environment is represented as a 2xk matrix where the
first row show the probability of a group 1 member in each unit and the
second row shows the number of individuals in the unit. If the environment
is passed in binomial form, set the \code{in_form} to \code{TRUE} argument to
signal it.

Finally, notice that the independent-units framework only has
k free parameters, since one probability per unit is sufficient to fully specify
it -- the numbers of individuals in each unit are not considered free in the
framework. For this reason, half of the derivatives under the
independent-units framework will result to be 0. The 0s are substantially
place-holders.
}
\section{Functions}{
\itemize{
\item \code{d_d_ind_ig}: D Index, independent-groups framework.

\item \code{d_gini_ind_ig}: Gini Index, independent-groups framework. As usual for
the Gini index, this is more complex than other derivatives under this
framework.

\item \code{d_mutual_info_ind_ig}: Mutual Information Index, independent-groups framework.

\item \code{d_theil_ind_ig}: Theil Index, independent-groups framework.

\item \code{d_atkinson_ind_ig}: Atkinson Index, independent-groups framework. Supposes
that all weights are equal.

\item \code{d_isolation_ind_ig}: Isolation Index, independent-groups framework.

\item \code{d_isolationinv_ind_ig}: Derivative for the Isolation Index (independent-groups
framework) applied to the second row of the environment matrix.

\item \code{d_v_ind_ig}: V Index, independent-groups framework.

\item \code{d_d_ind_fm}: D Index, full multinomial framework.

\item \code{d_gini_ind_fm}: Gini Index, full multinomial framework. As usual for
the Gini index, this is more complex than other derivatives under this
framework.

\item \code{d_mutual_info_ind_fm}: Mutual Information Index, full multinomial framework.

\item \code{d_theil_ind_fm}: Theil Index, full multinomial framework.

\item \code{d_atkinson_ind_fm}: Atkinson Index, full multinomial framework. All weights
are supposed to be equal.

\item \code{d_isolation_ind_fm}: Isolation Index, full multinomial framework.

\item \code{d_isolationinv_ind_fm}: Derivative for the Isolation Index (full multinomial
framework) applied to the second row of the environment matrix.

\item \code{d_v_ind_fm}: V Index, full multinomial framework.

\item \code{d_d_ind_iu}: D Index, independent-units framework.

\item \code{d_gini_ind_iu}: Gini Index, independent-units framework. As usual for
the Gini index, this is more complex than other derivatives under this
framework.

\item \code{d_mutual_info_ind_iu}: Mutual Information Index, independent-units framework.

\item \code{d_theil_ind_iu}: D Index, independent-units framework.

\item \code{d_atkinson_ind_iu}: Atkinson Index, independent-units framework. Supposes
that all weights are equal.

\item \code{d_isolation_ind_iu}: Isolation Index, independent-units framework.

\item \code{d_isolationinv_ind_iu}: Derivative for the Isolation Index (independent-units
framework) applied to the second row of the environment matrix.

\item \code{d_v_ind_iu}: V Index, independent-units framework.
}}

\section{Naming convention}{

All derivative functions' names start with a "d_". Then, they
have the complete name of the index function they take a derivative of, for
example "d_ind". Finally they have an abbreviated name of the statistical
framework, specifically:
\itemize{
\item{_ig for the independent-groups framework}
\item{_fm for the full-multinomial framework}
\item{_iu for the independent-units framework}
}
For example, the function ``d_d_ind_fm'' calculates the derivative of the
D index under the multinomial framework. Notice, everything is lowercase.
}

\examples{
env <- matrix(c(1,2,3,4,5,6,7,8), nrow = 2, byrow = TRUE)
d_d_ind_ig(env)
d_gini_ind_ig(env)
env_bin_form <- env
env_bin_form[2,] <- colSums(env_bin_form)
env_bin_form[1,] <- env_bin_form[1,]/env_bin_form[2,]
env_bin_form[2,] <- env_bin_form[2,]/sum(env_bin_form[2,])
der_d <- d_d_ind_iu(env)
der_d_bin <-  d_d_ind_iu(env_bin_form, in_form = TRUE)
all.equal(der_d, der_d_bin)
}
