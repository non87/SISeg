# #CI for one index
# #mat: data in matrix form (tested with only two rows!)
# #index: which index are you interested in? "D", "gini", "MI" & "Atk" are supported
# #method: Bootstrap, Bootstrap with bias correction, Studentized bootstrap
# #Studentized bias-corrected bootstrap, Bayesian Bootstrap, Normal approximation
# #sampling: 'fm': full multinomial, 'dm' double multinomial, 'mb': multiple binomial
# #n_b: number of bootstraps. You must provide this parameter if you use any bootstrap method.


#' Confidence Interval for an Index
#'
#' The \code{index_ci} function estimates different kinds of confidence
#' intervals for the supported indices in the SISeg package, under a given
#' framework, starting from the input data. The function returns the confidence
#' intervals and point estimates in a list. It can calculate the confidence
#' intervals for more than one segregation index at the time.
#'
#' \code{index_ci} starts from data, in the form of a \code{matrix}, and
#' calculates the confidence interval and point estimates for segregation
#' indices on the data. There are two main parameters to be chosen.
#'
#' First, all estimations are always done under a framework, representing the
#' data generating process. The available frameworks are
#' \itemize{
#' \item{`Independent_groups` ('ig')}
#' \item{`Full_Multinomial` ('fm') -- corresponding to standard non-parametric
#'     bootstrap}
#' \item{`Independent_Units` ('iu')}
#' }
#' The user can use the full names or the acronym. Framework names are
#' case-sensitive. By default the function will pick the 'Full_Multinomial'
#' framework, corresponding to standard bootstrap. The other frameworks
#' represent stratified sampling by group ('ig') or by units ('iu') with
#' proportionate allocations. If none of the frameworks represent the actual
#' data generating process 'fm' is the safest choice and it is usually
#' conservative with respect to 'ig' and 'iu'.
#'
#' The second choice the user needs to make is about the \code{method} argument.
#' \code{index_ci} can use six different methods to construct confidence
#' intervals and point-estimates:
#' \describe{
#' \item{`Asymp`}{Asymptotic confidence interval based on the delta method.}
#' \item{`Boot`}{Bootstrap percentile confidence interval, with no bias
#'     correction.}
#' \item{`CBoot`}{Bootstrap percentile confidence interval, with bias
#'     correction.}
#' \item{`TBoot`}{Studentized Bootstrap percentile confidence interval.}
#' \item{`BayModel`}{Credible interval and point-estimate based on a simple
#'     conjugate Bayesian model.}
#' \item{`CBayBoot`}{Credible interval and point-estimate based on a simple
#'     conjugate Bayesian model of the index bias (as opposed to a model of the
#'     whole environment).}
#' }
#' Each confidence interval is linked to a specif point-estimate. The methods
#' `Asymp`, `Boot`, `TBoot` are based on the (likely positively-biased)
#' plug-in estimate of the segregation index on the data. `CBoot` applies a
#' bootstrap correction to the plug-in point-estimate in order to diminish the
#' bias. `BayModel` draws upon the estimated Bayesian model and
#' use as point-estimate the expected value of the index based on the model.
#' `CBayBoot` creates a Bayesian model of the index plug-in's bias; it creates
#' the point-estimate by correcting the plug-in estimate with the fit model and
#' use the model's percentile to create credible intervals.
#'
#' All the methods and frameworks will be substantially equivalent in large
#' samples, but in small and mid-size samples they can give very different
#' answers. In very small samples, `CBayBoot` appears to have better
#' performances. In all other cases, `CBoot` provides the best point estimate.
#' Finally, in large and mid-size samples, `TBoot` seems the confidence
#' interval with coverage closest to nominal. In all cases, `Asymp` and `Boot`
#' should be avoided.
#'
#' @section Warnings:
#'
#' In small samples, it is common to receive warnings about few failure in the
#' bootstrap or Bayesian estimations. The warnings are due to some simulated
#' environments having a row (or a column) with only zeros. In these cases, the
#' methods can fail in calculating the standard deviation or even the index
#' value. These cases will therefore be discarded in the inference procedure and
#' the function will produce warnings about this. However, this will be hardly
#' consequential for the inference procedure.
#'
#' @references
#'
#' Allen, Rebecca, Simon Burgess, Russell Davidson, and Frank Windmeijer. 2015.
#' “More Reliable Inference for the Dissimilarity Index of Segregation.”
#' The Econometrics Journal 18 (1): 40–66.
#'
#' Ransom, Michael R. 2000. “Sampling Distributions of Segregation Indexes.”
#' Sociological Methods & Research 28 (4): 454–75.
#'
#' @param env nxk \code{matrix} of \code{numeric}.
#' @param seg_index \code{string} or \code{vector}. Which segregation indices
#'     should be analyzed? The function works for any index available in SISeg
#'     and any index the user manually added.
#' @param n_b \code{numeric}. How many simulations should be returned? Used by
#'     the bootstrap and Bayesian methods.
#' @param framework \code{string}. One of `fm`, `ig`, `iu`.
#' @param method \code{string} or \code{vector} indicating which method(s)
#'     should be used. The methods can be `Asymp`, `Boot`, `CBoot`, `TBoot`,
#'     `BayModel`, `CBayBoot`.
#' @param confidence_level \code{numeric} between 0 and 1 setting the confidence
#'     level for the intervals.
#' @param prior \code{list} of priors to be passed to the Bayesian models. See
#'     \code{\link{r_ig}}, \code{\link{r_iu}} and \code{\link{r_fm}}.
#'
#'
#' @return \code{list} of \code{lists} containing the confidence intervals and
#'     point estimates. The \code{list} is ordered by method and then by
#'     index.
#'
#' @family ci
#'
#' @examples
#' env <- matrix(seq(8), nrow = 2, byrow = TRUE)
#' # confidence interval for a whole lot of indices and methods, just 1 line
#' ci_lot <- index_ci(env, n_b = 5000, seg_index = c("Theil", "D", "Gini",
#'     "Atkinson", "Mutual_info"), method =  c("Asymp", "CBoot", "CBayBoot",
#'     "TBoot"))
#' # D index, CBoot (notice the lower case)
#' ci_lot[["CBoot"]][["d"]]
#' # Atkinso Asymptotic
#' ci_lot[["Asymp"]][["atkinson"]]
#' # Produce warnings
#' set.seed(1233)
#' ci1 <- index_ci(env, n_b = 500, seg_index = c("Mutual_info"), method =
#'     c("CBoot", "CBayBoot", "TBoot"))
#' # Doesn't produce warnings
#' set.seed(1232)
#' ci2 <- index_ci(env, n_b = 500, seg_index = c("Mutual_info"), method =
#'     c("CBoot", "CBayBoot", "TBoot"))
#' # Very similar estimates
#' ci1$CBoot$mutual_info - ci2$CBoot$mutual_info
#' @export
index_ci <- function(env, seg_index = c("D", "Gini", "Mutual_info", "Atkinson",
                                        "Theil"),
            method = c("Asymp", "CBoot", "TBoot", "CBayBoot"), framework = "fm",
            n_b = 5000, confidence_level = 0.95, prior = NULL){
  # Check framework is well specified
  if (framework %in% .index_env$frmews){
    framework <- .index_env$frmews_abbr[framework]
  }
  if (!(framework %in% .index_env$frmews_abbr)){
    stop("Framework not among known frameworks.")
  }
  # Check methods are well specified
  good <- method %in% .index_env$methods
  if (sum(good) < length(good)){
    stop(paste0("The methods ", method[!(good)], " are not among known methods"))
  }
  # Check indices are well specified
  seg_index <- tolower(seg_index)
  good <- seg_index %in% tolower(.index_env$index)
  if (sum(good) < length(good)){
    stop(paste0("The indices ", seg_index[!(good)], " are not among known indices"))
  }
  CIs <- list()
  # Check if confidence level makes sense
  alfa <- 1 - confidence_level
  if ((alfa < 0) | (alfa > 1)){
    stop("confidence_level must be between 0 and 1")
  }
  N <- sum(env)
  #We start by selecting the appropriate functions for bootstrap simulation, SE calculation
  #For SE
  sigma_f <- .index_env$sigmas[[framework]]
  ders_f <- .index_env$ders[[framework]]
  #For bootstrap simulation
  boot_f  <- .index_env$boots[[framework]]
  #Empirical indexes as calculated from the data
  centers <- unlist(lapply(seg_index, function(x) .index_env$indexes[[x]](env)))
  names(centers) <- seg_index
  #Here we calculate ws_var which express how the variance of shrinks with N getting bigger
  #This changes depending on the framework. Used to calculate the final variance
  if (framework == "iu"){
    #ws expresses how many individual each column has
    ws <- colSums(env)
    ws_var <- unlist(lapply(seq(length(ws)), function(x) rep(ws[x], 2)))
  } else {
    ws_var <- N
  }
  #Empirical Vcov matrix from the data and derivative
  Sigma_s <- sigma_f(env)
  #We cannot simply divide by sqrt(N) for mb. This solves it
  Sigma_s <- Sigma_s/ws_var
  lambdas <- lapply(seg_index, function(x) ders_f[[x]](env))
  v_est <- unlist(lapply(seq(length(lambdas)), function(x) as.numeric(
           lambdas[[x]] %*%  Sigma_s %*% t(t(lambdas[[x]])))))
  names(v_est) <- seg_index
  if ("Asymp" %in% method){
    CIs[["Asymp"]] <- lapply(seg_index, function(x) c(centers[x] + qnorm(alfa/2)*
                            sqrt(v_est[x]),  centers[x] - qnorm(alfa/2)*
                            sqrt(v_est[x]), centers[x]))
    names(CIs[["Asymp"]]) <- seg_index
  }
  if ( ("CBoot" %in% method) | ("TBoot" %in% method)
       | ("Boot" %in% method)
       | ("BayModel" %in% method) | ("CBayBoot" %in% method) ){
    #Now we construct the BS CIs!
    #We do this with a simple T-bootstrap
    #Bootstrap of categorical variable = resampling under the same scheme
    if ( ("CBoot" %in% method) | ("Boot" %in% method) | ("TBoot" %in% method)) {
      resample <- boot_f(env, n_b, list_out = TRUE)
      resample_indexes <- list()
    }
    #Calculate variance for studentized Boot
    if (("TBoot" %in% method)  ){
      # Sigma_est_res <- lapply(seq(nrow(resample)), function(x)
      #      sigma_f(matrix(resample[x,], nrow = 2, byrow = T )))
      Sigma_est_res <- lapply(resample,  sigma_f)
      Sigma_est_res <- lapply(Sigma_est_res, function(x) x / ws_var )
      v_est_res <- list()
      Ts <- list()
    }
    #Calculate indexes and Ts on the boot resample
    if ( ("CBoot" %in% method) | ("Boot" %in% method) | ("TBoot" %in% method)) {
      for (i in seg_index){
        ts_na <- rep(FALSE, n_b)
        print(i)
        # Apply the indices to each re-sample
        resample_indexes[[i]] <-  unlist(lapply(resample,
                                                .index_env$indexes[[i]]))
        if ( ("TBoot" %in% method)  ){
          # Calculate the approximate s.d.for each resampled index
          lambdas_est_res <- lapply(resample,  ders_f[[i]] )
          v_est_res[[i]] <- unlist(lapply(seq(length(resample)), function(x)
            as.numeric(lambdas_est_res[[x]] %*%  Sigma_est_res[[x]] %*%
                         t(t(lambdas_est_res[[x]])))))
          remove(lambdas_est_res)
          # Calculate the HDI for the sampled t-distribution
          Ts_temp <- (resample_indexes[[i]] - centers[i])/sqrt(v_est_res[[i]])
          ts_na <- is.na(Ts_temp)
          Ts[[i]] <- HDInterval::hdi(Ts_temp, credMass = confidence_level)
        }
        #Rarely, there may be NA in the resample due to the fact that the
        #indices are undefined when either a column or a row sums to zero
        #The probability of this happening are >0 in some bootstrap schemes.
        #We signal this to the user
        n_na <- length(which((is.na(resample_indexes[[i]]))) | (ts_na ))
        actual_n <- n_b - n_na
        if (actual_n < n_b){
          warning(paste0("Due to NAs, the actual number of bootstrap samples for
                         the bootstrap methods of ", i ," is ", n_b-
                           n_na))
        }
      }}
    #Calculate Bootstrap for the bayesian boot resample
    if ( ("BayModel" %in% method) | ("CBayBoot" %in% method)) {
      resample_bay <- boot_f(env, n_b = n_b, bay = T, prior = prior,
                             list_out = TRUE)
      resample_indexes_bay <- list()
      estimated_indexes_bay <- list()
      for (i in seg_index){
        print(i)
        resample_indexes_bay[[i]] <-  unlist(lapply(resample_bay[[1]],
                                                .index_env$indexes[[i]]))
        estimated_indexes_bay[[i]] <- unlist(lapply(resample_bay[[2]],
                                                .index_env$indexes[[i]]))
        # Rarely, there may be NA in the resample due to the fact that the
        # indices are undefined when either a column or a row sums to zero
        # The probability of this happening are >0 in some bootstrap schemes.
        # We signal this to the user
        n_na <- sum(is.na(resample_indexes_bay[[i]]))
        actual_n <- n_b - n_na
        if (actual_n < n_b){
          warning(paste0("Due to NAs, the actual number of bootstrap samples for
                       the corrected Bayesian method is of ", i ," is ", n_b-
                           n_na))
        }
      }
    }
    #Calculate the simple percentile Bootstrap CI
    if ("Boot" %in% method){
      CIs[["Boot"]] <- lapply(index, function(x) c(hdi(resample_indexes[[x]],
                                      credMass = confidence_level), centers[x]))
      names(CIs[["Boot"]]) <- index
    }
    #Calculate the bias-corrected percentile Bootstrap CI
    if ("CBoot" %in% method){
      centers_correct <- unlist(lapply(seq(length(centers)), function(x)
        2*centers[x] - mean(resample_indexes[[x]], na.rm = TRUE) ))
      centers_correct <- unlist(lapply(centers_correct, function(x) max(0,
                                                                    min(1, x))))
      names(centers_correct) <- seg_index
      CIs[["CBoot"]] <- lapply(seg_index, function(x) c(HDInterval::hdi(2*
                        centers[x] - resample_indexes[[x]], credMass =
                        confidence_level), centers_correct[x]))
      names(CIs[["CBoot"]]) <- seg_index
    }
    #Calculate the studentized Bootstrap CI
    if ("TBoot" %in% method){
      CIs[["TBoot"]] <- lapply(seg_index, function(x) c(centers[x] -
                        Ts[[x]][2]*sqrt(v_est[x]),   centers[x] -
                          Ts[[x]][1]*sqrt(v_est[x]) , centers[x]))
      names(CIs[["TBoot"]]) <- seg_index
    }
    #Calculate the Bayesian Model Index
    if ("BayModel" %in% method){
      CIs[["BayModel"]] <- lapply(seg_index, function(x)
         c(HDInterval::hdi(estimated_indexes_bay[[x]],
         credMass = confidence_level), mean(estimated_indexes_bay[[x]],
                                            na.rm = TRUE)))
      names(CIs[["BayModel"]]) <- seg_index
    }
    #Calculate the Bayesian Bootstrap CI
    if ("CBayBoot" %in% method){
      centers_correct_bay <- unlist(lapply(seq(length(centers)), function(x)
        2*centers[x] - mean(resample_indexes_bay[[x]], na.rm = TRUE) ))
      centers_correct_bay <- unlist(lapply(centers_correct_bay,
                                          function(x) max(0, min(1, x))))
      names(centers_correct_bay) <- seg_index
      CIs[["CBayBoot"]] <- lapply(seg_index, function(x)
                     c(HDInterval::hdi(2*centers[x] - resample_indexes_bay[[x]],
                          credMass = confidence_level), centers_correct_bay[x]))
      names(CIs[["CBayBoot"]]) <- seg_index
    }
  }
  #Nicer names inside the lists
  for (met in method) {
    for (ind in seg_index) {
      CIs[[met]][[ind]] <-
         c(max(0, CIs[[met]][[ind]][1]),
           min(1, CIs[[met]][[ind]][2]), CIs[[met]][[ind]][3])
      names(CIs[[met]][[ind]]) <- c("Lower", "Upper", "Estimate")
      # Just in case of a very strange index
      CIs[[met]][[ind]]['Upper'] <- max(0, CIs[[met]][[ind]]['Upper'])
      CIs[[met]][[ind]]['Lower'] <- min(1, CIs[[met]][[ind]]['Lower'])
      #Due to extreme bias in the measurement, sometime, the CIs obtained are invalid.
      #This happens with TBoot and CTBoot and is due to the fact that the T-distribution obtained
      #May be entirely positive (again, due to bias). When this happens we signal that
      if (CIs[[met]][[ind]]["Lower"] >= CIs[[met]][[ind]]['Upper']){
        CIs[[met]][[ind]] <- c(NA, NA, CIs[[met]][[ind]][3])
        print("The estimation produced an invalid CI: ")
        print(met)
        print(ind)
      }
    }
  }
  return(CIs)
}

#CI for differences between indexes
#Notice, that the function supposes that the matrices have the same number of columns and in general
#The same columns. It doesn't make a lot of sense to compare tables with different categories
#The function also supposes that the tables have the same sampling scheme. This is not strictly necessary
#env1 and env2 are the data table to be compared. Method regards how to calculate the CI
#index is a list of indexes to calculate based on the data

#' Confidence Interval for the Difference of Two Indices
#'
#' The \code{index_difference_ci} function estimates different kinds of
#' confidence intervals for the difference between indices
#' in two different datasets, under a given framework. The function returns the
#' confidence intervals and point estimates in a list. It can calculate the
#' confidence intervals for more than one segregation index at the time. It
#' assumes that the two data creation processes were independent.
#'
#' \code{index_difference_ci} starts from data, in the form of two
#' \code{matrices}, and calculates the confidence interval and point estimates
#' for difference of segregation indices in the two datasets. There are two main
#' parameters to be chosen. It assumes that the two data creation processes
#' were independent.
#'
#' First, all estimations are always done under a framework, representing the
#' data generating process. The available frameworks are
#' \itemize{
#' \item{`Independent_groups` ('ig')}
#' \item{`Full_Multinomial` ('fm') -- corresponding to standard non-parametric
#'     bootstrap}
#' \item{`Independent_Units` ('iu')}
#' }
#' The user can use the full names or the acronym. Framework names are
#' case-sensitive. By default the function will pick the 'Full_Multinomial'
#' framework, corresponding to standard bootstrap. The other frameworks
#' represent stratified sampling by group ('ig') or by units ('iu') with
#' proportionate allocations. If none of the frameworks represent the actual
#' data generating process 'fm' is the safest choice and it is usually
#' conservative with respect to 'ig' and 'iu'.
#'
#' The second choice the user needs to make is about the \code{method} argument.
#' \code{index_difference_ci} can use six different methods to construct
#' confidence intervals and point-estimates:
#' \describe{
#' \item{`Asymp`}{Asymptotic confidence interval based on the delta method.}
#' \item{`Boot`}{Bootstrap percentile confidence interval, with no bias
#'     correction.}
#' \item{`CBoot`}{Bootstrap percentile confidence interval, with bias
#'     correction.}
#' \item{`TBoot`}{Studentized Bootstrap percentile confidence interval.}
#' \item{`BayModel`}{Credible interval and point-estimate based on two simple
#'     conjugate Bayesian models.}
#' \item{`CBayBoot`}{Credible interval and point-estimate based on two simple
#'     conjugate Bayesian models of the index bias (as opposed to a model of the
#'     whole environment) in the two environments.}
#' }
#' Each confidence interval is linked to a specif point-estimate. The methods
#' `Asymp`, `Boot`, `TBoot` are based on the plug-in estimate of the segregation
#' index on the data. `CBoot` applies a
#' bootstrap correction to the plug-in point-estimate in order to diminish the
#' bias. `BayModel` draws upon two estimated Bayesian model (one per dataset)
#' and use as point-estimate the expected value of the indices' difference based
#' on the models.
#' `CBayBoot` creates two Bayesian models (one per environment) of the index
#' plug-in's bias; it creates the point-estimate of the difference between the
#' the bias-corrected plug-in estimates and use the models' percentile to create
#' credible intervals.
#'
#' All the methods and frameworks will be substantially equivalent in large
#' samples, but in small and mid-size samples they can give very different
#' answers. Based on what we know about the estimation of indices, the
#' reccomandation is to use `CBayBoot` in small samples.
#' In all other cases, `CBoot` provides the best point estimate.
#' Finally, in large and mid-size samples, `TBoot` seems the confidence
#' interval with coverage closest to nominal. In all cases, `Asymp` and `Boot`
#' should be avoided.
#'
#' @section Warnings:
#'
#' In small samples, it is common to receive warnings about few failure in the
#' bootstrap or Bayesian estimations. The warnings are due to some simulated
#' environments having a row (or a column) with only zeros. In these cases, the
#' methods can fail in calculating the standard deviation or even the index
#' value. These cases will therefore be discarded in the inference procedure and
#' the function will produce warnings about this. However, this will be hardly
#' consequential for the final results.
#'
#' @references
#'
#' Allen, Rebecca, Simon Burgess, Russell Davidson, and Frank Windmeijer. 2015.
#' “More Reliable Inference for the Dissimilarity Index of Segregation.”
#' The Econometrics Journal 18 (1): 40–66.
#'
#' @param env1 nxk1 \code{matrix} of \code{numeric}.
#' @param env2 nxk1 \code{matrix} of \code{numeric}.
#' @param seg_index \code{string} or \code{vector}. Which segregation indices
#'     should be analyzed? The function works for any index available in SISeg
#'     and any index the user manually added.
#' @param n_b \code{numeric}. How many simulations should be returned? Used by
#'     the bootstrap and Bayesian methods.
#' @param framework \code{string}. One of `fm`, `ig`, `iu`.
#' @param method \code{string} or \code{vector} indicating which method(s)
#'     should be used. The methods can be `Asymp`, `Boot`, `CBoot`, `TBoot`,
#'     `BayModel`, `CBayBoot`.
#' @param confidence_level \code{numeric} between 0 and 1 setting the confidence
#'     level for the intervals.
#' @param prior1 \code{list} of priors to be passed to the Bayesian models for
#'      of the \code{env1}. See \code{\link{r_ig}}, \code{\link{r_iu}} and
#'      \code{\link{r_fm}}.
#' @param prior2 \code{list} of priors to be passed to the Bayesian models for
#'      of the \code{env2}. See \code{\link{r_ig}}, \code{\link{r_iu}} and
#'      \code{\link{r_fm}}.
#'
#'
#' @return \code{list} of \code{lists} containing the confidence intervals and
#'     point estimates. The \code{list} is ordered by method and then by
#'     index.
#'
#' @family ci
#'
#' @examples
#' env <- matrix(seq(8), nrow = 2, byrow = TRUE)
#' env2 <- matrix(seq(9,16), nrow = 2, byrow = TRUE)
#' # confidence interval for a whole lot of indices and methods, just 1 line
#' ci_lot <- index_difference_ci(env, env2, n_b = 5000, seg_index = c("Theil",
#'     "D", "Gini", "Atkinson", "Mutual_info"), method =  c("Asymp", "CBoot",
#'     "CBayBoot", "TBoot"))
#' # D index, CBoot (notice the lower case)
#' ci_lot[["CBoot"]][["d"]]
#' # Atkinso Asymptotic
#' ci_lot[["Asymp"]][["atkinson"]]
#' # Produce warnings
#' set.seed(1233)
#' ci1 <- index_difference_ci(env, env2, n_b = 250, seg_index =
#'     c("Mutual_info"), method = c("CBoot", "CBayBoot", "TBoot"))
#' # Doesn't produce warnings
#' set.seed(1230)
#' ci2 <- index_difference_ci(env, env2, n_b = 250, seg_index =
#'     c("Mutual_info"), method = c("CBoot", "CBayBoot", "TBoot"))
#' # Very similar estimates
#' ci1$CBoot$mutual_info - ci2$CBoot$mutual_info
#' @export
index_difference_ci <- function(env1, env2, seg_index = c("D", "Gini",
                       "Mutual_info", "Atkinson", "Theil"),
                       method = c("Asymp", "CBoot", "TBoot", "CBayBoot"),
                       framework = "fm", n_b = 5000, confidence_level = 0.95,
                       prior1 = NULL, prior2 = NULL){
  # Check framework is well specified
  if (framework %in% .index_env$frmews){
    framework <- .index_env$frmews_abbr[framework]
  }
  if (!(framework %in% .index_env$frmews_abbr)){
    stop("Framework not among known frameworks.")
  }
  # Check methods are well specified
  good <- method %in% .index_env$methods
  if (sum(good) < length(good)){
    stop(paste0("The methods ", method[!(good)], " are not among known methods"))
  }
  # Check indices are well specified
  seg_index <- tolower(seg_index)
  good <- seg_index %in% tolower(.index_env$index)
  if (sum(good) < length(good)){
    stop(paste0("The indices ", seg_index[!(good)], " are not among known indices"))
  }
  CIs <- list()
  alfa <- 1 - confidence_level
  if ((alfa <= 0) | (alfa >= 1)){
    stop("confidence_level must be between 0 and 1")
  }
  N1 <- sum(env1)
  N2 <- sum(env2)
  #We start by selecting the appropriate functions for bootstrap simulation, SE
  #calculation. For SE
  sigma_f <-.index_env$sigmas[[framework]]
  ders_f <- .index_env$ders[[framework]]
  #For bootstrap simulation
  boot_f  <- .index_env$boots[[framework]]
  #Empirical indexes as calculated from the data
  centers1 <- unlist(lapply(seg_index, function(x) .index_env$indexes[[x]]
                            (env1) ))
  names(centers1) <- seg_index
  centers2 <- unlist(lapply(seg_index, function(x) .index_env$indexes[[x]]
                            (env2) ))
  names(centers2) <- seg_index
  centers_diff <- centers1 - centers2
  names(centers_diff) <- seg_index
  # Here we calculate ws_var which express how the variance of shrinks with N
  # getting bigger
  # This changes depending on the framework scheme. Used to calculate the final
  # variance
  if (framework == "iu"){
    #ws expresses how many individual each column has
    ws1 <- colSums(env1)
    ws2 <- colSums(env2)
    ws_var1 <- unlist(lapply(seq(length(ws1)), function(x) rep(ws1[x], 2)))
    ws_var2 <- unlist(lapply(seq(length(ws2)), function(x) rep(ws2[x], 2)))
  } else {
    ws_var1 <- N1
    ws_var2 <- N2
    # ws_var1 <- rep(N1, k)
    # ws_var2 <- rep(N2, k)
  }
  #Empirical Vcov matrix from the data and derivative
  Sigma_s1 <- sigma_f(env1)
  Sigma_s2 <- sigma_f(env2)
  #We cannot simply divide by sqrt(N) at the end for mb. This solves it
  Sigma_s1 <- Sigma_s1/ws_var1
  Sigma_s2 <- Sigma_s2/ws_var2
  lambdas1 <- lapply(seg_index, function(x) ders_f[[x]](env1))
  lambdas2 <- lapply(seg_index, function(x) ders_f[[x]](env2))
  v_est1 <- unlist(lapply(seq(length(lambdas1)), function(x)
               as.numeric(lambdas1[[x]] %*%  Sigma_s1 %*% t(t(lambdas1[[x]])))))
  v_est2 <- unlist(lapply(seq(length(lambdas2)), function(x)
               as.numeric(lambdas2[[x]] %*%  Sigma_s2 %*% t(t(lambdas2[[x]])))))
  names(v_est1) <- seg_index
  names(v_est2) <- seg_index
  if ("Asymp" %in% method){
    CIs[["Asymp"]] <- lapply(seg_index, function(x) c((centers1[x] - centers2[x])
                    + qnorm(alfa/2)*sqrt(v_est1[x] + v_est2[x]),
                    (centers1[x] - centers2[x]) -
                    qnorm(alfa/2)*sqrt(v_est1[x] + v_est2[x]), centers_diff[x]))
    names(CIs[["Asymp"]]) <- seg_index
  }
  if ( ("CBoot" %in% method) | ("TBoot" %in% method)
       | ("Boot" %in% method)
       | ("BayModel" %in% method) | ("CBayBoot" %in% method) ){
    if (!(n_b)){
      stop("You need to specify the n_b parameter to use bootstrap CIs")
    }
    #Now we construct the BS CIs!
    #Bootstrap of categorical variable = resampling under the same scheme
    if ( ("CBoot" %in% method) | ("Boot" %in% method) |
         ("TCBoot" %in% method) | ("TBoot" %in% method)  ) {
      resample1 <- boot_f(env1, n_b, list_out = TRUE)
      resample2 <- boot_f(env2, n_b, list_out = TRUE)
      resample_indexes1 <- list()
      resample_indexes2 <- list()
    }
    #Calculate variance for studentized Boot
    if ( ("TCBoot" %in% method) | ("TBoot" %in% method)  ){
      Sigma_est_res1 <- lapply(resample1, sigma_f)
      Sigma_est_res2 <- lapply(resample2, sigma_f)
      Sigma_est_res1 <- lapply(seq(length(resample1)),
                               function(x) Sigma_est_res1[[x]] / ws_var1 )
      Sigma_est_res2 <- lapply(seq(length(resample2)),
                               function(x) Sigma_est_res2[[x]] / ws_var2 )
      v_est_res1 <- list()
      v_est_res2 <- list()
      Ts <- list()
    }
    #Calculate indexes and Ts on the boot resample
    if ( ("CBoot" %in% method) | ("Boot" %in% method) |
         ("TBoot" %in% method)  ) {
      for (i in seg_index){
        ts_na <- rep(FALSE, n_b)
        print(i)
        resample_indexes1[[i]] <-  unlist(lapply(resample1,
                                                 .index_env$indexes[[i]]))
        resample_indexes2[[i]] <-  unlist(lapply(resample2,
                                                 .index_env$indexes[[i]]))
        if ("TBoot" %in% method){
          lambdas_est_res1 <- lapply(resample1,  ders_f[[i]] )
          lambdas_est_res2 <- lapply(resample2,  ders_f[[i]] )
          v_est_res1[[i]] <- unlist(lapply(seq(length(resample1)), function(x)
            as.numeric(lambdas_est_res1[[x]] %*%  Sigma_est_res1[[x]] %*%
                         t(t(lambdas_est_res1[[x]])))))
          v_est_res2[[i]] <- unlist(lapply(seq(length(resample2)), function(x)
             as.numeric(lambdas_est_res2[[x]] %*% Sigma_est_res2[[x]] %*%
                          t(t(lambdas_est_res2[[x]])))))
          remove(lambdas_est_res1, lambdas_est_res2)
          Ts_temp <- ((resample_indexes1[[i]] - resample_indexes2[[i]]) -
                      ( centers1[i] -centers2[i] ))/
            (sqrt(v_est_res1[[i]]+v_est_res2[[i]]))
          Ts[[i]] <- HDInterval::hdi(Ts_temp, confidence_level)
          ts_na <- is.na(Ts_temp)
        }
        # Rarely, there may be NA in the resample due to the fact that the
        # indices are undefined when either a column or a row sums to zero
        # The probability of this happening are >0 in some bootstrap schemes.
        # We signal this to the user. All indices are equally affected by this
        n_na <- sum( (is.na(resample_indexes[[i]])) | (ts_na) )
        actual_n <- n_b - n_na
        if (actual_n < n_b){
          warning(paste0("Due to NAs, the actual number of bootstrap samples for
                         the bootstrap methods of ", i ," is ", n_b-
                           n_na))
        }
      }}
    #Calculate Bootstrap for the bayesian boot resample
    if ( ("BayModel" %in% method) | ("CBayBoot" %in% method)) {
      resample_bay1 <- boot_f(env1, n_b, bay = TRUE,
                              prior = prior1, list_out = TRUE)
      resample_bay2 <- boot_f(env2, n_b, bay = TRUE,
                              prior = prior2, list_out = TRUE)
      resample_indexes_bay1 <- list()
      resample_indexes_bay2 <- list()
      estimated_indexes_bay1 <- list()
      estimated_indexes_bay2 <- list()
      resample_diff_bay <- list()
      estimated_diff_bay <- list()
      for (i in seg_index){
        print(i)
        resample_indexes_bay1[[i]] <-  unlist(lapply( resample_bay1$resample,
                                      .index_env$indexes[[i]]))
        resample_indexes_bay2[[i]] <-  unlist(lapply( resample_bay2$resample,
                                                      .index_env$indexes[[i]]))
        estimated_indexes_bay1[[i]] <- unlist(lapply(  resample_bay1$ps,
                                                       .index_env$indexes[[i]]))
        estimated_indexes_bay2[[i]] <- unlist(lapply(  resample_bay2$ps,
                                                       .index_env$indexes[[i]]))
        resample_diff_bay[[i]] <- resample_indexes_bay1[[i]] -
                                          resample_indexes_bay2[[i]]
        estimated_diff_bay[[i]] <- estimated_indexes_bay1[[i]] -
                                          estimated_indexes_bay2[[i]]
      }
      #Rarely, there may be NA in the resample due to the fact that the indices
      #are undefined when either a column or a row sums to zero
      #The probability of this happening are >0 in some bootstrap schemes.
      # We signal this to the user. All indices are equally affected by this
      n_na <- sum((is.na(resample_indexes_bay1[[i]])) |
                       (is.na(resample_indexes_bay2[[i]])))
      actual_n <- n_b - n_na
      if (actual_n < n_b){
        warning(paste0("Due to NAs, the actual number of bootstrap samples for
                       the corrected Bayesian method is of ", i ," is ", n_b-
                         n_na))
      }
    }
    #Calculate the simple percentile Bootstrap CI
    if ("Boot" %in% method){
      CIs[["Boot"]] <- lapply(seg_index, function(x)
        c(HDInterval::hdi(resample_indexes1[[x]] - resample_indexes2[[x]],
                          credMass = confidence_level), centers_diff[x]))
      names(CIs[["Boot"]]) <- seg_index
    }
    #Calculate the bias-corrected percentile Bootstrap CI
    if ("CBoot" %in% method){
      centers_correct1 <- unlist(lapply(seg_index, function(x) 2*centers1[x] -
                                  mean(resample_indexes1[[x]], na.rm = TRUE) ))
      names(centers_correct1) <- seg_index
      centers_correct2 <- unlist(lapply(seg_index, function(x) 2*centers2[x] -
                                  mean(resample_indexes2[[x]], na.rm = TRUE) ))
      names(centers_correct2) <- seg_index
      centers_correct_diff <- centers_correct1 -centers_correct2
      centers_correct_diff <- unlist(lapply(centers_correct_diff,
                                            function(x) max(-1, min(1, x))))
      names(centers_correct_diff) <- seg_index
      CIs[["CBoot"]] <- lapply(seg_index, function(x)
        c( HDInterval::hdi(2*centers_diff[x] - (resample_indexes1[[x]] -
                   resample_indexes2[[x]]), credMass = confidence_level),
                                                centers_correct_diff[x]))
      names(CIs[["CBoot"]]) <- seg_index
    }
    #Calculate the studentized Bootstrap CI
    if ("TBoot" %in% method){
      CIs[["TBoot"]] <- lapply(seg_index, function(x)
        c(centers_diff[x]  - Ts[[x]][2]*sqrt(v_est1[x] + v_est2[x]),
          centers_diff[x] - Ts[[x]][1]*sqrt(v_est1[x] + v_est2[x]),
          centers_diff[x]) )
      names(CIs[["TBoot"]]) <- seg_index
    }
    #Calculate the Bayesian Model Index
    if (("BayModel" %in% method)){
      CIs[["BayModel"]] <- lapply(seg_index, function(x)
                          c(HDInterval::hdi(estimated_diff_bay[[x]],
                          credMass = confidence_level),
                          mean(estimated_diff_bay[[x]])))
      names(CIs[["BayModel"]]) <- seg_index
    }
    #Calculate the Bayesian Bootstrap CI
    if ("CBayBoot" %in% method){
      centers_correct_bay1 <- unlist(lapply(seg_index, function(x) 2*centers1[x]
                             - mean(resample_indexes_bay1[[x]], na.rm = TRUE) ))
      names(centers_correct_bay1) <- seg_index
      centers_correct_bay2 <- unlist(lapply(seg_index, function(x) 2*centers2[x] -
                              mean(resample_indexes_bay2[[x]], na.rm = TRUE) ))
      names(centers_correct_bay2) <- seg_index
      centers_correct_diff <- centers_correct_bay1 -centers_correct_bay2
      centers_correct_diff <- unlist(lapply(centers_correct_diff,
                                            function(x) max(-1, min(1, x))))
      names(centers_correct_diff) <- seg_index
      CIs[["CBayBoot"]] <- lapply(seg_index, function(x)
        c(HDInterval::hdi(2*centers_diff[x] - resample_diff_bay[[x]],
                        credMass = confidence_level), centers_correct_diff[x]))
      names(CIs[["CBayBoot"]]) <- seg_index
    }
  }
  #Nicer names inside the lists
  for (met in method){
    for (ind in seg_index){
      CIs[[met]][[ind]] <- c(max(-1, CIs[[met]][[ind]][1]),
                            min(1, CIs[[met]][[ind]][2]), CIs[[met]][[ind]][3])
      # Due to extreme bias in the measurement, sometime, the CIs obtained are
      # invalid. This happens TBoot and CTBoot and is due to the fact that the
      # T-distribution obtained
      # May be entirely positive (again, due to bias). When this happens we
      # signal that
      if (CIs[[met]][[ind]][1] == CIs[[met]][[ind]][2]){
        CIs[[met]][[ind]] <- c(NA, NA, CIs[[met]][[ind]][3])
        print("The estimation produced an invalid CI: ")
        print(paste(met, ind))
      }
      names(CIs[[met]][[ind]]) <- c("Lower", "Upper", "Estimate")
    }
  }
  return(CIs)
}
