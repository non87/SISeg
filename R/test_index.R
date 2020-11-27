#### TEST THE INFERENCE METHODS THROUGH ALMOST ASYMPTOTIC CASE


#' Test Inference on New Indices
#'
#' The user can define new segregation indices and SISeg will perform
#' inference on those automatically. After using the function
#' \code{\link{add_index}} to add a new index to the supported indices in
#' SISeg, \code{test_index} will perform automatic checks that the new
#' index and its related functions work as expected. The tests are based on
#' asymptotic theory and, for this reason, it will take time to run a full
#' battery of tests.
#'
#' The inference techniques in SISeg are general and can be applied to any
#' (derivable) index. The user can specify a new index and then use the
#' inference functions in SISeg to get inference for the newly specified index.
#' The function \code{\link{add_index}} lets the user define new indices and
#' add them to SISeg list of supported indices. \code{test_index} tests that
#' the inference on the added index works as it should. It provides three
#' batteries of tests:
#' \describe{
#' \item{`asymp'}{Tests about the asymptotic normality of the index. Use these
#'      to test that the derivative functions for the new index are correct.}
#' \item{estimate}{Tests to check that the confidence interval for the index
#'      in a dataset works as expected. These tests check the coverage of the
#'      bootstrap and asymptotic confidence intervals.}
#' \item{difference}{Tests the coverage of the confidence intervals for the
#'      difference between indices in two environments. By default the
#'      function will not run these tests. The argument is that if the
#'      estimate tests succeed, the difference tests will also succeed.}
#' }
#'
#' All tests are based on the inference performances on a very convenient case.
#' That is, the tests are perfomed by checking the behavior of the inferential
#' methods on a simple environment with 350,000 individuals and only four
#' columns. In this environment, asymptotic normality is reached for the
#' indices that SISeg natively supports and it should be reached for any
#' well-behaved function. Therefore, the tests consist in checking that
#' asymptotic normality holds for this environment and that the different
#' techniques reach the nominal coverage.
#'
#' The user should be ware that running the entire set of tests will take few
#' hours, with the asymptotic tests only taking few minutes.
#'
#' @param seg_index A \code{string} with the name of the index to be tested.
#' @param n_b \code{numeric} How many simulations in the tests? The number of
#'     simulations in the estimate and difference tests will be capped to 1,000
#'     if n_b is bigger than 1,000.
#' @param asymp \code{logical} Should the asymptotic tests be run?
#' @param estimate \code{logical} Should the estimate tests be run?
#' @param difference \code{logical} Should the difference tests be run?
#' @param graph \code{logical} Should the sampling and asymptotic distribution
#'      be compared in a plot? Ignored if \code{asymp} is FALSE.
#' @param framework \code{vector} or `all`. If vector, it should contain the
#'      names of the framework to be tested.
#'
#' @return \code{list} containing the results of the tests. If \code{asymp} is
#'      \code{TRUE}, it will contain the results of the Kolgomorov-Smirnov test
#'      against normality. A successful test will return a large p-value. If
#'      \code{asymp} is \code{TRUE} and \code{graph} is also \code{TRUE},
#'      the function will also produce plots
#'      showing the sampling distribution for the tested index and the
#'      asymptotic normal distribution for the index -- the two distribution
#'      should be very similar. If \code{estimate} is \code{TRUE}, the
#'      \code{list} will contain the coverage of the different methods for the
#'      point-estimate. If \code{framework} is \code{TRUE}, the \code{list}
#'      will contain the coverage of the different methods for the difference.
#'
#' @examples
#' # We can use the function to test indices that are already in SISeg
#' # Test the Theil index
#' # Asymptotic test are quick
#' test_index("Theil", estimate = FALSE, asymp = TRUE)
#' # Bootstrap tests are not. We only do 10 simulations
#' test_index("Theil", n_b = 10, asymp = FALSE, estimate = TRUE)
#' # The actual tests take a long time.
#' \donttest{test_index("Theil")}
#'
#' @export
test_index <- function(seg_index, n_b = 100000, asymp = TRUE, estimate = TRUE,
                       difference = FALSE, framework = "all", graph = TRUE){
  final_results <- list()
  # Create an explicit vector of frameworks
  if ((sum( framework == "all")/ length(framework)) == 1){
    framework <- c("ig", "fm", "iu")
    }
  # Collect functions to be tested
  seg_index <- tolower(seg_index)
  index_f <- .index_env$indexes[[seg_index]]
  # Check we have everything in order
  if (is.null(index_f)){stop("The segregation index is not among supported
                             indices")}
  if (asymp) {
    der_f <- lapply(framework, function(x) .index_env$ders[[x]][[seg_index]])
    good <- unlist(lapply(der_f, is.null))
    if (sum(good)){stop("Not all derivatives found")}
    names(der_f) <- framework
  }
  # Boot functions
  boot_f <- lapply(framework, function(x) .index_env$boots[[x]])
  names(boot_f) <- framework
  # Sigma functions
  sigma_f <- lapply(framework, function(x) .index_env$sigmas[[x]])
  names(sigma_f) <- framework
  #For drawing
  x <- seq(-1, 1, 0.00001)
  n_f <- 150000
  n_m <- 200000
  N <- n_f + n_m
  x <- seq(-1, 1, 0.00001)
  #These are arbitrary probability vector
  P_1 =  c(0.7, 0.2, 0.1)
  P_2 =  c(0.2, 0.5, 0.3)
  #Variance under this sampling scheme
  real_matrix <- matrix(c(n_m/N* P_1, n_f/N* P_2), byrow = TRUE, nrow = 2)
  center_real <-  index_f(real_matrix)
  if (asymp){
    final_results[['asymp']] <- list()
    for (frm in framework){
      print(paste0("Asymptotic Test for ", frm))
      Sigma <- sigma_f[[frm]](real_matrix*N)
      if (frm == "iu"){
        #ws expresses how many individual each column has
        ws <- colSums(real_matrix*N)
        ws <- rep(ws, each = 2)
      } else {
        ws <- N
      }
      # Just for iu we need to divide
      Sigma <- Sigma/ws
      lambda <- der_f[[frm]](real_matrix*N)
      v <- (lambda %*% Sigma %*% t(t(lambda)))[1]
      m <- boot_f[[frm]](real_matrix*N, n_b = n_b, list_out = TRUE)
      #For HT we save the actual values in a vector
      resample_ind <- unlist( lapply(m, index_f))
      y <- stats::dnorm(x, center_real, sqrt(v))
      final_results[['asymp']][[frm]] <-
                   stats::ks.test((resample_ind-center_real)/(sqrt(v)), "pnorm")
      if (graph){
        title = paste(seg_index, frm)
        graphics::hist(resample_ind, breaks = 100, freq = F, main = title,
                       xlab = seg_index)
        graphics::lines(x,y, col = "red")
      }
    }
  }
  if (difference){
    #Asymptotic CI. Should work well in this case
    as_CIs <- list()
    #Bootstrap CI, the real test
    T_CIs <- list()
    #Bootstrap CI, corrected
    corr_CIs <- list()
    #Number of bootstrap
    for (frm in framework){
      #We only test part of the simulated data. n_b * 2 because we need to
      # construct double the samples
      m <- boot_f[[frm]](real_matrix*N, n_b = min(n_b, 1000)*2, list_out = TRUE)
      n_diff <- min(n_b, 1000)
      m1 <- m[1:n_diff]
      m2 <- m[(n_diff+1):(2*n_diff)]
      n_b_2 <- 1000
      # Calculate the different CIs for each simulated sample
      all_cis <- mapply(function(x,y) index_difference_ci(x, y,
                        seg_index = seg_index, n_b = n_b_2,
                        framework = frm, method = c("Asymp", "TBoot", "CBoot")),
                        x = m1, y = m2)
      # mapply returns the results in a strange shape...which is actually rather
      # convenient to Calculate the coverage of the CIs
      as_CIs[[frm]] <- unlist(lapply(all_cis["Asymp",],
                             function(x) (0 >= x[[1]][1]) &  (0 <= x[[1]][2])))
      as_CIs[[frm]] <- sum(as_CIs[[frm]])/min(n_b, 1000)
      T_CIs[[frm]] <- unlist(lapply(all_cis["TBoot",],
                             function(x) (0 >= x[[1]][1]) &  (0 <= x[[1]][2])))
      T_CIs[[frm]] <- sum(T_CIs[[frm]])/min(n_b, 1000)
      corr_CIs[[frm]] <- unlist(lapply(all_cis["CBoot",],
                             function(x) (0 >= x[[1]][1]) &  (0 <= x[[1]][2])))
      corr_CIs[[frm]] <- sum(corr_CIs[[frm]])/min(n_b, 1000)
    }
    final_results[['difference']][['estimate_asymp']] <- as_CIs
    final_results[['difference']][['boot_T_asymp']] <- T_CIs
    final_results[['difference']][['boot_C_asymp']] <- corr_CIs
  }
  if (estimate){
    final_results[['estimate']] <- list()
    #Asymptotic CI. Should work well in this case
    as_CIs <- list()
    #Bootstrap CI, the real test
    T_CIs <- list()
    #Bootstrap CI, corrected
    corr_CIs <- list()
    #Number of bootstrap
    for (frm in framework){
      print(frm)
      #We only test part of the simulated data
      m <- boot_f[[frm]](real_matrix*N, n_b = min(n_b, 1000), list_out = TRUE)
      n_b_2 <- 1000
      # Calculate the different CIs for each simulated sample
      all_cis <- lapply(m, function(x) index_ci(x, seg_index = seg_index,
                                                n_b = n_b_2, framework = frm,
                                                method = c("Asymp", "TBoot", "CBoot")))
      # Calculate the coverage of the CIs
      as_CIs[[frm]] <- unlist(lapply(all_cis, function(x)
        (center_real >= x[["Asymp"]][[seg_index]][1]) &
          (center_real <= x[["Asymp"]][[seg_index]][2])))
      as_CIs[[frm]] <- sum(as_CIs[[frm]])/min(n_b, 1000)
      T_CIs[[frm]] <- unlist(lapply(all_cis, function(x)
        (center_real >= x[["TBoot"]][[seg_index]][1]) &
          (center_real <= x[["TBoot"]][[seg_index]][2])))
      T_CIs[[frm]] <- sum(T_CIs[[frm]])/min(n_b, 1000)
      corr_CIs[[frm]] <- unlist(lapply(all_cis, function(x)
        (center_real >= x[["CBoot"]][[seg_index]][1]) &
          (center_real <= x[["CBoot"]][[seg_index]][2])))
      corr_CIs[[frm]] <- sum(corr_CIs[[frm]])/min(n_b, 1000)
    }
    final_results[['estimate']][['estimate_asymp']] <- as_CIs
    final_results[['estimate']][['boot_T_asymp']] <- T_CIs
    final_results[['estimate']][['boot_C_asymp']] <- corr_CIs
  }
  return(final_results)
}





