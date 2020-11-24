#CI for one index
#mat: data in matrix form (tested with only two rows!)
#index: which index are you interested in? "D", "gini", "MI" & "Atk" are supported
#method: Bootstrap, Bootstrap with bias correction, Studentized bootstrap
#Studentized bias-corrected bootstrap, Bayesian Bootstrap, Normal approximation
#sampling: 'fm': full multinomial, 'dm' double multinomial, 'mb': multiple binomial
#n_b: number of bootstraps. You must provide this parameter if you use any bootstrap method.
index_ci <- function(mat, index = c("D", "gini", "MI", "Atk", "norm_MI"),
            method = c("Norm", "CBoot", "TBoot", "CBayBoot"), framework = "fm",
            n_b = 5000, confidence_level = 0.95, prior = NULL){
  if (framework %in% .index_env$frmews){
    framework <- .index_env$frmews_abbr[framework]
  }
  if (!(framework %in% .index_env$frmews_abbr)){
    stop("Framework not among known frameworks.")
  }
  CIs <- list()
  alfa <- 1 - confidence_level
  k <- ncol(mat)
  N <- sum(mat)
  #We start by selecting the appropriate functions for bootstrap simulation, SE calculation
  #For SE
  sigma_f <- .index_env$sigmas[[framework]]
  ders_f <- .index_env$ders[[framework]]
  #For bootstrap simulation
  boot_f  <- .index_env$boots[[framework]]
  #Empirical indexes as calculated from the data
  centers <- unlist(lapply(index, function(x) .index_env$indexes[[x]](mat) ))
  names(centers) <- index
  #Here we calculate ws_var which express how the variance of shrinks with N getting bigger
  #This changes depending on the sampling scheme. Used to calculate the final variance
  if (sampling == "iu"){
    #ws expresses how many individual each column has
    ws <- colSums(mat)
    ws_var <- unlist(lapply(seq(length(ws)), function(x) rep(ws[x], 2)))
  } else {
    ws_var <- rep(N, k)
  }
  #Empirical Vcov matrix from the data and derivative
  Sigma_s <- sigma_f(mat)
  #We cannot simply divide by sqrt(N) for mb. This solves it
  Sigma_s <- Sigma_s/ws_var
  lambdas <- lapply(index, function(x) ders_f[[x]](mat))
  v_est <- unlist(lapply(seq(length(lambdas)), function(x) as.numeric(
           lambdas[[x]] %*%  Sigma_s %*% t(t(lambdas[[x]])))))
  names(v_est) <- index
  if ("Norm" %in% method){
    CIs[["Norm"]] <- lapply(index, function(x) c(centers[x] + qnorm(alfa/2)*
                            sqrt(v_est[x]),  centers[x] - qnorm(alfa/2)*
                            sqrt(v_est[x]), centers[x]))
    names(CIs[["Norm"]]) <- index
  }
  if ( ("CBoot" %in% method) | ("TBoot" %in% method)
       | ("Boot" %in% method)
       | ("BayModel" %in% method) | ("CBayBoot" %in% method) ){
    #Now we construct the BS CIs!
    #We do this with a simple T-bootstrap
    #Bootstrap of categorical variable = resampling under the same scheme
    if ( ("CBoot" %in% method) | ("Boot" %in% method) |
         ("TCBoot" %in% method) | ("TBoot" %in% method)  ) {
      resample <- boot_f(mat, n_b)
      resample_indexes <- list()
    }
    #Calculate variance for studentized Boot
    if ( ("TCBoot" %in% method) | ("TBoot" %in% method)  ){
      Sigma_est_res <- lapply(seq(nrow(resample)), function(x)
           sigma_f(matrix(resample[x,], nrow = 2, byrow = T )))
      Sigma_est_res <- lapply(seq(nrow(resample)),
                              function(x) Sigma_est_res[[x]] / ws_var )
      v_est_res <- list()
      Ts <- list()
    }
    #Calculate indexes and Ts on the boot resample
    if ( ("CBoot" %in% method) | ("Boot" %in% method) |
        ("TBoot" %in% method)  ) {
      for (i in index){
        print(i)
        # Apply the indices to each re-sample
        resample_indexes[[i]] <-  unlist(lapply(  seq(nrow(resample)),function(x)
          .index_env$indexes[[i]](matrix(c(resample[x,]), nrow = 2, byrow = T))
          ))
        if ( ("TBoot" %in% method)  ){
          # Calculate the approximate s.d.for each resampled index
          lambdas_est_res <- lapply(  seq(nrow(resample)), function(x)
            ders_f[[i]](matrix(resample[x,], nrow = 2, byrow = T)) )
          v_est_res[[i]] <- unlist(lapply(seq(nrow(resample)), function(x)
            as.numeric(lambdas_est_res[[x]] %*%  Sigma_est_res[[x]] %*%
                         t(t(lambdas_est_res[[x]])))))
          remove(lambdas_est_res)
          # Calculate the HDI for the sampled t-distribution
          Ts_temp <- (resample_indexes[[i]] - centers[i])/sqrt(v_est_res[[i]])
          Ts[[i]] <- HDInterval::hdi(Ts_temp, credMass = confidence_level)
        }
        #Rarely, there may be NA in the resample due to the fact that the
        #indices are undefined when either a column or a row sums to zero
        #The probability of this happening are >0 in some bootstrap schemes.
        #We signal this to the user
        n_na <- length(which(is.na(resample_indexes[[i]])))
        actual_n <- n_b - n_na
        if (actual_n < n_b){
          warning(paste0("Due to NAs, the actual number of bootstrap sample for
                         the studentized bootstrap method is ", n_b- n_na))
        }
      }}
    #Calculate Bootstrap for the bayesian boot resample
    if ( ("BayModel" %in% method) | ("CBayBoot" %in% method)) {
      resample_bay <- boot_f(mat, n_b, bay = T, prior)
      resample_indexes_bay <- list()
      estimated_indexes_bay <- list()
      for (i in index){
        print(i)
        resample_indexes_bay[[i]] <-  unlist(lapply(  seq(n_b),function(x)
          .index_env$indexes[[i]](matrix(c(resample_bay$resample[x,]), nrow = 2,
                                         byrow = T))))
        estimated_indexes_bay[[i]] <- unlist(lapply(  seq(n_b),function(x)
                                      .index_env$indexes[[i]](matrix(
                                      c(resample_bay$ps[x,]), nrow = 2,
                                      byrow = T))))
        # Rarely, there may be NA in the resample due to the fact that the
        # indices are undefined when either a column or a row sums to zero
        # The probability of this happening are >0 in some bootstrap schemes.
        # We signal this to the user
        n_na <- length(which(is.na(resample_indexes_bay[[i]])))
        actual_n <- n_b - n_na
        if (actual_n < n_b){
          print(paste0("Due to NAs, the actual number of bootstrap sample for the corrected Bayesian method is ", n_b- n_na))
        }
      }
    }
    #Calculate the simple percentile Bootstrap CI
    if ("Boot" %in% method){
      CIs[["Boot"]] <- lapply(index, function(x) c(hdi(resample_indexes[[x]], credMass = confidence_level), centers[x]))
      names(CIs[["Boot"]]) <- index
    }
    #Calculate the bias-corrected percentile Bootstrap CI
    if ("CBoot" %in% method){
      centers_correct <- unlist(lapply(seq(length(centers)), function(x) 2*centers[x] - mean(resample_indexes[[x]], na.rm = TRUE) ))
      centers_correct <- unlist(lapply(centers_correct, function(x) max(0, min(1, x))))
      names(centers_correct) <- index
      CIs[["CBoot"]] <- lapply(index, function(x) c(hdi(2*centers[x] - resample_indexes[[x]], credMass = confidence_level), centers_correct[x]))
      names(CIs[["CBoot"]]) <- index
    }
    #Calculate the studentized Bootstrap CI
    if ("TBoot" %in% method){
      CIs[["TBoot"]] <- lapply(index, function(x) c(centers[x] - Ts[[x]][2]*sqrt(v_est[x]),   centers[x] - Ts[[x]][1]*sqrt(v_est[x]) , centers[x]))
      names(CIs[["TBoot"]]) <- index
    }
    #Calculate the bias-corrected studentized Bootstrap CI
    if ("TCBoot" %in% method){
      #This is left here if I want to pick it up in the future
      #For each re-sample, we calculate the correction
      # re_resample <- lapply(seq(nrow(resample)),function(x) t(rbind(rmultinom(sub_n_b, n_m, resample[x,c(1:3)]/sum(resample[x,c(1:3)])), rmultinom(sub_n_b, n_f, resample[x,c(4:6)]/sum(resample[x,c(4:6)])))) )
      # resample_correction <- c()
      # for (i in seq(length(indexes))){
      #   for (j in seq(nrow(resample))){
      #     resample_correction[length(resample_correction)+1] <- mean(unlist(lapply(seq(sub_n_b), function(x) indexes[[i]](matrix(re_resample[[j]][x,], nrow = 2, byrow = T))   )))
      # }
      # }
      # resample_correction <- matrix(resample_correction, nrow = n_b, byrow = F)
      centers_correct <- unlist(lapply(seq(length(centers)), function(x) 2*centers[x] - mean(resample_indexes[[x]], na.rm = FALSE) ))
      centers_correct <- unlist(lapply(centers_correct, function(x) max(0, min(1, x))))
      names(centers_correct) <- index
      CIs[["TCBoot"]] <- lapply(index, function(x) c(centers_correct[x] - Ts[[x]][2]*sqrt(v_est[x]),   centers_correct[x] - Ts[[x]][1]*sqrt(v_est[x]) , centers_correct[x]))
      names(CIs[["TCBoot"]]) <- index
    }
    #Calculate the Bayesian Model Index
    if (("BayModel" %in% method)){
      CIs[["BayModel"]] <- lapply(index, function(x) c(hdi(estimated_indexes_bay[[x]], credMass = confidence_level), mean(estimated_indexes_bay[[x]], na.rm = TRUE)))
      names(CIs[["BayModel"]]) <- index
    }
    #Calculate the Bayesian Bootstrap CI
    if (("CBayBoot" %in% method)){
      centers_correct_bay <- unlist(lapply(seq(length(centers)), function(x) 2*centers[x] - mean(resample_indexes_bay[[x]], na.rm = TRUE) ))
      centers_correct_bay <- unlist(lapply(centers_correct_bay, function(x) max(0, min(1, x))))
      names(centers_correct_bay) <- index
      CIs[["CBayBoot"]] <- lapply(index, function(x) c(hdi(2*centers[x] - resample_indexes_bay[[x]], credMass = confidence_level), centers_correct_bay[x]))
      names(CIs[["CBayBoot"]]) <- index
    }
  }
  #Nicer names inside the lists
  for (met in method) {
    for (ind in index) {
      CIs[[met]][[ind]] <- c(max(0, CIs[[met]][[ind]][1]), min(1, CIs[[met]][[ind]][2]), CIs[[met]][[ind]][3])
      #Due to extreme bias in the measurement, sometime, the CIs obtained are invalid.
      #This happens with TBoot and CTBoot and is due to the fact that the T-distribution obtained
      #May be entirely positive (again, due to bias). When this happens we signal that
      if (CIs[[met]][[ind]][1] >= CIs[[met]][[ind]][2]){
        CIs[[met]][[ind]] <- c(NA, NA, CIs[[met]][[ind]][3])
        print("The estimation produced an invalid CI: ")
        print(met)
        print(ind)
      }
      names(CIs[[met]][[ind]]) <- c("Lower", "Upper", "Estimate")
    }
  }
  return(CIs)
}

#CI for differences between indexes
#Notice, that the function supposes that the matrices have the same number of columns and in general
#The same columns. It doesn't make a lot of sense to compare tables with different categories
#The function also supposes that the tables have the same sampling scheme. This is not strictly necessary
#mat1 and mat2 are the data table to be compared. Method regards how to calculate the CI
#index is a list of indexes to calculate based on the data
index_difference_ci <- function(mat1, mat2, index = c("D", "gini", "MI", "Atk", "norm_MI"), method = c("Norm", "Boot", "CBoot", "TBoot", "TCBoot", "BayModel", "CBayBoot"), sampling = "fm", n_b = FALSE, confidence_level = 0.95){
  CIs <- list()
  alfa <- 1 - confidence_level
  k <- ncol(mat1)
  N1 <- sum(mat1)
  N2 <- sum(mat2)
  #We start by selecting the appropriate functions for bootstrap simulation, SE calculation
  #For SE
  sigma_f <- sigmas[[sampling]]
  ders_f <- ders[[sampling]]
  #For bootstrap simulation
  boot_f  <- boots[[sampling]]
  #Empirical indexes as calculated from the data  indexes[[x]](mat)
  centers1 <- unlist(lapply(index, function(x) indexes[[x]](mat1) ))
  names(centers1) <- index
  centers2 <- unlist(lapply(index, function(x) indexes[[x]](mat2) ))
  names(centers2) <- index
  centers_diff <- centers1 - centers2
  names(centers_diff) <- index
  #Here we calculate ws_var which express how the variance of shrinks with N getting bigger
  #This changes depending on the sampling scheme. Used to calculate the final variance
  if (sampling == "mb"){
    #ws expresses how many individual each column has
    ws1 <- colSums(mat1)
    ws2 <- colSums(mat2)
    ws_var1 <- unlist(lapply(seq(length(ws1)), function(x) rep(ws1[x], 2)))
    ws_var2 <- unlist(lapply(seq(length(ws2)), function(x) rep(ws2[x], 2)))
  } else {
    ws_var1 <- rep(N1, k)
    ws_var2 <- rep(N2, k)
  }
  #Empirical Vcov matrix from the data and derivative
  Sigma_s1 <- sigma_f(mat1)
  Sigma_s2 <- sigma_f(mat2)
  #We cannot simply divide by sqrt(N) at the end for mb. This solves it
  Sigma_s1 <- Sigma_s1/ws_var1
  Sigma_s2 <- Sigma_s2/ws_var2
  lambdas1 <- lapply(index, function(x) ders_f[[x]](mat1))
  lambdas2 <- lapply(index, function(x) ders_f[[x]](mat2))
  v_est1 <- unlist(lapply(seq(length(lambdas1)), function(x) as.numeric(lambdas1[[x]] %*%  Sigma_s1 %*% t(t(lambdas1[[x]])))))
  v_est2 <- unlist(lapply(seq(length(lambdas2)), function(x) as.numeric(lambdas2[[x]] %*%  Sigma_s2 %*% t(t(lambdas2[[x]])))))
  names(v_est1) <- index
  names(v_est2) <- index
  if ("Norm" %in% method){
    CIs[["Norm"]] <- lapply(index, function(x) c((centers1[x] - centers2[x]) + qnorm(alfa/2)*sqrt(v_est1[x] + v_est2[x]),  (centers1[x] - centers2[x]) - qnorm(alfa/2)*sqrt(v_est1[x] + v_est2[x]), centers_diff[x]))
    names(CIs[["Norm"]]) <- index
  }
  if ( ("CBoot" %in% method) | ("TBoot" %in% method)
       | ("TCBoot" %in% method) | ("Boot" %in% method)
       | ("BayModel" %in% method) | ("CBayBoot" %in% method) ){
    if (!(n_b)){
      print("You need to specify the n_b parameter to use bootstrap CIs")
    }
    #Now we construct the BS CIs!
    #Bootstrap of categorical variable = resampling under the same scheme
    if ( ("CBoot" %in% method) | ("Boot" %in% method) |
         ("TCBoot" %in% method) | ("TBoot" %in% method)  ) {
      resample1 <- boot_f(mat1, n_b)
      resample2 <- boot_f(mat2, n_b)
      resample_indexes1 <- list()
      resample_indexes2 <- list()
    }
    #Calculate variance for studentized Boot
    if ( ("TCBoot" %in% method) | ("TBoot" %in% method)  ){
      Sigma_est_res1 <- lapply(seq(nrow(resample1)), function(x) sigma_f(matrix(resample1[x,], nrow = 2, byrow = T )))
      Sigma_est_res2 <- lapply(seq(nrow(resample2)), function(x) sigma_f(matrix(resample2[x,], nrow = 2, byrow = T )))
      Sigma_est_res1 <- lapply(seq(nrow(resample1)), function(x) Sigma_est_res1[[x]] / ws_var1 )
      Sigma_est_res2 <- lapply(seq(nrow(resample2)), function(x) Sigma_est_res2[[x]] / ws_var2 )
      v_est_res1 <- list()
      v_est_res2 <- list()
      Ts <- list()
    }
    #Calculate indexes and Ts on the boot resample
    if ( ("CBoot" %in% method) | ("Boot" %in% method) |
         ("TCBoot" %in% method) | ("TBoot" %in% method)  ) {
      for (i in index){
        print(i)
        resample_indexes1[[i]] <-  unlist(lapply(  seq(nrow(resample1)),function(x) indexes[[i]](matrix(c(resample1[x,]), nrow = 2, byrow = T)) ))
        resample_indexes2[[i]] <-  unlist(lapply(  seq(nrow(resample2)),function(x) indexes[[i]](matrix(c(resample2[x,]), nrow = 2, byrow = T)) ))
        if ( ("TCBoot" %in% method) | ("TBoot" %in% method)  ){
          lambdas_est_res1 <- lapply(  seq(nrow(resample1)), function(x) ders_f[[i]](matrix(resample1[x,], nrow = 2, byrow = T)) )
          lambdas_est_res2 <- lapply(  seq(nrow(resample2)), function(x) ders_f[[i]](matrix(resample2[x,], nrow = 2, byrow = T)) )
          v_est_res1[[i]] <- unlist(lapply(seq(nrow(resample1)), function(x) as.numeric(lambdas_est_res1[[x]] %*%  Sigma_est_res1[[x]] %*% t(t(lambdas_est_res1[[x]])))))
          v_est_res2[[i]] <- unlist(lapply(seq(nrow(resample2)), function(x) as.numeric(lambdas_est_res2[[x]] %*%  Sigma_est_res2[[x]] %*% t(t(lambdas_est_res2[[x]])))))
          remove(lambdas_est_res1, lambdas_est_res2)
          Ts_temp <- ((resample_indexes1[[i]] - resample_indexes2[[i]]) - ( centers1[i] -centers2[i] ))/(sqrt(v_est_res1[[i]]+v_est_res2[[i]]))
          Ts[[i]] <- hdi(Ts_temp, confidence_level)
        }
        #Rarely, there may be NA in the resample due to the fact that the indices are undefined when either a column or a row sums to zero
        #The probability of this happening are >0 in some bootstrap schemes. We signal this to the user. All indices are equally affected by this
        n_na <- length(which(is.na(resample_indexes1[[i]])))
        actual_n <- n_b - n_na
        if (actual_n < n_b){
          print(paste0("Due to NAs, the actual number of bootstrap sample for the first table for the studentized bootstrap method is ", n_b- n_na))
        }
        n_na <- length(which(is.na(resample_indexes2[[i]])))
        actual_n <- n_b - n_na
        if (actual_n < n_b){
          print(paste0("Due to NAs, the actual number of bootstrap sample for the second table for the studentized bootstrap method is ", n_b- n_na))
        }
      }}
    #Calculate Bootstrap for the bayesian boot resample
    if ( ("BayModel" %in% method) | ("CBayBoot" %in% method)) {
      resample_bay1 <- boot_f(mat1, n_b, bay = T)
      resample_bay2 <- boot_f(mat2, n_b, bay = T)
      resample_indexes_bay1 <- list()
      resample_indexes_bay2 <- list()
      estimated_indexes_bay1 <- list()
      estimated_indexes_bay2 <- list()
      resample_diff_bay <- list()
      estimated_diff_bay <- list()
      for (i in index){
        print(i)
        resample_indexes_bay1[[i]] <-  unlist(lapply(  seq(n_b),function(x) indexes[[i]](matrix(c(resample_bay1$resample[x,]), nrow = 2, byrow = T)) ))
        resample_indexes_bay2[[i]] <-  unlist(lapply(  seq(n_b),function(x) indexes[[i]](matrix(c(resample_bay2$resample[x,]), nrow = 2, byrow = T)) ))
        estimated_indexes_bay1[[i]] <- unlist(lapply(  seq(n_b),function(x) indexes[[i]](matrix(c(resample_bay1$ps[x,]), nrow = 2, byrow = T)) ))
        estimated_indexes_bay2[[i]] <- unlist(lapply(  seq(n_b),function(x) indexes[[i]](matrix(c(resample_bay2$ps[x,]), nrow = 2, byrow = T)) ))
        resample_diff_bay[[i]] <- resample_indexes_bay1[[i]] - resample_indexes_bay2[[i]]
        estimated_diff_bay[[i]] <- estimated_indexes_bay1[[i]] - estimated_indexes_bay2[[i]]
      }
      #Rarely, there may be NA in the resample due to the fact that the indices are undefined when either a column or a row sums to zero
      #The probability of this happening are >0 in some bootstrap schemes. We signal this to the user. All indices are equally affected by this
      n_na <- length(which(is.na(resample_indexes_bay1[[i]])))
      actual_n <- n_b - n_na
      if (actual_n < n_b){
        print(paste0("Due to NAs, the actual number of bootstrap sample for the first table for the Bayesian corrected method is ", n_b- n_na))
      }
      n_na <- length(which(is.na(resample_indexes_bay2[[i]])))
      actual_n <- n_b - n_na
      if (actual_n < n_b){
        print(paste0("Due to NAs, the actual number of bootstrap sample for the second table for the Bayesian corrected method is ", n_b- n_na))
      }
    }
    #Calculate the simple percentile Bootstrap CI
    if ("Boot" %in% method){
      CIs[["Boot"]] <- lapply(index, function(x) c(hdi(resample_indexes1[[x]] - resample_indexes2[[x]], credMass = confidence_level), centers_diff[x]))
      names(CIs[["Boot"]]) <- index
    }
    #Calculate the bias-corrected percentile Bootstrap CI
    if ("CBoot" %in% method){
      centers_correct1 <- unlist(lapply(index, function(x) 2*centers1[x] - mean(resample_indexes1[[x]], na.rm = TRUE) ))
      names(centers_correct1) <- index
      centers_correct2 <- unlist(lapply(index, function(x) 2*centers2[x] - mean(resample_indexes2[[x]], na.rm = TRUE) ))
      names(centers_correct2) <- index
      centers_correct_diff <- centers_correct1 -centers_correct2
      centers_correct_diff <- unlist(lapply(centers_correct_diff, function(x) max(-1, min(1, x))))
      names(centers_correct_diff) <- index
      CIs[["CBoot"]] <- lapply(index, function(x) c(hdi(2*centers_diff[x] - (resample_indexes1[[x]] - resample_indexes2[[x]]), credMass = confidence_level), centers_correct_diff[x]))
      names(CIs[["CBoot"]]) <- index
    }
    #Calculate the studentized Bootstrap CI
    if ("TBoot" %in% method){
      CIs[["TBoot"]] <- lapply(index, function(x) c(centers_diff[x]  - Ts[[x]][2]*sqrt(v_est1[x] + v_est2[x]), centers_diff[x] - Ts[[x]][1]*sqrt(v_est1[x] + v_est2[x]), centers_diff[x]) )
      names(CIs[["TBoot"]]) <- index
    }
    #Calculate the bias-corrected studentized Bootstrap CI
    if ("TCBoot" %in% method){
      centers_correct1 <- unlist(lapply(index, function(x) 2*centers1[x] - mean(resample_indexes1[[x]], na.rm = TRUE) ))
      names(centers_correct1) <- index
      centers_correct2 <- unlist(lapply(index, function(x) 2*centers2[x] - mean(resample_indexes2[[x]], na.rm = TRUE) ))
      names(centers_correct2) <- index
      centers_correct_diff <- centers_correct1 -centers_correct2
      centers_correct_diff <- unlist(lapply(centers_correct_diff, function(x) max(-1, min(1, x))))
      names(centers_correct_diff) <- index
      CIs[["TCBoot"]] <- lapply(index, function(x) c(centers_diff[x]  - Ts[[x]][2]*sqrt(v_est1[x] + v_est2[x]),  centers_correct_diff[x] - Ts[[x]][1]*sqrt(v_est1[x] + v_est2[x]), centers_correct_diff[x]) )
      names(CIs[["TCBoot"]]) <- index
    }
    #Calculate the Bayesian Model Index
    if (("BayModel" %in% method)){
      CIs[["BayModel"]] <- lapply(index, function(x) c(hdi(estimated_diff_bay[[x]], credMass = confidence_level), mean(estimated_diff_bay[[x]])))
      names(CIs[["BayModel"]]) <- index
    }
    #Calculate the Bayesian Bootstrap CI
    if (("CBayBoot" %in% method)){
      centers_correct_bay1 <- unlist(lapply(index, function(x) 2*centers1[x] - mean(resample_indexes_bay1[[x]], na.rm = TRUE) ))
      names(centers_correct_bay1) <- index
      centers_correct_bay2 <- unlist(lapply(index, function(x) 2*centers2[x] - mean(resample_indexes_bay2[[x]], na.rm = TRUE) ))
      names(centers_correct_bay2) <- index
      centers_correct_diff <- centers_correct_bay1 -centers_correct_bay2
      centers_correct_diff <- unlist(lapply(centers_correct_diff, function(x) max(-1, min(1, x))))
      names(centers_correct_diff) <- index
      CIs[["CBayBoot"]] <- lapply(index, function(x) c(hdi(2*centers_diff[x] - resample_diff_bay[[x]], credMass = confidence_level), centers_correct_diff[x]))
      names(CIs[["CBayBoot"]]) <- index
    }
  }
  #Nicer names inside the lists
  for (met in method){
    for (ind in index){
      CIs[[met]][[ind]] <- c(max(-1, CIs[[met]][[ind]][1]), min(1, CIs[[met]][[ind]][2]), CIs[[met]][[ind]][3])
      #Due to extreme bias in the measurement, sometime, the CIs obtained are invalid.
      #This happens TBoot and CTBoot and is due to the fact that the T-distribution obtained
      #May be entirely positive (again, due to bias). When this happens we signal that
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
