#For HT we need to have all the derivatives of the various
#indexes under the various sampling schemes. This amounts to number-indices
# times 3 (frameworks) different derivatives
##########DERIVATIVE FOR INDEPENDENT-GROUPS FRAMEWORK##########

d_d_ind_dm <- function(mat){
  mat <- mat/rowSums(mat)
  der <- c(sign(mat[1,] - mat[2,]), -sign(mat[1,] - mat[2,]))
  return(1/2*der)
}

d_gini_ind_dm <- function(mat){
  # p <- rowSums(mat)[1]/sum(mat)
  # norm <- 1/2*1/(p*(1-p))
  mat <- mat/rowSums(mat)
  abs_matrix <- t(t(rep(1, ncol(mat)))) %*% t(mat[1,]) * mat[2,] %*% t(rep(1, ncol(mat)))
  abs_matrix <- t(abs_matrix) - abs_matrix
  sign_matrix <- sign(abs_matrix)
  # abs_matrix <- abs(abs_matrix)
  # A <- sum(abs_matrix)
  val_mat_1 <- matrix(rep(mat[2,],ncol(mat)), byrow = T,  nrow = ncol(mat))
  val_mat_1 <- val_mat_1 * sign_matrix
  val_mat_2 <- matrix(rep(mat[1,],ncol(mat)), byrow = T,  nrow = ncol(mat))
  #Notice the minus
  val_mat_2 <- -val_mat_2 * sign_matrix
  val_mats <- list(val_mat_1, val_mat_2)
  ders <- rep(0, ncol(mat)*2)
  for (i in seq(2*ncol(mat))){
    sel = 1
    der_n <- i
    if (i > ncol(mat)){
      i <- i - ncol(mat)
      sel <- 2}
    ders[der_n] <- sum(val_mats[[sel]][i,])
  }
  return(ders)
}

#Derivative of Mutual Info for double multinomial
d_mutual_info_ind_dm <- function(mat){
  k <- ncol(mat)
  p <- rowSums(mat)/sum(mat)
  mat <- mat/rowSums(mat)
  mat_ext <- unlist(lapply(seq(nrow(mat)), function(x) mat[x, ]))
  #To avoid problems we add a little probability to 0 cells
  MN <- min(mat_ext[mat_ext!=0])
  mat_ext <- ifelse(mat_ext !=0, mat_ext, MN/length(mat_ext) )
  der <- unlist(lapply(seq(k*nrow(mat)),function(x) p[ceiling(x/k)] * log(mat_ext[x]/(p %*% mat[, x-k*(ceiling(x/k)-1)])[1]) ))
  return(der)
}

#Derivative of Mutual Info for full multinomial
d_theil_ind_dm <- function(mat){
  p <- rowSums(mat)
  ent <- entropy(p)
  der <- 1/ent*d_MI_dm(mat)
  return(der)
}
#Works


#Supposes there are 2 rows and not more!
d_atkinson_ind_dm<- function(mat){
  mat <- mat/rowSums(mat)
  k <- ncol(mat)
  mat_ext <- unlist(lapply(seq(nrow(mat)), function(x) mat[x, ]))
  #To avoid problems we add a little probability to 0 cells
  MN <- min(mat_ext[mat_ext!=0])
  mat_ext <- ifelse(mat_ext !=0, mat_ext, MN/length(mat_ext) )
  der <- -0.5*unlist(lapply(seq(k*nrow(mat)), function(x) sqrt(prod(mat[, x-k*(ceiling(x/k)-1)])/(mat_ext[x]^2)) ))
  return(der)
}

# Derivative of xPx for double Multinomial
d_isolation_ind_dm <- function(mat){
  p <- rowSums(mat)
  c_sm <- colSums(mat)
  mat <- mat/rowSums(mat)
  # derivative for the first group
  der <- p[1]*mat[1,]/c_sm^2 * (mat[1,]*p[1] + 2*mat[2,]*p[2])
  # derivative for the second group
  der <- c(der, -p[1]*p[2]*(mat[1,]/c_sm)^2  )
  return(der)
}

# Derivative of xPx_inv for double Multinomial
d_isolationInv_ind_dm <- function(mat){
  mat <-matrix(c(mat[2,], mat[1,]), nrow = 2, byrow = TRUE)
  return(d_xpx_dm(mat))
}



# Derivative of V for double Multinomial
d_v_ind_dm <- function(mat){
  p <- rowSums(mat)[1]/sum(mat)
  # It is very convenient to use d_xpx_dm
  der <-  d_xpx_dm(mat)/(1-p)
  return(der)
}


##########DERIVATIVE FOR FULL MULTINOMIAL##########

d_d_ind_fm <- function(mat){
  k <- ncol(mat)
  sums <- rowSums(mat)/(sum(mat))
  mat <- mat/rowSums(mat)
  signs <- matrix(c(sign(mat[1,] - mat[2,]), -sign(mat[1,] - mat[2,])), nrow = 2, byrow = TRUE)
  addendo <- signs %*% t(mat)
  der <- 1/2 * unlist(lapply(seq(k*nrow(mat)), function(x) 1/sums[ceiling(x/k)] *
                               (signs[ceiling(x/k), x-k*((ceiling(x/k))-1)]  -
                                  addendo[ceiling(x/k), ceiling(x/k)]) ))
  return(der)
}

d_gini_ind_fm <- function(mat){
  p <- rowSums(mat)[1]/sum(mat)
  norm <- 1/2*1/(p*(1-p))
  mat <- mat/sum(mat)
  abs_matrix <- t(t(rep(1, ncol(mat)))) %*% t(mat[1,]) * mat[2,] %*% t(rep(1, ncol(mat)))
  abs_matrix <- t(abs_matrix) - abs_matrix
  sign_matrix <- sign(abs_matrix)
  abs_matrix <- abs(abs_matrix)
  A <- sum(abs_matrix)
  val_mat_1 <- matrix(rep(mat[2,],ncol(mat)), byrow = T,  nrow = ncol(mat))
  val_mat_1 <- val_mat_1 * sign_matrix
  val_mat_2 <- matrix(rep(mat[1,],ncol(mat)), byrow = T,  nrow = ncol(mat))
  #Notice the minus
  val_mat_2 <- -val_mat_2 * sign_matrix
  val_mats <- list(val_mat_1, val_mat_2)
  ders <- rep(0, ncol(mat)*2)
  for (i in seq(2*ncol(mat))){
    sel = 1
    der_n <- i
    if (i > ncol(mat)){
      i <- i - ncol(mat)
      sel <- 2}
    ders[der_n] <- sum(val_mats[[sel]][i,])
  }
  addend <- c(rep(-A/p, ncol(mat)), rep(-A/(1-p), ncol(mat)))
  norm * (addend+ 2*ders)
}

#Derivative of Mutual Info for full multinomial
d_mutual_info_ind_fm <- function(mat){
  k <- ncol(mat)
  mat <- mat/sum(mat)
  mat_ext <- unlist(lapply(seq(nrow(mat)), function(x) mat[x, ]))
  #To avoid problems we add a little probability to 0 cells
  MN <- min(mat_ext[mat_ext!=0])
  mat_ext <- ifelse(mat_ext !=0, mat_ext, MN/length(mat_ext) )
  # mat_ext <- unlist(lapply(seq(k*nrow(mat)), function(x) mat[ceiling(x/k), x-k*(ceiling(x/k)-1)]))
  p <- rowSums(mat)
  m1 <- colSums(mat)
  # Assumes 2 rows
  d_ent_m1 <- rep(-(log(m1) + 1), nrow(mat))
  der <- d_ent_m1 - unlist(lapply(seq(k*nrow(mat)),function(x) log(p[ceiling(x/k)]) -log(mat_ext[x])  ))
  return(der)
}

#Derivative of Normalized Mutual Info for full multinomial
d_theil_ind_fm <- function(mat){
  p <- rowSums(mat)
  ent <- entropy(p)
  mi <- mutual_info(mat)
  der <- d_MI_fm(mat)
  p <- unlist(lapply(seq(nrow(mat)), function(x) rep(p[x], ncol(mat))))
  der <- 1/ent * ( der + (1 + log(p))/ent * mi )
}
#Works

#Supposes there are 2 rows and not more!
d_atkinson_ind_fm<- function(mat){
  k <- ncol(mat)
  p <- rowSums(mat)/sum(mat)
  mat_ext <- unlist(lapply(seq(nrow(mat)), function(x) mat[x, ]))
  #To avoid problems we add a little probability to 0 cells
  MN <- min(mat_ext[mat_ext!=0])
  mat_ext <- ifelse(mat_ext !=0, mat_ext, MN/length(mat_ext) )
  norm <- 1/sqrt(prod(p))
  c_prod <- sum(unlist(lapply(seq(k), function(x) sqrt(prod(mat[,x]/sum(mat))))))
  # c_prod <- c_prod*norm
  # der <- unlist(lapply(seq(k*nrow(mat)),function(x) 1/p[ceiling(x/k)] ))
  der <- -0.5*norm*unlist(lapply(seq(k*nrow(mat)), function(x) sqrt(prod(mat[, x-k*(ceiling(x/k)-1)])/(mat_ext[x]^2)) - c_prod/p[ceiling(x/k)] ))
  return(der)
}

d_isolation_ind_fm <- function(mat){
  mat <- mat/sum(mat)
  p <- rowSums(mat)[1]
  # derivative for the first group
  addendo <- -1/p^2*sum(mat[1, ]^2/colSums(mat))
  der <- 1/p*(1/colSums(mat))*mat[1,]*(
    2 - mat[1,]/ colSums(mat)
  )
  der <- addendo + der
  # Add the derivative for the second group (very different)
  der <- c(der, -1/p*(mat[1,]/colSums(mat))^2)
  return(der)
}

d_isolationInv_ind_fm <- function(mat){
  mat <- matrix(c(mat[2,], mat[1,]), nrow = 2, byrow = T)
  return(d_xpx_fm(mat))
}

# Supposes two rows
d_v_ind_fm <- function(mat){
  mat <- mat/sum(mat)
  p <- rowSums(mat)[1]
  # It is useful to use xpx
  der <- d_xpx_fm(mat)
  # Second group already done
  der <- 1/(1-p)*der
  # First group case
  der[1:ncol(mat)] <- V(mat)/(1-p) + (der[1:ncol(mat)] - 1/(1-p))
  return(der)
}


##########DERIVATIVE FOR INDEPENDENT-UNITS##########

#Each function has a 'in_form' parameter.
#This parameter check if the matrix is already in the same form as
#real_matrix is.

#Derivative for the D index mult-binomial version
d_D_ind_mb <- function(mat, in_form = F){
  if (!(in_form)){
    mat[2,] <- colSums(mat)
    mat[1,] <- mat[1,]/mat[2,]
    mat[2,] <- mat[2,]/sum(mat[2,])
  }
  p <- (mat[1,] %*% mat[2,])[1]
  norm <- 1/2 * 1/(p*(1-p))
  signs <- sign(mat[1,] - p)
  #There is a fixed part that does not need to be calculated for each p
  addendo <- -((1-2*p)/(p*(1-p)))*(mat[2,] %*% (abs(mat[1,] - p)))[1]
  addendo <- addendo - (mat[2,] %*%  signs)[1]
  der <- mat[2,]
  der <- norm * der * unlist(lapply( seq(ncol(mat)), function(x) signs[x] + addendo ))
  der <- unlist(lapply(seq(length(der)), function(x) c(der[x],0) ))
  return(der)
}

d_gini_ind_mb <- function(mat, in_form = F){
  if (!(in_form)){
    mat[2,] <- colSums(mat)
    mat[1,] <- mat[1,]/mat[2,]
    mat[2,] <- mat[2,]/sum(mat[2,])
  }
  p <- c(mat[1,]%*% mat[2,])
  norm <- (1/(p*(1-p)))
  abs_matrix <- t(t(rep(1, ncol(mat)))) %*% t(mat[1,])
  abs_matrix <- t(abs_matrix) - abs_matrix
  sign_matrix <- sign(abs_matrix)
  w_matrix <- matrix(rep(mat[2,],ncol(mat)), byrow = T,  nrow = ncol(mat))
  w_matrix <- sign_matrix*w_matrix
  ders <- rowSums(w_matrix)
  ders <- mat[2,]*norm*(ders - norm*1/2*(1-2*p)*(t(mat[2,]) %*% abs(abs_matrix) %*% t(t(mat[2,])))[1])
  #We put the derivatives (= 0) for the probabilities of the second group
  ders <- unlist(lapply(seq(length(ders)), function(x) c(ders[x],0) ))
  return(ders)
}

#Derivative of Mutual info for multiple binomial
#Supposes two rows
d_mutual_info_ind_mb <- function(mat, in_form = FALSE) {
  if (!(in_form)){
    mat[2,] <- colSums(mat)
    mat[1,] <- mat[1,]/mat[2,]
    mat[2,] <- mat[2,]/sum(mat[2,])
  }
  p <- (mat[1,]%*% mat[2,])[1]
  mat_ext <- mat[1,]
  #To avoid problems we add a little probability to 0 cells
  MN <- min(mat_ext[mat_ext!=0])
  mat_ext <- ifelse(mat_ext !=0, mat_ext, MN/length(mat_ext) )
  mat_ext <- ifelse(mat_ext !=1, mat_ext, 1- MN/length(mat_ext) )
  mat[1,] <- mat_ext
  der <- mat[2,]*(log(mat[1,]/(1-mat[1,])) -  log(p/(1-p)))
  der <- unlist(lapply(seq(length(der)), function(x) c(der[x],0) ))
  return(der)
}

#Derivative of Normalized Mutual info for multiple binomial
#Supposes two rows
d_theil_ind_mb <- function(mat, in_form = FALSE) {
  if (!(in_form)){
    mi <- mutual_info(mat)
    mat[2,] <- colSums(mat)
    mat[1,] <- mat[1,]/mat[2,]
    mat[2,] <- mat[2,]/sum(mat[2,])
  }
  p <- (mat[1,]%*% mat[2,])[1]
  p <- c(p, 1-p)
  ent <- entropy(p)
  #We need to calculate mi if the matrix was in_form from the beginning
  if (in_form){
    mi <- mutual_info(matrix(c(mat[1,]*mat[2,], (1-mat[1,])*mat[2,]), nrow =2, byrow = T))
  }
  der <- d_MI_mb(mat, in_form = TRUE)
  mat_2 <- unlist(lapply(seq(ncol(mat)), function(x) c(mat[2,x],0) ))
  der <- 1/ent * (der + mat_2/ent* log(p[1]/p[2])*mi)
  return(der)
}


#Derivative of the Atkinson index under a
#mult. binomial sampling scheme.
#Supposes there are 2 rows and not more!
d_atkinson_ind_mb<- function(mat, in_form = FALSE){
  if (!(in_form)){
    mat[2,] <- colSums(mat)
    mat[1,] <- mat[1,]/mat[2,]
    mat[2,] <- mat[2,]/sum(mat[2,])
  }
  k <- ncol(mat)
  p <- c(mat[1,] %*% mat[2,])
  norm <- 1/sqrt(p*(1-p))
  mat_ext <- mat[1,]
  #To avoid problems we add a little probability to 0 cells
  MN <- min(mat_ext[mat_ext!=0])
  mat_ext <- ifelse(mat_ext !=0, mat_ext, MN/length(mat_ext) )
  mat_ext <- ifelse(mat_ext !=1, mat_ext, 1- MN/length(mat_ext) )
  mat[1,] <- mat_ext
  #This is fixed for every term
  addendo <-norm^2*(1-2*p)*c(sqrt(mat[1,]*(1- mat[1,])) %*% mat[2,])
  der <- 1/2*norm*mat[2,]
  der <- der * (addendo -(1-2*mat[1,])/sqrt(mat[1,]*(1-mat[1,])))
  #We add 0s for the derivative of the 1-p terms
  der <- unlist(lapply(seq(length(der)), function(x) c(der[x],0) ))
  return(der)
}

#Derivative of xPx for multiple binomial
# Supposes two rows
d_isolation_ind_mb <- function(mat, in_form = FALSE){
  if (!(in_form)){
    mi <- mutual_info(mat)
    mat[2,] <- colSums(mat)
    mat[1,] <- mat[1,]/mat[2,]
    mat[2,] <- mat[2,]/sum(mat[2,])
  }
  p <- (mat[1,]%*% mat[2,])[1]
  der <- mat[2,]/p*(2*mat[1,] - xpx_bin(mat, in_form = TRUE))
  #We add 0s for the derivative of the 1-p terms
  der <- unlist(lapply(seq(length(der)), function(x) c(der[x],0) ))
  return(der)
}

#Derivative of xPx_inv for multiple binomial
# Supposes two rows
d_isolationInv_ind_mb <- function(mat, in_form = FALSE){
  if (!(in_form)){
    mi <- mutual_info(mat)
    mat[2,] <- colSums(mat)
    mat[1,] <- mat[1,]/mat[2,]
    mat[2,] <- mat[2,]/sum(mat[2,])
  }
  mat[1,] <- 1 - mat[1,]
  return(d_xpx_mb(mat, in_form=TRUE))
}



#Derivative of V for multiple binomial
# Supposes two rows
d_v_ind_mb <- function(mat, in_form = FALSE){
  if (!(in_form)){
    mi <- mutual_info(mat)
    mat[2,] <- colSums(mat)
    mat[1,] <- mat[1,]/mat[2,]
    mat[2,] <- mat[2,]/sum(mat[2,])
  }
  p <- (mat[1,]%*% mat[2,])[1]
  # convenient for the actual calculation
  p <- 1 - p
  addendo <- mat[2,]/p * V_bin(mat, in_form = TRUE)
  # Useful to use d_xpx_mb, but only non-zero values
  # Here we are making the assumption of only two rows`
  non_zero_d_xpx <- d_xpx_mb(mat, in_form = TRUE)[seq(1,ncol(mat)*2, 2)]
  der <- addendo + 1/(p) * (non_zero_d_xpx - mat[2,])
  #We add 0s for the derivative of the 1-p terms
  der <- unlist(lapply(seq(length(der)), function(x) c(der[x],0) ))
  return(der)
}
