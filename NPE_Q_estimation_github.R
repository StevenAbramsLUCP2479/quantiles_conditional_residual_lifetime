######################################################################
#  Nonparametric estimation of quantiles of the conditional residual # 
#  lifetime distribution                                             #
#                                                                    #
#  Version: GitHub                                                   #
#  Author: Steven Abrams                                             #
#  Last update: July, 2025                                           #
######################################################################

####################
# Default settings #
####################

options(warn = -1)

#####################################################
# 1. Empirical survival functions for the marginals #
#####################################################

S1n = Vectorize(function(t, t1, n){(1/(n+1))*sum(as.numeric(t1 > t))}, "t")
S2n = Vectorize(function(t, t2, n){(1/(n+1))*sum(as.numeric(t2 > t))}, "t")

####################################################################
# 2. Empirical cumulative distribution functions for the marginals #
####################################################################

F1n = Vectorize(function(t, t1, n){(1/(n+1))*sum(as.numeric(t1 <= t))}, "t")
F2n = Vectorize(function(t, t2, n){(1/(n+1))*sum(as.numeric(t2 <= t))}, "t")

#######################################################
# 3. Akritas-Van Keilegom estimator F(t1 | T2 <= t21) #
#######################################################

library(npcure)
Fn <- function(t, t1, t2, delta, bws_local = TRUE){
  Bfit = beran(x = t2, t = t1, d = delta, x0 = t, local = bws_local)
  
  ## Rows: t2 dimension; Columns: t1 dimension; (see algorithm)
  Fn_mat = matrix(nrow = length(t), ncol = length(t1));
  Sn_mat = matrix(nrow = length(t), ncol = length(t1));
  W_mat = matrix(nrow = length(t), ncol = length(t1));
  
  for (j in 1:length(t)){
    Sn_mat[j, ] = Bfit$S[[j]];
    Fn_mat[j, ] = 1 - Sn_mat[j, ];
    W_mat[j, ] = diff(c(0,Fn_mat[j, ]));
  }
  
  h_mat = Bfit$h

  return(list(t1_grid = Bfit$testim, t2_grid = sort(t2),
              Fn_mat = Fn_mat, Sn_mat = Sn_mat, 
              W_mat = W_mat, h_mat = h_mat, Bfit = Bfit))
}

F_t2n <- function(s1, s2, t1, t2, delta){
  ## Beran estimator with jump sizes W_j for j = 1, ..., n at 
  ## uncensored event times T_1j (otherwise W_j = 0)
  Beran_estim = Fn(t = t2, t1 = t1, t2 = t2, delta = delta)
  
  ## W_ij = W_j(T2i) (dimension: t2 x t1)
  W_mat = Beran_estim$W_mat
  
  ## Ordered T1 and T2 values
  ordered_t2 = sort(t2); ordered_t1 = sort(t1);
  
  F_t2n_tilde = matrix(nrow = length(s1), ncol = 1)
  F_t2n_dtilde = matrix(nrow = length(s1), ncol = 1) 
  W_tilde = matrix(nrow = length(s2), ncol = length(t1)) 
  W_dtilde = matrix(nrow = length(s2), ncol = length(t1)) 

  for (k in 1:length(s1)){
    ind_tilde = as.numeric(ordered_t2 <= s2[k])
    ind_dtilde = as.numeric(ordered_t2 > s2[k])
    ind_t1 = as.numeric(ordered_t1 <= s1[k])
    
    ## W_tilde and W_dtilde values (see algorithm)
    W_tilde[k,] = t(W_mat) %*% ind_tilde
    W_dtilde[k,] = t(W_mat) %*% ind_dtilde
    W_tilde[k,] = W_tilde[k,]/sum(ind_tilde)
    W_dtilde[k,] = W_dtilde[k,]/sum(ind_dtilde)
    
    F_t2n_tilde[k] = W_tilde[k,] %*% ind_t1  
    F_t2n_dtilde[k] = W_dtilde[k,] %*% ind_t1
  }
  
  return(list(F_t2n_tilde = F_t2n_tilde, F_t2n_dtilde = F_t2n_dtilde,
              W_tilde = W_tilde, W_dtilde = W_dtilde))
} 

###################################################################
# X. Estimation of the cumulative distribution function F(t1, t2) #
###################################################################

kmG_fun <- function(x, t, delta, left_lim = F){
  kmG <- survfit(Surv(t, delta) ~ 1, type = "kaplan-meier") 
  sm <- summary(kmG)
  indices <- findInterval(x, sm$time, left.open = left_lim)
  value <- c(1,sm$surv)[indices + 1]
  return(value)
}

#####################################
# X.1 Stute estimators of F(t1, t2) #
#####################################

stute_F_estimator_oc <- function(t1, t2, z1, z2, delta1){
  n <- length(z1)
  t <- z1
  delta <- 1 - delta1
  G <- kmG_fun(t, t, delta, left_lim = T)
  
  value <- (1/n)*sum(((delta1)/G)*as.numeric(z1 <= t1 & z2 <= t2));
  return(value)
}

F_stute_oc <- Vectorize(stute_F_estimator_oc, vectorize.args = c("t1"))

stute_F_estimator <- function(t1, t2, z1, z2, delta1, delta2){
  n <- length(z1)
  t <- pmax(z1, z2)
  delta <- 1 - delta1*delta2
  G <- kmG_fun(t, t, delta, left_lim = T)
  
  value <- (1/n)*sum(((delta1*delta2)/G)*as.numeric(z1 <= t1 & z2 <= t2));
  return(value)
}

F_stute <- Vectorize(stute_F_estimator, vectorize.args = c("t1"))

########################################
# X.2 Lin-Ying estimators of F(t1, t2) #
########################################

linying_F_estimator_oc <- function(t1, t2, z1, z2, delta1){
  n <- length(z1)
  t <- z1
  delta <- 1 - delta1
  
  G_t1 <- pmax(pmin(kmG_fun(t1, t, delta),1),1e-12)
  
  value <- ((1/n)*sum(as.numeric(z2 <= t2)) - 
             ((1/n)*(1/G_t1)*sum(as.numeric(z1 > t1 & z2 <= t2))))
  return(value)
}

F_linying_oc <- Vectorize(linying_F_estimator_oc, vectorize.args = c("t1"))

linying_F_estimator <- function(t1, t2, z1, z2, delta1, delta2){
  n <- length(z1)
  t <- pmax(z1, z2)
  delta <- 1 - delta1*delta2
  
  G_t1t2 <- pmax(pmin(kmG_fun(pmax(t1, t2), t, delta),1),1e-12)
  G_t2 <- pmax(pmin(kmG_fun(t2, t, delta),1),1e-12)
  G_t1 <- pmax(pmin(kmG_fun(t1, t, delta),1),1e-12)
  
  value <- ((1/n)*(1/G_t1t2)*sum(as.numeric(z1 > t1 & z2 > t2))) - 
           ((1/n)*(1/G_t1)*sum(as.numeric(z1 > t1))) - 
           ((1/n)*(1/G_t2)*sum(as.numeric(z2 > t2))) + 1
  return(value)
}

F_linying <- Vectorize(linying_F_estimator, vectorize.args = c("t1"))

########################################################################
# 4. Estimation of the quantile function under one-component censoring #
########################################################################

##############################################################
# 4.1 Stute estimator of Fv(t1 | t21, t22) and F2v(t21, t22) #
##############################################################

stute_Fv_estimator_oc <- function(t1, t21, t22, z1, t2, delta1){
  n <- length(z1)
  t <- z1
  delta <- 1 - delta1
  G <- kmG_fun(t, t, delta, left_lim = T)
  
  value <- (1/n)*sum(((delta1)/G)*as.numeric(z1 <= t1 & t21 < t2 & t2 <= t22));
  return(value)
}

stute_F2v_estimator_oc <- function(t21, t22, z1, t2, delta1){
  n <- length(z1)
  t <- z1
  delta <- 1 - delta1
  G <- kmG_fun(t, t, delta, left_lim = T)
  
  value <- (1/n)*sum(((delta1)/G)*as.numeric(t21 < t2 & t2 <= t22));
  return(value)
} 

##################################################################
# 4.2 Lin and Ying estimator Fv(t1 | t21, t22) and F2v(t21, t22) #
##################################################################

linying_Fv_estimator_oc <- function(t1, t21, t22, z1, t2, delta1){
  n <- length(z1)
  t <- z1
  delta <- 1 - delta1
  
  G_t1 <- pmax(pmin(kmG_fun(t1, t, delta),1),1e-12)
  
  diff <- ((1/n)*sum(as.numeric(t2 > t21 & t2 <= t22)) - 
          ((1/n)*(1/G_t1)*sum(as.numeric(z1 > t1 & t21 < t2 & t2 <= t22))))
  value <- pmax(pmin(diff,1),0)
  return(value)
}

linying_F2v_estimator_oc <- function(t21, t22, z1, t2, delta1){
  n <- length(z1)  
  diff <- ((1/n)*sum(as.numeric(t21 < t2 & t2 <= t22)))
  value <- pmax(pmin(diff,1),0)
  return(value)
}

linying_Fv_estimator_oc2 <- function(t1, t21, t22, z1, t2, delta1){
  diff <- linying_F_estimator_oc(t1, t22, z1, z2 = t2, delta1) - 
          linying_F_estimator_oc(t1, t21, z1, z2 = t2, delta1);
  value <- pmax(pmin(diff,1),0)
  return(value)
}

##################################################
# 4.3 Empirical estimator for F2vtilde(t21, t22) #
##################################################

library(survival)
emp_F2v_estimator <- function(t21, t22, t2){
  n <- length(t2)
  
  F_t21 = F2n(t21, t2, n)
  F_t22 = F2n(t22, t2, n)
  
  value <- pmax(pmin(F_t22,1),0) - pmax(pmin(F_t21,1),0)
  return(value)
} 

#####################################################################
# 5. Estimation of the quantile function under univariate censoring #
#####################################################################

############################################################
# 5.1 IPW estimator of Fv(t1 | t21, t22) and F2v(t21, t22) #
############################################################

stute_Fv_estimator <- function(t1, t21, t22, z1, z2, delta1, delta2){
  n <- length(z1)
  t <- pmax(z1, z2)
  delta <- 1 - delta1*delta2
  G <- kmG_fun(t, t, delta, left_lim = T)
  
  value <- (1/n)*sum(((delta1*delta2)/G)*as.numeric(z1 <= t1 & t21 < z2 & z2 <= t22));
  return(value)
}

stute_F2v_estimator <- function(t21, t22, z1, z2, delta1, delta2){
  n <- length(z1)
  t <- pmax(z1, z2)
  delta <- 1 - delta1*delta2
  G <- kmG_fun(t, t, delta, left_lim = T)
  
  value <- (1/n)*sum(((delta1*delta2)/G)*as.numeric(t21 < z2 & z2 <= t22));
  return(value)
} 

original_stute_Fv_estimator <- function(t1, t21, t22, z1, z2, delta1, delta2, method = "km"){
  n <- length(z2)
  zz2 <- sort(z2)
  ddelta2 <- delta2[order(z2)]
  
  if (method == "direct"){
    prod_vec <- matrix(1, nrow = 1, ncol = n)
    for (i in 2:n){
      for (j in 1:(i-1)){
        prod_vec[i] = prod_vec[i]*(((n-j)/(n-j+1))^(ddelta2[j]))
      }
    }
    weights <- (ddelta2/(n - seq(1,n,1) + 1))*prod_vec
  }
  
  if (method == "km"){
    prod_vec <- matrix(1, nrow = 1, ncol = n)
    for (i in 2:n){
      for (j in 1:(i-1)){
        prod_vec[i] = prod_vec[i]*(((n-j)/(n-j+1))^(ddelta2[j]))
      }
    }
    
    km2 <- survfit(Surv(z2, delta2) ~ 1, type = "kaplan-meier") 
    km2_est <- summary(km2)$surv
    km_weights <- matrix(0, nrow = 1, ncol = n)
    km_weights[ddelta2 == 1] <- (1-(km2_est/c(1,km2_est[-length(km2_est)])))
    
    weights <- km_weights*prod_vec
  }
  
  value <- sum(weights*as.numeric(z1[order(z2)] <= t1 & t21 < sort(z2) & sort(z2) <= t22));
  return(value)
}

stute_Fv_estimator_vec <- function(t1, t21, t22, z1, z2, delta1, delta2){
  n <- length(z1)
  t <- pmax(z1, z2)
  delta <- 1 - delta1*delta2
  G <- kmG_fun(t, t, delta, left_lim = T)
  
  value <- matrix(nrow = 1, ncol = length(t1))
  for (j in 1:length(t1)){
    value[j] <- (1/n)*sum(((delta1*delta2)/G)*as.numeric(z1 <= t1[j] & t21 < z2 & z2 <= t22));
  }
  return(value)
}

##################################################################
# 5.2 Lin and Ying estimator Fv(t1 | t21, t22) and F2v(t21, t22) #
##################################################################

linying_Fv_estimator <- function(t1, t21, t22, z1, z2, delta1, delta2){
  n <- length(z1)
  t <- pmax(z1, z2)
  delta <- 1 - delta1*delta2
  
  G_t1t22 <- pmax(pmin(kmG_fun(pmax(t1, t22), t, delta),1),1e-12)
  G_t1t21 <- pmax(pmin(kmG_fun(pmax(t1, t21), t, delta),1),1e-12) 
  G_t22 <- pmax(pmin(kmG_fun(t22, t, delta),1),1e-12)
  G_t21 <- pmax(pmin(kmG_fun(t21, t, delta),1),1e-12)
  G_t1 <- pmax(pmin(kmG_fun(t1, t, delta),1),1e-12)
  
  diff <- ((1/n)*(1/G_t1t22)*sum(as.numeric(z1 > t1 & z2 > t22)) - 
           (1/n)*(1/G_t1t21)*sum(as.numeric(z1 > t1 & z2 > t21))) - 
          ((1/n)*(1/G_t22)*sum(as.numeric(z2 > t22)) - 
           (1/n)*(1/G_t21)*sum(as.numeric(z2 > t21)))
  
  value <- pmax(pmin(diff,1),0)
  return(value)
}

linying_F2v_estimator <- function(t21, t22, z1, z2, delta1, delta2){
  n <- length(z1)
  t <- pmax(z1, z2)
  delta <- 1 - delta1*delta2
  G_t1t22 <- pmax(pmin(kmG_fun(t22, t, delta),1),1e-12)
  G_t1t21 <- pmax(pmin(kmG_fun(t21, t, delta),1),1e-12)
  
  diff <- ((1/n)*(1/G_t1t21)*sum(as.numeric(z2 > t21))) - 
          ((1/n)*(1/G_t1t22)*sum(as.numeric(z2 > t22)))
  value <- pmax(pmin(diff,1),0)
  return(value)
} 

###########################################################
# 5.3 Kaplan-Meier based estimator for F2vtilde(t21, t22) #
###########################################################

library(survival)
km_F2v_estimator <- function(t21, t22, z1, z2, delta1, delta2){

  km <- survfit(Surv(z2, delta2) ~ 1, type = "kaplan-meier") 
  sm <- summary(km)
 
  indices_t22 <- findInterval(t22, sm$time, left.open = T)
  indices_t21 <- findInterval(t21, sm$time, left.open = T)
  F_t22 <- 1 - c(1,sm$surv)[indices_t22 + 1]
  F_t21 <- 1 - c(1,sm$surv)[indices_t21 + 1]
  
  value <- pmax(pmin(F_t22,1),0) - pmax(pmin(F_t21,1),0)
  return(value)
} 

#########################################################################
# 6. Estimators for Ftilde(t1 | t21, t22) under one-component censoring #
#########################################################################

###########################################################################
# 6.1 Estimators for Ftilde(t1 | t21, t22) based on Stute or LY estimator #
###########################################################################

Fn_tilde_stute_oc <- function(t1, t21, t22, z1, t2, delta1){
  num = stute_Fv_estimator_oc(t1, t21, t22, z1, t2, delta1)
  denom = pmax(pmin(stute_F2v_estimator_oc(t21, t22, z1, t2, delta1),1),1e-12)
  
  return(num/denom)
}

Fn_tilde_stute_oc <- Vectorize(Fn_tilde_stute_oc, vectorize.args = c("t1"))

Fn_tilde_linying_oc <- function(t1, t21, t22, z1, t2, delta1){
  num = linying_Fv_estimator_oc(t1, t21, t22, z1, t2, delta1)
  denom = pmax(pmin(linying_F2v_estimator_oc(t21, t22, z1, t2, delta1),1),1e-12) 
  
  return(num/denom)
}

Fn_tilde_linying_oc <- Vectorize(Fn_tilde_linying_oc, vectorize.args = c("t1"))

Fn_tilde_stute2_oc <- function(t1, t21, t22, z1, t2, delta1){
  num = stute_Fv_estimator_oc(t1, t21, t22, z1, t2, delta1)
  denom = pmax(pmin(emp_F2v_estimator(t21, t22, t2),1),num)
  
  return(num/denom)
}

Fn_tilde_stute2_oc <- Vectorize(Fn_tilde_stute2_oc, vectorize.args = c("t1"))

Fn_tilde_linying2_oc <- function(t1, t21, t22, z1, t2, delta1){
  num = linying_Fv_estimator_oc(t1, t21, t22, z1, t2, delta1)
  denom = pmax(pmin(emp_F2v_estimator(t21, t22, t2),1),num)

  return(num/denom)
}

Fn_tilde_linying2_oc <- Vectorize(Fn_tilde_linying2_oc, vectorize.args = c("t1"))

######################################################################
# 7. Estimators for Ftilde(t1 | t21, t22) under univariate censoring #
######################################################################

###################################################################################################
# 7.1. Estimators for Ftilde(t1 | t21, t22) based on inverse probability weighted or LY estimator #
###################################################################################################

Fn_tilde_stute <- function(t1, t21, t22, z1, z2, delta1, delta2){
  num = stute_Fv_estimator(t1, t21, t22, z1, z2, delta1, delta2)
  denom = pmax(pmin(stute_F2v_estimator(t21, t22, z1, z2, delta1, delta2),1),1e-12)
  
  return(num/denom)
}

Fn_tilde_stute <- Vectorize(Fn_tilde_stute, vectorize.args = c("t1"))

Fn_tilde_linying <- function(t1, t21, t22, z1, z2, delta1, delta2){
  num = linying_Fv_estimator(t1, t21, t22, z1, z2, delta1, delta2)
  denom = pmax(pmin(linying_F2v_estimator(t21, t22, z1, z2, delta1, delta2),1),1e-12) 
  
  return(num/denom)
}

Fn_tilde_linying <- Vectorize(Fn_tilde_linying, vectorize.args = c("t1"))

Fn_tilde_stute2 <- function(t1, t21, t22, z1, z2, delta1, delta2){
  num = stute_Fv_estimator(t1, t21, t22, z1, z2, delta1, delta2)
  denom = pmax(pmin(km_F2v_estimator(t21, t22, z1, z2, delta1, delta2),1),num)
  
  return(num/denom)
}

Fn_tilde_stute2 <- Vectorize(Fn_tilde_stute2, vectorize.args = c("t1"))

Fn_tilde_linying2 <- function(t1, t21, t22, z1, z2, delta1, delta2){
  num = linying_Fv_estimator(t1, t21, t22, z1, z2, delta1, delta2)
  denom = pmax(pmin(km_F2v_estimator(t21, t22, z1, z2, delta1, delta2),1),num)

  return(num/denom)
}

Fn_tilde_linying2 <- Vectorize(Fn_tilde_linying2, vectorize.args = c("t1"))

Fn_tilde_stute_all <- function(t1, t21, t22, z1, z2, delta1, delta2){
  num = stute_Fv_estimator(t1, t21, t22, z1, z2, delta1, delta2)
  denom = pmax(pmin(stute_F2v_estimator(t21, t22, z1, z2, delta1, delta2),1),1e-12)
  
  return(list(value = num/denom, num = num, denom = denom))
}

Fn_tilde_stute_all <- Vectorize(Fn_tilde_stute_all, vectorize.args = c("t1"))

Fn_tilde_linying_all <- function(t1, t21, t22, z1, z2, delta1, delta2){
  num = linying_Fv_estimator(t1, t21, t22, z1, z2, delta1, delta2)
  denom = pmax(pmin(linying_F2v_estimator(t21, t22, z1, z2, delta1, delta2),1),1e-12) 
  
  return(list(value = num/denom, num = num, denom = denom))
}

Fn_tilde_linying_all <- Vectorize(Fn_tilde_linying_all, vectorize.args = c("t1"))

Fn_tilde_stute2_all <- function(t1, t21, t22, z1, z2, delta1, delta2){
  num = stute_Fv_estimator(t1, t21, t22, z1, z2, delta1, delta2)
  denom = pmax(pmin(km_F2v_estimator(t21, t22, z1, z2, delta1, delta2),1),num)
  
  return(list(value = num/denom, num = num, denom = denom))
}

Fn_tilde_stute2_all <- Vectorize(Fn_tilde_stute2_all, vectorize.args = c("t1"))

Fn_tilde_linying2_all <- function(t1, t21, t22, z1, z2, delta1, delta2){
  num = linying_Fv_estimator(t1, t21, t22, z1, z2, delta1, delta2)
  denom = pmax(pmin(km_F2v_estimator(t21, t22, z1, z2, delta1, delta2),1),num)
  
  return(list(value = num/denom, num = num, denom = denom))
}

Fn_tilde_linying2_all <- Vectorize(Fn_tilde_linying2_all, vectorize.args = c("t1"))

########################################################
# 7.2. Vectorized estimators for Ftilde(t1 | t21, t22) #
########################################################

Fn_tilde_stute_vec <- function(t1, t21, t22, z1, z2, delta1, delta2){
	num = vector()
	denom = vector()
	
	for (k in 1:length(t1)){
		num[k] = stute_Fv_estimator(t1[k], t21[k], t22[k], z1, z2, delta1, delta2)
  	denom[k] = pmax(pmin(stute_F2v_estimator(t21[k], t22[k], z1, z2, delta1, delta2),1),1e-12)
	}

  return(num/denom)
}

Fn_tilde_linying_vec <- function(t1, t21, t22, z1, z2, delta1, delta2){
	num = vector()
	denom = vector()
	
	for (k in 1:length(t1)){
  		num[k] = linying_Fv_estimator(t1[k], t21[k], t22[k], z1, z2, delta1, delta2)
  		denom[k] = pmax(pmin(linying_F2v_estimator(t21[k], t22[k], z1, z2, delta1, delta2),1),1e-12) 
	}
  
  return(num/denom)
}

Fn_tilde_stute2_vec <- function(t1, t21, t22, z1, z2, delta1, delta2){
  num = stute_Fv_estimator_vec(t1, t21, t22, z1, z2, delta1, delta2)
  denom = pmax(pmin(km_F2v_estimator(t21, t22, z1, z2, delta1, delta2),1),num)
  
  return(num/denom)
}

######################################################################################
# 7.3. Vectorized estimators for Ftilde(t1 | t21, t22) using smoothing for F(t1, t2) #
######################################################################################

pavit_fun <- function(values){
  new_values = values
  for (k in 1:(length(values)-1)){
    if (new_values[k] > new_values[k+1]){
      new_values[k:(k+1)] = (new_values[k] + new_values[k+1])/2
      j <- k + 1
      for (l in (j - 1):1){
        if (new_values[l] > new_values[j]){
          pool <- sum(values[l:j])/length(l:j)
          new_values[l:j] <- pool
        }
      }
    }
  }
  return(new_values)
}

Fn_tilde_linying_vec2 <- function(t1, t21, t22, z1, z2, delta1, delta2, smoothing = F){
  num1 = vector()
  num2 = vector()
  denom = vector()
  
  num1 = F_linying(t1, t22, z1, z2, delta1, delta2)
  num2 = F_linying(t1, t21, z1, z2, delta1, delta2)
  
  if (smoothing == F){
    num = num1 - num2 
  }
  
  if (smoothing == T){
    num = pavit_fun(num1) - pavit_fun(num2) 
  }
  
  denom = pmax(pmin(km_F2v_estimator(t21, t22, z1, z2, delta1, delta2),1),num)
  
  return(num/denom)
}

##############################################################
# 8. Estimation of the quantile function Q(p | t1, t21, t22) #
##############################################################

#############################################################################################
# 8.1 Estimation of the quantile function Q(p | t1, t21, t22) under one-component censoring #
#############################################################################################

Q_tilde_oc <- function(p, tt1, t21, t22, z1, t2, delta1, method = "km", 
                       pavit = TRUE, show.plots = FALSE, correct = TRUE){
  
  z1_grid = c(0,sort(z1))
  
  if (method == "km"){
    v = Fn_tilde_stute2_oc(tt1, t21, t22, z1, t2, delta1)
    w = Fn_tilde_linying2_oc(tt1, t21, t22, z1, t2, delta1)  
    
    vgrid_vals = Fn_tilde_stute2_oc(z1_grid, t21, t22, z1, t2, delta1)
    wgrid_vals = Fn_tilde_linying2_oc(z1_grid, t21, t22, z1, t2, delta1)
    
    # Non-monotonicity of the Lin-Ying estimator
    linying_error_flag = 0
    if (max(wgrid_vals) == 1 & wgrid_vals[length(wgrid_vals)] < 1){
      linying_error_flag = 1}
    
    # Create a monotonically increasing function 
    if (pavit == TRUE){wgrid_vals = pavit_fun(wgrid_vals)}
  }
  
  if (method == "full"){
    v = Fn_tilde_stute_oc(tt1, t21, t22, z1, t2, delta1)
    w = Fn_tilde_linying_oc(tt1, t21, t22, z1, t2, delta1)  
    
    vgrid_vals = Fn_tilde_stute_oc(z1_grid, t21, t22, z1, t2, delta1)
    wgrid_vals = Fn_tilde_linying_oc(z1_grid, t21, t22, z1, t2, delta1)
    
    # Non-monotonicity of the Lin-Ying estimator
    linying_error_flag = 0
    if (max(wgrid_vals) == 1 & wgrid_vals[length(wgrid_vals)] < 1){
      linying_error_flag = 1}
    
    # Create a monotonically increasing function 
    if (pavit == TRUE){wgrid_vals = pavit_fun(wgrid_vals)}
  }

  # T1-values for which the Stute estimator of Fv makes a jump
  jumps_stute = delta1[order(z1)]*as.numeric(t21 <= t2[order(z1)] & t2[order(z1)] <= t22)

  # Maximum values for the Stute and Lin-Ying estimator
  max_Fn_stute = max(vgrid_vals);
  max_Fn_linying = max(wgrid_vals);
  
  # For which Z1-value is the maximum reached
  max_Fn_inv_stute = z1_grid[min(which(vgrid_vals == max_Fn_stute))];
  max_Fn_inv_linying = z1_grid[min(which(wgrid_vals == max_Fn_linying))]; 
  
  # Error debugging
  stute_error_flag = 0
  if (is.na(max_Fn_stute) == T){
    message("Stute estimator Fn_tilde takes NA value")}
  else if (max_Fn_stute > 1){
    stute_error_flag = 1
    message("Stute estimator Fn_tilde takes non-sensible result")}
  
  if (is.na(max_Fn_linying) == T){
    message("Lin-Ying estimator Fn_tilde takes NA value")}
  else if (max_Fn_linying > 1){
    message("Lin-Ying estimator Fn_tilde takes non-sensible result")} 
  
  # Additional constraints on Ftilde
  vgrid_vals = pmin(vgrid_vals,1);
  wgrid_vals = pmin(wgrid_vals,1);
  
  # Arguments of Fn_tilde_inv
  arg1 = p + (1-p)*v;
  arg2 = p + (1-p)*w; 

  # Solving the equation x* = Fn_tilde_inv(y)
  lower1 = matrix(NA, nrow = 1, ncol = length(arg1));
  upper1 = matrix(NA, nrow = 1, ncol = length(arg1));
  lower2 = matrix(NA, nrow = 1, ncol = length(arg2));
  upper2 = matrix(NA, nrow = 1, ncol = length(arg2));

  for (jj in 1:length(p)){
    lower1[jj] = max(which(vgrid_vals < as.numeric(arg1[jj])))
    lower2[jj] = max(which(wgrid_vals < as.numeric(arg2[jj])))
    
    upper1[jj] = ifelse(arg1[jj] <= max(vgrid_vals), min(which(vgrid_vals >= as.numeric(arg1[jj]))), NA)
    upper2[jj] = ifelse(arg2[jj] <= max(wgrid_vals), min(which(wgrid_vals >= as.numeric(arg2[jj]))), NA)
    
    if (is.na(upper1[jj]) == F){
      if (lower1[jj] > upper1[jj]){message(paste0("Stute estimator is non-monotone: p_index ",jj))}
    }
    if (is.na(upper2[jj]) == F){
      if (lower2[jj] > upper2[jj]){message(paste0("Lin-Ying estimator is non-monotone: p_index ",jj))}
    }
  }

  qt1_lower = -tt1 + pmax(tt1, z1_grid[lower1])
  qt1_upper = -tt1 + pmax(tt1, z1_grid[upper1])
  qt2_lower = -tt1 + pmax(tt1, z1_grid[lower2])
  qt2_upper = -tt1 + pmax(tt1, z1_grid[upper2])

  if (show.plots == TRUE){
    par(mfrow = c(2,2))
    plot(z1_grid, vgrid_vals, xlab = "Z1-values", ylab = "Fn_tilde"); 
    abline(h = arg1[1], col = 2, lwd = 2); 
    abline(v = z1_grid[lower1[1]], col = 4, lwd = 2); 
    abline(v = tt1, col = 5, lwd = 2);
    abline(v = tt1 + qt1_lower[1], col = 6, lwd = 2); 
    
    plot(z1_grid, vgrid_vals, xlab = "Z1-values", ylab = "Fn_tilde"); 
    abline(h = arg1[length(p)], col = 2, lwd = 2); 
    abline(v = z1_grid[lower1[length(p)]], col = 4, lwd = 2); 
    abline(v = tt1, col = 5, lwd = 2); 
    abline(v = tt1 + qt1_lower[length(p)], col = 6, lwd = 2); 
    
    plot(z1_grid, wgrid_vals, xlab = "Z1-values", ylab = "Fn_tilde"); 
    abline(h = arg2[1], col = 2, lwd = 2); 
    abline(v = z1_grid[lower2[1]], col = 4, lwd = 2); 
    abline(v = tt1, col = 5, lwd = 2);
    abline(v = tt1 + qt2_lower[1], col = 6, lwd = 2); 
    
    plot(z1_grid, wgrid_vals, xlab = "Z1-values", ylab = "Fn_tilde"); 
    abline(h = arg2[length(p)], col = 2, lwd = 2); 
    abline(v = z1_grid[lower2[length(p)]], col = 4, lwd = 2); 
    abline(v = tt1, col = 5, lwd = 2); 
    abline(v = tt1 + qt2_lower[length(p)], col = 6, lwd = 2); 
  }
  
  # Error checking
  if (sum(qt1_lower - qt1_upper > 0, na.rm = T) > 0){message("Natural ordering of limits incorrect")}
  if (sum(qt2_lower - qt2_upper > 0, na.rm = T) > 0){message("Natural ordering of limits incorrect")}
  
  upper_Fn_stute = ifelse(sum(is.na(qt1_upper)) > 0, p[min(which(is.na(qt1_upper) == T))], p[min(which(z1_grid[upper1] == max_Fn_inv_stute))])
  upper_Fn_linying = ifelse(sum(is.na(qt2_upper)) > 0, p[min(which(is.na(qt2_upper) == T))], p[min(which(z1_grid[upper2] == max_Fn_inv_linying))])
  
  if (correct == TRUE){
    condition1_ids = which(qt1_lower == max(qt1_lower))
    if (length(condition1_ids) > 1){
      qt1_lower[condition1_ids[-1]] = NA
    }
    condition2_ids = which(qt1_upper == max(qt1_upper))
    if (length(condition2_ids) > 1){
      qt1_upper[condition2_ids[-1]] = NA
    }
    condition3_ids = which(qt2_lower == max(qt2_lower))
    if (length(condition3_ids) > 1){
      qt2_lower[condition3_ids[-1]] = NA
    }
    condition4_ids = which(qt2_upper == max(qt2_upper))
    if (length(condition4_ids) > 1){
      qt2_upper[condition4_ids[-1]] = NA
    }
  }
  
  return(list(pvals = p,
              arg1 = arg1, arg2 = arg2,
              vgrid_vals = vgrid_vals,
              wgrid_vals = wgrid_vals,
              qt_tilde_lower = qt1_lower, qt_tilde_upper = qt1_upper, 
              qt_tilde_midpoint = (qt1_lower + qt1_upper)/2,
              qt_dtilde_lower = qt2_lower, qt_dtilde_upper = qt2_upper, 
              qt_dtilde_midpoint = (qt2_lower + qt2_upper)/2,
              l1 = lower1, u1 = upper1, l2 = lower2, u2 = upper2,
              max_Fn_stute = max_Fn_stute,
              max_Fn_linying = max_Fn_linying,
              max_Fn_inv_stute = max_Fn_inv_stute,
              max_Fn_inv_linying = max_Fn_inv_linying,
              max_Q_stute = max_Fn_inv_stute - tt1,
              max_Q_linying = max_Fn_inv_linying - tt1,
              upper_Fn_stute = upper_Fn_stute,
              upper_Fn_linying = upper_Fn_linying,
              stute_error_flag = stute_error_flag,
              linying_error_flag = linying_error_flag))
}

Q_fun_oc <- function(p, tt1, tt21, tt22, z1, t2, delta1, method = "km",
                     pavit = "TRUE", add_constraint = "TRUE"){
  if (length(tt1) != length(tt21) | length(tt1) != length(tt22)){
    paste0("tt1, tt21 and tt22 should have the same dimensions")};
  
  z1_grid = c(0,sort(z1))

  if (length(unique(tt21)) == 1 & length(unique(tt22)) == 1){
  	message("Unique limits t21 and t22 specified");
  	  
    if (method == "km"){
      v1 = Fn_tilde_stute2_oc(tt1, tt21[1], tt22[1], z1, t2, delta1)
      w1 = Fn_tilde_linying2_oc(tt1, tt21[1], tt22[1], z1, t2, delta1)

      v_mat = matrix(rep(v1, length(p)), nrow = length(p), ncol = length(tt1), byrow = T)
      w_mat = matrix(rep(w1, length(p)), nrow = length(p), ncol = length(tt1), byrow = T)  

      vz = Fn_tilde_stute2_oc(z1_grid, tt21[1], tt22[1], z1, t2, delta1)
      wz = Fn_tilde_linying2_oc(z1_grid, tt21[1], tt22[1], z1, t2, delta1)
    }
    
    if (method == "full"){
      v1 = Fn_tilde_stute_oc(tt1, tt21[1], tt22[1], z1, t2, delta1)
      w1 = Fn_tilde_linying_oc(tt1, tt21[1], tt22[1], z1, t2, delta1)
      
      v_mat = matrix(rep(v1, length(p)), nrow = length(p), ncol = length(tt1), byrow = T)
      w_mat = matrix(rep(w1, length(p)), nrow = length(p), ncol = length(tt1), byrow = T)  
      
      vz = Fn_tilde_stute_oc(z1_grid, tt21[1], tt22[1], z1, t2, delta1)
      wz = Fn_tilde_linying_oc(z1_grid, tt21[1], tt22[1], z1, t2, delta1)
    }
    
    # Additional constraint to ensure that Stute estimator jumps to one
    if (add_constraint == TRUE){
      vz[which(z1_grid == max(z1_grid))] = 1
    }

    vgrid_vals = matrix(rep(pmin(vz,1), length(tt1)), nrow = length(z1_grid), ncol = length(tt1), byrow = F)
      
    # Create a monotonically increasing function 
    if (pavit == TRUE){
      wgrid_vals = matrix(rep(pmin(pavit_fun(wz),1), length(tt1)), nrow = length(z1_grid), ncol = length(tt1), byrow = F)
    }    
    if (pavit == FALSE){
      wgrid_vals = matrix(rep(pmin(wz,1), length(tt1)), nrow = length(z1_grid), ncol = length(tt1), byrow = F)
    }  
    message("Preprocessing done");    
  }

  if (length(unique(tt21)) != 1 | length(unique(tt22)) != 1){
  	message("Different limits t21 and t22 specified");
  	  
  	v_mat = matrix(nrow = length(p), ncol = length(tt1))
  	w_mat = matrix(nrow = length(p), ncol = length(tt1))
  	  
  	vgrid_vals = matrix(nrow = length(z1_grid), ncol = length(tt1))
  	wgrid_vals = matrix(nrow = length(z1_grid), ncol = length(tt1))
  	  
  	for (col_id in 1:length(tt1)){
  	  if (method == "km"){
  	    v = Fn_tilde_stute2_oc(tt1[col_id], tt21[col_id], tt22[col_id], z1, t2, delta1)
  	    w = Fn_tilde_linying2_oc(tt1[col_id], tt21[col_id], tt22[col_id], z1, t2, delta1)
  	    
  	    v_mat[, col_id] = rep(v, length(p))
  	    w_mat[, col_id] = rep(w, length(p))  
  	    
  	    vgrid_vals[, col_id] = Fn_tilde_stute2_oc(z1_grid, tt21[col_id], tt22[col_id], z1, t2, delta1)
  	    
  	    if (pavit == TRUE){wgrid_vals[, col_id] = pavit_fun(Fn_tilde_linying2_oc(z1_grid, tt21[col_id], tt22[col_id], z1, t2, delta1))}
  	    if (pavit == FALSE){wgrid_vals[, col_id] = Fn_tilde_linying2_oc(z1_grid, tt21[col_id], tt22[col_id], z1, t2, delta1)}
  	  }
  	  if (method == "full"){
  	    v = Fn_tilde_stute_oc(tt1[col_id], tt21[col_id], tt22[col_id], z1, t2, delta1)
  	    w = Fn_tilde_linying_oc(tt1[col_id], tt21[col_id], tt22[col_id], z1, t2, delta1)
  	    
  	    v_mat[, col_id] = rep(v, length(p))
  	    w_mat[, col_id] = rep(w, length(p))  
  	    
  	    vgrid_vals[, col_id] = Fn_tilde_stute_oc(z1_grid, tt21[col_id], tt22[col_id], z1, t2, delta1)
  	    
  	    if (pavit == TRUE){wgrid_vals[, col_id] = pavit_fun(Fn_tilde_linying_oc(z1_grid, tt21[col_id], tt22[col_id], z1, t2, delta1))}
  	    if (pavit == FALSE){wgrid_vals[, col_id] = Fn_tilde_linying_oc(z1_grid, tt21[col_id], tt22[col_id], z1, t2, delta1)}
  	  }
  	 }
    message("Preprocessing done");    
  }
 
  # Create matrix of p-values for different tt1-values (with corresponding limits)  
  p_mat = matrix(p, nrow = length(p), ncol = length(tt1), byrow = F)
  
  # Arguments of Fn_tilde_inv: Dimension length(p) x length(tt1)
  arg1 = p_mat + (1-p_mat)*v_mat;
  arg2 = p_mat + (1-p_mat)*w_mat;

  # Solving the equation x* = Fn_tilde_inv(y)  
  lower1_mat = matrix(NA, nrow = nrow(arg1), ncol = ncol(arg1));
  upper1_mat = matrix(NA, nrow = nrow(arg1), ncol = ncol(arg1));
  lower2_mat = matrix(NA, nrow = nrow(arg2), ncol = ncol(arg2));
  upper2_mat = matrix(NA, nrow = nrow(arg2), ncol = ncol(arg2));
  
  qt1_lower_mat = matrix(NA, nrow = nrow(arg1), ncol = ncol(arg1));
  qt1_upper_mat = matrix(NA, nrow = nrow(arg1), ncol = ncol(arg1));
  qt2_lower_mat = matrix(NA, nrow = nrow(arg2), ncol = ncol(arg2));
  qt2_upper_mat = matrix(NA, nrow = nrow(arg2), ncol = ncol(arg2));
  
  for (row_id in 1:nrow(arg1)){
    for (col_id in 1:ncol(arg1)){
    	#message(paste0("Computation: p", row_id, " from ", nrow(arg1), " t1-value ", col_id, " from ", ncol(arg1)));
    	lower1_mat[row_id, col_id] = max(which(vgrid_vals[,col_id] < as.numeric(arg1[row_id, col_id])))
      lower2_mat[row_id, col_id] = max(which(wgrid_vals[,col_id] < as.numeric(arg2[row_id, col_id])))
    
      upper1_mat[row_id, col_id] = ifelse(arg1[row_id, col_id] <= max(vgrid_vals[,col_id]), min(which(vgrid_vals[,col_id] >= as.numeric(arg1[row_id, col_id]))), NA)
      upper2_mat[row_id, col_id] = ifelse(arg2[row_id, col_id] <= max(wgrid_vals[,col_id]), min(which(wgrid_vals[,col_id] >= as.numeric(arg2[row_id, col_id]))), NA)
    
    	if (is.na(upper1_mat[row_id, col_id]) == F){
      		if (lower1_mat[row_id, col_id] > upper1_mat[row_id, col_id]){message(paste0("Stute estimator is non-monotone: p_index ", row_id, col_id))}
    	}
    	if (is.na(upper2_mat[row_id, col_id]) == F){
      		if (lower2_mat[row_id, col_id] > upper2_mat[row_id, col_id]){message(paste0("Lin-Ying estimator is non-monotone: p_index ", row_id, col_id))}
    	}
  	}
  	
  	qt1_lower_mat[row_id, ] = -tt1 + pmax(tt1,z1_grid[lower1_mat[row_id, ]])
    qt1_upper_mat[row_id, ] = -tt1 + pmax(tt1,z1_grid[upper1_mat[row_id, ]])
    qt2_lower_mat[row_id, ] = -tt1 + pmax(tt1,z1_grid[lower2_mat[row_id, ]])
    qt2_upper_mat[row_id, ] = -tt1 + pmax(tt1,z1_grid[upper2_mat[row_id, ]])
  }

  return(list(pvals = p,
              arg1 = arg1, arg2 = arg2,
              vgrid_vals = vgrid_vals,
              wgrid_vals = wgrid_vals,
              qt_tilde_lower = qt1_lower_mat, qt_tilde_upper = qt1_upper_mat, 
              qt_tilde_midpoint = (qt1_lower_mat + qt1_upper_mat)/2,
              qt_dtilde_lower = qt2_lower_mat, qt_dtilde_upper = qt2_upper_mat, 
              qt_dtilde_midpoint = (qt2_lower_mat + qt2_upper_mat)/2,
              l1 = lower1_mat, u1 = upper1_mat, l2 = lower2_mat, u2 = upper2_mat))
}

###########################################################################################
# 8.2. Estimation of the quantile function Q(p | t1, t21, t22) under univariate censoring #
###########################################################################################

Q_tilde <- function(p, tt1, t21, t22, z1, z2, delta1, delta2, method = "km", 
                    pavit = TRUE, show.plots = FALSE, correct = TRUE){
  
  z1_grid = c(0,sort(z1))
  
  if (method == "km"){
    v = Fn_tilde_stute2(tt1, t21, t22, z1, z2, delta1, delta2)
    w = Fn_tilde_linying2(tt1, t21, t22, z1, z2, delta1, delta2)
    
    vgrid_vals = Fn_tilde_stute2(z1_grid, t21, t22, z1, z2, delta1, delta2)
    wgrid_vals = Fn_tilde_linying2(z1_grid, t21, t22, z1, z2, delta1, delta2)
    
    # Non-monotonicity of the Lin-Ying estimator
    linying_error_flag = 0
    if (max(wgrid_vals) == 1 & wgrid_vals[length(wgrid_vals)] < 1){
      linying_error_flag = 1}
    
    # Create a monotonically increasing function 
    if (pavit == TRUE){wgrid_vals = pavit_fun(wgrid_vals)}
  }

  if (method == "full"){
    v = Fn_tilde_stute(tt1, t21, t22, z1, z2, delta1, delta2)
    w = Fn_tilde_linying(tt1, t21, t22, z1, z2, delta1, delta2)  

    vgrid_vals = Fn_tilde_stute(z1_grid, t21, t22, z1, z2, delta1, delta2)
    wgrid_vals = Fn_tilde_linying(z1_grid, t21, t22, z1, z2, delta1, delta2)
    
    # Non-monotonicity of the Lin-Ying estimator
    linying_error_flag = 0
    if (max(wgrid_vals) == 1 & wgrid_vals[length(wgrid_vals)] < 1){
      linying_error_flag = 1}

    # Create a monotonically increasing function 
    if (pavit == TRUE){wgrid_vals = pavit_fun(wgrid_vals)}
  }

  # T1-values for which the Stute estimator of Fv makes a jump
  jumps_stute = delta1[order(z1)]*delta2[order(z1)]*as.numeric(t21 < z2[order(z1)] & z2[order(z1)] <= t22)

  # Maximum values for the Stute and Lin-Ying estimator
  max_Fn_stute = max(vgrid_vals);
  max_Fn_linying = max(wgrid_vals);
  
  # For which Z1-value is the maximum reached
  max_Fn_inv_stute = z1_grid[min(which(vgrid_vals == max_Fn_stute))];
  max_Fn_inv_linying = z1_grid[min(which(wgrid_vals == max_Fn_linying))]; 
  
  # Error debugging
  stute_error_flag = 0
  if (is.na(max_Fn_stute) == T){
    message("Stute estimator Fn_tilde takes NA value")}
  else if (max_Fn_stute > 1){
    stute_error_flag = 1
    message("Stute estimator Fn_tilde takes non-sensible result")}
  
  if (is.na(max_Fn_linying) == T){
    message("Lin-Ying estimator Fn_tilde takes NA value")}
  else if (max_Fn_linying > 1){
    message("Lin-Ying estimator Fn_tilde takes non-sensible result")} 
  
  # Additional constraints on Ftilde
  vgrid_vals = pmin(vgrid_vals,1);
  wgrid_vals = pmin(wgrid_vals,1);
  
  # Arguments of Fn_tilde_inv
  arg1 = p + (1-p)*v;
  arg2 = p + (1-p)*w; 

  # Solving the equation x* = Fn_tilde_inv(y)
  lower1 = matrix(NA, nrow = 1, ncol = length(arg1));
  upper1 = matrix(NA, nrow = 1, ncol = length(arg1));
  lower2 = matrix(NA, nrow = 1, ncol = length(arg2));
  upper2 = matrix(NA, nrow = 1, ncol = length(arg2));

  for (jj in 1:length(p)){
    lower1[jj] = max(which(vgrid_vals < as.numeric(arg1[jj])))
    lower2[jj] = max(which(wgrid_vals < as.numeric(arg2[jj])))
    
    upper1[jj] = ifelse(arg1[jj] <= max(vgrid_vals), min(which(vgrid_vals >= as.numeric(arg1[jj]))), NA)
    upper2[jj] = ifelse(arg2[jj] <= max(wgrid_vals), min(which(wgrid_vals >= as.numeric(arg2[jj]))), NA)
    
    if (is.na(upper1[jj]) == F){
      if (lower1[jj] > upper1[jj]){message(paste0("Stute estimator is non-monotone: p_index ",jj))}
    }
    if (is.na(upper2[jj]) == F){
      if (lower2[jj] > upper2[jj]){message(paste0("Lin-Ying estimator is non-monotone: p_index ",jj))}
    }
  }

  qt1_lower = -tt1 + pmax(tt1, z1_grid[lower1])
  qt1_upper = -tt1 + pmax(tt1, z1_grid[upper1])
  qt2_lower = -tt1 + pmax(tt1, z1_grid[lower2])
  qt2_upper = -tt1 + pmax(tt1, z1_grid[upper2])

  if (show.plots == TRUE){
    par(mfrow = c(2,2))
    plot(z1_grid, vgrid_vals, xlab = "Z1-values", ylab = "Fn_tilde"); 
    abline(h = arg1[1], col = 2, lwd = 2); 
    abline(v = z1_grid[lower1[1]], col = 4, lwd = 2); 
    abline(v = tt1, col = 5, lwd = 2);
    abline(v = tt1 + qt1_lower[1], col = 6, lwd = 2); 
    
    plot(z1_grid, vgrid_vals, xlab = "Z1-values", ylab = "Fn_tilde"); 
    abline(h = arg1[length(p)], col = 2, lwd = 2); 
    abline(v = z1_grid[lower1[length(p)]], col = 4, lwd = 2); 
    abline(v = tt1, col = 5, lwd = 2); 
    abline(v = tt1 + qt1_lower[length(p)], col = 6, lwd = 2); 
    
    plot(z1_grid, wgrid_vals, xlab = "Z1-values", ylab = "Fn_tilde"); 
    abline(h = arg2[1], col = 2, lwd = 2); 
    abline(v = z1_grid[lower2[1]], col = 4, lwd = 2); 
    abline(v = tt1, col = 5, lwd = 2);
    abline(v = tt1 + qt2_lower[1], col = 6, lwd = 2); 
    
    plot(z1_grid, wgrid_vals, xlab = "Z1-values", ylab = "Fn_tilde"); 
    abline(h = arg2[length(p)], col = 2, lwd = 2); 
    abline(v = z1_grid[lower2[length(p)]], col = 4, lwd = 2); 
    abline(v = tt1, col = 5, lwd = 2); 
    abline(v = tt1 + qt2_lower[length(p)], col = 6, lwd = 2); 
  }
  
  # Error checking
  if (sum(qt1_lower - qt1_upper > 0, na.rm = T) > 0){message("Natural ordering of limits incorrect")}
  if (sum(qt2_lower - qt2_upper > 0, na.rm = T) > 0){message("Natural ordering of limits incorrect")}
  
  upper_Fn_stute = ifelse(sum(is.na(qt1_upper)) > 0, p[min(which(is.na(qt1_upper) == T))], p[min(which(z1_grid[upper1] == max_Fn_inv_stute))])
  upper_Fn_linying = ifelse(sum(is.na(qt2_upper)) > 0, p[min(which(is.na(qt2_upper) == T))], p[min(which(z1_grid[upper2] == max_Fn_inv_linying))])
  
  if (correct == TRUE){
    condition1_ids = which(qt1_lower == max(qt1_lower))
    if (length(condition1_ids) > 1){
      qt1_lower[condition1_ids[-1]] = NA
    }
    condition2_ids = which(qt1_upper == max(qt1_upper))
    if (length(condition2_ids) > 1){
      qt1_upper[condition2_ids[-1]] = NA
    }
    condition3_ids = which(qt2_lower == max(qt2_lower))
    if (length(condition3_ids) > 1){
      qt2_lower[condition3_ids[-1]] = NA
    }
    condition4_ids = which(qt2_upper == max(qt2_upper))
    if (length(condition4_ids) > 1){
      qt2_upper[condition4_ids[-1]] = NA
    }
  }
  
  return(list(pvals = p,
              arg1 = arg1, arg2 = arg2,
              vgrid_vals = vgrid_vals,
              wgrid_vals = wgrid_vals,
              qt_tilde_lower = qt1_lower, qt_tilde_upper = qt1_upper, 
              qt_tilde_midpoint = (qt1_lower + qt1_upper)/2,
              qt_dtilde_lower = qt2_lower, qt_dtilde_upper = qt2_upper, 
              qt_dtilde_midpoint = (qt2_lower + qt2_upper)/2,
              l1 = lower1, u1 = upper1, l2 = lower2, u2 = upper2,
              max_Fn_stute = max_Fn_stute,
              max_Fn_linying = max_Fn_linying,
              max_Fn_inv_stute = max_Fn_inv_stute,
              max_Fn_inv_linying = max_Fn_inv_linying,
              max_Q_stute = max_Fn_inv_stute - tt1,
              max_Q_linying = max_Fn_inv_linying - tt1,
              upper_Fn_stute = upper_Fn_stute,
              upper_Fn_linying = upper_Fn_linying,
              stute_error_flag = stute_error_flag,
              linying_error_flag = linying_error_flag))
}

Q_fun <- function(p, tt1, tt21, tt22, z1, z2, delta1, delta2, method = "km",
                  pavit = "TRUE", add_constraint = "TRUE"){
  if (length(tt1) != length(tt21) | length(tt1) != length(tt22)){paste0("tt1, tt21 and tt22 should have the same dimensions")};
  
  z1_grid = c(0,sort(z1))
  
  if (length(unique(tt21)) == 1 & length(unique(tt22)) == 1){
    message("Unique limits t21 and t22 specified");
    
    if (method == "km"){
      v1 = Fn_tilde_stute2(tt1, tt21[1], tt22[1], z1, z2, delta1, delta2)
      w1 = Fn_tilde_linying2(tt1, tt21[1], tt22[1], z1, z2, delta1, delta2)
      
      v_mat = matrix(rep(v1, length(p)), nrow = length(p), ncol = length(tt1), byrow = T)
      w_mat = matrix(rep(w1, length(p)), nrow = length(p), ncol = length(tt1), byrow = T)  
      
      vz = Fn_tilde_stute2(z1_grid, tt21[1], tt22[1], z1, z2, delta1, delta2)
      wz = Fn_tilde_linying2(z1_grid, tt21[1], tt22[1], z1, z2, delta1, delta2)
    }
    
    if (method == "full"){
      v1 = Fn_tilde_stute(tt1, tt21[1], tt22[1], z1, z2, delta1, delta2)
      w1 = Fn_tilde_linying(tt1, tt21[1], tt22[1], z1, z2, delta1, delta2)
      
      v_mat = matrix(rep(v1, length(p)), nrow = length(p), ncol = length(tt1), byrow = T)
      w_mat = matrix(rep(w1, length(p)), nrow = length(p), ncol = length(tt1), byrow = T)  
      
      vz = Fn_tilde_stute(z1_grid, tt21[1], tt22[1], z1, z2, delta1, delta2)
      wz = Fn_tilde_linying(z1_grid, tt21[1], tt22[1], z1, z2, delta1, delta2)
    }
    
    # Additional constraint to ensure that Stute estimator jumps to one
    if (add_constraint == TRUE){
      vz[which(z1_grid == max(z1_grid))] = 1
    }
    
    vgrid_vals = matrix(rep(pmin(vz,1), length(tt1)), nrow = length(z1_grid), ncol = length(tt1), byrow = F)
    
    # Create a monotonically increasing function 
    if (pavit == TRUE){
      wgrid_vals = matrix(rep(pmin(pavit_fun(wz),1), length(tt1)), nrow = length(z1_grid), ncol = length(tt1), byrow = F)
    }    
    if (pavit == FALSE){
      wgrid_vals = matrix(rep(pmin(wz,1), length(tt1)), nrow = length(z1_grid), ncol = length(tt1), byrow = F)
    }  
    message("Preprocessing done");    
  }
  
  if (length(unique(tt21)) != 1 | length(unique(tt22)) != 1){
    message("Different limits t21 and t22 specified");
    
    v_mat = matrix(nrow = length(p), ncol = length(tt1))
    w_mat = matrix(nrow = length(p), ncol = length(tt1))
    
    vgrid_vals = matrix(nrow = length(z1_grid), ncol = length(tt1))
    wgrid_vals = matrix(nrow = length(z1_grid), ncol = length(tt1))
    
    for (col_id in 1:length(tt1)){
      if (method == "full"){
        v = Fn_tilde_stute(tt1[col_id], tt21[col_id], tt22[col_id], z1, z2, delta1, delta2)
        w = Fn_tilde_linying(tt1[col_id], tt21[col_id], tt22[col_id], z1, z2, delta1, delta2)
        
        v_mat[, col_id] = rep(v, length(p))
        w_mat[, col_id] = rep(w, length(p))  
        
        vgrid_vals[, col_id] = Fn_tilde_stute(z1_grid, tt21[col_id], tt22[col_id], z1, z2, delta1, delta2)
        
        if (pavit == TRUE){wgrid_vals[, col_id] = pavit_fun(Fn_tilde_linying(z1_grid, tt21[col_id], tt22[col_id], z1, z2, delta1, delta2))}
        if (pavit == FALSE){wgrid_vals[, col_id] = Fn_tilde_linying(z1_grid, tt21[col_id], tt22[col_id], z1, z2, delta1, delta2)}
      }
      if (method == "km"){
        v = Fn_tilde_stute2(tt1[col_id], tt21[col_id], tt22[col_id], z1, z2, delta1, delta2)
        w = Fn_tilde_linying2(tt1[col_id], tt21[col_id], tt22[col_id], z1, z2, delta1, delta2)
        
        v_mat[, col_id] = rep(v, length(p))
        w_mat[, col_id] = rep(w, length(p))  
        
        vgrid_vals[, col_id] = Fn_tilde_stute2(z1_grid, tt21[col_id], tt22[col_id], z1, z2, delta1, delta2)
        
        if (pavit == TRUE){wgrid_vals[, col_id] = pavit_fun(Fn_tilde_linying2(z1_grid, tt21[col_id], tt22[col_id], z1, z2, delta1, delta2))}
        if (pavit == FALSE){wgrid_vals[, col_id] = Fn_tilde_linying2(z1_grid, tt21[col_id], tt22[col_id], z1, z2, delta1, delta2)}
      }
    }
    message("Preprocessing done");    
  }
  
  # Create matrix of p-values for different tt1-values (with corresponding limits)  
  p_mat = matrix(p, nrow = length(p), ncol = length(tt1), byrow = F)
  
  # Arguments of Fn_tilde_inv: Dimension length(p) x length(tt1)
  arg1 = p_mat + (1-p_mat)*v_mat;
  arg2 = p_mat + (1-p_mat)*w_mat;
  
  # Solving the equation x* = Fn_tilde_inv(y)  
  lower1_mat = matrix(NA, nrow = nrow(arg1), ncol = ncol(arg1));
  upper1_mat = matrix(NA, nrow = nrow(arg1), ncol = ncol(arg1));
  lower2_mat = matrix(NA, nrow = nrow(arg2), ncol = ncol(arg2));
  upper2_mat = matrix(NA, nrow = nrow(arg2), ncol = ncol(arg2));
  
  qt1_lower_mat = matrix(NA, nrow = nrow(arg1), ncol = ncol(arg1));
  qt1_upper_mat = matrix(NA, nrow = nrow(arg1), ncol = ncol(arg1));
  qt2_lower_mat = matrix(NA, nrow = nrow(arg2), ncol = ncol(arg2));
  qt2_upper_mat = matrix(NA, nrow = nrow(arg2), ncol = ncol(arg2));
  
  for (row_id in 1:nrow(arg1)){
    for (col_id in 1:ncol(arg1)){
      #message(paste0("Computation: p", row_id, " from ", nrow(arg1), " t1-value ", col_id, " from ", ncol(arg1)));
      lower1_mat[row_id, col_id] = max(which(vgrid_vals[,col_id] < as.numeric(arg1[row_id, col_id])))
      lower2_mat[row_id, col_id] = max(which(wgrid_vals[,col_id] < as.numeric(arg2[row_id, col_id])))
      
      upper1_mat[row_id, col_id] = ifelse(arg1[row_id, col_id] <= max(vgrid_vals[,col_id]), min(which(vgrid_vals[,col_id] >= as.numeric(arg1[row_id, col_id]))), NA)
      upper2_mat[row_id, col_id] = ifelse(arg2[row_id, col_id] <= max(wgrid_vals[,col_id]), min(which(wgrid_vals[,col_id] >= as.numeric(arg2[row_id, col_id]))), NA)
      
      if (is.na(upper1_mat[row_id, col_id]) == F){
        if (lower1_mat[row_id, col_id] > upper1_mat[row_id, col_id]){message(paste0("Stute estimator is non-monotone: p_index ", row_id, col_id))}
      }
      if (is.na(upper2_mat[row_id, col_id]) == F){
        if (lower2_mat[row_id, col_id] > upper2_mat[row_id, col_id]){message(paste0("Lin-Ying estimator is non-monotone: p_index ", row_id, col_id))}
      }
    }
    
    qt1_lower_mat[row_id, ] = -tt1 + pmax(tt1,ifelse(length(z1_grid[lower1_mat[row_id, ]]) != ncol(arg1),
    												 rep(NA, ncol(arg1)), z1_grid[lower1_mat[row_id, ]]))
    qt1_upper_mat[row_id, ] = -tt1 + pmax(tt1,ifelse(length(z1_grid[upper1_mat[row_id, ]]) != ncol(arg1),
    												 rep(NA, ncol(arg1)), z1_grid[upper1_mat[row_id, ]]))
    qt2_lower_mat[row_id, ] = -tt1 + pmax(tt1,ifelse(length(z1_grid[lower2_mat[row_id, ]]) != ncol(arg1),
    												 rep(NA, ncol(arg1)), z1_grid[lower2_mat[row_id, ]]))
    qt2_upper_mat[row_id, ] = -tt1 + pmax(tt1,ifelse(length(z1_grid[upper2_mat[row_id, ]]) != ncol(arg1),
    												 rep(NA, ncol(arg1)), z1_grid[upper2_mat[row_id, ]]))
  }
  
  return(list(pvals = p,
              arg1 = arg1, arg2 = arg2,
              vgrid_vals = vgrid_vals,
              wgrid_vals = wgrid_vals,
              qt_tilde_lower = qt1_lower_mat, qt_tilde_upper = qt1_upper_mat, 
              qt_tilde_midpoint = (qt1_lower_mat + qt1_upper_mat)/2,
              qt_dtilde_lower = qt2_lower_mat, qt_dtilde_upper = qt2_upper_mat, 
              qt_dtilde_midpoint = (qt2_lower_mat + qt2_upper_mat)/2,
              l1 = lower1_mat, u1 = upper1_mat, l2 = lower2_mat, u2 = upper2_mat))
}

Q_fun_stute <- function(p, tt1, tt21, tt22, z1, z2, delta1, delta2, method = "km",
                        pavit = "TRUE", add_constraint = "TRUE"){
  if (length(tt1) != length(tt21) | length(tt1) != length(tt22)){paste0("tt1, tt21 and tt22 should have the same dimensions")};
  
  z1_grid = c(0,sort(z1))
  
  if (length(unique(tt21)) == 1 & length(unique(tt22)) == 1){
    message("Unique limits t21 and t22 specified");
    
    if (method == "km"){
      v1 = Fn_tilde_stute2_vec(tt1, tt21[1], tt22[1], z1, z2, delta1, delta2)
      
      v_mat = matrix(rep(v1, length(p)), nrow = length(p), ncol = length(tt1), byrow = T)
      
      vz = Fn_tilde_stute2_vec(z1_grid, tt21[1], tt22[1], z1, z2, delta1, delta2)
    }
    
    if (method == "full"){
      v1 = Fn_tilde_stute_vec(tt1, tt21[1], tt22[1], z1, z2, delta1, delta2)
      
      v_mat = matrix(rep(v1, length(p)), nrow = length(p), ncol = length(tt1), byrow = T)
      
      vz = Fn_tilde_stute_vec(z1_grid, tt21[1], tt22[1], z1, z2, delta1, delta2)
    }
    
    # Additional constraint to ensure that Stute estimator jumps to one
    if (add_constraint == TRUE){
      vz[which(z1_grid == max(z1_grid))] = 1
    }
    
    vgrid_vals = matrix(rep(pmin(vz,1), length(tt1)), nrow = length(z1_grid), ncol = length(tt1), byrow = F)
    message("Preprocessing done");    
  }
  
  if (length(unique(tt21)) != 1 | length(unique(tt22)) != 1){
    message("Different limits t21 and t22 specified");
    
    v_mat = matrix(nrow = length(p), ncol = length(tt1))
    
    vgrid_vals = matrix(nrow = length(z1_grid), ncol = length(tt1))
    
    for (col_id in 1:length(tt1)){
      if (method == "full"){
        v = Fn_tilde_stute_vec(tt1[col_id], tt21[col_id], tt22[col_id], z1, z2, delta1, delta2)
        
        v_mat[, col_id] = rep(v, length(p))
        
        vgrid_vals[, col_id] = Fn_tilde_stute_vec(z1_grid, tt21[col_id], tt22[col_id], z1, z2, delta1, delta2)
      }
      if (method == "km"){
        v = Fn_tilde_stute2_vec(tt1[col_id], tt21[col_id], tt22[col_id], z1, z2, delta1, delta2)
        
        v_mat[, col_id] = rep(v, length(p))
        
        vgrid_vals[, col_id] = Fn_tilde_stute2_vec(z1_grid, tt21[col_id], tt22[col_id], z1, z2, delta1, delta2)
      }
    }
    message("Preprocessing done");    
  }
  
  # Create matrix of p-values for different tt1-values (with corresponding limits)  
  p_mat = matrix(p, nrow = length(p), ncol = length(tt1), byrow = F)
  
  # Arguments of Fn_tilde_inv: Dimension length(p) x length(tt1)
  arg1 = p_mat + (1-p_mat)*v_mat;
  
  # Solving the equation x* = Fn_tilde_inv(y)  
  lower1_mat = matrix(NA, nrow = nrow(arg1), ncol = ncol(arg1));
  upper1_mat = matrix(NA, nrow = nrow(arg1), ncol = ncol(arg1));
  
  qt1_lower_mat = matrix(NA, nrow = nrow(arg1), ncol = ncol(arg1));
  qt1_upper_mat = matrix(NA, nrow = nrow(arg1), ncol = ncol(arg1));
  
  for (row_id in 1:nrow(arg1)){
    for (col_id in 1:ncol(arg1)){
      #message(paste0("Computation: p", row_id, " from ", nrow(arg1), " t1-value ", col_id, " from ", ncol(arg1)));
      lower1_mat[row_id, col_id] = max(which(vgrid_vals[,col_id] < as.numeric(arg1[row_id, col_id])))
      upper1_mat[row_id, col_id] = ifelse(arg1[row_id, col_id] <= max(vgrid_vals[,col_id]), min(which(vgrid_vals[,col_id] >= as.numeric(arg1[row_id, col_id]))), NA)
      
      if (is.na(upper1_mat[row_id, col_id]) == F){
        if (lower1_mat[row_id, col_id] > upper1_mat[row_id, col_id]){message(paste0("Stute estimator is non-monotone: p_index ", row_id, col_id))}
      }
    }
    
    qt1_lower_mat[row_id, ] = -tt1 + pmax(tt1,ifelse(length(z1_grid[lower1_mat[row_id, ]]) != ncol(arg1),
                                                     rep(NA, ncol(arg1)), z1_grid[lower1_mat[row_id, ]]))
    qt1_upper_mat[row_id, ] = -tt1 + pmax(tt1,ifelse(length(z1_grid[upper1_mat[row_id, ]]) != ncol(arg1),
                                                     rep(NA, ncol(arg1)), z1_grid[upper1_mat[row_id, ]]))
  }
  
  return(list(pvals = p,
              arg1 = arg1,
              vgrid_vals = vgrid_vals,
              qt_tilde_lower = qt1_lower_mat, qt_tilde_upper = qt1_upper_mat, 
              qt_tilde_midpoint = (qt1_lower_mat + qt1_upper_mat)/2,
              l1 = lower1_mat, u1 = upper1_mat))
}