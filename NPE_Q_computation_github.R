######################################################################
#  Nonparametric estimation of quantiles of the conditional residual # 
#  lifetime distribution                                             #
#                                                                    #
#  Version: GitHub                                                   #
#  Author: Steven Abrams                                             #
#  Last update: July, 2025                                           #
######################################################################

rm(list=ls())

###############################
# Loading important functions #
###############################

source("./NPE_Q_estimation_github.R")

##############################################
# Application of estimation of the quantiles #
##############################################

# Simulated data example 
#------------------------
set.seed(1207)
random_samples <- rCopula(copula=claytonCopula(0.5, dim = 2), n=500)
u = random_samples[,1]; v = random_samples[,2];
t1 = (-log(1-u)/0.5)**(1/1.5); 
t2 = (-log(1-v)/0.5)**(1/1.5);    

w = runif(500, 0, 1)
ct = (-log(w)/(0.15))**(1/1.5);   
delta1 = as.numeric(t1 <= ct);
z1 = pmin(t1, ct)

# One-component censoring
#-------------------------
Qres = Q_tilde_oc(p = seq(0.1,0.9,0.1), 
                  tt1 = 0.5, 
                  t21 = quantile(t2)[2], 
                  t22 = quantile(t2)[3], 
                  z1 = z1, 
                  t2 = t2, 
                  delta1 = delta1, 
                  method = "km")

par(mfrow = c(1,2))
plot(seq(0.1,0.9,0.1), Qres$qt_tilde_upper, type = "s", lwd = 2, col = 1, xlab = "p",
     ylab = "Estimated quantiles", ylim = c(0,2.5))
lines(seq(0.1,0.9,0.1), Qres$qt_tilde_upper_ll, col = 1, lwd = 2, type = "s", lty = 2)
lines(seq(0.1,0.9,0.1), Qres$qt_tilde_upper_ul, col = 1, lwd = 2, type = "s", lty = 2)

plot(seq(0.1,0.9,0.1), Qres$qt_tilde_lower, type = "s", lwd = 2, col = 1, xlab = "p",
     ylab = "Estimated quantiles", ylim = c(0,2.5))
lines(seq(0.1,0.9,0.1), Qres$qt_tilde_lower_ll, col = 1, lwd = 2, type = "s", lty = 2)
lines(seq(0.1,0.9,0.1), Qres$qt_tilde_lower_ul, col = 1, lwd = 2, type = "s", lty = 2)
