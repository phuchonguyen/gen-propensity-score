rm(list=ls())
library(tidyverse)
source("utils.R")

# Gaussian linear DRF
#   correctly specified treatment assignment model
set.seed(2)
n <- 2000
X <- rnorm(n, 0.5, 1)
T <- rnorm(n, X + X^2, 1)
Yt <- function(t, x) {
  rnorm(1, 2*x + t, 1)
}
df <- data.frame(Y=sapply(1:n, function(i) Yt(T[i], X[i])), X=X, T=T)
ts <- qnorm(seq(0.025, 0.975, length.out = 11), 0.5, sqrt(0.5))
y0 <- sapply(ts, function(t) (1+t))
d0 <- diff(y0, 1)/diff(ts, 1)

# Treatment assignment model for PF and GPS
ptx <- lm(T~poly(X, 2), data=df)
rtx <- function(t, x) {
  dnorm(t, mean = predict(ptx, list(X=x))[[1]], sd = sigma(ptx))
}
R <- sapply(1:n, function(i) rtx(T[i], X[i]))
df <- df %>%
  mutate(theta = ptx$fitted.values) %>%  # unique parametrization of PF
  mutate(R = R)

# Response model
## Naive model: correctly specified outcome model
est_naive <- function(df) {
  naivelm <- lm(Y~X*T, data=df)
  y_naive <- sapply(ts, function(t) naive(t=t, X=df$X, ymod=naivelm))
  dif_naive <- diff(y_naive, 1) / diff(ts, 1)
  return(list(y=y_naive, dif=dif_naive))
}

## Nonparametric/flexible response model for GPS
est_gps <- function(df, ptx, rtx) {
  eytr <- mgcv::gam(Y ~ s(T, R, k=27), data=df, method="REML", sp=0.001)  # TODO: tune sp and k
  y_gps <- sapply(ts, function(t) gps(t=t, rtx=rtx, X=df$X, eytr=eytr))
  dif_gps <- diff(y_gps, 1) / diff(ts, 1)
  return(list(y=y_gps, dif=dif_gps))
}

## IPW of GPS
est_ipw <- function(df, ptx, rtx, h=0.01) {
  y_ipw <- sapply(ts, function(t) ipw(t=t, X=df$X, rtx=rtx, Y=df$Y, T=df$T, h=h))
  dif_ipw <- diff(y_ipw, 1) / diff(ts, 1)
  return(list(y=y_ipw, dif=dif_ipw))
}

## Smooth coefficient model for PF
est_pf <- function(df, ptx) {
  eytt <- mgcv::gam(Y ~ s(T, theta, k=6), data=df, method="REML")  # TODO: tune sp and k
  y_pf <- sapply(ts, function(t) pf(t=t, eytt=eytt, theta=df$theta))
  dif_pf <- diff(y_pf, 1) / diff(ts, 1)
  return(list(y=y_pf, dif=dif_pf))
}

# Estimation
## Naive model
temp <- est_naive(df)
y0_naive <- temp$y
d0_naive <- temp$dif

## Nonparametric/flexible response model for GPS
temp <- est_gps(df, ptx, rtx)
y0_gps <- temp$y
d0_gps <- temp$dif

## IPW of GPS
temp <- est_ipw(df, ptx, rtx, 1/sqrt(n))
y0_ipw <- temp$y
d0_ipw <- temp$dif

## Smooth coefficient model for PF
temp <- est_pf(df, ptx)
y0_pf <- temp$y
d0_pf <- temp$dif

save(y0_naive, y0_gps, y0_ipw, y0_pf, d0_naive, d0_pf, d0_gps, d0_ipw, 
     ts, y0, d0,
     file = "data/sim2-observed.RData")

# Boostrap Variance
s <- 50
Y_naive <- Y_gps <- Y_ipw <- Y_pf <- matrix(NA, s, length(ts))
D_naive <- D_gps <- D_ipw <- D_pf <- matrix(NA, s, length(ts)-1)
for (i in 1:s) {
  idx <- sample(1:n, size=n, replace=TRUE)
  df_cur <- df[idx, ]
  
  # Treatment assignment model for PF and GPS
  ptx_cur <- lm(T~poly(X, 2), data=df_cur)
  rtx_cur <- function(t, x) {
    dnorm(t, mean = predict(ptx_cur, list(X=x))[[1]], sd = sigma(ptx_cur))
  }
  R <- sapply(1:n, function(i) rtx_cur(df_cur$T[i], df_cur$X[i]))
  df_cur <- df_cur %>%
    mutate(theta = ptx_cur$fitted.values) %>%  # unique parametrization of PF
    mutate(R = R)
  
  # Response model
  ## Naive model
  temp <- est_naive(df_cur)
  Y_naive[i, ] <- temp$y
  D_naive[i, ] <- temp$dif
  
  ## Nonparametric/flexible response model for GPS
  temp <- est_gps(df_cur, ptx_cur, rtx_cur)
  Y_gps[i, ] <- temp$y
  D_gps[i, ] <- temp$dif
  
  ## IPW of GPS
  temp <- est_ipw(df_cur, ptx_cur, rtx_cur, h=1/sqrt(n))
  Y_ipw[i, ] <- temp$y
  D_ipw[i, ] <- temp$dif
  
  ## Smooth coefficient model for PF
  temp <- est_pf(df_cur, ptx_cur)
  Y_pf[i, ] <- temp$y
  D_pf[i, ] <- temp$dif
}
save(Y_naive, Y_ipw, Y_gps, Y_pf, D_naive, D_gps, D_ipw, D_pf, 
     file="data/sim2-bootstrap.RData" )
