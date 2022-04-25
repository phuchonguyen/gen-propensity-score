rm(list=ls())
library(tidyverse)

# Simulation 1: Zero effect
#   correctly specified treatment assignment and flexible response models
set.seed(123)
n <- 2000
X <- rnorm(n, 0.5, sqrt(0.25))
T <- rnorm(n, X, sqrt(0.25))
Y <- function(t) {
  rnorm(n, 10*X, 1)
}
df <- data.frame(Y=Y(1), X=X, T=T)
ts <- qnorm(seq(0.05, 0.95, length.out = 11), 0.5, sqrt(0.5))
y0 <- sapply(ts, function(t) 10*0.5)
d0 <- diff(y0, 1)/diff(ts, 1)

# Treatment assignment model for PF and GPS
ptx <- lm(T~X, data=df)
rtx <- function(t, X) {
  dnorm(t, mean = predict(ptx, list(X=X))[[1]], sd = sigma(ptx))
}

# Response model
## Naive model
est_naive <- function(df) {
  naivelm <- lm(Y~X*T, data=df)
  naive <- function(t) {
    yt <- predict(naivelm, list(X=df$X, T=rep(t, length(df$X))))
    return(mean(yt))
  }
  #dnaive <- derivative(naive, ts)
  y_naive <- sapply(ts, naive)
  dif_naive <- diff(y_naive, 1)/ diff(ts, 1)
  return(list(y=y_naive, dif=dif_naive))
}

## Nonparametric/flexible response model for GPS
est_gps <- function(df, ptx, rtx) {
  df <- df %>%
    mutate(R = rtx(T, X))
  eytr <- mgcv::gam(Y ~ s(T, R), data=df, method="REML")
  gps <- function(t) {
    r <- sapply(df$X, rtx, t=t)
    yt <- predict(eytr, list(R=r, T=rep(t, length(df$X))))
    return(mean(yt))
  }
  y_gps <- sapply(ts, gps)
  dif_gps <- diff(y_gps, 1)/ diff(ts, 1)
  return(list(y=y_gps, dif=dif_gps))
}

## IPW of GPS
est_ipw <- function(df, ptx, rtx, h=0.01) {
  # ipw <- function(t) {
  #   r <- sapply(df$X, rtx, t=t)
  #   w <- dnorm((df$X - t) / h) / r
  #   w <- w / sum(w)
  #   yt <- sum(df$Y * w)
  #   return(yt)
  # }
  Sj <- function(t, j, w) {
    sum(w * (df$T - t)^j)
  }
  Dj <- function(t, j, w) {
    sum(df$Y * (df$T - t)^j * w)
  }
  ipw <- function(t) {
    r <- sapply(df$X, rtx, t=t)
    w <- dnorm((df$T - t) / h) / r
    res <- (Dj(t, 0, w) * Sj(t, 2, w) - Dj(t, 1, w) * Sj(t, 1, w)) / (Sj(t, 0, w) * Sj(t, 2, w) - Sj(t, 1, w)^2)
    return(res)
  }
  y_ipw <- sapply(ts, ipw)
  dif_ipw <- diff(y_ipw, 1)/ diff(ts, 1)
  return(list(y=y_ipw, dif=dif_ipw))
}

## Smooth coefficient model for PF
est_pf <- function(df, ptx) {
  df <- df %>%
    mutate(theta = ptx$fitted.values)  # unique parametrization of PF
  eytt <- mgcv::gam(Y ~ s(T, theta), data=df, method="REML")
  pf <- function(t) {
    yt <- predict(eytt, list(theta=df$theta, T=rep(t, length(X))))
    return(mean(yt))
  }
  #dpf <- derivative(pf, ts)
  y_pf <- sapply(ts, pf)
  dif_pf <- diff(y_pf, 1)/ diff(ts, 1)
  return(list(y=y_pf, dif=dif_pf))
}

# Response model
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

save(y0_naive, y0_gps, y0_ipw, y0_pf, d0_naive, d0_pf, d0_gps, d0_ipw, ts, d0,
     file = "data/sim1-observed.RData")

## Incrementatl causal effect
## Interpretation: what is the change in outcome if all observed treatment increases
##   by a small amount
#### Using correct model
#ice_cor <- coef(lm(Y~X+T, data = df))[['T']]
#### Using a flexible model
 
# Boostrap Variance
s <- 50
Y_naive <- Y_gps <- Y_ipw <- Y_pf <- matrix(NA, s, length(ts))
D_naive <- D_gps <- D_ipw <- D_pf <- matrix(NA, s, length(ts)-1)
for (i in 1:s) {
  idx <- sample(1:n, size=n, replace=TRUE)
  df_cur <- df[idx, ]
  
  # Treatment assignment model for PF and GPS
  ptx_cur <- lm(T~X, data=df_cur)
  rtx_cur <- function(t, X) {
    dnorm(t, mean = predict(ptx_cur, list(X=X))[[1]], sd = sigma(ptx_cur))
  }
  
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
  temp <- est_ipw(df_cur, ptx_cur, rtx_cur)
  Y_ipw[i, ] <- temp$y
  D_ipw[i, ] <- temp$dif
  
  ## Smooth coefficient model for PF
  temp <- est_pf(df_cur, ptx_cur)
  Y_pf[i, ] <- temp$y
  D_pf[i, ] <- temp$dif
}
save(Y_naive, Y_ipw, Y_gps, Y_pf, D_naive, D_gps, D_ipw, D_pf, 
     file="data/sim1-bootstrap.RData" )
