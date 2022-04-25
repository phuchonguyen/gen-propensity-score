# E(Y(t)) based on outcome modeling
naive <- function(t, X, ymod) {
  if (!is.matrix(X)) X <- as.matrix(X, ncol=1)
  yt <- predict(ymod, data.frame(X, T=rep(t, nrow(X))))
  return(mean(yt))
}

# E(Y(t)) based on SCM(GPS)
gps <- function(t, rtx, X, eytr) {
  if (!is.matrix(X)) X <- as.matrix(X, ncol=1)
  r <- apply(X, 1, function(x) rtx(t=t, x=x))
  yt <- predict(eytr, data.frame(R=r, T=rep(t, nrow(X))))
  return(mean(yt))
}

# E(Y(t)) based on kernel regression weighted by GPS IPW
ipw_naive <- function(t, X, rtx, Y, T, h=NULL) {
  if (!is.matrix(X)) X <- as.matrix(X, ncol=1)
  if (is.null(h)) h <- 1/sqrt(length(X))
  r <- apply(X, 1, function(x) rtx(t=t, x=x))
  w <- dnorm((T - t) / h) / r
  w <- w / sum(w)
  yt <- sum(Y * w)
  return(yt)
}

ipw <- function(t, X, rtx, Y, T, h=NULL) {
  if (!is.matrix(X)) X <- as.matrix(X, ncol=1)
  if (is.null(h)) h <- 1/sqrt(length(X))
  r <- apply(X, 1, function(x) rtx(t=t, x=x))
  w <- dnorm((T - t) / h) / r
  Sj <- function(j) {
    sum(w * (T - t)^j)
  }
  Dj <- function(j) {
    sum(Y * (T - t)^j * w)
  }
  res <- (Dj(0) * Sj(2) - Dj(1) * Sj(1)) / (Sj(0) * Sj(2) - Sj(1)^2)
  return(res)
}

# E(Y(t)) based on SCM(PF)
pf <- function(t, eytt, theta) {
  yt <- predict(eytt, data.frame(theta=theta, T=rep(t, length(theta))))
  return(mean(yt))
}

# get derivatives of m at ts
derivative <- function(m, ts) {
  d <- rep(NA, length(ts))
  d[1] <- (m(ts[2]) - m(ts[1])) / (ts[2]-ts[1])
  for (i in 2:length(ts)) {
    if (i == length(ts)) {
      d[i] <- (m(ts[i]) - m(ts[i-1])) / (ts[i]-ts[i-1])
    } else {
      d[i] <- 0.5*((m(ts[i+1]) - m(ts[i])) / (ts[i+1]-ts[i]) + (m(ts[i]) - m(ts[i-1])) / (ts[i]-ts[i-1]))
    }
  }
  return(d)
}