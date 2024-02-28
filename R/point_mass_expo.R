
#from ebnm for point mass expo
wpost_exp <- function(x, s, w, a) {
  if (w == 0) {
    return(rep(0, length(x)))
  }

  if (w == 1) {
    return(rep(1, length(x)))
  }

  lf <- dnorm(x, 0, s, log = TRUE)
  lg <- log(a) + s^2 * a^2 / 2 - a * x + pnorm(x / s - s * a, log.p = TRUE)
  wpost <- w / (w + (1 - w) * exp(lf - lg))

  return(wpost)
}





pe_summres_untransformed <- function(x, s, w, a, mu, output) {
  x <- x - mu

  wpost <- wpost_exp(x, s, w, a)

  post <- list()

  if (result_in_output(output)) {
    post$mean  <- wpost * my_etruncnorm(0, Inf, x - s^2 * a, s)
    post$mean2 <- wpost * my_e2truncnorm(0, Inf, x - s^2 * a, s)
    post$mean2 <- pmax(post$mean2, post$mean^2)

    if (any(is.infinite(s))) {
      post$mean[is.infinite(s)]  <- w / a
      post$mean2[is.infinite(s)] <- 2 * w / a^2
    }
    post$sd <- sqrt(pmax(0, post$mean2 - post$mean^2))

    post$mean2 <- post$mean2 + mu^2 + 2 * mu * post$mean
    post$mean  <- post$mean + mu
  }

  if ("lfsr" %in% output) {
    post$lfsr <- 1 - wpost
    if (any(is.infinite(s))) {
      post$lfsr[is.infinite(s)] <- 1 - w
    }
  }

  return(post)
}
