#' @param q q(mu, var, alpha)
#' @param R LD matrix-- assumes diag(R) = rep(1, p)
#' @param tau0 prior effect variance
#' @param prior_logit p-vector with prior log odds for gamma = 1
rssvb <- function(zhat, q, R, tau0, prior_logit){
  # unpack
  mu <- q$mu
  var <- q$var
  alpha <- q$alpha

  p <- length(zhat)
  psi <- (R %*% (mu * alpha))[,1] # prediction
  for(i in 1:p){
    # remove effect of this variable
    psi <- psi - R[i,] * (mu[i]*alpha[i])

    # compute q(beta | gamma = 1)
    nu <- zhat[i] - psi[i]
    tau <- 1 + tau0
    mu[i] <- nu/tau
    var[i] <- 1/tau

    # logit <- zhat[i] * mu[i]
    #   - 0.5 * (psi[i] * mu[i] +  mu[i]^2 + var[i])
    #   -0.5 * tau0 * (mu[i]^2 + var[i]) + prior_logit[i]
    logit <- 0.5 * (mu[i]^2/var[i] + log(var[i]) + log(tau0)) + prior_logit[i]
    alpha[i] <- 1/(1 + exp(-logit))

    alpha[i]
    psi <- psi + R[i,] * (mu[i]*alpha[i])
  }
  return(list(mu=mu, var=var, alpha=alpha))
}

zhat <- rnorm(100)
zhat[10] <- zhat[10] + 5

p <- 50
X <- logisticsusie:::sim_X(p = p, length_scale = 5)
R <- cor(X)
z <- rep(0, p)
z[10] <- 5
zhat <- (R %*% z)[,1] + mvtnorm::rmvnorm(1, sigma=R)[1,]

q = list(
  mu = rep(0, p),
  var = rep(1, p),
  alpha = rep(1/p, p)
)
tau0 = 0.1
prior_logit <- rep(-3, p)

q <- rssvb(zhat, q, R, tau0, prior_logit)
