


init.mococomo <- function (data, max_class, mult = 2,upper=FALSE, nullweight, ...)
  UseMethod("init.mococomo")

#' @rdname init.mococomo
#'
#'
#' @method init.mococomo normal
#'
#' @export init.mococomo.normal
#'
#' @export
#'
init.mococomo.normal <- function(data, max_class, mult = 2,upper=FALSE,nullweight,...) {
  if(missing( nullweight)){
    nullweight <- 2.3
  }
  data$X2 <- data$X^2

  # TODO check input data has X, Z, betahat, se, etc.
  #### We can force the user to input data as being from a certain  class generarted by a set_data  function
  # w
  # scale mixture of normals w/ zero mean-- TODO: make this initializable


  #WIlliam to check here
  scales <- autoselect.mococomo(data)
  # scales <- cumprod(c(1., rep(sqrt(2), 5)))
  f_list <- purrr::map(scales, ~ normal_component(mu = 0, var = .x^2))

  K <- length(scales)
  p <- ncol(data$X)
  n <- nrow(data$X)

  # initialize posterior assignment
  Y <- matrix(rep(1, K * n) / K, nrow = n)
  N <- .expected_trials(Y)

  data$y <- Y
  data$N <- N
  if(is.null(data$Z))
  { data$Z <- matrix(1, nrow= nrow(data$X), ncol=1) }

  # it's multinational regression + some other stuff
  fit                 <- init.mnsusie(data, from_init_moco=TRUE)
  fit$f_list          <- f_list
  fit$post_assignment <- Y
  fit$N               <- N
  fit$model            <-"normal"
  fit$K               <- length(fit$f_list)
  fit$data            <- data
  fit$nullweight      <- nullweight
  class(fit) <- "mococomo_normal"
  # only need to do this once when component probabilities are fixed
  fit$data_loglik <- compute_data_loglikelihood(fit,data)

  return(fit)
}


#' @rdname init.mococomo
#'
#'
#' @method init.mococomo beta
#'
#' @export init.mococomo.beta
#'
#' @export
#'
init.mococomo.beta <- function(data, max_class, mult = 2,upper=FALSE,nullweight,...) {

  if(missing( nullweight)){
    nullweight <- 10
  }
  # different set up depending on fitting p values or fitting z-scores
  #need to polish the upper limit for alpha and beta, now set as 100
  if(upper){ #mixture to fit component close to 1
    alpha   <- seq(1.1,100,length.out=max_class)
    upper_l <- list( alpha=alpha, beta=1)
  }else{
    alpha   = 1
    upper_l = NULL
  }
  #mixture to fit component close to 0
  beta    <- c(1,seq(1.1,100,length.out=max_class)) #contains null component!!!
  lower_l <-list( alpha=1, beta=beta)
  f_list  <- beta_component(alpha = alpha,
                            beta  = beta
  )
  K      <- length(upper_l[[1]])+length(lower_l[[2]])
  p      <- ncol(data$X)
  n      <- nrow(data$X)

  # initialize posterior assignment
  Y <- matrix(rep(1, K * n) / K, nrow = n)
  N <- .expected_trials(Y)
  data$y <- Y
  data$N <- N

  data$X2 <- data$X^2
  if(is.null(data$Z))
  { data$Z <- matrix(1, nrow= nrow(data$X), ncol=1) }

  # it's multinational regression + some other stuff
  fit                 <- init.mnsusie(data, from_init_moco=TRUE)
  fit$f_list          <- f_list
  fit$upper_l         <- upper_l
  fit$lower_l         <- lower_l
  fit$post_assignment <- Y
  fit$N               <- N
  fit$model            <- "beta"
  fit$K               <- K
  fit$data            <- data
  fit$nullweight      <- nullweight
  class(fit)          <- "mococomo_beta"
  # only need to do this once when component probabilities are fixed
  fit$data_loglik     <- compute_data_loglikelihood(fit,data)

  return(fit)

}
#' Compute Expected Assignment Log Likelihood
#' @return N x K matrix each row has E[p(z=k | \beta, \omega)] for k = 1,..., K
compute_assignment_loglikelihood.mococomo <- function(fit, normalize = TRUE) {
  K <- fit$K
  Xb <- do.call(cbind, lapply( 1: length(fit$logreg_list),
                               function (k)  compute_Xb.binsusie(fit,k=k)
  )
  )

  Xbcum <- do.call(rbind, apply(Xb, 1, cumsum, simplify = F))
  kln2 <- seq(K - 1) * log(2)

  # second moments only appear in normalizing constant
  C <- 0
  if (normalize) {
    Xb2 <- do.call(cbind, lapply( 1: length(fit$logreg_list),
                                  function (k)  compute_Xb2.binsusie(fit,k=k)
    )
    )  # N x K - 1
    omega <- do.call(cbind, lapply( 1: length(fit$logreg_list),
                                    function (k)  compute_omega(fit ,
                                                                k =k)
    )
    )  # N x K - 1
    C <- -0.5 * rowSums(Xb2 * omega)
  }

  logp <- Xb - 0.5 * Xbcum + C
  logp <- do.call(rbind, apply(logp, 1, function(x) x - kln2, simplify = F))

  logpK <- -(K - 1) * log(2) - 0.5 * Xbcum[, K - 1] + C
  logp <- cbind(logp, logpK)
  return(logp)
}



compute_assignment_jj_bound.mococomo <- function(fit) {
  K <- fit$K
  Xb <- do.call(cbind,
                lapply(1:length(fit$logreg_list),
                       function( k) compute_Xb.binsusie(fit,k)
                )
  ) # N x K-1

  Xi <- do.call(cbind, purrr::map(fit$logreg_list, ~ purrr::pluck(.x, "params", "xi"))) # N x K-1

  f <- function(xi, xb) {
    tmp <- cumsum(log(sigmoid(xi)) - 0.5 * xi - 0.5 * xb) + xb
    tmpK <- sum(log(sigmoid(xi)) - 0.5 * xi - 0.5 * xb )
    jj <- c(tmp, tmpK)
    return(jj)
  }
  jj <- do.call(rbind, purrr::map(seq(nrow(Xb)), ~ f(Xi[.x, ], Xb[.x, ])))

  # JJ bound, should be faster with matrix operations?
  # tmp <- rowCumSum(log(sigmoid(xi)) - 0.5 * xi - 0.5 * Xb)
  # jjK <- tmp[, K-1]
  # jj <- cbind(Xb + tmp, jjK)
  return(jj)
}

#' Compute Posterior Assignment Probabilities
#' For each data point return posterior assignment probabilities
#' @param fit a MoCoCoMo fit object
#' @return an n x K matrix of log posterior probabilities for each data point
compute_posterior_assignment <- function(fit, log = TRUE) {
  data_loglik <- fit$data_loglik
  assignment_loglik <- compute_assignment_jj_bound.mococomo(fit)  # n x L
  assignment_loglik[, 1] <- assignment_loglik[, 1] + fit$nullweight

  # normalize
  res <- do.call(
    rbind,
    apply(data_loglik + assignment_loglik, 1, function(x) x - logSumExp(x), simplify = F)
  )
  if (!log) {
    res <- exp(res)
  } # exponentiate if log=FALSE
  return(res)
}


compute_elbo3.mococomo <- function(fit) {
  post_assignment <- fit$post_assignment
  data_loglik <- fit$data_loglik
  jj <- compute_assignment_jj_bound.mococomo(fit)

  kl_susie <- Reduce("+", lapply(1:length(fit$logreg_list),
                                 function(k) compute_kl.binsusie(fit,k)
  )
  )

  assignment_entropy <- sum(apply(post_assignment, 1, categorical_entropy))

  elbo <- sum(post_assignment * (jj + data_loglik)) + assignment_entropy - kl_susie
  return(elbo)
}


compute_elbo2.mococomo<-  function(fit) {
  post_assignment <- fit$post_assignment
  assignment_loglik <- compute_assignment_loglikelihood.mococomo(fit)
  ll <- sum(post_assignment * (fit$data_loglik + assignment_loglik))
  assignment_entropy <- sum(apply(post_assignment, 1, categorical_entropy))

  kl_susie <- Reduce("+", lapply( 1:length(fit$logreg_list),
                                  function(k) compute_kl.binsusie(fit,k)))  # sum(purrr::map_dbl(fit$logreg_list, compute_kl.binsusie))
  kl_omega <-  Reduce("+", lapply( 1:length(fit$logreg_list),
                                   function(k) sum(pg_kl(fit$data$N[,k],
                                                         fit$logreg_list[[k]]$params$xi)#sum(purrr::map_dbl(fit$logreg_list, function(x) sum(pg_kl(x$data$N, x$params$xi))))
                                   )
  )
  )


  kl <- kl_susie + kl_omega

  return(ll + assignment_entropy - kl)

}

compute_elbo.mococomo <- function(fit) {


  K <- fit$K
  N <- .expected_trials(fit$post_assignment)

  # MAKE SURE susie has access to the right posteriors!
  # and update omega so jj bound is tight
  for (k in seq(K - 1)) {
    fit$logreg_list[[k]]$data$y <- fit$post_assignment[, k]
    fit$logreg_list[[k]]$data$N <- N[, k]
    fit$logreg_list[[k]]$params$xi <- update_xi.binsusie(fit ,k=k)
    fit$logreg_list[[k]]$params$tau <- compute_tau(fit ,

                                                   k =k)
  }


  # TODO: do we want to just store post_assignment after each update?
  post_assignment <- fit$post_assignment
  data_loglik <- fit$data_loglik

  ll <- sum(post_assignment * data_loglik)
  assignment_entropy <- sum(apply(post_assignment, 1, categorical_entropy))
  mnsusie_elbo <- compute_elbo.mnsusie(fit, from_init_moco = TRUE)
  elbo <- ll + assignment_entropy + mnsusie_elbo

  return(elbo)
}


#' Expected Binomial Trials
#' Compute expected  number of binomial trials
#' MoCoCoMo fits K-1 binomial regressions, this is the E-step for N
#' @param Y an N x K matrix of assignment probabilities
.expected_trials <- function(Y) {
  cumY <- do.call(rbind, apply(Y, 1, cumsum, simplify = F))
  N <- 1. - cumY + Y #TODO add penalty here
}


iter.mococomo <- function(fit, update_assignment = T, update_logreg = T) {
  K <- fit$K

  # updates posterior assignments
  if (update_assignment) {
    fit$post_assignment <- compute_posterior_assignment(fit, log = F) #TODO add penalty
    fit$N <- .expected_trials(fit$post_assignment)
  }

  # pass assignments to logreg
  for (k in seq(K - 1)) {
    fit$logreg_list[[k]]$data$y <- fit$post_assignment[, k]
    fit$logreg_list[[k]]$data$N <- fit$N[, k]
  }

  if (update_logreg) {
    fit <- iter.mnsusie(fit, from_init_moco = TRUE)
  }

  return(fit)
}





#' @title Compute individual posterior variance from marginal normal mean model
#' @description internal function to compute posterior mean and sds


t_ind_var.mococomo <- function(fit, i) {
  do.call(
    c,
    lapply(
      1:length(fit$f_list),
      function(k) {
        1 / ((1 / fit$data$se[i]^2) + (1 / fit$f_list[[k]]$var))
      }
    )
  )
}


#' @title Compute individual posterior first and second moment
#' @description Compute individual posterior first and second moment
#'
#' # TODO: currently only for center prior
#' @param fit a mococomo object
#' @param  t_ind_var output of t_ind_var (using same mocomo object!)
#' @param i individual of interest
#' @exemple
#' t_post_var <-   do.call(rbind,
#'                        lapply( 1:length(fit$data$betahat),
#'                                function(i)t_ind_var.mococomo(fit, i)
#'                        )
#' )
#'
#'
#' post_beta <-     do.call(c, lapply( 1: length(fit$data$betahat), function(i)cal_ind_postmean(fit, t_post_var,i,) ))

cal_ind_moment12 <- function(fit, t_post_var, i) {
  temp <- do.call(
    c,
    lapply(
      1:length(fit$f_list),
      function(k) {
        (t_post_var[i, k] / (fit$data$se[i]^2) )*
          (fit$data$betahat[i])
      }
    )
  )

  ind_mean <- sum(fit$post_assignment[i, ] * temp)
  ind_second_moment <- sum(fit$post_assignment[i, ] * (t_post_var[i, ] + temp^2))

  ind_var <- ind_second_moment - ind_mean^2

  return(list(
    mean = ind_mean,
    sd = sqrt(ind_var)
  ))
}





#' @title Compute individual posterior mean and sd under mococomo model
#' @description Compute individual posterior mean and sd under mococomo model
#'
#' # TODO: allow using new observation from another data set (e.g. testing)
#' @param fit a mococomo object
#' @export
#' @example
#' see \link{\code{fit.mococomo}}

post_mean_sd.mococomo <- function(fit) {
  t_post_var <- do.call(
    rbind,
    lapply(
      1:length(fit$data$betahat),
      function(i) t_ind_var.mococomo(fit, i)
    )
  )
  out <- do.call(
    rbind,
    lapply(
      1:length(fit$data$betahat),
      function(i) cal_ind_moment12(fit, t_post_var, i)
    )
  )
  out <- data.frame(
    mean = do.call(c, out[, 1]),
    sd = do.call(c, out[, 2])
  ) # could be skip for speed



  return(out)
}





# adapted from autoselect.mxisqp
# try to select a default range for the sigmaa values
# that should be used, based on the values of betahat and sebetahat
# mode is the location about which inference is going to be centered
# gridmult is the multiplier by which the sds differ across the grid
#'@export
autoselect.mococomo <- function(data, max_class, mult = 2) {
  betahat <- data$betahat
  sebetahat <- data$se


  sigmaamin <- min(sebetahat) / 10 # so that the minimum is small compared with measurement precision
  if (all(betahat^2 <= sebetahat^2)) {
    sigmaamax <- 8 * sigmaamin # to deal with the occassional odd case where this could happen; 8 is arbitrary
  } else {
    sigmaamax <- 2 * sqrt(max(betahat^2 - sebetahat^2)) # this computes a rough largest value you'd want to use, based on idea that sigmaamax^2 + sebetahat^2 should be at least betahat^2
  }

  if (mult == 0) {
    return(c(0, sigmaamax / 2))
  } else {
    npoint <- ceiling(log2(sigmaamax / sigmaamin) / log2(mult))
    out <- c(0,mult^((-npoint):0) * sigmaamax)#  mult^((-npoint):0) * sigmaamax
    if (!missing(max_class)) {
      if (length(out) > max_class) {
        out <- seq(min(out), max(out), length.out = max_class)
      }
    }
    return(out)
  }
}



get_KL.mococomo <- function(fit) {
  kl_susie <- sum(purrr::map_dbl(fit$logreg_list, compute_kl.binsusie))
  kl_omega <- sum(purrr::map_dbl(fit$logreg_list, function(x) sum(pg_kl(x$data$N, x$params$xi))))
  kl <- kl_susie + kl_omega
  return(kl)
}


#' @title Compute individual fdr value  mococomo model with centered normal mixture
#' @descriptionCompute individual fdr value  mococomo model with centered normal mixture
#'
#' @param fit a mococomo object
#' @export
get_fdr <- function(fit) {
  tt1 <- fit$post_assignment[, 1] * dnorm(fit$data$betahat, mean = 0, sd = fit$data$se)
  tt2 <- Reduce("+", lapply(
    2:ncol(fit$post_assignment),
    function(k) {
      fit$post_assignment[, k] * dnorm(fit$data$betahat,
                                       mean = 0,
                                       sd = sqrt(fit$data$se^2 + fit$f_list[[k]]$var)
      )
    }
  ))
  out <- tt1 / (tt1 + tt2)

  return(out)
}
