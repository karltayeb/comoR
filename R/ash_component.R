# Component distribution for ASH

# ASH Component ----
ash_null <- function(){
  list(fitted_g=NULL)
}

#' Ash component distribution
#'
#' Initialize an ash component distribution
#' when the distribution is updated, it initialized the new ash fit with the estimated `g` from the last iterations
#' this has the effect of only updating the mixture weights, but no the mixture distribution.
#'
#' By default we use a normal mixture distribution for g (note: no point mass at 0)
#'
#' @param ash_args a list of params that get passed to `ashr::ash`
#' @export
ash_component <- function(ash_args = list()){
  # default arguments should be reasonable
  # no point mass because it is handled by the null componenent
  default_ash_args <- list(
    pointmass = F,
    optmethod ='mixSQP',
    mixcompdist = 'normal',
    outputlevel = 2
  )
  ash_args <- modifyList(default_ash_args, ash_args)
  f <- list(ash = ash_null(), ash_args = ash_args)
  class(f) <- c("ash_component", "component_distribution")
  return(f)
}

ll_obs_normalmix <- function(betahat, sebetahat, sds, pi){
  matrixStats::logSumExp(dnorm(betahat, 0, sqrt(sebetahat^2 + sds^2), log=T) + log(pi))
}

#' @export
convolved_logpdf.ash_component <- function(dist, betahat, se){
  # note: this assumes mixcompdist = normal
  # if you use outputlevel >= 3 in ash you get normalized likelihoods which could be more general
  # but you still need a consistent way to un-normalize back to log-likelihoods
  if(!is.null(dist$ash$fitted_g)){
    if(!class(dist$ash$fitted_g == 'normalmix')){
      stop('convolved_logpdf.ash_component only implemented for normalmix')
    }
    sd <- dist$ash$fitted_g$sd
    pi <- dist$ash$fitted_g$pi
    ll <- purrr::map2_dbl(betahat, sebetahat, ~ll_obs_normalmix(.x, .y, sd, pi))
  } else{
    # case: we haven't called ash yet
    # return loglik under null -> weights = 1/2, this will sort out next iteration...
    ll <- dnorm(betahat, sd = sd, log=T)
  }
  return(ll)
}

#' @export
update_params.ash_component <- function(dist, betahat, se, weights) {
  args <- c(
    list(betahat = betahat, sebetahat = se, weights=weights),
    dist$ash_args
  )
  # initialize with g = fitted g from last iteration--
  # keeps grid of mixcompdist the same across iterations
  # but lets us reestimate pi at each iteration
  args$g <- dist$ash$fitted_g
  dist$ash <- rlang::exec(ashr::ash, !!!args)
  return(dist)
}
