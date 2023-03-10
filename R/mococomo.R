# Implements covariate moderated ASH "MOre COmponents COvariate MOderated"

#' @title Function implementation the mococomo mode
#' @details Function implementation the mococomo mode
#'
#' @param data an object of class data_mococomo  see \link{\code{set_data_mococomo}}
#' @param modeltype of model currently supported (normal and beta )
#' @param maxiter numeric, maximum numerous of iteration set to 100 by defaults
#' @param tol tolerance in term of change in ELBO value for stopping criterion
#' @param upper, logical, set to FALSE by default. Specific to beta distribution.
#'  If true use a to set of mixture for fitting both end of the of the distribution as in the ZAP paper by Leung and Sunn
#' @parma nullweight  numeric value for penalizing likelihood at point mass 0/null component (should be larger than  1, 1 corresponds to no penalty , 2 corresponds to considering 1 individual "being null" and so on)
#' (usefull in small sample size)
#' @export
#' @example
#' #Simulate data under the mococomo model
#' sim  <- sim_twococomo()
#' #preparing the data
#' data <- set_data_mococomo(betahat = sim$betahat,
#'                                se = sim$se ,
#'                                 X = sim$X)
#' #fit mococomo model
#' fit <- fit.mococomo(data, maxiter=20)
#' plot(fit$elbo)
#' .monotone(fit$elbo)
#'
#' #get posterior quantities
#' est<- post_mean_sd.mococomo (fit)
#' head(est)
#'  plot( est$mean, data$betahat)
#'
#' #comparison with ash
#'
#' t_ash <- ash(sim $betahat, sim $se, mixcompdist = "normal")
#' post_mean_ash <- t_ash$result$PosteriorMean
#' plot(est$mean, post_mean_ash)
#' # TODO make a more convincing example
#'
#'  sim  <- logisticsusie:::sim_mococomo_beta(n=100)
#'#preparing the data
#'data <- set_data_mococomo(p = sim$p,
#'                          X = sim$X)
#'
#' fit <- fit.mococomo(data, maxiter=20)


# TODO modulate L and decreasing number of CS if obviously dummy cs
fit.mococomo <- function(data,
                         model     = "normal",
                         maxiter   = 100,
                         tol       = 1e-3,
                         max_class = 10,
                         mult      = 2,
                         upper     = FALSE,
                         nullweight ) {

  if("data_mococomo"%!in% class(data))
  {stop("Please provide object of class data_mococomo")}
  if(missing( nullweight)){
    nullweight <- 4
  }

  fit <- init.mococomo(data       = data,
                       model      = model,
                       max_class  = max_class,
                       mult       = mult,
                       upper      = upper,
                       nullweight = nullweight
                       )

  fit$elbo <- compute_elbo.mococomo(fit)
  for (i in 1:maxiter) {

    fit <- iter.mococomo(fit,
                         update_assignment =is.even(i),
                         update_logreg = is.odd(i))

    fit$elbo <- c(fit$elbo, compute_elbo.mococomo(fit))

    # print(paste('asgn:', is.even(i), 'logreg:', is.odd(i), 'elbo: ', tail(fit$elbo, 1)))
    if (.converged(fit, tol)) {
      break
    }
  }
  return(fit)
}


