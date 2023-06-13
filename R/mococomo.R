# Implements covariate moderated ASH "MOre COmponents COvariate MOderated"

#' @title Function implementation the mococomo mode
#' @details Function implementation the mococomo mode
#'
#' @param betahat the estimated coefficient
#' @param se, the  corresponding standard error (set to 1 if not provided)
#' @param X covariate of interest
#' @param  Z additional covariate for adjustement
#' @param mnreg_type  character that specify the type of regression method used for
#' latent states default 'mult_reg'
#' @param max_class maximum number of class
#' @param scales  vector of positive value to define the prior variance manually.
#'  If not specified otherwise done following the ash procedure
#' @param nullweight penalty for the null component
#' @param tol stopping criterion
#' @param max_iter maximum number of iteration
#' @export
#'
#'effect_var <- 3
#'
#'
#'N=5000
#'x1 <- rnorm(N,sd=3)
#'beta0=-2
#'beta1=1
#'samp_prob <- 1/(1 +exp(-(beta0+beta1*x1)))
#'P=20
#'mix <- c()
#'betahat <- c()
#'betatrue <- c()
#'
#'X <- cbind( x1, matrix(rnorm(P*N), ncol=P))
#'
#'se <-   rchisq(N, df=1 )
#'#se <- rep(1,N)
#'for ( i in 1:N){
#'  mix <-c(mix, sample(c(0,1), size=1, prob = c(1- samp_prob[i], samp_prob[i])))
#'  betatrue <- c(betatrue,  mix[i] *rnorm(1,sd=effect_var))
#'  betahat <- c( betahat ,   betatrue[i]+rnorm(1,sd=se[i] ) )
#'}
#'
#'#p <- runif(N)
#'X <- cbind( x1, matrix(rnorm(P*N), ncol=P))
#'plot( x1,betahat, col=mix+1)
#'plot( x1,betatrue)
#'
#'fit <- mococomo(betahat=betahat, se=se ,X=X,max_iter=5, nullweight = 0.1)
#'df1 <- data.frame( betahat =betahat,
#'                   x =x1, col=fit$post_assignment[,1])
#'P1 <- ggplot(df1 , aes( y=betahat ,x =x1, col=fit$post_assignment[,1]))+
#'            geom_point()+
#'            ggtitle(("Posterior assignment for\n null comp  mococomo"))
#'
#'res_ash <- ash(betahat , se)
#'df3 <- data.frame( betahat =betahat,
#'                   x =x1, col=res_ash$result$lfdr)
#'P3 <- ggplot(df3 , aes( y=betahat ,x =x , col=col))+geom_point()+ggtitle(("Posterior assignment for\n null comp ash"))
#'library(gridExtra)
#'grid.arrange(P1,P3,ncol=2)
#'plot(  betahat,fit$result$mean)
#'abline(a=0,b=1)
#'points(betahat,res_ash$result$PosteriorMean, col="green")
#'plot(  betatrue, fit$result$mean)
#'abline(a=0,b=1)
#'points(betahat,res_ash$result$PosteriorMean, col="green")
#'#Lower RMSE for mococomo
#'sqrt(sum( (fit$result$mean  -betatrue )^2))
#'sqrt(sum( (res_ash$result$PosteriorMean-betatrue )^2))



mococomo <- function( betahat,
                      se,
                      X,
                      Z,
                      mnreg_type='mult_reg',
                      max_class=20,
                      scales=NULL,
                      nullweight=.1,
                      tol= 1e-3,
                      max_iter=100){

  if (missing( Z)){
    Z <-  rep( 1, length(betahat))
  }
  if( missing( se)){
    se=rep( 1, length(betahat))
  }

  data <- prep_data_como2(betahat=betahat,
                          se=se,
                          X=X ,
                          Z=Z)


  fit <- data_initialize_como(data       = data,
                              max_class  = max_class,
                              scales     = scales,
                              mnreg_type = mnreg_type ,
                              nullweight = nullweight )

  #fit <-logisticsusie:::fit_model.default(fit, data, max_iter = max_iter)
  for(i in 1:max_iter){
    fit <- update_model(x=fit, data =data)


    # Check convergence, assumes fit is tracking elbo
    if(abs(diff(tail(fit$elbo, 2))) < tol){
      message('converged')
      break
    }
  }
  fit$result  <-    post_mean_sd.mococomo(fit,data)
  return(fit)
}










initialize_como <- function(scales, n, p, p2, mu0=0, var0=1, nullweight=0, mnreg_type='constant'){
  # initialize multinomial susie-- but could be any multinomial regression
  K <- length(scales)


    mnreg <- initialize_mnreg (mnreg_type = mnreg_type,
                               K          = K,
                               n          = n,
                               p          = p)

  #mn_reg <- logisticsusie:::initialize_sbmn_susie(K, n, p, p2, L, mu0, var0)

  # initialize_scales
  f_list <- purrr::map(scales, ~ normal_component(mu = 0, var = .x^2))

  fit <- list(
    mnreg = mnreg, # multinomial regression function. takes X, returns pi
    f_list = f_list, # component distributions
    nullweight = nullweight, # penalty promoting the first component,
    K = K,
    elbo = -Inf
  )
  class(fit) <- c('como')
  return(fit)
}

#' Use data to autoselect scales
data_initialize_como <- function(data, max_class, scales=NULL, mu0=0, var0=1, nullweight=0, mnreg_type='constant') {
  como_check_data(data)

  if(is.null(scales)){
    scales <- autoselect_scales(data$betahat, data$se, max_class)
  }

  K <- length(scales) # K <= max_class
  p <- ncol(data$X)
  n <- nrow(data$X)
  p2 <- ncol(data$Z)

  fit <- initialize_como(scales=scales,
                         n=n,
                         p=p,
                         p2=p2,
                         mu0=mu0,
                         var0=var0,
                         nullweight=nullweight,
                         mnreg_type=mnreg_type)
  return(fit)
}

#' @export
update_model.como <- function(x, data, update_assignment = T, update_logreg=T, fit_prior_variance=F, track_elbo=T){
  fit <- x
  K <- fit$K

  # pre-compute data likelihood, if we haven't already
  # TODO: use digest::digest to hash the scales and recompute if scales change?
  if(is.null(fit$data_loglik)){
    fit$data_loglik <- compute_data_loglikelihood(fit, data)
  }

  # updates posterior assignments probabilities
  # these are the response variable for mn_regression
  if (update_assignment) {
    fit$post_assignment <- compute_posterior_assignment(fit = fit, data = data)
    #data$Y <- fit$post_assignment
    #data$Nk <- logisticsusie:::stick_breaking_N(data$Y)
  }

  if (update_logreg) {
    fit$mnreg <- update_prior(mnreg = fit$mnreg,
                              resps = fit$post_assignment,
                              data  = data)
  }

  if (track_elbo){
    fit$elbo <- c(fit$elbo, compute_elbo(fit, data))
  }
  return(fit)
}

#' @export
compute_elbo.como <- function(fit, data) {
  # E[log p(beta | y)] -- expected data likelihood
  ll <- sum(fit$post_assignment * fit$data_loglik)

  # Entropy term # E[-log q(y)]
  assignment_entropy <- sum(apply(fit$post_assignment, 1, logisticsusie:::categorical_entropy))

  # E[log p(y | X, theta)] - KL[q(theta) || p(theta)] SuSiE ELBO
  elbo <- compute_elbo(fit$mnreg, fit$post_assignment, data)

  # put it all together
  elbo <- ll - assignment_entropy + elbo
  return(elbo)
}


