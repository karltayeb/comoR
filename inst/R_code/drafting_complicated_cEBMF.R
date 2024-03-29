rm(list=ls())
library(comoR)
library(nnet)
library(comoR)
simdata <- sim_cfactor()
 devtools::load_all(".")

X <- simdata$X_l
Y <- simdata$X_f


cebnm_L <- function( x,s ){

  Z <- matrix( 1, nrow=length(x), ncol=1)

  param_como  = list(max_class=5,mnreg_type="mult_reg")
  param_nnet  =list( size=1, decay=1,MaxNWts = 10000)

  data <- prep_data_como2(betahat=x,
                          se=s, X=X,
                          Z =Z )
  fit <-  rlang::exec( "data_initialize_como", !!! param_como ,
                       data= data,
                       param_nnet=  param_nnet) # initialize the model from the data
  fit  <- fit_model(  fit,  data, max_iter = 5 )


  est    <- post_mean_sd  (fit,data)


  g <- ashr::normalmix(rep(1/length(fit$f_list),length(fit$f_list)),
                       rep( 0, length(fit$f_list)),
                       do.call(c, lapply( 1: length(fit$f_list)  ,
                                          function(k) {sqrt(fit$f_list [[k]]$var) } )
                              )
                       )


  out <- list( data= data.frame(x=data$betahat,
                                s=data$se),
               posterior = data.frame(mean= est$mean,
                                      sd= est$sd
               ) ,
               fitted_g =  g,
               log_likelihood=sum(compute_data_loglikelihood(fit, data) * (fit$post_assignment))

  )

  return( out)

}
cebnm_F <- function( x,s ){


  Z <- matrix( 1, nrow=length(x), ncol=1)
  param_como  = list(max_class=5,mnreg_type="mult_reg")
  param_nnet  =list( size=1, decay=1,MaxNWts = 10000)


  data <- prep_data_como2(betahat=x,
                          se=s, X=X,
                          Z =Z )
  fit <-  rlang::exec( "data_initialize_como", !!! param_como ,
                       data= data,
                       param_nnet=  param_nnet) # initialize the model from the data
  fit  <- fit_model(  fit,  data, max_iter = 5 )


  est    <- post_mean_sd  (fit,data)

  g <- ashr::normalmix(rep(1/length(fit$f_list),length(fit$f_list)),
                       rep( 0, length(fit$f_list)),
                       do.call(c, lapply( 1: length(fit$f_list)  ,
                                          function(k) {sqrt(fit$f_list [[k]]$var) } )
                       )
  )
  out <- list( data= data.frame(x=data$betahat,
                                s=data$se),
               posterior = data.frame(mean= est$mean,
                                      sd= est$sd
                                      ) ,
               fitted_g =  g,
               log_likelihood=sum(compute_data_loglikelihood(fit, data) * (fit$post_assignment))

          )

  return( out)

}



library(flashier)
fit_custom <- flash_init(simdata$Y_obs, var_type = 2) %>%
  flash_set_verbose(0) %>%
  flash_greedy(
    Kmax = 2,
    ebnm_fn = c(cebnm_L, cebnm_F)
  )
