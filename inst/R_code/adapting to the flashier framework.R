rm(list=ls())
library(comoR)
library(nnet)
library(comoR)
simdata <- comoR:::sim_cfactor()
 devtools::load_all(".")

X <- simdata$X_l
Y <- simdata$X_f




cebnm_L <- function( x,s,g_init=FALSE,fix_g=TRUE, output){


  if (length(x) == 3){ ### just to satisfy check of custom function
    return (ebnm_flat(x))
  }
  Z <- matrix( 1, nrow=length(x), ncol=1)

  param_como  = list(max_class=5,mnreg_type="mult_reg")
  param_nnet  =list( size=1, decay=1,MaxNWts = 10000)

  data <- comoR:::prep_data_como2 (betahat=x,
                          se=s, X=X,
                          Z =Z )
  fit <-  rlang::exec( "data_initialize_como", !!! param_como ,
                       data= data,
                       param_nnet=  param_nnet) # initialize the model from the data
  fit  <- comoR:::fit.como (  fit,  data, max_iter = 5 )


  est    <- comoR:::post_mean_sd  (fit,data)


  g <- ashr::normalmix(rep(1/length(fit$f_list),length(fit$f_list)),
                       rep( 0, length(fit$f_list)),
                       do.call(c, lapply( 1: length(fit$f_list)  ,
                                          function(k) {sqrt(fit$f_list [[k]]$var) } )
                              )
                       )


  out <- list( data= data.frame(x=data$betahat,
                                s=data$se),
               posterior = data.frame(mean= est$mean,
                                      second_moment= est$sd^2
               ) ,
               fitted_g =  g,
               log_likelihood=sum( comoR:::compute_data_loglikelihood(fit, data) * (fit$post_assignment))

  )

  return( out)

}
cebnm_F <- function( x,s,g_init,fix_g=TRUE, output){
  if (length(x) == 3){ ### just to satisfy check of custom function
    return (ebnm_flat(x))
  }

  Z <- matrix( 1, nrow=length(x), ncol=1)
  Z <- matrix( 1, nrow=length(x), ncol=1)

  param_como  = list(max_class=5,mnreg_type="mult_reg")
  param_nnet  =list( size=1, decay=1,MaxNWts = 10000)

  data <- comoR:::prep_data_como2 (betahat=x,
                                   se=s, X=Y,
                                   Z =Z )
  fit <-  rlang::exec( "data_initialize_como", !!! param_como ,
                       data= data,
                       param_nnet=  param_nnet) # initialize the model from the data
  fit  <- comoR:::fit.como (  fit,  data, max_iter = 5 )


  est    <- comoR:::post_mean_sd  (fit,data)


  g <- ashr::normalmix(rep(1/length(fit$f_list),length(fit$f_list)),
                       rep( 0, length(fit$f_list)),
                       do.call(c, lapply( 1: length(fit$f_list)  ,
                                          function(k) {sqrt(fit$f_list [[k]]$var) } )
                       )
  )


  out <- list( data= data.frame(x=data$betahat,
                                s=data$se),
               posterior = data.frame(mean= est$mean,
                                      second_moment= est$sd^2
               ) ,
               fitted_g =  g,
               log_likelihood=sum( comoR:::compute_data_loglikelihood(fit, data) * (fit$post_assignment))

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



X_l <- simdata$X_l
X_f <- simdata$X_f

res_cEBMF <- cEBMF  (Y=simdata$Y_obs,
        X_l=X_l,
        X_f=X_f,
        reg_method="logistic_susie",
        K=2, init_type = "udv_si",
        param_como = list(max_class=10,mnreg="mult_reg"),
        param_como2 = list(f1_dist = 'ash'),
        maxit_como =5 ,
        param_nnet= list(size=3, decay=1.2),
        maxit=5)


