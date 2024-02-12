x <-runif(1000)
y <-runif(1000)
X = cbind(x,y)
plot (x,y)

set.seed(1)
f <- matrix(NA, nrow = 2, ncol =200)
for ( i in 1:ncol (f)){

  t1<- sample (c(0,1), size=1)
  t2<- sample (c(0,1), size=1)

  f[1,i] <- t1*rnorm(n=1)
  f[2,i] <- t2*rnorm(n=1)

}
L <- matrix(NA, ncol=2, nrow=length(x))

for (i in 1:length(x)){

  if ( (x[i] <.5 & y[i] <.5 )|(x[i] >.5 & y[i] >.5 )){
    L[i,] <- c(1,0)
  }else{
    L[i,] <- c(0,1)
  }


}


plot ( x* c( L[,1]), y* c( L[,1]))

Z = L%*%f + matrix(rnorm(nrow(L)* ncol(f)), nrow = nrow(L))


library(flashier)
library(keras)
fit_default <- flash(Z, greedy_Kmax = 5)

Y = matrix(rnorm(ncol(f)*2), ncol=2)

cebnm_L <- function( x,s,g_init=FALSE,fix_g=TRUE, output){

  if (length(x) == 3){ ### just to satisfy check of custom function
    return (ebnm_flat(x))
  }
  Z <- matrix( 1, nrow=length(x), ncol=1)
  param_como = list(max_class= 10,
                    mnreg_type="keras")
  data <- comoR:::prep_data_como2 (betahat=x,
                                   se=s, X=X,
                                   Z =Z )

  # you need to retreive the actual number of mixture component in the model
  num_classes <- length( autoselect_scales(data$betahat, data$se,10))

  #define the nnet paramet using Keras syntax
  param_nnet =keras_model_sequential() %>%
    layer_dense(units = 64,
                activation = 'relu',
                input_shape = c(ncol(X))) %>%
    layer_dense(units = 64,
                activation = 'relu' ) %>%
    layer_dense(units = 64,
                activation = 'relu' ) %>%
    layer_dense(units = num_classes,
                activation = 'softmax')

  # run comoR
  fit  <- rlang::exec( "data_initialize_como", !!! param_como ,
                            data= data,
                            param_nnet= param_nnet) # initialize the model from the data
  fit <- comoR:::fit.como (  fit, data, max_iter = 3 )






  g <- ashr::normalmix(rep(1/length(fit$f_list),length(fit$f_list)),
                       rep( 0, length(fit$f_list)),
                       do.call(c, lapply( 1: length(fit$f_list) ,
                                          function(k) {sqrt(fit$f_list [[k]]$var) } )
                       )
  )

  out <- list( data= data.frame(x=data$betahat,
                                s=data$se),
               posterior = data.frame(mean= est$mean,
                                      second_moment=(est$sd^2+est$mean^2)
               ) ,
               fitted_g = g,
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

  param_como = list(max_class=10,mnreg_type='constant')
  param_nnet =list( )

  data <- comoR:::prep_data_como2 (betahat=x,
                                   se=s, X=Y,
                                   Z =Z )
  fit <- rlang::exec( "data_initialize_como", !!! param_como ,
                      data= data ) # initialize the model from the data
  fit <- comoR:::fit.como ( fit, data, max_iter = 5 )

  est <- comoR:::post_mean_sd (fit,data)

  g <- ashr::normalmix(rep(1/length(fit$f_list),length(fit$f_list)),
                       rep( 0, length(fit$f_list)),
                       do.call(c, lapply( 1: length(fit$f_list) ,
                                          function(k) {sqrt(fit$f_list [[k]]$var) } )
                       )
  )

  out <- list( data= data.frame(x=data$betahat,
                                s=data$se),
               posterior = data.frame(mean= est$mean,
                                      second_moment= (est$sd^2+est$mean^2)
               ) ,
               fitted_g = g,
               log_likelihood=sum( comoR:::compute_data_loglikelihood(fit, data) * (fit$post_assignment))

  )
  return( out)

}

library(flashier)
fit_custom <- flash_init(Z, var_type = 2) %>%
  flash_set_verbose(0) %>%
  flash_greedy(
    Kmax = 2,
    ebnm_fn = c(cebnm_L, ebnm_point_laplace),maxiter = 2
  )


cor (c(fitted(fit_default )) ,c(L%*%f))
cor (c(fitted(fit_custom )) ,c(L%*%f))


plot(fitted(fit_default ) ,L%*%f )
points(fitted(fit_custom ) ,L%*%f  , col="lightgreen")
cor (c(fitted(fit_default )) ,c(L%*%f))
cor (c(fitted(fit_custom )) ,c(L%*%f))



plot(fit_default$L_pm,x)
