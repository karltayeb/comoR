
#### Simulation cEBMF----
sim_func_cEBMF <- function( N=200, # number of row
                            L=100, #number of columns
                            K=2, #number of factor
                            P1=200, # number of cov for row /loadings
                            P2=200, # number of cov for col /factors
                            beta0=-2,
                            beta1=2,
                            noise_level= 3,
                            max_iter_cEBMF=20,
                            max_iter_como=20,
                            max_class=10,
                            seed
){

library(softImpute)
  library(susieR)
  library(mvtnorm)
  data(N3finemapping)
  attach(N3finemapping)
library(comoR)

  if( missing( seed)){
    set.seed(rpois(lambda = 100,n=1))
  }else{
    set.seed(seed)
  }
  L_l <-  sample (1:20,size=1)
  L_f <-  sample (1:20,size=1)





  X_l <-   rmvnorm(N,sigma=cov(N3finemapping$X[1:100,1:P1]))
  X_f <-    rmvnorm(L,sigma=cov(N3finemapping$X[1:100,1:P2]))


  true_pos_l <- sample( 1:P1, size=L_l, replace=FALSE)
  true_pos_f <- sample( 1:P2, size=L_f, replace=FALSE)




  true_l  <- list()
  true_f  <- list()

  for( k in 1:K){



    samp_prob <- 1/(1 +exp(-(beta0+beta1*X_l[,k])))
    lk <- c()
    mix <- c()
    for ( i in 1:N){
      mix <-c(mix, sample(c(0,1), size=1, prob = c(1- samp_prob[i], samp_prob[i])))
      lk <- c(lk, mix[i]*rnorm(1,sd=3))
    }


    samp_prob <- 1/(1 +exp(-(beta0+beta1*X_f[,k])))
    fk <- c()
    mix <- c()
    for ( j in 1:L){

      mix <-c(mix, sample(c(0,1), size=1, prob = c(1- samp_prob[j], samp_prob[j])))
      fk <- c(fk, mix[j]*rnorm(1,sd=3))
    }

    true_l[[k]] <- lk
    true_f[[k]] <- fk

  }


  Y_true <- Reduce("+", lapply( 1:K, function(k) true_l[[k]]%*%t(true_f[[k]])
  )
  )

  Y_obs <- Y_true+ matrix( rnorm(N*L, sd= noise_level), ncol=L)


  res <- cEBMF  (Y=Y_obs,
                 X_l=X_l,
                 X_f=X_f,
                 reg_method="logistic_susie",
                 K=K, init_type = "udv_si",
                 param_como = list(max_class=max_class,mnreg="mult_reg"),
                 maxit_como =max_iter_como ,
                 param_nnet= list(size=3, decay=1.2),
                 maxit=max_iter_cEBMF)



  f <- flashier::flash(Y_obs)


  rmse_cEBMF   <- sqrt(mean( (Y_true-Y_est)^2))
  rmse_flash   <-  sqrt(mean( (Y_true- fitted(f))^2))
  rmse         <- c(rmse_cEBMF, rmse_flash)
  names(rmse ) <- c("rmse_cEBMF", "rmse_flash")


  par <-  c( N,
             L ,
             K ,
             P1 ,
             P2 ,
             beta0 ,
             beta1 ,
             noise_level ,
             max_iter_cEBMF ,
             max_iter_como
  )


  out <- list(rmse      = rmse,
              par       = par,
              flash.obj = f,
              cEBMF.obj = cEBMF.obj
  )
  return( out)
}


tt <-sim_func_cEBMF ()

