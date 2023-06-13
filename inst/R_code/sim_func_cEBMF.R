sim_func_cEBMF <- function( N=200, # number of row
                            L=100, #number of columns
                            K=2, #number of factor
                            P1=20, # number of cov for row /loadings
                            P2=20, # number of cov for col /factors
                            beta0=-2,
                            beta1=2,
                            noise_level= 3,
                            max_iter_cEBMF=20,
                            max_iter_mococomo=20
){

  X_l <-   matrix(rnorm(P1*N, sd=3), ncol=P1)
  X_f <-   matrix(rnorm(P2*L, sd=3), ncol=P2)

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


  cEBMF.obj <- init_cEBMF (Y=Y_obs,
                           X_l=X_l,
                           X_f=X_f,
                           K=K, init_type = "udv_si")






  for ( o in 1:max_iter_cEBMF){
    for ( k in 1:K){
      Rk <- cal_partial_residuals.cEBMF(cEBMF.obj,k)
      l_k <- cal_expected_loading( cEBMF.obj, Rk,k)

      t_fit <- mococomo(betahat    = l_k$l_i_hat,
                        se         = l_k$s_i,
                        X          = X,
                        max_iter   = 20,
                        nullweight = 0.1)


      cEBMF.obj$loading[,k]  <-  t_fit$result$mean
      cEBMF.obj$loading2[,k] <-  t_fit$result$sd^2+ t_fit$result$mean^2

      #factor update

      f_k <- cal_expected_factor( cEBMF.obj, Rk,k)
      t_fit <- mococomo(betahat    = f_k$f_j_hat,
                        se         = f_k$s_j,
                        X          = Y,
                        max_iter   = max_iter_mococomo,
                        nullweight = 0.1)


      cEBMF.obj$factor[,k]  <-   t_fit$result$mean
      cEBMF.obj$factor2[,k] <-   t_fit$result$sd^2+  t_fit$result$mean^2


      cEBMF.obj<- update_tau.cEBMF (cEBMF.obj )

      Y_est <- cEBMF.obj$loading[,1]%*%t(cEBMF.obj$factor[,1])
      plot( Y_est, Y_obs )
      abline(a=0,b=1)
    }

  }

  Y_est <- Reduce("+", lapply(1:ncol(cEBMF.obj$factor), function(k)
                                                      cEBMF.obj$loading[,k]%*%t(cEBMF.obj$factor[,k])
                              )
                  )
  f <- flash(Y_obs)


  rmse_cEBMF   <- sqrt(sum( (Y_true-Y_est)^2))
  rmse_flash   <-  sqrt(sum( (Y_true- fitted(f))^2))
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
            max_iter_mococomo
     )


 out <- list(rmse      = rmse,
             par       = par,
             flash.obj = f,
             cEBMF.obj = cEBMF.obj
             )

}
