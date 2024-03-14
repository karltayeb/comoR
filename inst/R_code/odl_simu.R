devtools::load_all(".")
 N=2000
                            L=100  #number of columns
                            K=2  #number of factor
                            P1=20  # number of cov for row /loadings
                            P2=20 # number of cov for col /factors
                            beta0=-2
                            beta1=2
                            noise_level= 3
                            max_iter_cEBMF=20
                            max_iter_como=20
                            max_class=10
                            seed=seed+10

  set.seed(seed)
  library(softImpute)
  library(susieR)
  library(mvtnorm)
  data(N3finemapping)
  attach(N3finemapping)
  library(comoR)


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



  res_nnet <- cEBMF  (Y=Y_obs,
                      X_l=X_l,
                      X_f=X_f,
                      reg_method="nnet",
                      K=K, init_type = "udv_si",
                      param_como = list(max_class=max_class,mnreg="mult_reg"),
                      maxit_como =max_iter_como ,
                      param_nnet= list(size=3, decay=1.2),
                      maxit=max_iter_cEBMF)

  f <- flashier::flash(Y_obs)
  Y_est_nnet <- Reduce("+", lapply( 1:res_nnet$K, function(k) res_nnet $loading[,k] %*%t(res_nnet $factor[,k] ) ))



  library(irlba)
  library(PMA)
  ssvd_res = ssvd(Y_obs, k=3)
  svd_res  = svd(Y_obs)

  rmse = function(mat1 , mat2){

    squared_diff <- (mat1-  mat2)^2

    # Compute mean of squared differences
    mean_squared_diff <- mean(squared_diff)

    # Compute RMSE
    rmse <- sqrt(mean_squared_diff)
    return (rmse)
  }


  rmse(Y_true, svd_res$u%*%diag(svd_res$d)%*%t(svd_res$v))


  cv.out <- PMD.cv(Y_obs, type="standard", sumabss=seq(0.1, 0.6, len=20))
  PMD_res <- PMD(Y_obs,
                 type="standard",
                 sumabs=cv.out$bestsumabs,
                 K=3, v=cv.out$v.init
  )



  rmse_cEBMF_nnet   <-  rmse(Y_true , Y_est_nnet )
  rmse_flash        <-  rmse(Y_true ,fitted(f))

  rmse_PMD         <- rmse(Y_true, PMD_res$u%*%diag(PMD_res$d)%*%t(PMD_res$v))
  rmse_svd         <- rmse(Y_true, svd_res$u%*%diag(svd_res$d)%*%t(svd_res$v))
  rmse_ssvd        <- rmse(Y_true, ssvd_res$u%*%ssvd_res$d%*%t(ssvd_res$v))
  rmse_out         <- c( rmse_cEBMF_nnet ,rmse_flash  , rmse_PMD, rmse_svd, rmse_ssvd)

  rmse_out
