sim_func_single_effect <- function( N=1000,
                                    beta0=-2,
                                    beta1=1,
                                    var_cov_effect=3,
                                    effect_var=3,
                                    noise_level=1,
                                    P=20,
                                    se_type="random",
                                    df_se=2,
                                    max_iter=5

){

  x1 <- rnorm(N,sd=var_cov_effect)
  samp_prob <- 1/(1 +exp(-(beta0+beta1*x1)))
  betahat <- c()
  betatrue <- c()
  mix <- c()
  X <- cbind( x1, matrix(rnorm(P*N), ncol=P))

  for ( i in 1:N){
    mix <-c(mix, sample(c(0,1), size=1, prob = c(1- samp_prob[i], samp_prob[i])))
    betatrue <- c( betatrue, mix[i] *rnorm(1,sd=effect_var))
    betahat <- c( betahat ,   betatrue[i]+rnorm(1,sd=noise_level ) )
  }
  if(se_type=="constant"){
    se <- rep(1, length(betahat))

  }
  if( se_type=="random"){
    se <-   rchisq(N, df=df_se)
  }


  X <- cbind( x1, matrix(rnorm(P*N), ncol=P))
  res <- mococomo(betahat = betahat,se=se, X,max_iter=max_iter)
  tt <- ash(betahat, se)


  #Power
  power_ash <-  sum( mix[which( tt$result$lfdr<0.05)])/length(which( tt$result$lfdr<0.05))
  power_mco <- sum( mix[which(res $post_assignment[,1]<0.05)])/length(which(res $post_assignment[,1]<0.05))
  #T1 error
  T1_ash <-  length(which( mix[which( tt$result$lfdr<0.05)]==0))/length(which( tt$result$lfdr<0.05))
  T1_mco <-length(which( mix[which(res $post_assignment[,1]<0.05)]==0))/length(which(res $post_assignment[,1]<0.05))

  rmse_mco <-  sqrt(sum( (res$result$mean -  betatrue)^2))

  rmse_ash <- sqrt(sum( (( tt$result$PosteriorMean -  betatrue)^2)))


  out <- c( rmse_mco,
            rmse_ash,
            power_mco,
            power_ash,
            T1_mco,
            T1_ash,
            N ,
            beta0  ,
            beta1  ,
            effect_var ,
            noise_level ,
            var_cov_effect ,
            P ,
            ifelse(se_type=="random",1,0) ,
            df_se,
            max_iter
  )
  names (out) <-  c("rmse_mco",
                    "rmse_ash",
                    "power_mco",
                    "power_ash",
                    "T1_mco",
                    "T1_ash",
                    "N" ,
                    "beta0"  ,
                    "beta1"  ,
                    "effect_var" ,
                    "noise_level" ,
                    "var_cov_effect" ,
                    "P" ,
                    "se_type"  ,
                    "df_se",
                    "max_iter"
  )

  return(out)


}


