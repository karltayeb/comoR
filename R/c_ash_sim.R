#'@export
c_ash_sim <- function( N=1000,
                       beta0=-2,
                       beta1=2,
                       var_cov_effect=3,
                       effect_var=3,
                       noise_level=1,
                       nullweight=.1,
                       P=20,
                       se_type="constant",
                       df_se=2,
                       max_iter=5,
                       dist="normal",
                       extended=FALSE

){
  L <- sample(1:20, size=1)
  library(susieR)
  library(mvtnorm)
  data(N3finemapping)
  attach(N3finemapping)
  x1 <- rmvnorm(N,sigma=cov(N3finemapping$X[1:100,]))

  true_pos <- sample( 1:ncol(x1), L)
  lin_pred <- rep(0,N)


    for ( l in 1:L){
      lin_pred <- lin_pred+beta1*x1 [ ,true_pos[[l]]]
    }

  lin_pred <- lin_pred+beta0

  samp_prob <- 1/(1 +exp(-(lin_pred)))
  betahat <- c()
  betatrue <- c()
  mix <- c()






  if( dist=="spiky"){

    samp_fun <- function(){
      id <- sample(size=1, 1:4)
      if(id==1){
        out <-  rnorm( 1, sd=0.25)
      }
      if(id==2){
        out <-  rnorm( 1, sd=0.5)
      }
      if(id==3){
        out <-  rnorm( 1, sd=1)
      }
      if(id==4){
        out <-  rnorm( 1, sd=2)
      }
      return( out)
    }
    for ( i in 1:N){
      mix <-c(mix, sample(c(0,1), size=1, prob = c(1- samp_prob[i], samp_prob[i])))
      betatrue <- c( betatrue, mix[i] *samp_fun())
      betahat <- c( betahat ,   betatrue[i]+rnorm(1,sd=noise_level ) )
    }
    if(se_type=="constant"){
      se <- rep(1, length(betahat))

    }
    if( se_type=="random"){
      se <-   rchisq(N, df=df_se)
    }


    tt1 <-  ( 1- samp_prob)  *dnorm(betahat, sd=noise_level)
    tt2 <- list( )
    for ( i in 1: length(betahat)){
      tt2[[i]]<-   (   samp_prob[i])*0.25*

        dnorm(  betahat[i] ,
                mean = 0,sd= sqrt(noise_level^2+ c(0.25,0.5,1,2)^2))
    }

    tt2 <- as.matrix(do.call(rbind,tt2) )
    true_lfdr <-  tt1/ ( tt1+apply(tt2,1, sum))

  }

  if( dist=="near_normal"){
    samp_fun <- function(){
      id <- sample(size=1, 1:2, prob=c(2/3,1/3))
      if(id==1){
        out <-  rnorm( 1, sd=1)
      }
      if(id==2){
        out <-  rnorm( 1, sd=2)
      }

      return( out)
    }

    for ( i in 1:N){
      mix <-c(mix, sample(c(0,1), size=1, prob = c(1- samp_prob[i], samp_prob[i])))
      betatrue <- c( betatrue, mix[i] *samp_fun())
      betahat <- c( betahat ,   betatrue[i]+rnorm(1,sd=noise_level ) )
    }
    if(se_type=="constant"){
      se <- rep(1, length(betahat))

    }
    if( se_type=="random"){
      se <-   rchisq(N, df=df_se)
    }
    tt1 <-  ( 1- samp_prob)  *dnorm(betahat, sd=noise_level)
    tt2 <- list( )
    for ( i in 1: length(betahat)){
      tt2[[i]]<-   (  samp_prob[i])*( (2/3)* dnorm(  betahat[i] ,
                                                     mean = 0,sd= sqrt(noise_level^2+ 1^2))+
                                        (1/3)* dnorm(  betahat[i] ,
                                                       mean = 0,sd= sqrt(noise_level^2+ 2^2)))


    }

    tt2 <- as.matrix(do.call(rbind,tt2) )
    true_lfdr <-  tt1/ ( tt1+apply(tt2,1, sum))

  }

  if( dist=="normal"){
    for ( i in 1:N){
      mix <-c(mix, sample(c(0,1), size=1, prob = c(1- samp_prob[i], samp_prob[i])))
      betatrue <- c( betatrue, mix[i] *rnorm(1,sd=1))
      betahat <- c( betahat ,   betatrue[i]+rnorm(1,sd=noise_level ) )
    }
    if(se_type=="constant"){
      se <- rep(1, length(betahat))

    }
    if( se_type=="random"){
      se <-   rchisq(N, df=df_se)
    }

    tt1 <-  ( 1- samp_prob)  *dnorm(betahat, sd=noise_level)
    tt2 <- list( )
    for ( i in 1: length(betahat)){
      tt2[[i]]<-   (  samp_prob[i])*(   dnorm(  betahat[i] ,
                                                mean = 0,sd= sqrt(noise_level^2+ 1^2))
      )


    }

    tt2 <-  (do.call(c,tt2) )
    true_lfdr <-  tt1/ ( tt1+tt2)
  }

  if( dist=="flattop"){

    samp_fun <- function(){
      id <- sample(size=1, 1:7 )
      if(id==1){
        out <-  rnorm( 1,mean=-1.5, sd=5)
      }
      if(id==2){
        out <-  rnorm( 1,mean=-1, sd=5)
      }
      if(id==3){
        out <-  rnorm( 1,  mean=-.5, sd=5)
      }
      if(id==4){
        out <-  rnorm( 1,   sd=5)
      }
      if(id==5){
        out <-  rnorm( 1,  mean= .5, sd=5)
      }
      if(id==6){
        out <-  rnorm( 1, mean= 1, sd=5)
      }
      if(id==7){
        out <-  rnorm( 1,  mean= 1.5, sd=5)
      }
      return( out)
    }

    for ( i in 1:N){
      mix <-c(mix, sample(c(0,1), size=1,
                          prob = c(1- samp_prob[i], samp_prob[i])))
      betatrue <- c( betatrue, mix[i] *samp_fun())
      betahat <- c( betahat ,   betatrue[i]+rnorm(1,sd=noise_level ) )
    }
    if(se_type=="constant"){
      se <- rep(1, length(betahat))

    }
    if( se_type=="random"){
      se <-   rchisq(N, df=df_se)
    }
    tt1 <-  ( 1- samp_prob)  *dnorm(betahat, sd=noise_level)
    tt2 <- list( )
    for ( i in 1: length(betahat)){
      tt2[[i]]<-   (  samp_prob[i])*(1/7)*(dnorm(  betahat[i] ,
                                                   mean = c(-1.5,-1,-.5,0,.5,1,1.5),
                                                   sd= sqrt(noise_level^2+ 5^2)) )


    }

    tt2 <-  (do.call(rbind,tt2) )
    true_lfdr <-  tt1/ ( tt1+apply(tt2,1, sum))

  }

  if( dist=="skew"){



    samp_fun <- function(){
      id <- sample(size=1, 1:4, prob = c(1/4,1/4,1/3,1/6) )
      if(id==1){
        out <-  rnorm( 1,mean=-2, sd=2)
      }
      if(id==2){
        out <-  rnorm( 1,mean=-1, sd=1.5)
      }
      if(id==3){
        out <-  rnorm( 1,  mean=0, sd=1)
      }
      if(id==4){
        out <-  rnorm( 1,mean=1,   sd=1)
      }

      return( out)
    }

    for ( i in 1:N){
      mix <-c(mix, sample(c(0,1), size=1, prob = c(1- samp_prob[i], samp_prob[i])))
      betatrue <- c( betatrue, mix[i] *samp_fun())
      betahat <- c( betahat ,   betatrue[i]+rnorm(1,sd=noise_level ) )
    }
    if(se_type=="constant"){
      se <- rep(1, length(betahat))

    }
    if( se_type=="random"){
      se <-   rchisq(N, df=df_se)
    }
    tt1 <-  ( 1- samp_prob)  *dnorm(betahat, sd=noise_level)
    tt2 <- list( )
    for ( i in 1: length(betahat)){
      tt2[[i]]<-   (  samp_prob[i])*(  (1/4)* dnorm(betahat[i], mean=-2, sd= sqrt(noise_level^2+ 2^2))+
                                         (1/4)* dnorm(betahat[i], mean=-1, sd= sqrt(noise_level^2+ 1.5^2))+
                                         (1/3)* dnorm(betahat[i], mean=0, sd= sqrt(noise_level^2+ 1^2))+
                                         (1/6)* dnorm(betahat[i], mean=1, sd= sqrt(noise_level^2+ 1^2))


      )



    }

    tt2 <-  (do.call(c,tt2) )
    true_lfdr <-  tt1/ ( tt1+tt2)

  }

  if( dist=="big-normal"){
    for ( i in 1:N){
      mix <-c(mix, sample(c(0,1), size=1, prob = c(1- samp_prob[i], samp_prob[i])))
      betatrue <- c( betatrue, mix[i] *rnorm(1,sd=4))
      betahat <- c( betahat ,   betatrue[i]+rnorm(1,sd=noise_level ) )
    }
    if(se_type=="constant"){
      se <- rep(1, length(betahat))

    }
    if( se_type=="random"){
      se <-   rchisq(N, df=df_se)
    }

    tt1 <-  ( 1- samp_prob)  *dnorm(betahat, sd=noise_level)
    tt2 <- list( )
    for ( i in 1: length(betahat)){
      tt2[[i]]<-   (  samp_prob[i])*(   dnorm(  betahat[i] ,
                                                mean = 0,sd= sqrt(noise_level^2+ 4^2))
      )


    }

    tt2 <-  (do.call(c,tt2) )
    true_lfdr <-  tt1/ ( tt1+tt2)

  }
  if( dist=="bimodal"){
    samp_fun <- function(){
      id <- sample(size=1 ,1:2   )
      if(id==1){
        out <-  rnorm( 1,mean=-2, sd=1)
      }
      if(id==2){
        out <-  rnorm( 1,mean=2, sd=1)
      }


      return( out)
    }


    for ( i in 1:N){
      mix <-c(mix, sample(c(0,1), size=1, prob = c(1- samp_prob[i], samp_prob[i])))
      betatrue <- c( betatrue, mix[i] *rnorm(1,sd=4))
      betahat <- c( betahat ,   betatrue[i]+rnorm(1,sd=noise_level ) )
    }
    if(se_type=="constant"){
      se <- rep(1, length(betahat))

    }
    if( se_type=="random"){
      se <-   rchisq(N, df=df_se)
    }
    tt1 <-  ( 1- samp_prob)  *dnorm(betahat, sd=noise_level)
    tt2 <- list( )
    for ( i in 1: length(betahat)){
      tt2[[i]]<-   (  samp_prob[i])*( 0.5*dnorm(  betahat[i] , mean=-2, sd=sqrt(noise_level^2+ 1^2))  +
                                        0.5*dnorm(  betahat[i] , mean=2, sd=sqrt(noise_level^2+ 1^2))    )



    }

    tt2 <-  (do.call(c,tt2) )
    true_lfdr <-  tt1/ ( tt1+tt2)


  }
  X=x1
  res <- cFDR(betahat = betahat,se=se, X=X,max_iter=max_iter, nullweight = nullweight)
  tt <- ashr::ash(betahat, se)


  #Power
  if (length(which( tt$result$lfdr<0.05))==0){
    power_ash <- 0
    T1_ash <-  0
  }else{
    power_ash <-  sum( mix[which( tt$result$lfdr<0.05)])/length(which( tt$result$lfdr<0.05))
    T1_ash <-  length(which( mix[which( tt$result$lfdr<0.05)]==0))/length(which( tt$result$lfdr<0.05))

  }
  if (length(which(res$result$lfdr<0.05))==0){
    power_mco <- 0
    T1_mco <-  0
  }else{
    power_mco <- sum( mix[which(res$result$lfdr <0.05)])/length(which(res$result$lfdr<0.05))
    T1_mco <-length(which( mix[which(res$result$lfdr<0.05)]==0))/length(which(res$result$lfdr<0.05))

  }
  #T1 error

  rmse_mco <-  sqrt(mean( (res$result$PosteriorMean -  betatrue)^2))

  rmse_ash <- sqrt(mean( (( tt$result$PosteriorMean -  betatrue)^2)))

  is_dummy = res$is_dummy
  n_cs =  length(res$cs)  #number of CS
  n_effect <-   length(which(true_pos%in% do.call(c,

                                    lapply(1:length(res$cs), function(k)
                                      res$cs[[k]]$cs
                                      )
                                    )
               )
         )
 #number of effect found
  nfalse_effect <-  Reduce("+",sapply(1:length(res$cs), function(k)
    ifelse( length(which(true_pos%in%res$cs[[k]]$cs ))==0, 1,0)
   )
  )

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
            max_iter,
            dist,
            is_dummy,
            n_cs,
            n_effect,
            nfalse_effect
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
                    "max_iter",
                    "dist",
                    "is_dummy",
                    "n_cs",
                    "n_effect",
                    "nfalse_effect"
  )


  if( extended ){


    lfdr_est = data.frame(lfdr_ash  =tt$result$lfdr,
                          lfdr_cash =res$result$lfdr,
                          true_lfdr = true_lfdr)

    out <- list( summary=out,
                 lfdr_est =lfdr_est
    )
  }
  return(out)


}


