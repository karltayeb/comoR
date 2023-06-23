#'@export
c_ash_sim <- function( N=1000,
                                    beta0=-2,
                                    beta1=1,
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

  x1 <- rnorm(N,sd=var_cov_effect)
  samp_prob <- 1/(1 +exp(-(beta0+beta1*x1)))
  betahat <- c()
  betatrue <- c()
  mix <- c()
  X <- cbind( x1, matrix(rnorm(P*N), ncol=P))





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

  X <- cbind( x1, matrix(rnorm(P*N), ncol=P))
  res <- mococomo(betahat = betahat,se=se, X,max_iter=max_iter, nullweight = nullweight)
  tt <- ash(betahat, se)


  #Power
  power_ash <-  sum( mix[which( tt$result$lfdr<0.05)])/length(which( tt$result$lfdr<0.05))
  power_mco <- sum( mix[which(res $post_assignment[,1]<0.05)])/length(which(res $post_assignment[,1]<0.05))
  #T1 error
  T1_ash <-  length(which( mix[which( tt$result$lfdr<0.05)]==0))/length(which( tt$result$lfdr<0.05))
  T1_mco <-length(which( mix[which(res $post_assignment[,1]<0.05)]==0))/length(which(res $post_assignment[,1]<0.05))

  rmse_mco <-  sqrt(mean( (res$result$mean -  betatrue)^2))

  rmse_ash <- sqrt(mean( (( tt$result$PosteriorMean -  betatrue)^2)))


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
            dist
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
                    "dist"
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


