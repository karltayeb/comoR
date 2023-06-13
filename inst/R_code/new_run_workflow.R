rm(list=ls())
library(nnet)
library(ggplot2)
devtools::load_all(".")

effect_var <- 3


N=1000
x1 <- rnorm(N,sd=3)
beta0=-2
beta1=1
samp_prob <- 1/(1 +exp(-(beta0+beta1*x1)))
P=20
mix <- c()
betahat <- c()
betatrue <- c()

X <- cbind( x1, matrix(rnorm(P*N), ncol=P))

se <-   rchisq(N, df=1 )
#se <- rep(1,N)
for ( i in 1:N){
  mix <-c(mix, sample(c(0,1), size=1, prob = c(1- samp_prob[i], samp_prob[i])))
  betatrue <- c(betatrue,  mix[i] *rnorm(1,sd=effect_var))
  betahat <- c( betahat ,   betatrue[i]+rnorm(1,sd=se[i] ) )
}

#p <- runif(N)
X <- cbind( x1, matrix(rnorm(P*N), ncol=P))
plot( x1,betahat, col=mix+1)
plot( x1,betatrue)


plot( x1,betahat/se, col=mix+1)
plot( betatrue,betahat/se, col=mix+1)


mix# if 0 correspond to beta=0


data <- prep_data_como2(betahat=betahat, se=se ,
                        X=X ,
                        Z= rep( 1, length(betahat)))


fit <- data_initialize_como(data=data,max_class =  10, scales = seq(from=0, to=10, length.out=10)
                            ,mnreg_type='mult_reg' ,
                            nullweight=1 )

fit$data_loglik <- compute_data_loglikelihood(fit, data)
fit$post_assignment <- compute_posterior_assignment(fit = fit, data = data)
image(fit$data_loglik)


plot( fit$post_assignment[,1],se)
plot( fit$post_assignment[,2],se)






df1 <- data.frame( betahat =betahat,
                   x =x1, col=fit$post_assignment[,1])
 P1 <- ggplot(df1 , aes( y=betahat ,x =x1, col=fit$post_assignment[,1]))+geom_point()
 ggplot(df1 , aes( y=betahat ,x =x, col=col))+geom_point()

fit <- data_initialize_como(data=data,max_class =  10, scales = seq(from=0, to=20, length.out=25)
                            ,mnreg_type='mult_reg' ,
                            nullweight=0.1 )

#fit <- fit_model(fit, data, max_iter =20)
max_iter=20
tol=1e-3
for(i in 1:max_iter){
  fit <- update_model(x=fit, data =data)

  plot( betahat,  fit$post_assignment[,1],col=mix+1 )
  points(  betahat, fit$post_assignment[,10],col=mix+3 )

  points(betahat,  fit$post_assignment[,15], col=mix+5 )
  points( betahat, fit$post_assignment[,20] , col=mix+7 )

  # Check convergence, assumes fit is tracking elbo
  if(abs(diff(tail(fit$elbo, 2))) < tol){
    message('converged')
    break
  }
}
library(ggplot2)
df2 <- data.frame( betahat =betahat,
                   x =x1, col=fit$post_assignment[,1])
P2 <- ggplot(df2 , aes( y=betahat ,x =x , col=col))+geom_point()

P2
 library(gridExtra)






library(ashr)
tt <- ash(betahat, se)
df3 <- data.frame( betahat =betahat,
                   x =x1, col=tt$result$lfdr)
P3 <- ggplot(df3 , aes( y=betahat ,x =x , col=col))+geom_point()

P3

plot( tt$result$lfdr, fit $post_assignment[,1], col=mix+1)


resb <- mococomo(betahat = betahat,se=se,X=X, max_iter = 5)

res <- post_mean_sd.mococomo(fit,data)

plot(  betahat,res$mean)
abline(a=0,b=1)
points(betahat, tt$result$PosteriorMean, col="green")
points(betahat, resb$result$mean, col="red")


plot(  resb$result$mean,res$mean)
abline(a=0,b=1)


plot(  betatrue,res$mean)
abline(a=0,b=1)
points(betatrue, tt$result$PosteriorMean, col="green")



#Power
sum( mix[which( tt$result$lfdr<0.05)])/length(which( tt$result$lfdr<0.05))
sum( mix[which(fit $post_assignment[,1]<0.05)])/length(which(fit $post_assignment[,1]<0.05))
#T1 error
length(which( mix[which( tt$result$lfdr<0.05)]==0))/length(which( tt$result$lfdr<0.05))
length(which( mix[which(fit $post_assignment[,1]<0.05)]==0))/length(which(fit $post_assignment[,1]<0.05))


library(ggplot2)

sqrt(sum( (res$mean -  betatrue)^2))

sqrt(sum( (( tt$result$PosteriorMean -  betatrue)^2)))
