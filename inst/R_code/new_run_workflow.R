rm(list=ls())
library(nnet)
devtools::load_all(".")
N=10000
x1 <- rnorm(N,sd=3)
beta0=-2
beta1=1
samp_prob <- 1/(1 +exp(-(beta0+beta1*x1)))

P=20
mix <- c()
betahat <- c()
for ( i in 1:N){
  mix <-c(mix, sample(c(0,1), size=1, prob = c(1- samp_prob[i], samp_prob[i])))
  betahat <- c( betahat , ifelse( mix[i]==1, rnorm(1,sd=1 ), rnorm(1,sd=3)))
}
#p <- runif(N)
X <- cbind( x1, matrix(rnorm(P*N), ncol=P))
plot( x1,betahat)
mix# if 0 correspond to beta=0


data <- prep_data_como2(betahat=betahat, se=rep( 1, length(betahat)),
                        X=X ,
                        Z= rep( 1, length(betahat)))


fit <- data_initialize_como(data=data,max_class =  10, scales = seq(from=0, to=10, length.out=10)
                            ,mnreg_type='mult_reg' ,
                            nullweight=1 )

fit$data_loglik <- compute_data_loglikelihood(fit, data)
fit$post_assignment <- compute_posterior_assignment(fit = fit, data = data)

df1 <- data.frame( betahat =betahat,
                   x =x1, col=fit$post_assignment[,1])
 P1 <- ggplot(df1 , aes( y=betahat ,x =x1, col=fit$post_assignment[,1]))+geom_point()
 ggplot(df1 , aes( y=betahat ,x =x, col=col))+geom_point()

fit <- data_initialize_como(data=data,max_class =  10, scales = seq(from=0, to=20, length.out=25)
                            ,mnreg_type='mult_reg' ,
                            nullweight=0.1 )

#fit <- fit_model(fit, data, max_iter =20)
max_iter=5
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

P1
P2
 library(gridExtra)






library(ashr)
tt <- ash(betahat, rep( 1, length(betahat)))
plot( tt$result$lfdr, fit $post_assignment[,1])



res <- post_mean_sd.mococomo(fit,data)

plot(  betahat,res$mean)
abline(a=0,b=1)
points(betahat, tt$result$PosteriorMean, col="green")




sum( mix[which( tt$result$lfdr<0.05)])/length(which( tt$result$lfdr<0.05))
sum( mix[which(fit $post_assignment[,1]<0.05)])/length(which(fit $post_assignment[,1]<0.05))


library(ggplot2)

get_all_cs( fit$logreg_list[[1]] )

 #get posterior quantities
 est<- post_mean_sd.mococomo (fit)
 head(est)
  plot( est$mean, data$betahat)

 #comparison with ash

 t_ash <- ash(sim $betahat, sim $se, mixcompdist = "normal")
 post_mean_ash <- t_ash$result$PosteriorMean
 plot(est$mean, post_mean_ash)
 # TODO make a more convincing example

  sim  <- logisticsusie:::sim_mococomo_beta(n=100)
#preparing the data
data <- set_data_mococomo(p = sim$p,
                          X = sim$X)

 fit <- fit.mococomo(data, maxiter=20)





 tt<- compute_posterior_assignment(fit, log = F)
boxplot(tt[,1]~as.factor(mix))
