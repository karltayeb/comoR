#rm(list = ls())
set.seed(1)
library(logisticsusie)  #Simulate data under the mococomo model
library(comoR)

#### Loadings sampling -----
N=200
x1 <- rnorm(N,sd=3)
beta0=-2
beta1=2
samp_prob <- 1/(1 +exp(-(beta0+beta1*x1)))
P=20
mix <- c()
betahat <- c()
betatrue <- c()

X <- cbind( x1, matrix(rnorm(P*N), ncol=P))
se <- rep(1, length(betahat))


for ( i in 1:N){
  mix <-c(mix, sample(c(0,1), size=1, prob = c(1- samp_prob[i], samp_prob[i])))
  betatrue <- c(betatrue, mix[i]*rnorm(1,sd=3))
  betahat <- c( betahat ,   betatrue[i]+rnorm(1,sd=1 ) )
}
beta11 <- betahat
betatrue11 <- betatrue

#p <- runif(N)
#preparing the data

data11 <- prep_data_como2(betahat =  betahat,
                           se = se ,
                           X  =  X,
                          Z= rep( 1, length(se)))


samp_prob <- 1/(1 +exp(-(beta0+beta1*X[,2])))
mix <- c()
betahat <- c()
betatrue <- c()

for ( i in 1:N){
  mix <-c(mix, sample(c(0,1), size=1, prob = c(1- samp_prob[i], samp_prob[i])))
  betatrue <- c(betatrue, mix[i]*rnorm(1,sd=3))
  betahat <- c( betahat ,   betatrue[i]+rnorm(1,sd=1 ) )
}
betahat12 <- betahat
betatrue12 <- betatrue



data12 <- prep_data_como2(betahat =  betahat,
                          se = se ,
                          X  =  X,
                          Z= rep( 1, length(se)))


#### factors sampling -----
N=100
P=20
y1 <- rnorm(N,sd=3)

Y <- cbind( x1, matrix(rnorm(P*N), ncol=P))
beta0=-2
beta1=1
samp_prob <- 1/(1 +exp(-(beta0+beta1*y1)))

mix <- c()
betahat <- c()
betatrue <- c()
for ( i in 1:N){
  mix <-c(mix, sample(c(0,1), size=1, prob = c(1- samp_prob[i], samp_prob[i])))
  betatrue <- c(betatrue, mix[i]*rnorm(1,sd=3))
  betahat <- c( betahat ,   betatrue[i]+rnorm(1,sd=1 ) )
}

betatrue21 <- betatrue
betahat21 <- betahat

se <- rep(1, length(betahat))

#preparing the data

#p <- runif(N)
data21 <- prep_data_como2(betahat =  betahat,
                          se = se ,
                          X  =  Y,
                          Z= rep( 1, length(se)))


samp_prob <- 1/(1 +exp(-(beta0+beta1*Y[,2])))

mix <- c()
betahat <- c()
betatrue <- c()
for ( i in 1:N){
  mix <-c(mix, sample(c(0,1), size=1, prob = c(1- samp_prob[i], samp_prob[i])))
  betatrue <- c(betatrue, mix[i]*rnorm(1,sd=3))
  betahat <- c( betahat ,   betatrue[i]+rnorm(1,sd=1 ) )
}

betatrue22 <- betatrue
betahat22 <- betahat

se <- rep(1, length(betahat))

#preparing the data

#p <- runif(N)
data22 <- prep_data_como2(betahat =  betahat,
                          se = se ,
                          X  =  Y,
                          Z= rep( 1, length(se)))


Y_true <- betatrue11%*%t(betatrue21)+ betatrue12%*%t(betatrue22)
Y <- betatrue11%*%t(betatrue21) +betatrue12%*%t(betatrue22)+ matrix(rnorm(100*100,sd=4), ncol=100)

plot( Y,Y_true)
X_l =data11$X

X_f =data22$X
image(Y)
image(Y_true)



K=2
dim(Y)
library(softImpute)
cEBMF.obj <- init_cEBMF (Y, X_l,X_f,K=2, init_type = "udv_si")
plot(cEBMF.obj$loading[,1] ,data11$betahat)



for ( o in 1:4){
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
                      max_iter   = 20,
                      nullweight = 0.1)


    cEBMF.obj$factor[,k]  <-   t_fit$result$mean
    cEBMF.obj$factor2[,k] <-   t_fit$result$sd^2+  t_fit$result$mean^2


    cEBMF.obj<- update_tau.cEBMF (cEBMF.obj )

    Y_est <- cEBMF.obj$loading[,1]%*%t(cEBMF.obj$factor[,1])
    plot( Y_est, Y )
    abline(a=0,b=1)
  }

}


Y_est <- cEBMF.obj$loading[,1]%*%t(cEBMF.obj$factor[,1])+cEBMF.obj$loading[,2]%*%t(cEBMF.obj$factor[,2])
library(flashier)
f <- flash(Y)
plot( Y_est, Y )
plot( Y_true, Y )
plot( Y_true, Y_est )
points(Y_true, fitted(f), col="green")

sqrt(sum( (Y_true-Y_est)^2))
sqrt(sum( (Y_true- fitted(f))^2))





plot(cEBMF.obj$loading[,2], betatrue12 )


plot(cEBMF.obj$loading[,1], betatrue11 )

plot(cEBMF.obj$factor[,1], betatrue21 )
plot(cEBMF.obj$factor[,2], betatrue22 )


