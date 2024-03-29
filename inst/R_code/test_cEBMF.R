#rm(list = ls())
devtools::load_all(".")
library(logisticsusie)  #Simulate data under the como model
library(comoR)

library(softImpute)
sim11  <- sim_twococomo(n=1000)#contains all the info


#preparing the data



data11 <-set_data_como(betahat = sim11$betahat,
                       se = sim11$se ,
                       X  = sim11$X) # prepare the data

 fit <- data_initialize_como(data=data11, max_class=5, scales = c(0, 1, 5, 10)) # initialize the model from the data
fit <- fit_model(fit, data11, max_iter = 10)




sim12  <- sim_twococomo(n=1000)

data12 <-set_data_como(betahat = sim12$betahat,
                           se = sim12$se ,
                           X  = sim12$X)
sim21  <- sim_twococomo(20)
#preparing the data
data21 <-set_data_como(betahat = sim21$betahat,
                           se = sim21$se ,
                           X  = sim21$X)
sim22  <- sim_twococomo(20)
data22 <-set_data_como(betahat = sim22$betahat,
                           se = sim22$se ,
                           X  = sim22$X)

Y_true <- sim11$beta%*%t(sim21$beta)+ sim12$beta%*%t(sim22$beta)
Y <- data11$betahat%*%t(data21$betahat) +data12$betahat%*%t(data22$betahat)

plot( Y,Y_true)
X_l =data11$X

X_f =data22$X
image(Y)
image(Y_true)



K=2
dim(Y)
cEBMF.obj <- init_cEBMF (Y, X_l,X_f,K=2, init_type = "udv_si",
                         param_como  = list(max_class=10,mnreg_type="mult_reg"),
                         param_nnet  =list( size=1, decay=1),
                         param_como2 = list(),
                         param_susie =  list(L=5),
                         maxit_como  = 10)
plot(cEBMF.obj$loading[,1] ,data11$betahat)



for ( o in 1:20){
  for ( k in 1:K){
    Rk <- cal_partial_residuals.cEBMF(cEBMF.obj,k)
    l_k <- cal_expected_loading( cEBMF.obj, Rk,k)
    t_data <- set_data_como(betahat = l_k$l_i_hat,
                                se      = l_k$s_i,
                                X       = cEBMF.obj$X_l )

    t_fit <- data_initialize_como(t_data, max_class=5, scales = c(0, 1, 5, 10)) # initialize the model from the data
    t_fit <- fit_model( t_fit, t_data, max_iter = 10)


    fitted_loading <- post_mean_sd.como (fit= t_fit, data=t_data )
    cEBMF.obj$loading[,k] <-  fitted_loading$mean
    cEBMF.obj$loading2[,k] <- fitted_loading$sd^2+ fitted_loading$mean^2

    #factor update

    f_k <- cal_expected_factor( cEBMF.obj, Rk,k)
    t_data <- set_data_como(betahat    = f_k$f_j_hat,
                                se      = f_k$s_j,
                                X       = cEBMF.obj$X_f )
    t_fit <- data_initialize_como(t_data, max_class=5, scales = c(0, 1, 5, 10)) # initialize the model from the data
    t_fit <- fit_model( t_fit, t_data, max_iter = 10)


    fitted_factor <- post_mean_sd.como (fit= t_fit, data=t_data )
    cEBMF.obj$factor[,k] <-  fitted_factor $mean
    cEBMF.obj$factor2[,k] <-  fitted_factor$sd^2+ fitted_factor$mean^2


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
cor( c(Y_true), c(Y_est))
cor( c(Y_true), c(Y))


cor( c(Y_true), c(fitted(f)))







####Hand run update----

library(flashr)
library(softImpute)
ftrue = matrix(rnorm(200), ncol=2)
ltrue = matrix(rnorm(40), ncol=2)
ltrue[1:10, 1] = 0 # set up some sparsity
ltrue[11:20, 2] = 0
Y = ltrue %*% t(ftrue) + rnorm(2000) # set up a simulated matrix
f = flash(Y, K=3 ,#ebnm_fn= 'ebnm_ash' ,
          var_type = "constant")
ldf = f$ldf
Y_true <- ltrue %*% t(ftrue)
X_l <- matrix(rnorm(nrow(Y)*10), nrow= nrow(Y))
X_f <- matrix(rnorm(ncol(Y)*10), nrow= ncol(Y))

cEBMF.obj <- cEBMF (Y, X_l,X_f,K=3, init_type = "udv_si")

cEBMF.obj$elbo


cEBMF.obj <- init_cEBMF (Y, X_l,X_f,K=3, init_type = "udv_si")

for (i in 1:10) {
  cEBMF.obj <- cEBMF_iter  (cEBMF.obj)
  # print(i)

  print(cbind(c(cEBMF.obj$KL_f),c(cEBMF.obj$KL_l)))
  print(mean(cEBMF.obj$tau))
  print(cEBMF.obj$elbo)

}


library(flashier)
Y_est <- Reduce("+", lapply( 1:cEBMF.obj$K, function(k) cEBMF.obj$loading[,k]%*%t(cEBMF.obj$factor[,k]) ))
f = f <- flashier::flash(Y_obs)
plot( cEBMF.obj$elbo)
plot( Y_est,Y)
points( Y_est,Y_true, col="green")

plot( fitted(f),Y_true)
points( Y_est,Y_true, col="green")
plot( Y_est,fitted(f))

library(nnet)
cEBMF.fit <- cEBMF (Y, X_l,X_f,K=1, init_type = "udv_si")
cEBMF.fit$elbo

cEBMF.obj <- init_cEBMF (Y, X_l,X_f,K=1, init_type = "udv_si")

for (i in 1:3) {
  cEBMF.obj <- cEBMF_iter  (cEBMF.obj)
  # print(i)

  print(cbind(c(cEBMF.obj$KL_f),c(cEBMF.obj$KL_l)))
  print(mean(cEBMF.obj$tau))
  print(cEBMF.obj$elbo)

}
Y_est <- Reduce("+", lapply( 1:cEBMF.obj$K, function(k) cEBMF.obj$loading %*%t(cEBMF.obj$factor ) ))
Y_fit <- Reduce("+", lapply( 1:cEBMF.fit $K, function(k) cEBMF.fit $loading %*%t(cEBMF.fit $factor ) ))
plot( Y_est, cEBMF.fit$Y_fit)
abline(a=0,b=1)
 f <- flashier::flash(Y )
plot( cEBMF.obj$elbo)
plot( Y_est,Y)
points( Y_est,Y_true, col="green")

plot( fitted(f),Y_true)
points(Y_true, Y_est, col="green")
plot( Y_est,fitted(f))
