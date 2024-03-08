



rm(list=ls())
library(comoR)
library(nnet)
library(keras3)
library(comoR)
library(tensorflow)
set.seed(1)
y <-runif(2000, max=2)
X = cbind(y)
xtrue = 0*y
for (i in 1:length(xtrue)){

  if ( (    y[i] <.5 & y[i] >.2 ) | (y[i] >1.5 & y[i] <1.8) ){
    xtrue[i] = 0
  }else{
    xtrue[i]= rnorm(1,sd=0.5+ 2*abs(sin( pi*y[i])))
  }


}
plot (y,xtrue)
x = xtrue + rnorm(length(xtrue), sd=1)
plot (y,x)
s= rep(1,length(x))
Z <- matrix( 1, nrow=length(x), ncol=1)
#start by defining the como parameter using mixture of exponential priors and a neural net regressor
param_como = list(max_class= 10,
                  mnreg_type="keras",
                  prior ='mix_norm'# "mix_exp"
)
data <- comoR:::como_prep_data (betahat=x,
                                se=s, X=X,
                                Z =Z )


num_classes <- length( autoselect_scales_mix_exp(data$betahat, data$se,10))

#define the nnet paramet using Keras syntax
param_nnet =keras_model_sequential() %>%
  layer_dense(units = 64, activation = 'relu', input_shape = c(1))%>%
  layer_dense(units = 64, activation = 'relu' )%>%
  layer_dense(units = 64, activation = 'relu' )%>%
  layer_dense(units = 64, activation = 'relu' )%>%
 layer_dense(units = 10,#important to have the same number of units as the number of classes
              activation = 'softmax')


reticulate::py_last_error()





fit_como  <- rlang::exec( "data_initialize_como", !!! param_como ,
                          data= data,
                          param_nnet= param_nnet) # initialize the model from the data
fit =fit_como

#fit$data_loglik <- compute_data_loglikelihood(fit, data)
#hist(fit$data_loglik)

fit_como <- comoR:::fit.como ( fit_como, data, max_iter = 40 )
fit =fit_como

est <- comoR:::post_mean_sd (fit,data)

par(mfrow=c(1,1))
plot(est$mean, xtrue ,   col =ifelse(xtrue==0, 1,2))



lol <- ashr::ash(x, s )

plot(lol$result$PosteriorMean, xtrue,  col =ifelse(xtrue==0, 3,4))
cor (lol$result$PosteriorMean, xtrue)


plot(lol$result$PosteriorMean,
     est$mean,
     col =ifelse(xtrue==0, 3,4),
     main='comparison of ebnm vs cebnm',xlab='ebnm',ylab='cebnm')

rmse = function(x,y){
  sqrt(mean (x-y)^2)
}
rmse(lol$result$PosteriorMean, xtrue)

rmse(est$mean, xtrue )
