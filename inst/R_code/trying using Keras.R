
library(keras)
library(comoR)
library(nnet)
devtools::load_all(".")
set.seed(1)
y <-runif(10000)
X = cbind(y)
xtrue = 0*y
for (i in 1:length(xtrue)){

  if ( (  y[i] <.5 ) ){
    xtrue[i] = rnorm(1)
  }else{
    xtrue[i] = 0
  }


}
plot (y,xtrue)
x = xtrue + rnorm(length(xtrue), sd=0.5)
plot (y,x)
s= rep(0.5,length(x))
Z <- matrix( 1, nrow=length(x), ncol=1)

param_como = list(max_class=10,mnreg_type="mult_reg")
param_nnet =list( size=3, decay=1,MaxNWts = 10000)

data <- comoR:::prep_data_como2 (betahat=x,
                                 se=s, X=X,
                                 Z =Z )
fit <- rlang::exec( "data_initialize_como", !!! param_como ,
                    data= data,
                    param_nnet= param_nnet) # initialize the model from the data


fit$data_loglik <- compute_data_loglikelihood(fit, data)
fit$data_loglik






est <- comoR:::post_mean_sd (fit,data)




x_train=y
y_train=fit$data_loglik
num_classes = ncol(y_train)
model <- keras_model_sequential() %>%
  layer_dense(units = 64, activation = 'relu', input_shape = c(1)) %>%
  layer_dense(units = 64, activation = 'relu' ) %>%

  layer_dense(units = num_classes, activation = 'softmax')


custom_loss <- function( y_true , y_pred  ) {
  # Define your custom loss computation here
  # For example, mean squared error
    tt  = 0
    for (i in 1: nrow(y_true))
    {
      tt <- log ( sum( exp( y_true[i,])*y_pred[i,]))
    }

  mse <-    -tt
  return(mse)
}

blank_model <-  model %>% compile(
  loss = custom_loss,
  optimizer = 'adam',
  metrics = c('accuracy')
)

model1 <- blank_model

history <-blank_model %>% fit(
  x_train, y_train,
  epochs = 100,
  batch_size = 100
)

tt = predict(model, x_train)
image (tt)
plot(tt[,1], y, col =ifelse(y<0.5, 1,2))

model2 <- blank_model

x_train=y[-1]
y_train=fit$data_loglik[-1,]
history <- blank_model %>% fit(
  x_train, y_train,
  epochs = 100,
  batch_size = 100
)



update_model.como
function(fit, data, update_assignment = T, update_logreg=T, fit_prior_variance=F, track_elbo=T){
  K <- fit$K

  # pre-compute data likelihood, if we haven't already
  # TODO: use digest::digest to hash the scales and recompute if scales change?
  if(is.null(fit$data_loglik)){
    fit$data_loglik <- compute_data_loglikelihood(fit, data)
  }

  # updates posterior assignments probabilities
  # these are the response variable for mn_regression
  if (update_assignment) {
    fit$post_assignment <- compute_posterior_assignment(fit, data)
    #data$Y <- fit$post_assignment
    #data$Nk <- logisticsusie:::stick_breaking_N(data$Y)
  }

  if (update_logreg) {
    fit$mnreg <- update_prior(fit$mnreg, resps=fit$post_assignment, data=data)
  }

  if (track_elbo){
    fit$elbo <- c(fit$elbo, compute_elbo(fit, data))
  }
  return(fit)
}
