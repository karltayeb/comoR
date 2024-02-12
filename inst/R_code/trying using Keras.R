rm(list=ls())
library(keras)
library(comoR)
library(nnet)
devtools::load_all(".")
set.seed(1)
y <-runif(2000)
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
x = xtrue + rnorm(length(xtrue), sd=1)
plot (y,x)
s= rep(1,length(x))
Z <- matrix( 1, nrow=length(x), ncol=1)

param_como = list(max_class=10,mnreg_type="mult_reg")
param_nnet =list( size=3, decay=1,MaxNWts = 10000)

data <- comoR:::prep_data_como2 (betahat=x,
                                 se=s, X=X,
                                 Z =Z )
fit_como  <- rlang::exec( "data_initialize_como", !!! param_como ,
                    data= data,
                    param_nnet= param_nnet) # initialize the model from the data


fit_como $data_loglik <- compute_data_loglikelihood(fit_como , data)
fit_como $data_loglik


yb = y

x_train=yb
y_train=fit_como$data_loglik
num_classes = ncol(y_train)
model <- keras_model_sequential() %>%
  layer_dense(units = 64, activation = 'relu', input_shape = c(1)) %>%
  layer_dense(units = 64, activation = 'relu' ) %>%
  layer_dense(units = 64, activation = 'relu' ) %>%
  layer_dense(units = num_classes, activation = 'softmax')


custom_loss <- function(y_true, y_pred) {
  tt <- 0
  for (i in 1:nrow(y_true)) {
    tt <- tt + log(sum(exp(y_true[i,]) * y_pred[i,]))
  }
  mse <- -tt
  return(mse)
}

model1 <- clone_model(model )

model2 <- clone_model(model )

model1 <-  model1 %>% compile(
  loss = custom_loss,
  optimizer = 'adam',
  metrics = c('accuracy')
)


history <-model1 %>% fit(
  x_train, y_train,
  epochs = 40,
  batch_size = 100
)

tt1 = predict(model1, x_train)
image (tt1)
plot(tt1[,1], y, col =ifelse(y<0.5, 1,2))



model2 <-  model2 %>% compile(
  loss = custom_loss,
  optimizer = 'adam',
  metrics = c('accuracy')
)




x_train2=y[-c(1,2)]
y_train2=fit_como $data_loglik[-c(1,2),]

history <-model2  %>% fit(
  x_train2, y_train2,
  epochs = 100,
  batch_size = 111
)


tt1 = predict(model2, x_train2)

plot(tt1[,1], y[-c(1,2)], col =ifelse(y[-c(1,2)]<0.5, 1,2))

