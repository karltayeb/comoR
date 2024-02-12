rm(list=ls())
library(keras)

y = rnorm(1000)
x= rnorm (1000)
num_classes=10
mat <- matrix( rnorm (1000*num_classes), ncol=num_classes)


x_train =x
y_train =mat
length(x_train)
dim(y_train)

model <- keras_model_sequential() %>%
  layer_dense(units = 64, activation = 'relu', input_shape = c(1)) %>%
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

blank_model <-  model %>% compile(
  loss = custom_loss,
  optimizer = 'adam',
  metrics = c('accuracy')
)



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


model2 <- keras_model_sequential() %>%
  layer_dense(units = 64, activation = 'relu', input_shape = c(1)) %>%
  layer_dense(units = 64, activation = 'relu' ) %>%
    layer_dense(units = num_classes, activation = 'softmax')


  model2 <-  model2 %>% compile(
    loss = custom_loss,
    optimizer = 'adam',
    metrics = c('accuracy')
  )




x_train2 =x[-1]
y_train2 =mat[-1,]


history <-model2 %>% fit(
  x_train2, y_train2,
  epochs = 40,
  batch_size = 111
)
