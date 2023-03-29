
example_fit <- function(){
  sim <- logisticsusie::sim_ser()
  betahat <- rnorm(length(sim$y))
  betahat[sim$y==1] <- rnorm(sum(sim$y == 1), sd=10)
  se <- rep(1, length(sim$y))

  data <- como2_prep_data(betahat, se, X, Z)
  fit <- data_initialize_como2(data, L=5)

  fit <- update_model(fit, data, estimate_f1 = F)
  fit <- fit_model(fit, data, estimate_f1=T)
  fit$f1$var
  table(fit$logits > 0, sim$y)
}
