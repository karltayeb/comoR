example_fit <- function(){
  sim <- logisticsusie::sim_ser()
  betahat <- rnorm(length(sim$y))
  betahat[sim$y==1] <- rnorm(sum(sim$y == 1), sd=10)
  se <- rep(1, length(sim$y))

  data <- prep_data_como2(betahat, se, sim$X, sim$Z)
  fit <- data_initialize_como(data, 5, scales = c(0, 1, 5, 10))
  fit <- fit_model(fit, data, max_iter = 1000)

  assignments <- purrr::map_int(1:nrow(fit$post_assignment), ~which.max(fit$post_assignment[.x,]))
  table(sim$y, assignments)
  fit <- update_model(fit, data, fit_prior_variance=F)
  table(fit$logits > 0, sim$y)
}


example_fit_como <- function(){
  sim <- logisticsusie::sim_ser()
  betahat <- rnorm(length(sim$y))
  betahat[sim$y==1] <- rnorm(sum(sim$y == 1), sd=5)
  se <- rep(1, length(sim$y))

  data <- prep_data_como2(betahat, se, sim$X, sim$Z)
  fit <- data_initialize_como(data, 5, scales = c(0, 1, 5, 10))

  fit <- update_model(fit, data, fit_prior_variance=F)
  fit <- fit_model(fit, data)
  table(fit$logits > 0, sim$y)
}
