
como2_linear_susie <- function(){
  sim <- logisticsusie::sim_ser()
  betahat <- rnorm(length(sim$y))
  betahat[sim$y==1] <- rnorm(sum(sim$y == 1), sd=2)
  se <- rep(1, length(sim$y))

  data <- prep_data_como2(betahat, se, sim$X, sim$Z)
  fit <- data_initialize_como2(data, f1_params = list(mu=0, var = 4), logreg='linear_susie')
  fit <- fit_model(fit, data)
}


william_example <- function(){
  set.seed(2)
  P=20
  N <- 5000
  beta0 <- 0
  beta1 <-0.41
  x1 <- rnorm(N)
  samp_prob <- 1/(1 +exp(-(beta0+beta1*x1)))

  mix <- c()
  obs <- c()
  se <- runif(N)

  for (i in 1:N){
    mix <-c(mix, sample(c(0,1), size=1, prob = c(1- samp_prob[i], samp_prob[i])))
    obs <- c(obs, ifelse( mix[i]==1, rnorm(1, sd=sqrt(2+se[i])),rnorm(1, sd=se[i]) ))
  }
  X <- cbind( x1, matrix(rnorm(P*N), ncol=P))

  boxplot(obs~mix)

  data <- prep_data_como2(obs, se, X, Z = matrix(rep(1, N), nrow=N))
  fit <- data_initialize_como2(data, L=5)

  fit <- update_model(fit, data, estimate_f1 = F)
  fit <- fit_model(fit, data, estimate_f1=T)

  res <- cFDR( betahat=obs,
               se=se,
               X = X,
               n_sim= 1000, nullweight = 2.3,
               outputlevel = 2 )

  res$cs

  #or via mococomo

  model <- "normal"
  data <- set_data_mococomo(betahat = obs,
                            X = X,
                            se= se)
  fit <- fit.mococomo(data, maxiter=20)
  cs <- lapply( 1:length(fit$logreg_list),
                function(k)
                  get_all_cs(fit$logreg_list[[k]])
  )
  cs
}
