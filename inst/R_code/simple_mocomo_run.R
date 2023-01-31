set.seed(1)

sim  <- sim_mococomo(n=1000)
#preparing the data
data <- set_data_mococomo(betahat  = sim$betahat,
                          X = sim$X,
                          se= sim$se)
#fit mococomo model
maxiter   = 100
tol       = 1e-3
max_class = 10
mult      = 2


#working init
fit <- init.mococomo(data,
                     max_class = max_class,
                     mult   = mult,
                     upper     = upper
)
str(fit$logreg_list)
  compute_elbo.mococomo(fit)


 #works with two distribution (left and right)
  tfit <- fit.mococomo  (data )

  tfit$elbo
  tfit$f_list


res <-   cFDR( betahat = data$betahat,
        se      = data$se,
        X       = data$X)

tt <- ash(data$betahat ,data$se ,outputlevel=3, mixcompdist = "normal")

plot( res$result$betahat, res$result$PosteriorMean)
points(tt$result$betahat,tt$result$PosteriorMean, col="green")

plot(res$result$PosteriorMean,tt$result$PosteriorMean)

plot(res$full_obj$post_assignment[,5],
     tt$result$lfdr)
abline(a=0,b=1)

param_var <- Reduce("c", lapply(1: length(res$full_obj$f_list), function(k)res$full_obj$f_list[[k]]$var))
tt$fitted_g$sd
sqrt(param_var)






tt <- ash(data$betahat ,data$se ,outputlevel=3, mixcompdist = "normal")
tt$fitted_g$sd
post_asign <- list()
for ( i in 1:length(data$betahat)){

  vec_t <-  tt$fitted_g$pi *dnorm( fit$data$betahat[i],
                                   mean=0,
                                   sd=   sqrt( fit$data$se[i]^2+ tt$fitted_g$sd^2))

  temp  <- (tt$fitted_g$pi *dnorm( fit$data$betahat[i],
                                   mean=0,
                                   sd=  sqrt( fit$data$se[i]^2+ tt$fitted_g$sd^2)))/ (sum(vec_t))
  post_asign[[i]] <-  temp
}
post_asign <- do.call(rbind, post_asign)
plot( post_asign[,1], tt$result$lfdr)


cor(post_asign)
cor(res$full_obj$post_assignment)
