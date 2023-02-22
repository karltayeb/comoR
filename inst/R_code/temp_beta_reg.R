set.seed(1)
library(comoR)
 sim  <- logisticsusie:::sim_mococomo_beta(n=1000)
 #preparing the data
 data <- set_data_mococomo(p = sim$p,
                            X = sim$X)
 #fit mococomo model
 maxiter   = 100
 tol       = 1e-3
 max_class = 10
 mult      = 2
 upper     = TRUE


#working init
 fit <- init.mococomo(data,
                      max_class = max_class,
                      model= "beta",
                      mult   = mult,
                      upper     = upper,
                      nullweight = 1000
                      )

 str(fit)

 fit$N
 #working elbo
fit$elbo <- compute_elbo.mococomo(fit)

fit$elbo <- compute_elbo2.mococomo(fit)

fit$elbo <- compute_elbo3.mococomo(fit)

##working exemple
for (i in 1:maxiter) {
  fit <- iter.mococomo(fit, is.even(i), is.odd(i))
  fit$elbo <- c(fit$elbo, compute_elbo3.mococomo(fit))

  # print(paste('asgn:', is.even(i), 'logreg:', is.odd(i), 'elbo: ', tail(fit$elbo, 1)))
  if (.converged(fit, tol)) {
    break
  }
}

#works with two distribution (left and right)
tfit <- fit.mococomo  (data , nullweight = 1)

tfit$elbo
tfit$f_list
tfit$N
#works with two distribution (left and right)
tfit <- fit.mococomo  (data,model="beta",
                       upper=TRUE)

## We should see almost no weight on the upper distribution

 tfit$elbo
 tfit$f_list

 tfit$post_assignment



 #works with a signal left distribution
tfit <- fit.mococomo  (data,
                       model="beta",
                       upper=FALSE)
plot(data$p, tfit$post_assignment[,1]  )



tfit <- fit.mococomo  (data,
                       model="beta",
                       upper=TRUE, nullweight = 4 , tol=0.1)
plot(data$p, tfit$post_assignment[,1]  )


plot(tfit$elbo)



plot(log10(diff(tfit$elbo)))

tfit$f_list

str(tfit)




res <-   cFDR( pvalue =    data$p,
               X       = data$X)
res$result
plot(res$result$lfdr, res$result$p  )

sim  <- sim_mococomo(n=1000)
#preparing the data
data <- set_data_mococomo(zscore  = sim$betahat/sim$se,
                          X = sim$X  )

res <-   cFDR( pvalue =    data$p,
               X       = data$X)
res$result
plot(res$result$lfdr, res$result$p  )
