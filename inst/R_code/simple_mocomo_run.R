set.seed(1)

sim  <- sim_mococomo(n=100)
#preparing the data
data <- set_data_mococomo(betahat  = sim$betahat,
                          X = sim$X,
                          se= sim$se)
#fit mococomo model
maxiter   = 100
tol       = 1e-3
max_class = 10
mult      = 2
upper     = TRUE


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


  cFDR( betahat = data$betahat,
        se      = data$se,
        X       = data$X)
