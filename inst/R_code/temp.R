set.seed(1)
library(logisticsusie)
sim  <- sim_mococomo(n=10000)
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
fit$elbo <- compute_elbo3.mococomo(fit)

i=1
backfit=TRUE
tt <- proc.time()
fit <- iter.mococomo(fit, is.even(i), is.odd(i))
tt -proc.time()
  #str(fit)

  cs <- lapply( 1:length(fit$logreg_list),
                function(k)
                  get_all_cs(fit$logreg_list[[k]])
  )


#  str( fit$logreg_list[[1]])
 fit$logreg_list[[1]]$hypers$L <-1


 fit$logreg_list[[1]]$params$alpha <- fit$logreg_list[[1]]$params$alpha[1:2,]
 fit$logreg_list[[1]]$params$mu    <- fit$logreg_list[[1]]$params$mu   [1:2,]
 fit$logreg_list[[1]]$params$var   <- fit$logreg_list[[1]]$params$var  [1:2,]
 fit$logreg_list[[1]]$params$delta <- matrix(fit$logreg_list[[1]]$params$delta[1:2,],ncol = 1)
i=2
tt <- proc.time()
fit <- iter.mococomo(fit, is.even(i), is.odd(i))
tt -proc.time()
 # iter mococomo ------
#str( fit$logreg_list[[1]])
