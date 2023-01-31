set.seed(1)
library(logisticsusie)
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
fit$elbo <- compute_elbo3.mococomo(fit)

i=1


  fit <- iter.mococomo(fit, is.even(i), is.odd(i))
  str(fit)
i=2

fit <- iter.mococomo(fit, is.even(i), is.odd(i))
 # iter mococomo ------


  K <- fit$K

  # updates posterior assignments

    fit$post_assignment <- compute_posterior_assignment(fit, log = F)
    fit$N <- .expected_trials(fit$post_assignment)


  # pass assignments to logreg
  for (k in seq(K - 1)) {
    fit$logreg_list[[k]]$data$y <- fit$post_assignment[, k]
    fit$logreg_list[[k]]$data$N <- fit$N[, k]
  }

    fit <- iter.mnsusie(fit, from_init_moco = TRUE)
    ## iter.mnsusie ----
      fit_intercept = T
      fit_prior_variance = T
      from_init_moco = TRUE



        K <- length(fit$logreg_list) + 1

        logreg_list <- list()
        for (k in seq(K - 1)) {
          logreg <- fit$logreg_list[[k]]

          # update logreg
          fit$logreg_list[[k]] <- iter.binsusie(fit=fit,
                                                k=k,
                                                fit_intercept = fit_intercept,
                                                fit_prior_variance = fit_prior_variance
          )

        }
        ### iter.binsusie
        fit_intercept = TRUE
        fit_prior_variance = TRUE
        fit_xi = TRUE
        fit_alpha = TRUE



            # update b
            fit$logreg_list[[k]] <- update_b.binsusie(fit,
                                                      k                  = k,
                                                      fit_intercept      = fit_intercept,
                                                      fit_prior_variance = fit_prior_variance,
                                                      fit_alpha          = fit_alpha
            )
                 ####update_b.binsusie
                 fit_intercept = TRUE
                 fit_prior_variance = TRUE
                 fit_alpha = TRUE
                 update_idx = NULL

                  update_idx <- seq(fit$logreg_list[[k]]$hypers$L)

                shift <- compute_Xb.binsusie(fit, k = k)
                for (l in update_idx) {
                  # remove current effect estimate
                  shift <- Matrix::drop(shift - compute_Xb.binser(fit, k = k, idx = l))

                  # update SER
                  post_l <- update_b.binser(fit, k = k, idx = l, shift = shift)
                     #### update_b.binser
                  function(fit, idx = NULL, kidx = NULL, shift = 0
                    nu <- Matrix::drop(compute_nu.binser(fit, idx, kidx, shift))
                    tau0 <- 1 / .get_var0(fit, idx)
                    tau <- tau0 + fit$params$tau
                    post <- .update_b_ser(nu, tau, .get_pi(fit, idx))

                  fit$params$mu[l, ] <- post_l$mu
                  fit$params$var[l, ] <- post_l$var

                  if (fit_alpha) {
                    fit$params$alpha[l, ] <- post_l$alpha
                  }

                  # update intercept/fixed effect covariates
                  if (fit_intercept) {
                    fit$params$delta[l, ] <- update_delta.binser(fit, k = k, idx = l, shift = shift)
                  }

                  # update prior_variance
                  if (fit_prior_variance) {
                    fit$hypers$prior_variance[l] <- update_prior_variance.binser(fit, k = k, idx = l)
                  }

                  # add current effect estimate
                  shift <- shift + compute_Xb.binser(fit, k = k, idx = l)
                }




              fit$logreg_list[[k]]$params$xi <- update_xi.binsusie(fit, k)
              fit$logreg_list[[k]]$params$tau <- compute_tau(fit, k)






