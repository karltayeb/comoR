
prep_out_FDR_wrapper <- function(fit){
  if(dist=="normal"){
    est    <- post_mean_sd.mococomo (fit)
    lfdr   <- get_lfdr_mixnorm(fit)
    qvalue <- cal_qvalue(lfdr)
    resdf <- data.frame(betahat       = fit$data$betahat,
                        se            = fit$data$se,
                        lfdr          = lfdr,
                        qvalue        = qvalue,
                        PosteriorMean = est[,1],
                        PosteriorSD   = est[,2]
                        )
  }


  ### TODO: remove dummy CS
  cs <- lapply( 1:length(fit$logreg_list),
                function(k)
                 get_all_cs(fit$logreg_list[[k]])
                )
  fitted_effect <-  lapply(1:length(fit$logreg_list[[k]]),
                          function(k)
                            fit$logreg_list[[k]]$params$alpha*fit$logreg_list[[k]]$params$mu
                         )

  out <- list(result =resdf,
              cs=cs,
              fitted_effect)
}



get_lfdr_mixnorm <- function(fit){
  tt1 <-  fit$post_assignment[,1]* dnorm( fit$data$betahat, mean=0, sd= fit$data$se )
  tt2 <-    Reduce("+",lapply( 2: ncol(fit$post_assignment),
                               function(k)  fit$post_assignment[,k]*
                                            dnorm( fit$data$betahat,
                                                   mean=0,
                                                   sd= sqrt( fit$data$se^2+ fit$f_list[[k]]$var ) )))
  out <- tt1/(tt1+tt2)

  return(out)
}

cal_qvalue <- function(lfdr)
{

  torder <- order(lfdr)
  qvalue <- sapply( 1:length(lfdr), function(k )
                                    mean(lfdr[which(torder <= torder[k] )])
                    )
  return(qvalue)

}


assesor.beta <- function(fit){}
