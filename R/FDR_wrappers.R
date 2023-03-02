#' @param fit a fitted mococomo object
#' @param outputlevel	 Determines amount of output. outputlevel=1 provide simplest output (which should be enough for most applications)
#' outputlevel= also output the entire mococomo fitted object
#'
#'
 # TODO work on presentation of the output
prep_out_FDR_wrapper <- function(fit, outputlevel=1,n_sim ){
  if(fit$model=="normal"){
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

  if(fit$model=="beta"){

    # case where mixture of beta with increasing and decreasing value

      lfdr   <- get_lfdr_beta(fit)
      qvalue <- cal_qvalue(lfdr)
      FDR <- assesor.mococomo_beta(fit,n_sim=n_sim )
      if("p"%in% names(fit$data)){
        resdf <- data.frame(p       = fit$data$p,
                            lfdr    = lfdr,
                            qvalue  = qvalue,
                            FDR     = FDR
        )
      }
      if("zscore" %in% names(fit$data) ){
        resdf <- data.frame(zscore  = fit$data$zscore,
                            lfdr    = lfdr,
                            qvalue  = qvalue,
                            FDR     = FDR
        )
      }



  }


  ### TODO: remove dummy CS
  cs <- lapply( 1:length(fit$logreg_list),
                function(k)
                 get_all_cs(fit$logreg_list[[k]])
                )
  fitted_effect <-  lapply(1:length(fit$logreg_list),
                          function(k)
                            fit$logreg_list[[k]]$params$alpha*fit$logreg_list[[k]]$params$mu
                         )
  if(outputlevel==1){
    out <- list(result =resdf,
                cs=cs,
                fitted_effect=fitted_effect)
  }
  if(outputlevel==2){
    out <- list(result =resdf,
                cs=cs,
                fitted_effect=fitted_effect,
                full_obj = fit)
  }
  return(out)

}



get_lfdr_mixnorm <- function(fit){
  tt1 <-  fit$post_assignment[,1]
  out <- tt1

  return(out)
}


get_lfdr_beta <- function(fit){


  tt1 <-  fit$post_assignment[,1]
  out <- tt1
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


assesor.mococomo_beta <- function(fit,n_sim){

up <-  fit$post_assignment[, 1]
low <-  up + apply(fit$post_assignment[,-1]*compute_data_loglikelihood(fit, fit$data)[,-1],1,sum)

obs_assesor <- up/low #observed assessor



#simualted null
sim_null <- list()
tdata <- fit$data


temp_f <- function(i){
  tdata$p  <-  rep(runif(1), (nrow(fit$data$X)))
  low_sim <-  up + apply(fit$post_assignment[,-1]*exp(compute_data_loglikelihood(fit,tdata)[,-1]),1,sum)
  return(up/low_sim)
}





 H0_ind <- do.call(cbind, lapply(1:n_sim, function(i) temp_f(1)))




for( i in 1:nrow(H0_ind)){
  H0_ind[i,] <-  H0_ind[i,] [order( H0_ind[i,] )]
}



BC_est <- c()
 for( i in 1:length(obs_assesor)){
   if( length(which(H0_ind[i,]< obs_assesor[i]))==0)
   {
     obs_prob <- 1
   }else{
     obs_prob <-  ( length(which(H0_ind[i,]< obs_assesor[i]))+1)/length(H0_ind[i,])
   }
   BC_est <-  c(BC_est, obs_prob)

}




 # cutpt_max_BC <- max( blfdr_vec[FDR_est<= alpha])  # this can be -Inf
# rej_index_blfdr_BC  <- which(  blfdr_vec <=    cutpt_max_BC ) # this can be empty

return(BC_est)
}


