#' @param fit a fitted mococomo object
#' @param outputlevel	 Determines amount of output. outputlevel=1 provide simplest output (which should be enough for most applications)
#' outputlevel= also output the entire mococomo fitted object
prep_out_FDR_wrapper <- function(fit, outputlevel=1){
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

      lfdr   <- fit$post_assignment[,1]
      qvalue <- cal_qvalue(lfdr)
      if("p"%in% names(fit$data)){
        resdf <- data.frame(p       = fit$data$p,
                            lfdr    = lfdr,
                            qvalue  = qvalue
        )
      }
      if("zscore" %in% names(fit$data) ){
        resdf <- data.frame(zscore       = fit$data$zscore,
                            lfdr    = lfdr,
                            qvalue  = qvalue
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
                fitted_effect)
  }
  if(outputlevel==2){
    out <- list(result =resdf,
                cs=cs,
                fitted_effect,
                full_obj = fit)
  }
  return(out)

}



get_lfdr_mixnorm <- function(fit){
  tt1 <-  fit$post_assignment[,1]* dnorm( fit$data$betahat, mean=0, sd= fit$data$se )
  tt2 <-    Reduce("+",
            lapply( 2: ncol(fit$post_assignment),
                   function(k)
                   fit$post_assignment[,k]*
                   dnorm( fit$data$betahat,
                          mean=0,
                          sd= sqrt( fit$data$se^2+ fit$f_list[[k]]$var )
                          )
                   )
            )
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


assesor.mococomo_beta <- function(fit,n_sim){

up <- 1- apply(fit$post_assignment[,-1],1,sum)
low <-  up + apply(fit$post_assignment[,-1]*compute_data_loglikelihood(fit, fit$data)[,-1],1,sum)

obs_assesor <- up/low #observed assessor



#simualted null
sim_null <- list()
tdata <- fit$data
for (i in 1:n_sim){#TODO speed up Monte Carlo null dist

  tdata$p  <-  runif(nrow(fit$data$X))
  low_sim <-  up + apply(fit$post_assignment[,-1]*compute_data_loglikelihood(fit,tdata)[,-1],1,sum)
  sim_null[[i]] <-up/low_sim#observed assessor

}
 H0_ind <- do.call(cbind, sim_null) #individual null distribution

for( i in 1:nrow(H0_ind)){
  H0_ind[i,] <-  H0_ind[i,] [order( H0_ind[i,] )]
}


 T_m <- rep(NA, length(obs_assesor))

 for ( i in 1:length(obs_assesor)){
   if( length(which(H0_ind[i,]< obs_assesor[i]))==0)
   {
     obs_prob <- 1
   }else{
     obs_prob <-  length(which(H0_ind[i,]< obs_assesor[i]))/length(H0_ind[i,])
   }

   T_m[i] <- as.numeric(quantile(H0_ind[i,], probs = (1-obs_prob)))
 }


#code from zap R package, from Dennis Leung
 bottom <- sapply(X = obs_assesor,
                  FUN = function(cutpt, vec){max(1, sum( vec <= cutpt ))},
                  vec = obs_assesor)

 top_BC <-  1 +  sapply(X = obs_assesor,
                         FUN = function(cutpt, vec){ sum( vec <= cutpt )},
                         vec = T_m)

 FDR_est <- ( top_BC /bottom)
# cutpt_max_BC <- max( blfdr_vec[FDR_est<= alpha])  # this can be -Inf
# rej_index_blfdr_BC  <- which(  blfdr_vec <=    cutpt_max_BC ) # this can be empty

return( out)
}


