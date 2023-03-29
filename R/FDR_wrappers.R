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

  up <- exp(compute_assignment_jj_bound.mococomo(fit))[,1]
  tt <-exp(compute_assignment_jj_bound.mococomo(fit))
  tt2 <- 1- apply(tt,1,sum)
  low <-  up + apply(cbind(tt[,-1] )*exp(compute_data_loglikelihood(fit, fit$data))[,-1],1,sum)


obs_assesor <- up/low
sim_null <- list()
tdata <- fit$data


temp_f <- function(i){
  tdata$p  <-  rep(runif(1), (nrow(fit$data$X)))
  low_sim <-   up + apply(cbind(tt[,-1] )*exp(compute_data_loglikelihood(fit,tdata))[,-1],1,sum)
  return(up/low_sim)
}





 H0_ind <- do.call(cbind, lapply(1:n_sim, function(i) temp_f(1)))




for( i in 1:nrow(H0_ind)){
  H0_ind[i,] <-  H0_ind[i,] [order( H0_ind[i,] )]
}




for( i in 1:nrow(H0_ind)){
  H0_ind[i,] <-  H0_ind[i,] [order( H0_ind[i,] )]
}


T_m <- rep(NA, length(obs_assesor))

for ( i in 1:length(obs_assesor)){
  if( length(which(H0_ind[i,]< obs_assesor[i]))==0)
  {

  }else{
    obs_prob <-  length(which(H0_ind[i,]< obs_assesor[i]))/length(H0_ind[i,])
  }

  T_m[i] <- as.numeric(quantile(H0_ind[i,], probs = ( 1-obs_prob)))
}


#code from zap R package, from Dennis Leung
bottom <- sapply(X = obs_assesor,
                 FUN = function(cutpt, vec){max(1, sum( vec <= cutpt ))},
                 vec = obs_assesor)

top_BC <-  1 +  sapply(X = obs_assesor,
                       FUN = function(cutpt, vec){ sum( vec <= cutpt )},
                       vec = T_m)

FDR_est <- ( top_BC /(bottom ))



FDR_est_final <- do.call(c,lapply(1:length(obs_assesor), function(i)
                                                 min(FDR_est[which(obs_assesor>=obs_assesor[i])])
                                 )
                          )

#FDP_est <- FDR_est
#plot( log10(obs_assesor),FDP_est, xlim=c(-2.72,-2.4))
#points(log10(obs_assesor),FDR_est_final , col="red")

return(FDR_est_final )

}

