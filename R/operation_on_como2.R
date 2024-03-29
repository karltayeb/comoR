#' @title Compute individual posterior mean and sd under como model
#' @description Compute individual posterior mean and sd under como model
#'
#' # TODO: allow using new observation from another data set (e.g. testing)
#' @param fit a como object
#' @export
#' @example
#' see \link{\code{fit.como}}
post_mean_sd.como2 <- function(fit,data) {


  if( inherits(fit$f1,"ash_component")){

    tt <- compute_post_assignment (fit, data)

    fit$post_assignment <- cbind( 1-tt,tt)
    ash_data <-  ashr::set_data(data$betahat,data$se)

    post_mean <- tt*ashr::postmean(ashr::get_fitted_g(fit$f1$ash),ash_data )
    post_sd <- tt*ashr::postsd(ashr::get_fitted_g(fit$f1$ash),ash_data )
    out <- data.frame(
      mean =post_mean,
      sd = post_sd
    ) #
  }else{
    t_post_var <- do.call(
      rbind,
      lapply(
        1:length( data$betahat),
        function(i) t_ind_var.como2(fit=fit,
                                    data=data,
                                    i=i)
      )
    )
    t_post_var[,1] <- 0
    tt <- compute_post_assignment (fit, data)

    fit$post_assignment <- cbind( 1-tt,tt)
    out <- do.call(
      rbind,
      lapply(
        1:length(data$betahat),
        function(i) cal_ind_moment12(fit=fit,
                                     data=data,
                                     t_post_var=t_post_var,
                                     i=i)
      )
    )
    out <- data.frame(
      mean = do.call(c, out[, 1]),
      sd = do.call(c, out[, 2])
    ) # could be skip for speed

  }




  return(out)
}





#' @title Compute individual posterior variance from marginal normal mean model
#' @description internal function to compute posterior mean and sds
t_ind_var.como2 <- function(fit, data, i) {


   out <-  c (  data$se[i]^2 ,   1 / ((1 /  data$se[i]^2) + (1 /fit$f1$var))

    )

}
#' @title Compute individual fdr value   como2 model with centered normal mixture
#' @descriptionCompute individual fdr value  como2 model with centered normal mixture
#'
#' @param fit a como object
#' @export
get_lfdr.como2 <- function(fit,data) {






  return( 1-compute_post_assignment (fit, data))
}



cal_purity.como2  <- function(fit,data){

  l_cs <- fit$logreg$logistic_ibss$cs
  tt <- list()
  for (k in 1:length(l_cs)){
    if(l_cs[[k]]$size==1 ){
      tt[[k]] <- 1
    }else{

      x <-abs( cor(as.matrix(data$X[,unlist(l_cs[[k]]$cs   ) ])))


      tt[[k]] <-  min( x[col(x) != row(x)])
    }
  }
  return( do.call( c,  tt ))
}
