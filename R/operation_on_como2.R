#' @title Compute individual posterior mean and sd under como model
#' @description Compute individual posterior mean and sd under como model
#'
#' # TODO: allow using new observation from another data set (e.g. testing)
#' @param fit a como object
#' @export
#' @example
#' see \link{\code{fit.como}}
post_mean_sd.como2 <- function(fit,data) {


  t_post_var <- do.call(
    rbind,
    lapply(
      1:length( data$betahat),
      function(i) t_ind_var.como2(fit=fit,
                                 data=data,
                                 i=i)
    )
  )

  tt <- compute_post_assignment (fit, data)

  fit$post_assignment <- cbind(  tt,1- tt)
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

  tt <- compute_post_assignment (fit, data)

  fit$post_assignment <- cbind(  1-tt, tt)

  tt1 <- fit$post_assignment[, 1] * dnorm(data$betahat, mean = 0, sd = data$se)
  tt2 <- fit$post_assignment[, 2] * dnorm( data$betahat,
                                       mean = 0,
                                       sd = sqrt( data$se^2 +fit$f1$var)
      )


  out <- tt1 / (tt1 + tt2)

  return(out)
}



cal_purity.como2  <- function(fit,data){

  l_cs <- fit$logreg$logistic_ibss$cs
  tt <- list()
  for (k in 1:length(l_cs)){
    if(l_cs[[k]]$size==1 ){
      tt[[k]] <- 1
    }else{
      x <-abs( cor(as.matrix(X[,unlist(l_cs[[k]]$cs   ) ])))


      tt[[k]] <-  min( x[col(x) != row(x)])
    }
  }
  return( do.call( c,  tt ))
}
