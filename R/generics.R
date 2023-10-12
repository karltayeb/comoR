#' @importFrom logisticsusie fit_model
#' @export
logisticsusie::fit_model

#' @importFrom logisticsusie update_model
#' @export
logisticsusie::update_model

#' @importFrom logisticsusie compute_elbo
#' @export
logisticsusie::compute_elbo





#' @title Compute post first and second moment from como and como2 object
#
#' @export
post_mean_sd <- function (fit,data,...)
  UseMethod("post_mean_sd")



#' @title Compute individual posterior variance from marginal normal mean model
#' @description internal function to compute posterior mean and sds
t_ind_var <- function (fit, data, i,...)
  UseMethod("t_ind_var")


#' @title Compute individual fdr value como/como2 model with centered normal mixture
#' @description Compute individual fdr value como/como2 model with centered normal mixture
#'
get_fdr <- function (fit,data, i,...)
  UseMethod("get_fdr")

#' @title Compute individual lfdr value como/como2 model with centered normal mixture
#' @description Compute individual lfdr value como/como2 model with centered normal mixture
#'
get_lfdr <- function (fit,data, i,...)
  UseMethod("get_lfdr")




#' @title Cal purity ouput como2
#' @description  Cal purity ouput como2
#'
cal_purity <- function (fit,data, ...)
  UseMethod("cal_purity")
