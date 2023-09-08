#' adapted from autoselect.mxisqp
#' try to select a default range for the sigmaa values
#' that should be used, based on the values of betahat and sebetahat
#' mode is the location about which inference is going to be centered
#' gridmult is the multiplier by which the sds differ across the grid
#' @export
autoselect_scales <- function(betahat, sebetahat, max_class, mult = 2) {
  sigmaamin <- min(sebetahat) / 10 # so that the minimum is small compared with measurement precision
  if (all(betahat^2 <= sebetahat^2)) {
    sigmaamax <- 8 * sigmaamin # to deal with the occassional odd case where this could happen; 8 is arbitrary
  } else {
    sigmaamax <- 2 * sqrt(max(betahat^2 - sebetahat^2)) # this computes a rough largest value you'd want to use, based on idea that sigmaamax^2 + sebetahat^2 should be at least betahat^2
  }

  if (mult == 0) {
    return(c(0, sigmaamax / 2))
  } else {
    npoint <- ceiling(log2(sigmaamax / sigmaamin) / log2(mult))
    out <- c(0,mult^((-npoint):0) * sigmaamax)#  mult^((-npoint):0) * sigmaamax
    if (!missing(max_class)) {
      if (length(out) > max_class) {
        out <- seq(min(out), max(out), length.out = max_class)
      }
    }
    return(out)
  }
}

#' Prepare data for multinomial regression
#' @param betahat vector of effect estimates
#' @param se vector of standard errors for effect estimates
#' @param X matrix of features, each row an observation
#' @param Z matrix of fixed-effect covariates (including intercept), each row an observation
#' @param center boolean to center `X`
#' @param scale boolean to scale `X` to unit variance
como_prep_data <- function(betahat, se, X, Z=NULL, center = TRUE, scale = FALSE) {
  # center and scale data
  logisticsusie:::.check_X(X) # throws errors if something is wrong
  X <- Matrix::Matrix(scale(X, center = center, scale = scale))

  # set Z to intercept covariate if not provided
  n <- nrow(X)
  if (is.null(Z)) {
    Z <- matrix(rep(1, n), nrow = n)
  }

  # Make data list
  # TODO: store X means and standard errors
  data <- list(betahat = betahat, se = se, X = X, X2 = X^2, Z = Z)
  return(data)
}

como_check_data <- function(data){
  # check for betahat, se, X, X2, Z, etc
  TRUE
}
