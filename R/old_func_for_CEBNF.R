set_data_como <- function(betahat, se, p,zscore, X, ...) {
  if (!is.numeric(X)) {
    stop("X should be numercial vector")
  }
  if (!is.matrix(X)) {
    stop("X should be a matrix")
  }


  if (!missing(betahat)) {
    if (!is.numeric(betahat)) {
      stop("betahat should be numercial vector")
    }
    if (!is.numeric(se)) {
      stop("se should be numercial vector")
    }
    if (!(sum(c(length(se) == length(betahat), length(se) == nrow(X))) == 2)) {
      stop(" The number of lines in X should be equal to the number of entries in Betahat and se ")
    }
    dat <- list(
      betahat = betahat,
      se = se,
      X = X
    )
    class(dat) <- c("normal", "data_como")
  }


  if (!missing(p)) {
    if (!is.numeric(p)) {
      stop("p should be numercial vector")
    }
    if (!(length(p) == nrow(X)) == 1) {
      stop(" The number of lines in X should be equal to the number of entries in p ")
    }
    if( (min(p) <0)| max(p)>1){
      stop(" the provided p-values are not between 0 and 1")
    }
    dat <- list(
      p = p,
      se = rep(1, length(p)),
      X = X
    )
    class(dat) <- c("beta", "data_como")
  }
  if (!missing(zscore)) {
    if (!is.numeric(zscore)) {
      stop("zscore should be numercial vector")
    }
    if (!(length(zscore) == nrow(X)) == 1) {
      stop(" The number of lines in X should be equal to the number of entries in p ")
    }

    dat <- list(
      p = pnorm(zscore),
      se = rep(1, length(zscore)),
      X = X,
      zscore=zscore
    )
    class(dat) <- c("beta", "data_como")
  }

  return(dat)
}
