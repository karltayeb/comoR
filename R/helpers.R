
# Not in operator

"%!in%" <- function(x, y) !("%in%"(x, y))

#' Sigmoid
#' sigmoid function coverts log-odds to probability
#' @param x Log odds
#' @return Returns the probability
sigmoid <- function(x) {
  return(1 / (1 + exp(-x)))
}

#' @export
logodds <- function(p) {
  return(log(p) - log(1 - p))
}

logSumExp <- matrixStats::logSumExp

softmax <- function(x) {
  return(exp(x - logSumExp(x)))
}

.monotone <- function(v) {
  return(all(tail(v, -1) - head(v, -1) >= 0))
}

#' helper function gets difference in ELBO between iterations
.diff <- function(fit) {
  delta <- diff(tail(fit$elbo, 2))
  if (delta < 0) {
    warning("ELBO decreased")
  }
  return(delta)
}

#' check if ELBO is converged to some tolerance threshold
#' @export
.converged <- function(fit, tol = 1e-3) {
  is.converged <- F
  if ((length(fit$elbo) > 1) & abs(.diff(fit)) <= tol) {
    is.converged <- T
  }
  return(is.converged)
}


rowCumSum <- function(X) {
  do.call(rbind, apply(X, 1, cumsum, simplify = F))
}


#' Get Credible Sets
#' Wraps susieR::susie_get_cs
#' @export
binsusie_get_cs <- function(fit,
                            coverage = 0.95,
                            min_abs_corr = 0.5,
                            dedup = TRUE,
                            squared = FALSE,
                            check_symmetric = TRUE,
                            n_purity = 100) {
  res <- list(
    X = fit$data$X,
    alpha = fit$params$alpha,
    V = fit$hypers$prior_variance
  )
  cs <- susieR::susie_get_cs(res,
    X = fit$data$X,
    coverage = coverage,
    min_abs_corr = min_abs_corr,
    dedup = dedup,
    squared = squared,
    check_symmetric = check_symmetric,
    n_purity = n_purity
  )
  return(cs)
}

# implement normal and point mass component distributions

.clamp <- function(v, .upper = 100, .lower = -100) {
  return(pmax(.lower, pmin(v, .upper)))
}


is.odd <- function(x) {
  x %% 2 == 1
}
is.even <- function(x) {
  x %% 2 == 0
}


get_cs <- function(alpha, requested_coverage = 0.95) {
  rho <- order(-alpha)
  idx <- min(sum(cumsum(alpha[rho]) < requested_coverage) + 1, length(alpha))
  cs <- rho[1:idx]
  coverage <- sum(alpha[cs])
  return(list(cs = cs, prob = alpha[rho[1:idx]], size = idx, requested_coverage = requested_coverage, coverage = coverage))
}

get_all_cs <- function(fit, requested_coverage = 0.95) {
  sets <- purrr::map(1:fit$hypers$L, ~ get_cs(fit$params$alpha[.x, ], requested_coverage))
  names(sets) <- paste0("L", 1:fit$hypers$L)
  return(sets)
}


#' @title Preparing data for como fit
#' @details Preparing data for como fit for two type of input p-values or estimated regression coefficients
#' currently only support one type of entry, whether betahat or p-values
#'
#' @param betahat numeric, vector of regression coefficients  list of containing the following element see
#' @param se numeric, corresponding standard error
#' @param p numeric observed p-values
#' @param zscore numeric observed z-scores
#' @param maxiter numeric, maximum numerous of iteration set to 100 by defaults
#' @param tol tolerance in term of change in ELBO value for stopping criterion
#' @export
#' @example
#' see \link{\code{fit.como}}
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


get_all_cs2 <- function(alpha, requested_coverage = 0.95) {
  if (is.null(dim(alpha))) {
    alpha <- matrix(alpha, nrow = 1)
  }
  L <- dim(alpha)[1]
  sets <- purrr::map(1:L, ~ get_cs(alpha[.x, ], requested_coverage))
  names(sets) <- paste0("L", 1:L)
  return(sets)
}

#' Check if a CS covers a particular index
get_coverage <- function(cs, idx) {
  cs$covered <- idx %in% cs$cs
  cs$which_covered <- intersect(idx, cs$cs)
  names(cs$covered) <- idx
  return(cs)
}

#' Check if set of indices `idx` in a list of credible sets
get_all_coverage <- function(css, idx) {
  purrr::map(css, ~ get_coverage(.x, idx))
}

# convenient table of CSs from get_all_cs2
cs_tbl2 <- function(alpha) {
  get_all_cs2(alpha) %>%
    dplyr::tibble() %>%
    tidyr::unnest_wider(1) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(top_feaure = cs[1], top_feature_alpha = prob[1]) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(L = paste0("L", 1:dplyr::n())) # %>% unnest_longer(c(cs, prob))
}


#' Compute PIPs
get_pip <- function(alpha) {
  if (is.null(dim(alpha))) {
    alpha <- matrix(alpha, nrow = 1)
  }

  p <- dim(alpha)[2]
  pip <- purrr::map_dbl(1:p, ~ 1 - prod(1 - alpha[, .x]))
  return(pip)
}
#' set xi, update tau
set_xi <- function(fit, xi) {
  fit$params$xi <- xi
  fit$params$tau <- compute_tau(fit)
  return(fit)
}


# fit a como object
# cs an object computed using get_all_cs
# X matrix of covariate
# min.purity minimum purity level

which_dummy_cs_como <- function( cs, X,min.purity){
  dummy.cs <- list()
  for ( k in 1:length(cs)){
    dummy.cs[[k]] <- NULL
    for ( j in 1:length(cs[[k]]))
    {
      tt <- c()
      if( min(cor( X[,cs[[k]][[j]]$cs])) <  min.purity){
        #check if the purity of cs l is lower that min.purity
        tt<-  c( tt,j)
      }


    }
    if( length(tt)>0){
      dummy.cs[[k]] <- tt
    }

  }
  dummy.cs[[(k+1)]] <-NA #ensure that if everybody is null then the list exist with k+1 element

  return(dummy.cs)

}


# copied form; susiF.alpha
#l_cs a list of CS
#X the matrix of genotype

cal_purity_cFDR<- function(l_cs,X){
  tt <- list()
  for (k in 1:length(l_cs)){
    if(length(unlist(l_cs[[k]]$cs))==1 ){
      tt[[k]] <- 1
    }else{
      x <-abs( cor(X[,unlist(l_cs[[k]]$cs   ) ]))


      tt[[k]] <-  min( x[col(x) != row(x)])
    }
  }
  return( tt )
}



#from https://stackoverflow.com/questions/19465720/writing-a-function-to-calculate-divisors-in-r

divisors <- function(x){
  #  Vector of numberes to test against
  y <- seq_len(x)
  #  Modulo division. If remainder is 0 that number is a divisor of x so return it
  y[ x%%y == 0 ]
}
