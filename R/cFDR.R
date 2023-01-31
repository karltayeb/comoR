#'@title Main FDR function
#'@description Main function gathering all the routines for different FDR estimation procedures



#'@details if the user specifies (betahat, se) then cFDR will fit an covariate
#' moderated ash model. If the user specifies p-value, then cFDR fits a mixture of betas
#' which only contains distribution with decreasing density (mostly focusing on fitting the leftmost tails).
#' If the user specifies upper=TRUE when using the pvalue argument then cFDR fits a mixture of betas
#' which  contains distribution with decreasing density or distribution with increasing density (thus fitting both tail of the distribution).
#' If the user specifies the model essentially fit a extension of the ZAP procedure (projected zscore of 0-1) which
#' corresponds to fitting the same model as for pvalue argument with upper=TRUE;  (mixture of betas
#' which  contains distribution with decreasing density or distribution with increasing density)

cFDR <- function( betahat,
                  se,
                  pvalue,
                  zscore,
                  X,
                  coverage  = 0.95,
                  maxiter   = 100,
                  tol       = 1e-3,
                  max_class = 10,
                  mult      = 2,
                  upper     = FALSE
)
{
  if( !(missing(betahat))&missing(se)){
    stop("Please provide standard error when providing regression estimate or use zscore input")
  }
  if(   missing(betahat)&missing(pvalue)&missing(zscore)){
    stop("Please provide one of the following entry:betahat, pvalue or zscore")
  }


  if( !missing(betahat)& !missing( se)){
    dist <- "normal"
    data <- set_data_mococomo(betahat  = betahat,
                              X = X,
                              se=  se)
    #fit mococomo model
  }
  if(!missing(pvalue)){
    dist <- "beta"
    data <- set_data_mococomo(p = pvalue,
                              X = X)
  }

  if(!missing(zscore)){
    dist <- "beta"
    data <- set_data_mococomo(p = pnorm(zscore), #using U value as in the ZAP paper
                              X = X)
    upper =TRUE
  }

  fit =  fit.mococomo(data,
                      dist      = dist,
                      maxiter   = maxiter,
                      tol       = tol,
                      max_class = max_class,
                      mult      = mult,
                      upper     = upper)


  out <- prep_out_FDR_wrapper (fit)


  return( out)
}

