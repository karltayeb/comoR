#'@title Main FDR function
#'@description Main function gathering all the routines for different FDR estimation procedures
#' @param outputlevel	 Determines amount of output. outputlevel=1 provide simplest output (which should be enough for most applications)
#' outputlevel= also output the entire como fitted object
#'@details if the user specifies (betahat, se) then cFDR will fit an covariate
#' moderated ash model. If the user specifies p-value, then cFDR fits a mixture of betas
#' which only contains distribution with decreasing density (mostly focusing on fitting the leftmost tails).
#' If the user specifies upper=TRUE when using the pvalue argument then cFDR fits a mixture of betas
#' which  contains distribution with decreasing density or distribution with increasing density (thus fitting both tail of the distribution).
#' If the user specifies the model essentially fit a extension of the ZAP procedure (projected zscore of 0-1) which
#' corresponds to fitting the same model as for pvalue argument with upper=TRUE;  (mixture of betas
#' which  contains distribution with decreasing density or distribution with increasing density)
#' @export
#'
#'
#'
#'
#'
#'
#'
cFDR <- function( betahat,
                  se,
                  #  pvalue,
                  #  zscore,
                  X,
                  Z,
                  param_como2 = list(logreg='logistic_ibss',
                                     f1_dist = 'ash',
                                     logreg_params = list(L=5)
                                     ),
                  coverage  = 0.95,
                  maxiter   = 100,
                  tol       = 1e-3,
                  max_class = 10,
                  mult      = 2,
                  upper     = FALSE,
                  outputlevel = 1,
                  n_sim= 1000,
                  alpha,
                  nullweight=4,
                  verbose=TRUE,
                  max_iter_como =10,
                  min.purity=.3
)
{
    if( !(missing(betahat))&missing(se)){
    stop("Please provide standard error when providing regression estimate or use zscore input")
   }



  N <- length(betahat)
  if(missing(Z)){
  Z <- matrix(rep(1, N), nrow=N)
  }


    data <- prep_data_como2(betahat=betahat,
                            se=se, X=X,
                            Z =Z )

  if(verbose){
    print( "Fitting como model")
  }



  fit  <- rlang::exec(  "data_initialize_como2" , !!!param_como2, data= data)
  fit  <- fit_model( fit,  data, max_iter = max_iter_como, estimate_f1=T )



if(verbose){
  print( "Model fitting done, performing FDR computation")
}
  out <- prep_out_FDR_wrapper (fit=fit,
                               data= data,
                               outputlevel= outputlevel,
                               n_sim = n_sim,
                               min.purity=  min.purity)


  return( out)
}

