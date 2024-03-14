
#' @title covariate moderated Empirical Bayes Matrix Factorization
#' @description covariate moderated Empirical Bayes Matrix Factorization
#' @param Y numerical matrix size NxP
#' @param X_l matrix of size NxJ containning covariates affecting the factors
#' @param X_f matrix of size PxT containning covariates affecting the factors
#' @param reg_method specify underlying learner for within topic heterogeneity , two methods available 'nnet' or 'logistic_susie'. Default is nnet.
#' you can pass the nnet specification through the param_nnet argument and through the param_susie for logistic susie.
#' @param K numeric number of factors
#' @param type_noise specify which kind of noise structure is expected, currently three choices. Whether noise constant accross column ('column_wise'), constant 'constant' or constant across rown 'row_wise'
#' @param maxit maximum nuber of iterations
#' @param tol paramter for assessing convergence
#' @return a cEBMF object
#'
#' @export
#'

cEBMF <- function( Y,
                   X_l,
                   X_f,
                   reg_method="nnet",
                   K=1,
                   type_noise='constant',
                   init_type="udv_si",
                   maxit=100,
                   tol=1e-3 ,
                   param_como  = list(max_class=10,mnreg_type="mult_reg"),
                   param_nnet  =list( size=1, decay=1),
                   param_como2 = list(),
                   param_susie =  list(L=5),
                   maxit_como  = 10){

  if(reg_method %!in% c("nnet","logistic_susie")){
    stop("reg_method should be equal to nnet or logistic_susie ")
  }

  cEBMF.obj <- init_cEBMF (Y,
                           X_l,
                           X_f,
                           reg_method=reg_method,
                           K=K,
                           type_noise  = type_noise,
                           init_type   = init_type,
                           param_como  = param_como,
                           maxit_como  = maxit_como,
                           param_nnet  = param_nnet,
                           param_como2 = param_como2,
                           param_susie = param_susie

  )### Need to carry info about como obj

  for (i in 1:maxit) {
    cEBMF.obj <- cEBMF_iter  (cEBMF.obj)
    if (.converged(cEBMF.obj)) {
      break
    }
  }
  # plot( cEBMF.obj$elbo)
  cEBMF.obj <- out_prep.cEBMF(cEBMF.obj)
  return(cEBMF.obj)
}
