
reg_log = function(X,y,threshold = 1e-10, max_iter = 100, intercept=TRUE)

{
  if (intercept){
    X <- cbind(rep(1,nrow(X)),X)
  }


  calc_p = function(X,beta)
  {
    beta = as.vector(beta)
    return(exp(X%*%beta) / (1+ exp(X%*%beta)))
  }
  beta = rep(0,ncol(X))

  diff = 10000

  iter_count = 0


  while(diff > threshold )
  {
    #calculate probabilities using current estimate of beta
    p = as.vector(calc_p(X,beta))
    W =  diag(p*(1-p))

    #calculate the change in beta
    beta_change = solve(t(X)%*%W%*%X) %*% t(X)%*%(y - p)

    #update beta
    beta = beta + beta_change

    diff = sum(beta_change^2)

    #see if we've hit the maximum number of iterations
    iter_count = iter_count + 1
    if(iter_count > max_iter) {
      stop("Not converging")
    }
  }

  names(beta)[1] <- "Intercept"
  out <- list(coef     = beta,
              fitted_p = calc_p(X,beta))# first coef is intercept
              return(out)
}




#'X <- matrix(rnorm(30*20), ncol=20)
#'y <- sigmoid(0.3+ X[,1])

#'fit <- reg_log (X,y)
#'plot(y,fit$fitted_p)
