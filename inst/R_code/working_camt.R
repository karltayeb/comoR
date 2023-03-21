P=20
N <- 5000

beta1 <-0.41


beta0=-1
est_lfdr <- list()
est_FDR  <- list()


  x1 <- rnorm(N)
  samp_prob <- 1/(1 +exp(-(beta0+beta1*x1)))


  mix <- c()
  p <- c()
  for ( i in 1:N){

    #mix <- rep(0,N)
    mix <-c(mix, sample(c(0,1), size=1, prob = c(1- samp_prob[i], samp_prob[i])))
    p <- c(p, ifelse( mix[i]==1, rbeta(1,shape1=1,shape2=100 ), runif(1)))

  }
  #p <- runif(N)
  X <- cbind( x1, matrix(rnorm(P*N), ncol=P))





verbose=TRUE
if(verbose){
  print( "Fitting como model")
}


model <- "beta"
data <- set_data_mococomo(p = p ,
                          X = X)



verbose=TRUE
if(verbose){
  print( "Fitting como model")
}

max_class=10
tol=0.001
nullweight=2.3
maxiter=40
max_class=10
tol=0.001
nullweight=2.3
maxiter=40
fit =  fit.mococomo(data,
                    model      = model,
                    maxiter   = maxiter,
                    tol       = tol,
                    max_class = max_class,
                    mult      = mult,
                    upper     = FALSE,
                    nullweight = nullweight)

if(verbose){
  print( "Model fitting done, performing FDR computation")
}
#out <- prep_out_camt (fit=fit,
#                            outputlevel= outputlevel,
#                             n_sim = n_sim)







#low <-  1 + apply(cbind(tt[,-1] )*exp(compute_data_loglikelihood(fit, fit$data))[,-1],1,sum)
low <- 1 +apply( fit$post_assignment[,-1]*exp(compute_data_loglikelihood(fit, fit$data))[,-1],1,sum)
#obs_assesor <- up/low #observed assessor

#up <- 1- apply(fit$post_assignment[,-1],1,sum)
#low <-  up + apply(fit$post_assignment[,-1]*exp(compute_data_loglikelihood(fit, fit$data))[,-1],1,sum)
#p_i0 <-exp(compute_assignment_jj_bound.mococomo(fit))[,1]

p_i0 <-fit$post_assignment[,1]
qs <-apply( fit$post_assignment[,-1]*exp(compute_data_loglikelihood(fit, fit$data))[,-1],1,sum)*(1-p_i0)/(  (  p_i0))
#qs <-apply(cbind(tt[,-1] )*exp(compute_data_loglikelihood(fit, fit$data))[,-1],1,sum)*(1-p_i0)/(  (  p_i0))



ts <- 1 / (1 + qs)
out <- sort(ts,  index.return = TRUE)
ts <- out$x
index0 <- out$ix
rs <- (1-p_i0)/(  (  p_i0))[index0]
qs <- qs[index0]
pvals <- p
pvals <- pvals[index0]
p_i0 <- p_i0[index0]
ts.p <- 1/(low*p_i0/(  (1- p_i0)))

ts.p <- sort(c(ts,ts.p ), decreasing = TRUE)
m <- length(ts.p)
n <- length(p )
fdr <-  (1 + (m + 1) - match(ts, ts.p) - 1:n) / (1:n)
fdr
fdr <- rev(pmin(1, cummin(rev(fdr))))
index0 <- order(index0)
fdr = fdr[index0]
ts = ts[index0]
p_i0 = p_i0[index0]



plot( fdr, p, col=mix+1)
table(mix[which( fdr <0.05) ])  /(c(sum(mix[which( fdr <0.05) ])   , sum(mix)))
