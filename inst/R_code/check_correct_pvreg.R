set.seed(2)
P=20
N <- 5000

beta1 <-0.41


beta0=-1
est_lfdr <- list()
est_FDR  <- list()
for ( o in 1:100){

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



  res <-   cFDR( pvalue  =  p,
                 X       =  X,
                 n_sim= 1000, nullweight = 2.3,
                 outputlevel = 2 )






  plot(res$result$lfdr, res$result$p,col=ifelse(res$result$lfdr <0.05,1,2) , ylim = c(0,0.2),xlim=c(0,0.2) )

   points(res$result$FDR, res$result$p,col=ifelse(res$result$FDR <0.05,1,2)  )
   points(p.adjust(p,method = "BH"),res$result$p,col=ifelse(p.adjust(p,method = "BH") <0.05,1,2))

   abline(v=0.05)
   print(which(res$result$lfdr <0.05))
   print(which(res$result$FDR <0.05))
   est_lfdr[[o]] <-  table(mix[which(res$result$lfdr <0.05)])  /(c(sum(table(mix[which(res$result$lfdr <0.05)]) ) , sum(mix)))
   est_FDR [[o]] <-  table(mix[which(res$result$FDR <0.05)])  /(c(sum(table(mix[which(res$result$FDR <0.05)]) ) , sum(mix)))
  #print( apply( do.call(rbind, est_lfdr),2, mean))
}




fit <- res


cs <- list()
for( k in 1:length(fit$full_obj$logreg_list))
{
  tt <- get_all_cs(fit$full_obj$logreg_list[[k]])
  cs[[k]] <- list()
  for ( l in 1:length(tt)){
    cs[[k]] [[l]]  <-as.vector(unlist(tt[[l]]["cs"]))
  }
}

cs_como <- cs

cs_como


