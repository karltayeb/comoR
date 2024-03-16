library(softImpute)
library( comoR)
library(nnet)
library(ashr)
library(flashier)
library(irlba)
library(PMA)




if(file.exists(paste(getwd(),"/sim/local_res/check_cEBMF_1.RData",sep="" ))){
  load(paste(getwd(),"/sim/local_res/check_cEBMF_1.RData",sep="" ))

}else{
  res <-list()
}
for (o  in (length(res)+1):10000) {

  res[[o]] <-sim_func_cEBMF (noise_level= 0.1,
                             max_iter_cEBMF=50,
                             P1=200 , # number of cov for row /loadings
                             P2=200, # number of cov for col /factors
                             seed=o)

  save(res, file=paste(getwd(),"/sim/local_res/check_cEBMF_1.RData",sep="" ))
}

