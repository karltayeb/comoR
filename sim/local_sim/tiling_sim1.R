library(softImpute)
library( comoR)
library(nnet)
library(ashr)
library(flashier)
library(irlba)
library(PMA)




if(file.exists(paste(getwd(),"/sim/local_res/tiling_1.RData",sep="" ))){
  load(paste(getwd(),"/sim/local_res/tiling_1.RData",sep="" ))

}else{
  res <-list()
}
for (o  in (length(res)+1):10000) {

  res[[o]] <-comoR:::tiling_sim( noise_level= 1, seed=o
  )

  save(res, file=paste(getwd(),"/sim/local_res/tiling_1.RData",sep="" ))
}

