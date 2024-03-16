library(softImpute)
library( comoR)
library(nnet)
library(ashr)
library(flashier)
library(irlba)
library(PMA)




if(file.exists(paste(getwd(),"/sim/local_res/tiling_4.RData",sep="" ))){
  load(paste(getwd(),"/sim/local_res/tiling_4.RData",sep="" ))

}else{
  res <-list()
}
for (o  in (length(res)+1):10000) {

  res[[o]] <-comoR:::tiling_sim( noise_level= 5, seed=o+1
  )

  save(res, file=paste(getwd(),"/sim/local_res/tiling_4.RData",sep="" ))
}

