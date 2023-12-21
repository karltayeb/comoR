samp_simu <- function(dist, N=200000){
betatrue <- c()
  if( dist=="spiky"){

    samp_fun <- function(){
      id <- sample(size=1, 1:4)
      if(id==1){
        out <-  rnorm( 1, sd=sqrt(2)*0.25)
      }
      if(id==2){
        out <-  rnorm( 1, sd=sqrt(2)*0.5)
      }
      if(id==3){
        out <-  rnorm( 1, sd=sqrt(2)*1)
      }
      if(id==4){
        out <-  rnorm( 1, sd=sqrt(2)*2)
      }
      return( out)
    }
    for ( i in 1:N){

      betatrue <- c( betatrue,  samp_fun())
    }

  }

  if( dist=="near_normal"){
    samp_fun <- function(){
      id <- sample(size=1, 1:2, prob=c(2/3,1/3))
      if(id==1){
        out <-  rnorm( 1, sd=sqrt(2)*1)
      }
      if(id==2){
        out <-  rnorm( 1, sd=sqrt(2)*2)
      }

      return( out)
    }

    for ( i in 1:N){

      betatrue <- c( betatrue,  samp_fun())
    }

  }

  if( dist=="normal"){
    for ( i in 1:N){

      betatrue <- c( betatrue,  rnorm(1,sd=sqrt(2)*1))

    }

  }

  if( dist=="flattop"){

    samp_fun <- function(){
      id <- sample(size=1, 1:7 )
      if(id==1){
        out <-  rnorm( 1,mean=-1.5, sd=5)
      }
      if(id==2){
        out <-  rnorm( 1,mean=-1, sd=5)
      }
      if(id==3){
        out <-  rnorm( 1,  mean=-.5, sd=5)
      }
      if(id==4){
        out <-  rnorm( 1,   sd=5)
      }
      if(id==5){
        out <-  rnorm( 1,  mean= .5, sd=5)
      }
      if(id==6){
        out <-  rnorm( 1, mean= 1, sd=5)
      }
      if(id==7){
        out <-  rnorm( 1,  mean= 1.5, sd=5)
      }
      return( out)
    }

    for ( i in 1:N){

      betatrue <- c( betatrue,  samp_fun())

    }


  }

  if( dist=="skew"){



    samp_fun <- function(){
      id <- sample(size=1, 1:4, prob = c(1/4,1/4,1/3,1/6) )
      if(id==1){
        out <-  rnorm( 1,mean=-2, sd=sqrt(2)*2)
      }
      if(id==2){
        out <-  rnorm( 1,mean=-1, sd=sqrt(2)*1.5)
      }
      if(id==3){
        out <-  rnorm( 1,  mean=0, sd=sqrt(2)*1)
      }
      if(id==4){
        out <-  rnorm( 1,mean=1,   sd=sqrt(2)*1)
      }

      return( out)
    }
    for ( i in 1:N){
      betatrue <- c( betatrue,  samp_fun())
    }


  }

  if( dist=="big-normal"){
    for ( i in 1:N){
      betatrue <- c( betatrue, rnorm(1,sd=4))
    }


  }
  if( dist=="bimodal"){
    samp_fun <- function(){
      id <- sample(size=1 ,1:2   )
      if(id==1){
        out <-  rnorm( 1,mean=-2, sd=1)
      }
      if(id==2){
        out <-  rnorm( 1,mean=2, sd=1)
      }


      return( out)
    }


    for ( i in 1:N){

      betatrue <- c( betatrue, samp_fun())
    }



  }


  return(data.frame(x= betatrue,dist=rep( dist, length(betatrue))))
}


tt <- rbind( samp_simu(dist="spiky"),
             samp_simu(dist="near_normal"),
             samp_simu(dist="flattop"),
             samp_simu(dist="skew"),
             samp_simu(dist="big-normal"),
             samp_simu(dist="bimodal")


)


tt <-  data.frame(tt)
head(tt)
library(ggplot2)
ggplot( tt, aes(x) )+
  xlim(-10,10)+
  facet_wrap(.~dist , nrow=1)+
  theme_minimal()+
   geom_density(size=1.5
                )
