# Flower longevity

source("Functions.R")

simulation2 <- function(nbees,signalDelt,Rm) { # we could name some parameters we wanted to study in this line and then specify them in the initiation calls!
  
  #make a field of random flowers
  Field <- list()
  for(i in 1:500){
    Field[[i]] <- FlowerInit(signalPar = c(0,100), sugarPar = c(40,10), microbesP =c(0,0) )
  }
  for (i in 1:15){
    Field[[i]]$microbesP1 <- rnorm(1,40,10)
    Field[[i+15]]$microbesP2<-rnorm(1,40,10)
  }
  ti<- do.call(what = rbind,Field)
  signals<-ti$flowerSignal
  
  
  Learner<-list()
  for(i in 1:nbees){
    Learner[i] <- list(BeeInit(microbesPar = c(0,0),thresholdPar = c(median(signals)-rnorm(1,10,2),median(signals)+rnorm(1,10,2)))) #bee with 100 microbes
  }
  ti<- do.call(what = rbind,Field)
  ggplot(ti,aes(x=flowerSignal,y=sugar))+geom_point()+geom_smooth(method="lm")
  # 
  Microbe1 <- MicrobeInit(Rm = Rm ,signalDelt = signalDelt ,sugarDelt = 0,growthRate = 1.02)
  Microbe2 <- MicrobeInit(sugarDelt = 0,signalDelt = 0 , Rm = 0, growthRate = 1.02)
  
  ####COMBACK####
  for(i in 1:200){
    #first figure out which flower bees will go to by randomly sampling as many flower as there are bees
    # print(i) not needed for parallel
    encounters <- sample(1:length(Field), size = length(Learner), replace=F)
    encounterList <- vector("list",length(Field))
    encounterList[encounters] <- Learner
    
    test <- mapply(FUN=function(x,y){
      tmp <- daytimestep(flower=x,bee=y,microbe1 = Microbe1, microbe2 = Microbe2)
    },Field,encounterList,SIMPLIFY=F)
    
    Learner <- lapply(X = test[encounters], FUN = function(x){return(x$bee)})
    Field <- lapply(X=test, FUN = function(x){
      if(!is.null(x[["flower"]])){return(x[["flower"]])}
      else(return(x))  
    })
  }
  
  
  
  Field <- lapply(X=Field, FUN= function(x){x$cumulativeVisits<-cumsum(x$visit)
  return(x)})
  
  
  # ggplot(Learner[[1]],aes(x=time,y=visit))+geom_point()
  # ggplot(Learner[[1]],aes(x=time,y=thresholdF1))+geom_point()
  
  field<-do.call(rbind,Field)
  
  # ggplot(field,aes(x=time,y=necVol,group=ID)) + geom_line()
  # ggplot(field,aes(x=time,y=microbesP1,group=ID)) + geom_line()
  # ggplot(field,aes(x=time,y=microbeK,group=ID)) + geom_line()
  # ggplot(field,aes(x=time,y=signal,group=ID)) + geom_line()
  # ggplot(field,aes(x=time,y=cumulativeVisits,group=ID,color=signal)) + geom_line()
  # ggplot(Learner[[1]],aes(x=time,y=LearnValue,size=LearnFlowerSig))+geom_point()
  # 
  # 
  # pdf(paste("Plots/",nbees,"plot.pdf",sep=""),w=6,h=3)
  # print(
  #   ggarrange(ncol=2,nrow = 2,
  #             ggplot(Learner[[1]],aes(x=time,y=thresholdF1))+geom_point()+scale_y_continuous(limits=c(0,max(ti$signal))),
  #             ggplot(ti,aes(y=signal,x=sugar))+geom_point()+geom_smooth(method="lm")+scale_y_continuous(limits=c(0,max(ti$signal)))+geom_vline(xintercept=30),
  #             ggplot(Learner[[1]],aes(x=LearnValue,y=LearnFlowerSig))+geom_point()+scale_y_continuous(limits=c(0,max(ti$signal))),
  #             ggplot(Learner[[1]],aes(x=LearnRew,y=LearnFlowerSig))+geom_point()+scale_y_continuous(limits=c(0,max(ti$signal)))                             )
  # )
  # dev.off()
  return(list(
    bees=Learner,
    flowers=Field
  )
  )
  
}

