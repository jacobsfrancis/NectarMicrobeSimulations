# Agent based model of pollination and microbial dispersal

# This model is designed to test how microbes which change the signal of a flower
# and change its quality might shape bees decisions to visit flowers, and shape their 
# own dispersal.  It also will eventually include parameters for floral defenses that 
# change microbial growth and bee preference.

# There are two types of agents:
# Bees - Params:
# 1) their decision threshold (i.e. signal above which they will visit a flower) 
# 2) the quantity of microbes they are carying
# 3) the quantity of pollen they are carying (for plant perspective questions)
# 4) their innate preferences for nectar sugars (and alkaloids for later)
# Flowers - Params:
# 1) Their species (for plant perspective questions later)
# 2) Their signal this is a linear measure of a signal (TODO think of multiple signal axes?)       
# 3) their microbial carrying capacity
# 4) Their nectar chemistry (sugars and alkaloids) 
# 5) Their pollen load (for plant perspective questions later)
# 6) Their pollen deposition rate (for plant perspective questions later)

# Microbes - Params:

# 1) their growth rate - r in traditional logistic growth equations
# 2) their change to the signal, this is multiplied by the population number so it should be small
# 3) their change to the sugar, this is multiplied by the population number so it should be small
# 4) their change to the alkaloid, this is multiplied by the population number so it should be small

# At each time step microbes grow using a logarithmic growth function, then
# bees approach a flower (randomly) and decide whether to visit it
# If they visit it they remove nectar and microbes, deposit microbes, 
#   assess flower quality, and update their threshold.

# Some requirements - there must always be more flowers than pollinators!

# TODO
# add in a way for bees to make a decision based on 2 things... microbes and flower signal. 
# I think this should be the average of the difference between what they would accept and what the actual signal is. This should be standardized around zero
# different thresholds for each signal? then predict how good the flower will be based on signal, and then average them.

# think about how to integrate Amber's data into the decision model, maybe multiply the vuf``alue by 


#clear the work space (important for debugging code in production)

#rm(list=ls())
library(tidyverse)
library(ggplot2)
library(reshape2)
library(parallel)
library(emmeans)
library(ggpubr)
library(betareg)
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)

#Plotting parameters ####
my_theme <- theme_bw() +
  theme(legend.position="none" , 
        text = element_text(size = 12), 
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=.5)) 
my_colors<- function(){ 
  return(scale_color_manual(values = c(rgb(158,54,63,maxColorValue = 220), 
                                       rgb(85,142,168,maxColorValue = 220), 
                                       rgb(185,161,52,maxColorValue = 220), 
                                       rgb(70,70,70,maxColorValue = 220))
  )
  )
}




#Test the functions ####

##Gut check 1) do value functions work the way we think they should####

#check growth model

flowerT<-FlowerInit()
mic1 <- MicrobeInit()
mic2 <- MicrobeInit()
flowerT$microbesP1 <-20
flowerT$microbesP2 <-10

for (i in 1:50) {flowerT <- grow(microbe1 = mic1,microbe2 = mic2, flower = flowerT )}

ggplot(flowerT,aes(x=time)) + geom_line(aes(y=microbeK))


testMic <- MicrobeInit(Rm=0.002)
MicrobeValFun(mic1,mic2,flowerT)
flowerT$microbeValue<-NA
microbe1<-MicrobeInit(Rm = .001)
microbe2 <- MicrobeInit(Rm=.0)
for (i in 1:nrow(flowerT)){flowerT[i,"microbeValue"] <- (MicrobeValFun(microbe1,microbe2,flowerT[i,]))}

growthvalueplot<-ggplot(flowerT,aes(x=time)) + 
  geom_line(aes(y=microbesP1/necVol), size=1) + 
  geom_line(aes(y=microbesP2/necVol),color="dark orchid3", size=1) +
  geom_line(color="darkred",aes(y=microbeValue*5000), size=1) +
  geom_text(aes(x=5,y=4000,label="Microbial\nValue"),color="darkred")+
  geom_text(aes(x=50,y=4500,label="Microbe 1"))+
  geom_text(aes(x=50,y=800,label="Microbe 2"),color="darkorchid")+
  scale_y_continuous("Microbial Density (cells/μl)")+
  scale_x_continuous("Time")+
  geom_text(data=data.frame(time=rep(-2,7),y=c(0,1000,2000,2500,3000,4000,5000), label=c(-1,-.2,-.4,0.5,0.4,0.2,1)),aes(y=y,label=label),color="darkred")+
  my_theme
growthvalueplot

ggplot2::ggsave(filename = "testplot.pdf", 
                plot = growthvalueplot, 
                
                dpi = 1200, 
                width = 15,
                height = 10, 
                units = "cm")




testFlower <- FlowerInit( signalPar = c(30.1,100))
testBee <- BeeInit(microbesPar = c(0,200) , thresholdPar = c(0,30))
testMicrobe1 <- MicrobeInit(growthRate = 1.002,Rm = 0.02)
testMicrobe2 <- MicrobeInit(growthRate = 1.02,Rm = 0.02)

for (i in 1:100){
  daytimestep(bee = testBee,flower = testFlower, microbe1 = testMicrobe1, microbe2 = testMicrobe2) -> tmp
  testBee <- tmp[["bee"]]
  testFlower <- tmp[["flower"]]
}

plot <- ggplot(testFlower, aes(x=time,y=microbesP1))
plot+geom_line() +geom_line(aes(y=microbesP2),color="dark orchid")#+ geom_line(aes(x=time,y=microbeK),color="blue")

plot <- ggplot(testBee[1:100,], aes(x=time,y=microbe1))
plot+geom_line()+geom_line(aes(y=microbe2),color="dark orchid")

plot <- ggplot(testBee[1:100,], aes(x=time,y=visit))
plot+geom_point()

plot <- ggplot(testBee[2:100,], aes(x=time,y=LearnRew))
plot+geom_point()

plot <- ggplot(testBee[2:100,], aes(x=time,y=LearnFlowerSig))
plot+geom_point()


plot <- ggplot(testBee, aes(x=time,y=pollenLoadF1))
plot+geom_point()

plot <- ggplot(testFlower, aes(x=time,y=pollenF1))
plot+geom_point()

plot <- ggplot(testFlower, aes(x=time,y=microbesP1))
plot+geom_line() +geom_line(aes(y=microbesP2),color="dark orchid")#+ geom_line(aes(x=time,y=microbeK),color="blue")

plot <- ggplot(testFlower, aes(x=time,y=microbeK))
plot+geom_line()




#lets see if a bee learns####




# Question : how does the relative number of bees to flowers impact visitation
#  probabilities of flowers in the environment


#heres how to do it with FOREACH, but it doesnt really handle generating random numbers in the callwell...####
# didshelearn <- foreach(z = 1:40,
#                        beta = betas,
#                        .packages=c("ggplot2","ggpubr"),
#                        .combine='rbind') %dopar% {
#     Learner <- list(BeeInit(microbesPar = c(100,100))) #bee with 100 microbes
# 
#     Field <- list()
# 
#     for(i in 1:100){
#         Field[[i]] <- FlowerInit(signalPar = c(0,100),sugarPar = c(30,10))
#         Field[[i]]$sugar <-abs( Field[[i]]$signal*beta +rnorm(1,0,8))
#      }
# 
#     ti<- do.call(what = rbind,Field)
#     ggplot(ti,aes(x=signal,y=sugar))+geom_point()
# 
#     PassiveMicrobe <- MicrobeInit(sugarDelt = 0,signalDelt = 0)
# 
# 
# 
#     for(i in 1:100){
#           #first figure out which flower bees will go to by randomly sampling as many flower as there are bees
#           # print(i) not needed for parallel
#           encounters <- sample(1:length(Field), size = length(Learner), replace=F)
#           encounterList <- vector("list",length(Field))
#           encounterList[encounters] <- Learner
# 
#           test <- mapply(FUN=function(x,y){
#             tmp <- daytimestep(flower=x,bee=y,microbe=PassiveMicrobe)
#           },Field,encounterList,SIMPLIFY=F)
# 
#           Learner <- lapply(X = test[encounters], FUN = function(x){return(x$bee)})
#           Field <- lapply(X=test, FUN = function(x){
#             if(!is.null(x$flower)){return(x$flower)}
#             else(return(x))})
#     }
# 
# 
#     Field <- lapply(X=Field, FUN= function(x){x$cumulativeVisits<-cumsum(x$visit)
#       return(x)})
# 
# 
#     # ggplot(Learner[[1]],aes(x=time,y=visit))+geom_point()
#     # ggplot(Learner[[1]],aes(x=time,y=thresholdF1))+geom_point()
# 
#     field<-do.call(rbind,Field)
# 
#     # ggplot(field,aes(x=time,y=necVol,group=ID)) + geom_line()
#     # ggplot(field,aes(x=time,y=microbesP1,group=ID)) + geom_line()
#     # ggplot(field,aes(x=time,y=microbeK,group=ID)) + geom_line()
#     # ggplot(field,aes(x=time,y=signal,group=ID)) + geom_line()
#     # ggplot(field,aes(x=time,y=cumulativeVisits,group=ID,color=signal)) + geom_line()
#     # ggplot(Learner[[1]],aes(x=time,y=LearnValue,size=LearnFlowerSig))+geom_point()
#     #
# 
# 
#     # didshelearn[z,]<-c(NA,NA)
#     # didshelearn$learnedT[z] <- mean(tail(na.omit(Learner[[1]]$thresholdF1),30))
#     # didshelearn$realT[z] <- predict(lm(data=ti, signal~sugar),newdata = data.frame(sugar=30))
#     pdf(paste("Plots/",z,"plot.pdf",sep=""),w=6,h=3)
#     print(
#     ggarrange(
#         ggplot(Learner[[1]],aes(x=time,y=thresholdF1))+geom_point()+scale_y_continuous(limits=c(0,max(ti$signal))),
#         ggplot(ti,aes(y=signal,x=sugar))+geom_point()+geom_smooth(method="lm")+scale_y_continuous(limits=c(0,max(ti$signal)))+geom_vline(xintercept=30)
#       )
#     )
#     dev.off()
#     return(
#       data.frame(Threshold = mean(tail(na.omit(Learner[[1]]$thresholdF1),30)),
#                  TrueVal=predict(lm(data=ti, signal~sugar),newdata = data.frame(sugar=30))))
# 
# }


#instead lets do it with mcapply####

#step 1 designate a function 

simulation <- function(nbees) { # we could name some parameters we wanted to study in this line and then specify them in the initiation calls!
  
  #make a field of random flowers
  Field <- list()
  beta<-abs(rnorm(n=1,mean=4,sd=.6))
  signals<-NA
  for(i in 1:100){
    Field[[i]] <- FlowerInit(signalPar = c(0,100),sugarPar = c(30,10), microbesP =c(0,0) )
    Field[[i]]$flowerSignal <-abs( Field[[i]]$sugar*beta + rnorm(1,0,8) + 30)
    signals[i] <- Field[[i]]$flowerSignal
  }
  for(i in 1:12){
    Field[[i]]$microbesP1 <- runif(1,10,40)
    Field[[i]]$microbesP2 <- runif(1,10,40)
  }
  Learner<-list()
  for(i in 1:nbees){
    Learner[i] <- list(BeeInit(microbesPar = c(0,0),thresholdPar = c(median(signals)-rnorm(1,10,2),median(signals)+rnorm(1,10,2)))) #bee with 100 microbes who's initial threshold is somewhere near the median plant signal
  }
  ti<- do.call(what = rbind,Field)
  ggplot(ti,aes(x=flowerSignal,y=sugar))+geom_point()+geom_smooth(method="lm")
  
  PassiveMicrobe <- MicrobeInit(sugarDelt = 0,signalDelt = 0,growthRate = 1.02,Rm=0)
  
  
  
  for(i in 1:100){
    #first figure out which flower bees will go to by randomly sampling as many flower as there are bees
    # print(i) not needed for parallel
    encounters <- sample(1:length(Field), size = length(Learner), replace=F)
    encounterList <- vector("list",length(Field))
    encounterList[encounters] <- Learner
    
    test <- mapply(FUN=function(x,y){
      tmp <- daytimestep(flower=x,bee=y,microbe1 = PassiveMicrobe, microbe2 = PassiveMicrobe)
    },Field,encounterList,SIMPLIFY=F)
    
    Learner <- lapply(X = test[encounters], FUN = function(x){return(x$bee)})
    Field <- lapply(X=test, FUN = function(x){
      if(!is.null(x[["flower"]])){return(x[["flower"]])}
      else(return(x))  
      })
  }
                  
  Field <- lapply(X=Field, FUN= function(x){x$cumulativeVisits<-cumsum(x$visit)
  return(x)})
  
  # 
  # ggplot(Learner[[1]],aes(x=time,y=visit))+geom_point()
  # ggplot(Learner[[1]],aes(x=time,y=thresholdF1))+geom_point()
  
  field<-do.call(rbind,Field)
  learners <- do.call(rbind,Learner)
  
  learndata<- data.frame(learners %>% dplyr::group_by(ID) %>% summarise(
    Signal = mean(tail(LearnFlowerSig[!is.na(LearnFlowerSig)],3)),
    Value = mean(tail(LearnValue[!is.na(LearnFlowerSig)],3 ))
  ))
  trueval <- predict(lm(data=ti,flowerSignal~sugar),newdata=data.frame(sugar=30))
  learndata$TrueVal<-trueval
   # ggplot(field,aes(x=time,y=necVol,group=ID)) + geom_line()
   # ggplot(field,aes(x=time,y=microbesP1,group=ID)) + geom_line()
   # ggplot(field,aes(x=time,y=microbeK,group=ID)) + geom_line()
   # ggplot(field,aes(x=time,y=signal,group=ID)) + geom_line()
   # ggplot(field,aes(x=time,y=cumulativeVisits,group=ID,color=signal)) + geom_line()
   # ggplot(Learner[[1]],aes(x=time,y=LearnValue,size=LearnFlowerSig))+geom_point()


  # pdf(paste("Plots/",nbees,"plot.pdf",sep=""),w=6,h=3)
  # print(
  #   ggarrange(ncol=2,nrow = 2,
  #     ggplot(Learner[[1]],aes(x=time,y=thresholdF1))+geom_point()+scale_y_continuous(limits=c(0,max(ti$signal))),
  #     ggplot(ti,aes(y=signal,x=sugar))+geom_point()+geom_smooth(method="lm")+scale_y_continuous(limits=c(0,max(ti$signal)))+geom_vline(xintercept=30),
  #     ggplot(Learner[[1]],aes(x=LearnValue,y=LearnFlowerSig))+geom_point()+scale_y_continuous(limits=c(0,max(ti$signal))),
  #     ggplot(Learner[[1]],aes(x=LearnRew,y=LearnFlowerSig))+geom_point()+scale_y_continuous(limits=c(0,max(ti$signal)))                             )
  # )
  # dev.off()
  return(list(
    bees=Learner,
    flowers=Field,
    
    learndata= learndata
      
      
  )
  
  )
  
}

# lets do a simulation where we test how pollinator population impacts visit frequency
# 1- 5 bees with 10 reps per simulation
beePopParms <- rep(5,20)
simulation(nbees=1)
#

simulations  <- mclapply(X=beePopParms,FUN=function(x){simulation(nbees =x)}, mc.cores = 11)




bees<-data.frame(do.call(what = rbind,simulations[[1]]$bees))
bees$beePop <- 5
bees$simnum <- 1
for(i in 2:length(beePopParms)){
  df <- do.call(rbind,simulations[[i]]$bees)
  df$beePop <- beePopParms[i]
  df$simnum <-i
  bees<-rbind(bees,df)
}

ggplot(bees,aes(x=time,y=thresholdF1,group=ID,color=beePop))+
  geom_line(alpha=.2)+
  facet_grid(rows="beePop")



flowers<-do.call(rbind,simulations[[1]]$flowers)
flowers$beePop <- 5
flowers$simnum <-1
for(i in 2:length(beePopParms)){
  df <- do.call(rbind,simulations[[i]]$flowers)
  df$beePop <- beePopParms[i]
  df$simnum<-i
  flowers<-rbind(flowers,df)
}


ggplot(flowers,aes(x=time,y=cumulativeVisits,group=ID))+
  geom_line()

ggplot(flowers,aes(x=time,y=microbeK,group=ID))+
  geom_line()




Species1D <- ggplot(flowers,aes(x=time,y=microbesP1,group=ID))+ 
  geom_line(alpha=.2)+scale_y_continuous("Total Cells of Species 1",limits=c(0,75000))
Species2D <- ggplot(flowers,aes(x=time,y=microbesP1,group=ID))+ 
  geom_line(aes(y=microbesP2),color="dark orchid",alpha=.2)+scale_y_continuous("Total Cells of Species 2",limits=c(0,75000))

grid.arrange(Species1D,Species2D)


pdf("MicrobeTotal.pdf",w=8,h=4)
grid.arrange(Species1D,Species2D)
dev.off()

microbeDistribution <- ggplot(flowers[flowers$time==80,],aes(x=microbesP1))+
  geom_histogram(alpha=.5)+
  scale_x_continuous("Total Microbial Cells")+
  geom_histogram(aes(x=microbesP2),fill="darkorchid",alpha=.5)+scale_y_continuous("Number of Flowers")
pdf("MicrobePopulationHist.pdf",w=4,h=4)
microbeDistribution
dev.off()

competitiveexclusion <-
  ggplot(flowers,aes(x=time,y=(microbesP1-microbesP2)/(microbesP1+microbesP2+1),group=ID))+ 
  geom_line(alpha=.1)+
  scale_y_continuous("Microbe Population",breaks=c(-1,0,1),labels = c("All Microbe 2","Even or No Microbes","All Microbe 1"))
pdf("CompetitiveExclusion.pdf",w=12.5,h=6)
competitiveexclusion
dev.off()

ggplot(bees,aes(x=time,y=microbe1,group=ID))+ 
  geom_line() + geom_line(aes(y=microbe2),color="dark orchid")


SumFlowers <- flowers %>% group_by(time,simnum) %>%
  dplyr::summarise(
    nectarVol = mean(necVol),
    microbePop = mean(microbesP1),
    beePop=beePop[1],
    simnum=simnum[1],
    time=time[1])

NectarVolume <- ggplot(SumFlowers, aes(x=time,y=nectarVol,group=(simnum))) + geom_line(aes(color=simnum))
pdf("MeanNectarVolumeChange.pdf",w=5,h=5)
NectarVolume
dev.off()


finalstate <- flowers[flowers$time==max(flowers$time),]
finalstate$beePop <- as.factor(finalstate$beePop)

sugarLearningPlot <- ggplot(finalstate,aes(x=sugar,y=cumulativeVisits)) +
  geom_point(position=position_jitter(.3))+
  geom_smooth(method = "lm", formula = y ~ x + I(x^2),aes(color=as.factor(simnum)),alpha=.1)
sugarLearningPlot
pdf("SugarLearning.pdf",w=6,h=4)
sugarLearningPlot
dev.off()


learningPlotData <- simulations[[1]]$learndata
learningPlotData$simnum <- 1
for(i in 2:length(simulations)){
  new<-simulations[[i]]$learndata
  new$simnum <- i
  learningPlotData<-rbind(learningPlotData,new)
  
}

numlandings<-bees%>%dplyr::group_by(ID)%>% summarise(ID=ID[1],landings=sum(visit))
learningPlotData<-merge(numlandings,learningPlotData)

meanSignal <- flowers%>%group_by(simnum)%>%summarise(CommunitySignal=mean(flowerSignal))
meanThres <- bees%>%group_by(simnum,ID)%>%summarise(LearnedSignal=mean(tail(na.omit(LearnFlowerSig),3)))



learningPlotData<-merge(learningPlotData,meanSignal)

learningPlot <- ggplot(learningPlotData,aes(x=TrueVal,y=Signal))
didbeeslearn<-learningPlot + geom_point() +geom_abline(aes(intercept=0,slope=1)) + 
  geom_smooth(method="lm") +
  scale_x_continuous(limits=c(110,190),"True Signal Threshold\nIn Flowers")
didbeeslearn  
pdf("didbeeslearn.pdf",w=10,h=6)
didbeeslearn
dev.off()
ggplot(bees,aes(x=time,group=ID,y=LearnValue))+geom_line(alpha=.2)
ggplot(flowers[flowers$visit==1,],aes(x=sugar,y=cumulativeVisits))+geom_point()+geom_line(group=ID)
# lets look at some microbe parameters...####
#we are going to use a starting point of 5 pollinator per hundred plants


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

## Do the simulations####

RmParm <- c(-.002,0,.002)
SigDeltParm <- c(0,.002)      

microbe<-data.frame(Population = rep(times=5,x = seq(0,10000,length.out=100)),
                    Rm = rep(each=100,x=c(-.01,-.002,0,.002,.01)),times=5)

microbe$Value <- with(microbe, expr = exp((Rm+Rm*Population))/(1+exp((Rm+Rm*Population))))
RmPlot <- ggplot(microbe, aes(x=Population,y=Value,color=Rm))+geom_line(size=1.2)+facet_wrap(~Rm, nrow=1)+theme_bw()+
  scale_x_continuous("Microbe Population",n.breaks = 3)
RmPlot

pdf("Rm.pdf",w=8,h=4)
RmPlot
dev.off()

Parms1 <- data.frame(expand.grid(RmParm,SigDeltParm))
Parms <- data.frame(Rm=rep(Parms1$Var1,20),Sig=rep(Parms1$Var2,20))
ParmsList <- split(Parms,seq(nrow(Parms))) #need to feed mclapply a list

#simulation2(nbees=1,Rm=0,signalDelt = .001) #just here for trouble shooting...
cores<-detectCores()-2
microbeSensitivity <- mclapply(X=ParmsList,FUN =  function(x){simulation2(nbees=10 ,signalDelt = x$Sig, Rm = x$Rm)},mc.cores = cores) 
save(microbeSensitivity,list="microbeSensitivity",file=".rdata")




sensBees <- do.call(rbind,microbeSensitivity[[1]]$bees)
sensBees$Rm <- Parms[1,"Rm"]
sensBees$sigDelt <- Parms[1,"Sig"]
sensBees$simNum <- 1
for (i in 2:nrow(Parms)){
  d <-do.call(rbind,microbeSensitivity[[i]]$bee)
  d$Rm <- Parms[i,"Rm"]
  d$sigDelt <- Parms[i,"Sig"]
  d$simNum <- i
  sensBees<-rbind(sensBees,d)
}

sensFlowers <- do.call(rbind,microbeSensitivity[[1]]$flowers)
sensFlowers$Rm <- Parms[1,"Rm"]
sensFlowers$sigDelt <- Parms[1,"Sig"]
sensFlowers$simNum <- 1
for (i in 2:nrow(Parms)){
  d <-do.call(rbind,microbeSensitivity[[i]]$flowers)
  d$Rm<- Parms[i,"Rm"]
  d$sigDelt <- Parms[i,"Sig"]
  d$simNum <- i
  sensFlowers<-rbind(sensFlowers,d)
}

end <- sensFlowers[sensFlowers$time==max(sensFlowers$time),]

cumulativeviz <- ggplot(sensFlowers,aes(y=cumulativeVisits,x=time,group=ID))+
  geom_line(alpha=.01)+
  facet_grid(cols=vars(Rm))
cumulativeviz
pdf("cumulativevisits.pdf",w=8,h=4)
cumulativeviz
dev.off()

microbeCumulativeVisits <- ggplot(data=end,aes(x=microbesP1,y=cumulativeVisits))+
  geom_point(alpha=.05)+
  geom_smooth(aes(color=simNum),method=lm)+
  facet_grid(cols=vars(Rm),rows=vars(sigDelt))+
  scale_x_log10()
microbeCumulativeVisits

TastymicrobeCumulativeVisits <- ggplot(data=end[end$Rm==0.002,],aes(x=microbesP1,y=cumulativeVisits))+
  geom_point(alpha=.01)+
  geom_smooth(aes(color=simNum),method=lm)+
  scale_x_log10()+
  facet_grid(cols=vars(Rm),rows=vars(sigDelt))
  
TastymicrobeCumulativeVisits 

End <- end %>%group_by(simNum)%>%dplyr::summarise(
  AverageM1 = mean(microbesP1),
  AverageM2 = mean(microbesP2),
  AverageM1pul = mean(microbesP1)/tail(necVol,1),
  AverageM2pul = mean(microbesP2)/tail(necVol,1),
  NumM1 = length(microbesP1[microbesP1>0]),
  NumM2 = length(microbesP2[microbesP2>0]),
  NumM0 = length(microbesP1[(microbesP1+microbesP2)==0]),
  Num   = length(microbesP1),
  Rm=Rm[1],
  sigDelt=sigDelt[1])

endplot<-ggplot(data=End, aes(x=Rm))
DensPlot <- endplot + geom_point(aes(y=AverageM1),position=position_dodge2(.0001)) +
  geom_point(aes(y=AverageM2),color="dark orchid",position=position_dodge2(.0001))+
  geom_smooth(aes(y=AverageM1),method="lm",color="black")+
  geom_smooth(aes(y=AverageM2),method="lm",color="dark orchid")+
  facet_grid(~sigDelt)
DensPlot

DensPlot <- endplot + geom_point(aes(y=AverageM1pul),position=position_dodge2(.0001)) +
  geom_point(aes(y=AverageM2pul),color="dark orchid",position=position_dodge2(.0001))+
  geom_smooth(aes(y=AverageM1pul),method="lm",color="black")+
  geom_smooth(aes(y=AverageM2pul),method="lm",color="dark orchid")+
  
  facet_grid(~sigDelt)+
  scale_y_log10()
DensPlot
pdf("AverageDensity.pdf",w=8,h=4)
DensPlot
dev.off()

numberflowers <- endplot+
  geom_point(aes(y=NumM1,x=Rm-.0001),                   ,position=position_jitter(.0001))+
  geom_point(aes(y=NumM2, x=Rm+.0001),color="dark orchid",position=position_jitter(.0001))+
  geom_smooth(aes(y=NumM1),method="lm",color="black",se=F)+
  geom_smooth(aes(y=NumM2),method="lm",color="dark orchid",se=F)+
  facet_grid(~sigDelt)
numberflowers
pdf("numflw.pdf",w=8,h=4)
numberflowers
dev.off()

learningsignals <- ggplot(sensBees, aes(x=LearnMicrobeSig,y=LearnValue))+
  geom_smooth(method="lm",aes(color=ID),alpha=.01,se=F)+
  theme(legend.position="none")+
  facet_grid(cols=vars(Rm))
learningsignals
pdf("learningsignals.pdf",w=8,h=4)
learningsignals
dev.off()

learningsignals2 <- ggplot(sensBees, aes(x=LearnMicrobeSig,y=LearnValue))+
  stat_smooth(geom='line', aes(group=ID),alpha=0.1, se=FALSE, method="lm")+
  geom_smooth(method="lm",color="darkred",size=1.4)+
  geom_abline(intercept = .5, slope=0,color="grey",size=.8,linetype="dashed")+
  theme(legend.position="none")+
  
  facet_grid(cols=vars(Rm))
learningsignals2

pdf("learningplot2.pdf",w=10,h=6)
learningsignals2
dev.off()

TimeToCol <- sensFlowers %>% group_by(ID,simNum) %>% dplyr::summarise(
  MaxM1 = max(microbesP1),
  MaxM2 = max(microbesP2),
  MinM1 = min(time[microbesP1>0]),
  MinM2 = min(time[microbesP2>0]),
  Rm = Rm[1],
  sigDelt = sigDelt[1]
  )
TimeToCol$MinM1[TimeToCol$MinM1==Inf]<-200
TimeToCol$MinM2[TimeToCol$MinM2==Inf]<-200

MeltTimetoCol <- melt(TimeToCol,measure.vars =  c("MinM1","MinM2"))

timetocol <- ggplot(MeltTimetoCol, aes(x=Rm,color=variable,y=value))+
  geom_violin(aes(x=(Rm),group=interaction(Rm,variable),fill=variable),alpha=.2)+
  geom_point(alpha=.1,position=position_jitterdodge(jitter.widt=.0005,dodge.width=.001))+
  #geom_point(aes(y=MinM1),position=position_jitterdodge(.2,.2),alpha=.2)+
  #geom_point(aes(y=MinM2),,position=position_jitter(.0002),color="grey",alpha=.3)+
  geom_smooth(method="glm", method.args = list(family = "poisson"),se=T)+
  #geom_smooth(aes(y=MinM2),color="grey",method="lm")+
  facet_grid(~sigDelt)+
  scale_y_continuous("Minimum Time to First Colonization")+
  scale_x_continuous("Microbe Impact on Preference")+
  scale_color_manual(values=c("black","dark orchid"))+
  scale_fill_manual(values=c("black","dark orchid"))
  timetocol
pdf("timetocolonization.pdf",w=8,h=4)
timetocol
dev.off()
t1 <- lm(data=MeltTimetoCol[],value~variable*Rm)
summary(t1)
t<-glm(data=TimeToCol[],MinM1~Rm*sigDelt,family="poisson")
summary(t)

Parms1$Predict<-exp(predict(t,newdata=data.frame(Rm=Parms1$Var1,sigDelt=Parms1$Var2)))

FlowersVisited <- end[end$cumulativeVisits>0,"ID"]
sum(sensFlowers$ID%in%FlowersVisited)

ProbPlot<-ggplot(sensFlowers[sensFlowers$ID%in%FlowersVisited,],aes(y=visit,x=microbesP1))+
  geom_point(position=position_jitter(h=.2),alpha=.01)+
  geom_smooth(method="glm",method.args=list(family="binomial"))+
  facet_grid(~Rm)
#ProbPlot 



SigPlot<-ggplot(sensFlowers,aes(x=microbeSignal1,y=visit))+

  #geom_point(alpha=.1,position=position_jitter(h=.2))+
  geom_smooth(method="glm",method.args=list(family="binomial"))+
  facet_grid(~Rm)
pdf("signalPlot.pdf",w=12.5,h=12.5/2)
SigPlot
dev.off()
library(lme4)
SignalModel <- glm(data=sensFlowers[sample(nrow(sensFlowers),1000),],visit~microbeSignal1*as.factor(Rm)*sigDelt,family="binomial")
drop1(SignalModel,test="Chisq",.~.)

microbeSignal1<-seq(1,250,by=10)
Rm <- c(-0.002,0,0.002)
sigDelt<-c(0,0.02)
newdata<-data.frame(expand.grid(microbeSignal1=microbeSignal1,Rm=Rm,sigDelt=sigDelt))

logit2prob <- function(logit){
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
}

SignalModel <- glm(data=sensFlowers, visit~microbeSignal1*sigDelt,family="binomial")
drop1(SignalModel,test="Chisq")
newdata$predict <- logit2prob(predict(SignalModel,newdata=newdata))
predictionVisit <- ggplot(newdata,aes(x=microbeSignal1,y=predict,color=as.factor(Rm)))+
  geom_line()+
  facet_grid(~sigDelt)+
  scale_x_continuous("Microbial Signal")+
  scale_y_continuous("Predicted Probability of Visit")
predictionVisit

ggplot(sensBees, aes(x=time)) + 
  geom_line(aes(y=microbe1,group=ID),alpha=.2)+
  geom_line(aes(y=microbe2,group=ID),alpha=.2,color="dark orchid")+
  facet_grid(~Rm)

cumulativeVisits <- sensFlowers %>% group_by(ID) %>% dplyr::summarise(
  cumulativeVisits = cumsum(visit),
  time=time,
  Rm=Rm,
  sigDelt=sigDelt,
  microbe1Ratio =(microbesP1-microbesP2)/(microbesP1+microbesP2),
  pollenF1=pollenF1
  )

cumulativeVisits <- cumulativeVisits[cumulativeVisits$time==200,]

pdf("visitsbyMicrobe1.pdf",w=12.5,h=8)
ggplot(cumulativeVisits,aes(x=microbe1Ratio,y=cumulativeVisits,group=ID))+
  geom_point(aes(color=microbe1Ratio))+facet_grid(cols=vars(Rm),rows=vars(sigDelt))+geom_smooth(method = "lm",color="black")
dev.off()

ggplot(sensBees,aes(x=time,y=microbe1,group=ID))+ 
  geom_line(alpha=.02) #+ #geom_line(aes(y=microbe2),color="dark orchid",alpha=.02)+facet_grid(~Rm)

ggplot(cumulativeVisits,aes(x=microbe1Ratio,y=cumulativeVisits))+
  geom_point(aes(color=microbe1Ratio))+
  facet_grid(cols=vars(Rm),rows=vars(sigDelt)) +
  scale_y_log10()+
  geom_smooth(method="lm")


ggplot(cumulativeVisits,aes(x=microbe1Ratio,y=pollenF1)) +geom_point()+facet_grid(cols=vars(Rm),rows=vars(sigDelt))+geom_smooth(method="lm")

cumulativeBees <- sensBees%>%group_by(ID)%>%summarise(cumulativeVisits=cumsum(visit))
sensBees$cumulativeVisits<-cumulativeBees$cumulativeVisits
ggplot(sensBees,aes(x=time,y=cumulativeVisits,group=ID))+
  geom_line(aes())+
  facet_grid(rows=vars(sigDelt),cols=vars(Rm))

visitsbytime <- sensFlowers%>%group_by(time,Rm,sigDelt)%>%summarize(
  visits=sum(visit),
  Rm=Rm[1],
  sigDelt=sigDelt)

ggplot(visitsbytime,aes(x=time,y=visits))+geom_point()+facet_grid(cols=vars(Rm),rows=vars(sigDelt))+scale_y_log10()











####TRASH BELOW HERE 

onerun <-   function(sigDel, sugDel){
  
  # make a community of pollinators none of which have microbes! 
  #This is a list of bees. Each bee is a data frame
  nPol    <- 4   #number of indiv per species
  Pollinators <- list()
  for(i in 1:nPol){
    Pollinators[[i]] <- BeeInit(thresholdPar = c(0,1),
                                microbesPar = c(0,200),
                                preferencesPar = c(10,100)
    )
  }  
  
  # make a community of sterile flowers. this is a list of flowers and each is a data frame
  nFlow     <- 20
  Flowers <- list()
  for(i in 1:nFlow){
    Flowers[[i]] <- FlowerInit(signalPar = c(0,100))
  }
  
  # make a microbe with the parameters provided. we only have a single microbe right now
  ourMicrobe <- MicrobeInit(growthRate = 1.05,signalDelt = sigDel, sugarDelt = sugDel)
  
  for(i in 1:100){
    #first figure out which flowers bees will go to by randomly sampling as many flowers as there are bees
    encounters <- sample(1:length(Flowers), size = length(Pollinators), replace=F)
    encounterList <- vector("list",length(Flowers))
    encounterList[encounters] <- Pollinators
    
    test <- mapply(FUN=function(x,y){
      tmp <- daytimestep(flower=x,bee=y,microbe=ourMicrobe)
    },Flowers,encounterList,SIMPLIFY=F)
    
    Pollinators <- lapply(X = test[encounters], FUN = function(x){return(x$bee)})
    Flowers <- lapply(X=test, FUN = function(x){
      if(!is.null(x$flower)){return(x$flower)} 
      else(return(x))})
    
  }  
  
  #calculate cumulative visits!
  Flowers <- lapply(X=Flowers, FUN= function(x){x$cumulativeVisits<-cumsum(x$visit)
  x$signaldelt <- sigDel
  x$sugardelt <- sugDel
  return(x)})
  
  Pollinators <- lapply(X=Pollinators, FUN= function(x){x$cumulativeVisits<-cumsum(x$visit)
  x$signaldelt <- sigDel
  x$sugardelt <- sugDel
  return(x)})
  
  updatedFlowers<- do.call(rbind,Flowers)
  updatedBees <- do.call(rbind,Pollinators)
  return(list(flowers=updatedFlowers,bees=updatedBees))
}               
onerun(.1,.1)

#try parallel processing
detectCores()
Sys.time()
test<- mclapply(X = ParmsList, 
                FUN = function(x){onerun(sigDel = x[1,1],sugDel = x[1,2])},
                mc.cores = 11)
Sys.time()

flr<-test[[1]]$flowers
for(i in 2:length(test)){flr<-rbind(flr,test[[i]]$flowers)}

bbz<-test[[1]]$bees
for(i in 2:length(test)){bbz<-rbind(bbz,test[[i]]$bees)}
dev.new()      

flr$ID <- as.factor(flr$ID)

timetofirstvisit <-flr[flr$cumulativeVisits==1,]
timetofirstvisit <- ddply(.data = timetofirstvisit, .variables = .(ID),  function(x){x[x$time==min(x$time),]})
signalbytimetovis <- ggplot(timetofirstvisit, aes(x=signaldelt,y=time)) 
signalbytimetovis + geom_point(alpha=.2,position=position_jitter(width=.00002))
sugarbytimetovis <- ggplot(timetofirstvisit, aes(x=sugardelt,y=time)) 
sugarbytimetovis + geom_point(alpha=.2,position=position_jitter(width=.00002))

summary(timetofirstvisit)

firstVisSig <- glm(data=timetofirstvisit, time~signaldelt*sugardelt, family="poisson")
summary(firstVisSig)

summaryTtoV <- ddply(.data = timetofirstvisit, .variables = .(sugardelt,signaldelt), summarize, time=mean(time) )
ggplot(summaryTtoV, aes(x=signaldelt,y=time,color=(sugardelt),group=as.factor(sugardelt))) + geom_line()




names(Parms) <- c("signaldelt","sugardelt")
Parms$TtoFirst <- predict(firstVisSig,newdata = Parms,se.fit = T, type="response")$fit
Parms$TtoFirstSE <- predict(firstVisSig,newdata = Parms,se.fit = T,type="response")$se.fit
ggplot(Parms, aes(x=signaldelt,y=sugardelt,color=(TtoFirst),fill=TtoFirst)) +geom_tile()



ggplot(Parms[Parms$sugardelt==-0.010,], aes(x=signaldelt,y=TtoFirst)) + geom_point() +
  geom_errorbar(aes(ymin=TtoFirst-TtoFirstSE, ymax=TtoFirst+TtoFirstSE))

ggplot(Parms, aes(color=as.factor(sugardelt),x=signaldelt,y=TtoFirst)) + geom_point() 

maxvisits <- flr[flr$time==max(flr$time),]

maxvisitsplot <- ggplot(maxvisits, aes(x=signaldelt,y=cumulativeVisits, color=signaldelt))
maxvisitsplot + geom_point(alpha=.2,position=position_jitter(width=.0002))

summary(maxvisits)
table(maxvisits$cumulativeVisits)
signal <- glm(data=maxvisits, cumulativeVisits~signaldelt*sugardelt, family = "poisson")
summary(signal)    

ggplot(maxvisits, aes(x=signaldelt,y=cumulativeVisits, color=as.factor(sugardelt)))+
  geom_point()+
  facet_grid(vars(sugardelt))


Parms$Max <- predict(signal,newdata = Parms[,c("signaldelt","sugardelt")], type="response")
ggplot(Parms, aes(x=signaldelt,y=sugardelt,color=(Max),fill=Max)) +geom_tile()


pollenplot <- ggplot(maxvisits, aes(x=signaldelt,y=signaldelt, color=pollenF1))
pollenplot + geom_point(alpha=.2,position=position_jitter(width=.0002))



pollenmod <- lm(data=maxvisits, pollenF1~as.factor(signaldelt)*as.factor(sugardelt))
summary(pollenmod)

plotdata <- data.frame(emmeans(pollenmod, specs=c("signaldelt","sugardelt")))
ggplot(plotdata,aes(x=as.factor(signaldelt), color=as.factor(sugardelt),y=emmean)) +
  geom_point(position=position_dodge(.4)) +
  geom_line(aes(group=as.factor(sugardelt)),position=position_dodge(.4)) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), position=position_dodge(.4))+
  scale_y_continuous("Predicted pollen receipt")




Parms$Pollen <- predict(pollenmod,newdata = Parms[,c("signaldelt","sugardelt")], type="response")
ggplot(Parms, aes(x=signaldelt,y=Pollen,color=as.factor(sugardelt),fill=Pollen)) + geom_point()+geom_line(aes(group=as.factor(sugardelt)))




