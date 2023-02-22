#Critical functions ####

#Make agents####
BeeInit <- function(thresholdPar = c(0,100),        # we are going to specify the min and max bees have a threshold for each flower species for now just one 
                    microbesPar  = c(0,200),        # the microbial load for now just one species
                    pollenLoadPar = 0,       # the pollen load for now just one species
                    preferencesPar = c(30,.2)) # the initial floral preferences for sugar level and alkaloid
{
  
  time <- 1
  thresholdF1   <- runif(1,thresholdPar[1],thresholdPar[2])
  thresholdF1   <- runif(1,thresholdPar[1],thresholdPar[2])
  microbe1      <- runif(1,microbesPar[1],microbesPar[2])
  microbe2      <- runif(1,microbesPar[1],microbesPar[2])
  pollenLoadF1  <- pollenLoadPar[1]
  sugarPref     <- preferencesPar[1]
  alkaloidPref  <- preferencesPar[2]
  volPref       <- 5
  visit         <- 0
  ID <- paste(sample(size=10, x=letters),sep="",collapse="")
  LearnFlowerSig <- thresholdF1
  LearnMicrobeSig <- 0
  LearnRew <- 0
  LearnVol <- 0
  LearnValue <- 0
  # this is the value function of what they expect for their preferred flower, should be very close to 1...
  
  
  return( data.frame(
    time = time,
    thresholdF1   = thresholdF1,
    thresholdM1   = thresholdF1,
    microbe1      = microbe1,
    microbe2      = microbe2,
    microbe1avail = 0,
    microbe2avail = 0,
    pollenLoadF1  = pollenLoadF1,
    sugarPref     = sugarPref,
    alkaloidPref  = alkaloidPref,
    volPref       = volPref,
    visit         = visit,
    ID = ID,
    LearnFlowerSig = thresholdF1,  
    LearnMicrobeSig = 0,
    LearnRew = sugarPref,
    LearnVol = LearnVol,
    LearnValue = LearnValue,
    coef.Intercept = NA,
    coef.LearnMicrobeSig = NA,
    coef.LearnFlowerSig = NA,
    coef.phi=NA)  )                                           #these are the bees previous experiences.
}
#BeeInit()   # wont work without assess... it is defined later...        

# Initialize a flower
FlowerInit <- function(species ="flower1",             # flower species is coded as a character string
                       signalPar = c(0,100),              # this is a linear measure of a signal (TODO think of multiple signal axes?)       
                       microbesP = c(0,0),            #mean and SD of microbes starting
                       sugarPar = c(30,10),
                       alkaloidPar = c(.2,.01),  # mean and standard deviation for sugar then alkaloid
                       pollenLoadPar = 0,       # starting pollen load
                       pollenDepRate = 10,      # the amount of pollen deposited by a plant
                       maxNectar = 15,
                       refillRate = .025)             
{
  necVol <- runif(1,0,maxNectar)
  
  microbeSignal = 0
  return( data.frame(
    time = 1,
    species = species,
    flowerSignal = runif(1,min = signalPar[1], max = signalPar[2]),
    microbeSignal1 = 0, # will be updated at first growth
    microbeSignal2 = 0, # will be updated at first growth
    microbesP1 = abs(rnorm(1,microbesP[1],microbesP[2])),
    microbesP2 = abs(rnorm(1,microbesP[1],microbesP[2])),
    microbeK = 5000 * necVol, #900 cells per ul # https://bsapubs.onlinelibrary.wiley.com/doi/pdfdirect/10.1002/ajb2.1834 says higher but our data suggest otherwise?
    necVol = necVol,
    maxNectar = maxNectar,
    refillRate=refillRate,
    sugar = rnorm(1, mean = sugarPar[1], sd = sugarPar[2]),
    alkaloid = rnorm(1, mean= alkaloidPar[1],alkaloidPar[2]),
    pollenF1 = pollenLoadPar[1],
    pollenDepRate = pollenDepRate,
    visit=0,
    ID = paste(sample(size=10, x=letters),sep="",collapse="")
  )
  )
}
FlowerInit()
#Initialize a microbe
MicrobeInit <- function(species = "microbe1",
                        growthRate = 1.2,
                        signalDelt = .002,
                        Rm = -.001,
                        sugarDelt = .002,
                        alkaloidDelt = .001){
  return(data.frame (
    species = species,
    growthRate = growthRate,
    Rm = Rm,
    signalDelt = signalDelt,
    sugarDelt = sugarDelt,
    alkaloidDelt = alkaloidDelt
  ))
}



#Time step functions####
# microbial growth function

# #for debugging functions...
# microbe <- MicrobeInit()
# flower <- FlowerInit()
# bee <- BeeInit()

grow <- function(microbe1,microbe2,flower){
  # define some times
  prev <- max(flower$time)
  current <- prev + 1
  flower[current,]<-flower[prev,]
  flower$time[current]<- current
  
  # flowers secrete nectar if they are below max and we calculate capacity
  
  
  if(flower[prev,"necVol"] < flower[prev,"maxNectar"]){
    flower[current,"necVol"] <- flower[prev,"necVol"]+flower[prev,"refillRate"]}
  if(flower[current,"necVol"] > flower[current,"maxNectar"]){flower[current,"necVol"] <- flower[current,"maxNectar"]}
  flower[current,"microbeK"] <- flower[current,"necVol"] * 5000
  
  #first microbes grow using standard logistic growth 
  r1 <- microbe1$growthRate
  P1 <- flower[prev,"microbesP1"] # REMEMBER TO CHANGE IF WE ADD MORE MICROBES!
  r2 <- microbe2$growthRate
  P2 <- flower[prev,"microbesP2"]
  K <- flower[current,"microbeK"]
  flower[current,"microbesP1"] <- P1*(1+r1*(1-(P1+1.2*P2)/K)) #here the competition coefficients are equal to 1.2 so inter is a little higher than intra specific competition
  flower[current,"microbesP2"] <- P2*(1+r2*(1-(P2+1.2*P1)/K))
  
  if(flower[current,"microbesP1"]<0){flower[current,"microbesP1"]<-0}
  if(flower[current,"microbesP2"]<0){flower[current,"microbesP2"]<-0}
  
  # then microbes change the signal of the plant 
  flower[current,"flowerSignal"]<- flower$flowerSignal[1]  # could add faint nectar smell here
  flower[current,"microbeSignal1"]<-microbe1$signalDelt * flower[current,"microbesP1"]
  flower[current,"microbeSignal2"]<-microbe2$signalDelt * flower[current,"microbesP2"]
  
  #then microbes change the quality of the plant in this simulation the sugar content and
  #alkaloid content are static products that microbes act on
  flower[current,"sugar"] <- flower$sugar[1] + microbe1$sugarDelt* flower[current,"microbesP1"] +microbe2$sugarDelt* flower[current,"microbesP2"]
  flower[current,"alkaloid"] <- flower$alkaloid[1] + microbe1$alkaloidDelt * flower[current,"microbesP1"] + microbe2$alkaloidDelt * flower[current,"microbesP2"]
  flower$visit[current] <- 0 # this will get updated during the visit function
  #then return the changed flower
  return(flower)
}


# Value function
# calculate a reward function, here bees discount a flower if it has less than their 
# total capacity ie volPref very steeply. They have a linear increase in preference to
# as sugar increases from 0 - 45, then a steep decrease to 0 at 55
# they put twice as much weight on nectar sugar as they do volume

VolumeValFun <- function(bee,flower){
  # this gives us a value of ~0-1, it saturates towards 1 volume approaches the bees preference
  # the percentage before the expression in the exponent controls how fast that preference changes
  1/(1+exp(0.4*(-1*flower$necVol+bee$volPref-6)))
} 
VolumeValFun(bee=data.frame(volPref=20), flower=data.frame(necVol=8))

SugarValFun  <- function(bee,flower){
  #this is linear increase in preference up to a parameter, then steep drop to 0 preference parameter+10
  pref <-  bee$sugarPref[nrow(bee)]
  avail <- flower$sugar[nrow(flower)]
  
  if(avail <= pref & avail > 0)  {return(1/pref*avail)}
  if(avail > pref & avail <= pref+10) {return((-1/10)*(avail-(pref+10)))}
  if(avail > pref+10)             {return(0)}
  if(avail < 0 )               {return(0)}
}

SugarValFun(bee=data.frame(sugarPref=12), flower=data.frame(sugar=13))

MicrobeValFun <- function(microbe1,microbe2,flower){
  M1Dens<- with(tail(flower,1),microbesP1/necVol)
  M2Dens<- with(tail(flower,1),microbesP2/necVol)
  M1val <- exp((microbe1$Rm[1]+microbe1$Rm[1]*M1Dens))/(1+exp((microbe1$Rm[1]+microbe1$Rm[1]*M1Dens)))
  #M2val <- 1/(2+exp(microbe2$Rm[1]*(5000-M2Dens)))
  # return(mean(M1val,M2val))
  return(M1val)
}

#  
# testMic <- MicrobeInit(Rm=0.002)
# 
# MicrobeValFun(mic1,mic2,flowerT)
# flowerT$microbeValue<-NA
# microbe1<-MicrobeInit(Rm = .001)
# microbe2 <- MicrobeInit(Rm=.0)
# for (i in 1:nrow(flowerT)){flowerT[i,"microbeValue"] <- (MicrobeValFun(microbe1,microbe2,flowerT[i,]))}
# ggplot(flowerT,aes(x=time)) + geom_line(aes(y=microbesP1/necVol/5000)) + geom_line(aes(y=microbesP2/necVol/5000),color="dark orchid")+geom_line(color="red",aes(y=microbeValue))


assess <- function (bee,flower,microbe1,microbe2){
  #CHANGE THIS IF WE WANT BEES TO TAKE SOMETHING ELSE INTO ACCOUNT!
  # Syntax Value <- assess(bee=bee[currentTime,],flower=flower[currentTime,])
  
  volumeValue <- VolumeValFun(bee=bee[nrow(bee),], flower=flower[nrow(flower),])
  sugarValue  <- SugarValFun(bee=bee[nrow(bee),],flower=flower[nrow(flower),])
  micValue <- MicrobeValFun(microbe1,microbe2,flower)
  
  # they put twice as much weight on nectar sugar as they do volume (i.e. rep 2)
  return(mean(c(volumeValue,sugarValue,micValue)))
  #return(mean(c(rep(sugarValue,2))))} #took out volume
}




# decision function
decide <- function(bee, flower, microbe1, microbe2){
  # bee senses a flower, decides whether to visit based on its threshold
  # flowsp <- flower$species # may need this when we add more species
  currentTime <- max(flower$time)
  
  #build a regression for flower and microbe signal
  
  experience <- bee[bee$visit==1,] |> tail(30)
  NumberOfLandings <- nrow(experience)
  
  
  if(NumberOfLandings<=3){
    return (visit(bee=bee,flower=flower,microbe1 = microbe1, microbe2 = microbe2))
  }  
  
  # see https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1439-0310.2006.01174.x?casa_token=o_kiaZSxTLQAAAAA:ojoLikVyLYlcP81oB2JOf69ZH1pBiqheJSLgCEcKD2cZdL-qRTuXfnALXk5CE3EOo-sb4s3C741qCu4
  #step 2, bees build a linear model about the world based on what they know
  # and the information they just got. They weight previous experience higher if it was more 
  # recent (3X more), and their most recent experience is even higher weighted (2X highest previous)
  
  if(NumberOfLandings>3){
    #on 1/24 I commented out the 2 separate models and built a multivariate beta regression model
    # it cant handle all of one variable being the same number so we need to check if any of the data are rank deficient
    # 1) if both predictors are rank deficient the polliantor doesnt have enough info and decides to land to get more (optomistic :))
    if(var(experience$LearnMicrobeSig) ==0 & var(experience$LearnFlowerSig) ==0){
      return (visit(bee=bee,flower=flower,microbe1 = microbe1, microbe2 = microbe2))
    }
    
    # 2) If either predictor is rank deficient then construct a model without it
    predictors <- NA
    if(var(experience$LearnMicrobeSig) >0 & var(experience$LearnFlowerSig) ==0){predictors <- "LearnMicrobeSig"}
    if(var(experience$LearnMicrobeSig) ==0 & var(experience$LearnFlowerSig) > 0){predictors <- "LearnFlowerSig"}
    if(var(experience$LearnMicrobeSig) >0 & var(experience$LearnFlowerSig)  >0){predictors <- "LearnMicrobeSig+LearnFlowerSig"}
    
    
    experience$wgts <- c(seq(from=1,to=3,length.out=NumberOfLandings-1),6)
    
    #these were some processes to try and keep track of the coefficients from the previous
    #decide model and use them as a starting point for the optimizers the next time a bee had to decide.
    # startingCoef <- experience[nrow(experience),grep("coef",names(experience))]
    # names(startingCoef)<-sub("coef.","",names(startingCoef))
    # startingCoef <- startingCoef[,c("Intercept",predictors,"phi")]
    
    # ficients <- coef(decideMod)
    # 
    # names(ficients) <- sub(pattern = "\\(","",names(ficients)) 
    # names(ficients) <- sub(pattern = "\\)","",names(ficients))
    # 
    # bee[currentTime,paste("coef",names(ficients),sep=".")] <- ficients
    # 
    
    decideMod <- betareg(data=experience,paste("LearnValue~",predictors,sep=""), 
                         weights = wgts, control=betareg.control(phi=F,
                                                                 maxit=10000 ))
    
    prediction <- predict(decideMod,
                          newdata=data.frame(
                            LearnFlowerSig=flower$flowerSignal[currentTime],
                            LearnMicrobeSig=flower$microbeSignal1[currentTime]+flower$microbeSignal2[currentTime])
    )
    
    # decideflowerMod<-lm(data=experience, LearnValue~LearnFlowerSig, weights=c(seq(1,3,length.out=NumberOfLandings-1),6))
    # decidemicrobeMod<-lm(data=experience, LearnValue~LearnMicrobeSig, weights=c(seq(1,3,length.out=NumberOfLandings-1),6))
    # 
    # FlowPred <- predict(object=decideflowerMod,
    #          newdata=data.frame(LearnFlowerSig=flower$flowerSignal[currentTime]))
    # 
    # MicPred <- predict(object = decidemicrobeMod, 
    #          newdata=data.frame(LearnMicrobeSig=flower$microbeSignal1[currentTime]+flower$microbeSignal2[currentTime])) #COMEBACK HERE ON 1/18/23
    
    #Nigel Raine
    # second half of if statement makes a small chance (8%) that bees visit flowers that are below their threshold
    chance <- runif(1,0,100)   
    if(prediction >= .50 | chance<=10){ 
      return (visit(bee=bee,flower=flower,microbe1 = microbe1, microbe2 = microbe2))
    }
    if(prediction<.5){
      noChngBee <- bee[currentTime-1,]
      noChngBee$time <- noChngBee$time+1
      noChngBee$LearnRew<- NA
      noChngBee$LearnFlowerSig<-NA
      noChngBee$LearnVol <- NA
      bee[currentTime,] <- noChngBee
      flower$visit[currentTime] <- 0
      bee$visit[currentTime] <- 0
      
      # no need to change flower because we already updated it...
      return(list(bee=bee,flower=flower))
    }
  }
}



# visit function 
bee <- BeeInit()
visit <- function(bee,flower,microbe1,microbe2){
  
  #bee gets microbes, assess flower, updates threshold
  
  ## flowsp <- flower[["species"]]   # may need if we add more species of flowers
  currentTime <-  max(flower$time)
  prevBee <- bee[currentTime-1,]
  prevFlow <- flower[currentTime,]
  
  bee[currentTime,] <- bee[currentTime-1,]
  bee[currentTime,"time"]  <- currentTime
  bee[currentTime,"visit"] <- 1
  flower[currentTime,"visit"] <- 1
  
  
  bee[currentTime,"LearnRew"]    <- prevFlow$sugar
  bee[currentTime,"LearnVol"]    <- prevFlow$necVol
  bee[currentTime,"LearnValue"]  <- assess(bee=prevBee,flower=prevFlow,microbe1=microbe1, microbe2=microbe2)
  bee[currentTime,"LearnFlowerSig"]    <- prevFlow$flowerSignal
  bee[currentTime,"LearnMicrobeSig"]    <- prevFlow$microbeSignal1+prevFlow$microbeSignal2
  #transfer microbes 
  
  
  bee[currentTime,"microbe1"] <-  .95 * prevBee$microbe1 + .8 * prevFlow$microbesP1 # amount bee will get, notice they eat about 92% of microbes so not all are avail for dispersal
  flower[currentTime,"microbesP1"] <- .2 * prevFlow$microbesP1 + .05 * prevBee$microbe1avail  # amount bee will deposit notice that these decrease by the amount transfered in other steps (e.g. 1-val here)
  
  bee[currentTime,"microbe2"] <-  .95 * prevBee$microbe2 + .8 * prevFlow$microbesP2 # amount bee will get, notice they eat about 92% of microbes so not all are avail for dispersal
  flower[currentTime,"microbesP2"] <- .2 * prevFlow$microbesP2 + .05 * prevBee$microbe2avail  # amount bee will deposit notice that these decrease by the amount transfered in other steps (e.g. 1-val here)
  
  totalMicrobe <- bee[currentTime,"microbe1"] + bee[currentTime,"microbe2"]
  if(totalMicrobe>0){
    bee[currentTime,"microbe1avail"] <- bee[currentTime,"microbe1"]/totalMicrobe * 200 # here bees consume most of the microbes they pick up and only 200 total microbes stick to their proboscis and are avail for transfer (proportional to how many they consume)
    bee[currentTime,"microbe2avail"] <- bee[currentTime,"microbe2"]/totalMicrobe * 200
  }
  
  flower[currentTime,"necVol"] <- .2 * prevFlow$necVol # pollinator takes 80% of the nectar
  
  #Bees assess flower (todo add volume) and update threshold
  
  #check out https://www.sciencedirect.com/science/article/pii/S0896627319308402 for bayesian decision model
  #TODO make a function for how threshold moves i.e. bayesian updating
  
  #step 1 - look back at the bees experience. 
  #right now we will take up to the 20 most recent visits, but this could change
  #either based on literature or w/e
  
  experience <- bee[bee$visit==1,] |> tail(30)
  experience$wgts <- c(seq(from=1,to=3,length.out=nrow(experience)-1),6)
  
  
  # see https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1439-0310.2006.01174.x?casa_token=o_kiaZSxTLQAAAAA:ojoLikVyLYlcP81oB2JOf69ZH1pBiqheJSLgCEcKD2cZdL-qRTuXfnALXk5CE3EOo-sb4s3C741qCu4
  #step 2, bees build a linear model about the world based on what they know
  # and the information they just got. They weight previous experience higher if it was more 
  # recent (3X more), and their most recent experience is even higher weighted (2X highest previous)
  # if(nrow(experience)>3){
  #     flowerMod<-lm(data=experience, LearnFlowerSig~LearnValue, weights=c(seq(1,3,length.out=nrow(experience)-1),6))
  #     microbeMod<-lm(data=experience, LearnMicrobeSig~LearnValue, weights=c(seq(1,3,length.out=nrow(experience)-1),6))
  #     
  #     
  #     valueMod <- betareg(data=experience,LearnValue~LearnFlowerSig+LearnMicrobeSig,weights=wgts)
  #     # then we predict what signal they should use as their threshold based on the model
  #     
  #     FlowerSigPrediction <- predict(object = flowerMod, 
  #                                              newdata = data.frame(LearnValue=.5))
  #     
  #     MicrobeSigPrediction <- predict(object = microbeMod, 
  #                                    newdata = data.frame(LearnValue=.5))
  #     
  #     bee[currentTime,"thresholdF1"] <- FlowerSigPrediction
  #     bee[currentTime,"thresholdM1"] <- MicrobeSigPrediction
  # }
  #transfer pollen 
  bee[currentTime,"pollenLoadF1"] <- prevBee$pollenLoadF1 + prevFlow$pollenDepRate - .3 * prevBee$pollenLoadF1
  flower[currentTime,"pollenF1"] <- prevFlow$pollenF1     + .3 * prevBee$pollenLoadF1
  
  
  return(list(bee=bee,flower=flower))      
  
}



#
# time step 

daytimestep <- function(bee = NULL, flower=NULL, microbe1=NULL, microbe2=NULL){
  # microbes grow 
  
  
  flower1 <- grow(flower = flower, microbe1 = microbe1, microbe2=microbe2)
  
  # bees forage
  
  if(is.null(bee) == F){
    return(decide(bee = bee, flower = flower1, microbe1=microbe1, microbe2=microbe2))
  }
  
  if(is.null(bee)==T){
    return(flower1)
  }
  ###NOT FOR NOW 
  # flowers get closer to death
  # new flowers open
}




nighttimestep <- function(bee, flower){
  #microbes grow 
  #flowers get closer to death
  #new flowers open
}


####end functions