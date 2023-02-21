set.seed(612211323)

LearnMicrobeSig<-runif(100,0,1500)
LearnFlowerSig<-runif(100,0,100)
LearnValue<-runif(100,0,1)


decideMod <- betareg::betareg(LearnValue~LearnFlowerSig+LearnMicrobeSig, weights = c(seq(1,3,length.out=nrow(experience)-1),6))
summary(decideMod)


plotdata<-data.frame(LearnMicrobeSig,LearnFlowerSig,LearnValue)
ggplot(plotdata, aes(x=LearnMicrobeSig,y=LearnValue))+geom_point()


newMicrobeSig<-seq(0,1500,by=100)
newFlowerSig<-seq(0,100,by=10)

NewData <- expand_grid(LearnFlowerSig=newFlowerSig, LearnMicrobeSig=newMicrobeSig)

  NewData$Prediction <- predict(decideMod,newdata = NewData)

  ggplot(NewData,aes(x=LearnFlowerSig,y=Prediction))+geom_point(aes(size=LearnMicrobeSig))
  
  summary(decideMod)
  

bee<-BeeInit()
flower<-FlowerInit(
)
microbe1<-MicrobeInit()
microbe2<-MicrobeInit()

flower<-grow(microbe1 = microbe1, microbe2 = microbe2, flower= flower)

experience<-data.frame(LearnMicrobeSig,LearnFlowerSig,LearnValue)
if(var(experience$LearnMicrobeSig) ==0 & var(experience$LearnFlowerSig) ==0){
  return (visit(bee=bee,flower=flower,microbe1 = microbe1, microbe2 = microbe2))
}

# 2) If either predictor is rank deficient then construct a model without it
predictors <- NA
if(var(experience$LearnMicrobeSig) >0 & var(experience$LearnFlowerSig) ==0){predictors = "LearnMicrobeSig"}
if(var(experience$LearnMicrobeSig) ==0 & var(experience$LearnFlowerSig) > 0){predictors = "LearnFlowerSig"}
if(var(experience$LearnMicrobeSig) >0 & var(experience$LearnFlowerSig)  >0){predictors = "LearnMicrobeSig+LearnFlowerSig"}

decideMod <- betareg::betareg(data=experience,paste("LearnValue~",predictors,sep=""), weights = c(seq(1,3,length.out=nrow(experience)-1),6))

prediction <- predict(decideMod,
                      newdata=data.frame(
                        LearnFlowerSig=LearnFlowerSig,
                        LearnMicrobeSig=LearnMicrobeSig))
learners %>% dplyr::group_by(ID) %>% summarise(mean(tail(thresholdF1,3)))

ggplot(learners,aes(x=time,y=LearnFlowerSig))+geom_line(aes(group=ID))
ggplot(field,aes(x=time,y=flowerSignal))+geom_line(aes(group=ID))


x=seq(0,1,length.out=10)
y=seq(0,1,length.out=10)

parameterspace <- data.frame(expand.grid(mic1=x,mic2=y))
parameterspace <- parameterspace %>% rowwise %>% mutate(avg=mean(c(mic1,mic2)))

ggplot(data=parameterspace,aes(x=mic1,y=mic2,fill=avg)) +geom_tile()+geom_text(aes(label=round(avg,2)))


x<-seq(0,10000)
y<- 1/(2+exp(.002*(5000-x)))
z<- 1/(2+exp(-.002*(5000-x)))
a <- exp((.002+.002*x))/(1+exp((.002+.002*x)))
b <- exp(-.0002*(x))/(1+exp(-.0002*(x)))

ggplot(data=NULL,aes(x=x))+
  #geom_line(aes(y=z))+
  geom_line(aes(y=a),col="dark green")+
  geom_line(aes(y=b),color="dark orchid")


