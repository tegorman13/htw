---
title: Group Level Fits
date: last-modified
categories: [Modeling, ALM, R]
page-layout: full
code-fold: true
code-tools: true
execute: 
  warning: false
  eval: false
---

## test fold

### Prep data

```{r}
#| eval: true
#| code-fold: show
pacman::p_load(tidyverse,data.table,here)
options(dplyr.summarise.inform=FALSE)

# select first row for each id in d, then create histogram for nTrain
#  d  %>% group_by(id) %>% slice(1) %>% ggplot(aes(nTrain)) + geom_histogram() + facet_wrap(~condit)

d <- readRDS(here("data/dPrune-01-19-23.rds"))
dtest <- d %>% filter(expMode %in% c("test-Nf","test-train-nf")) %>% group_by(id,lowBound) %>% 
  mutate(nBand=n(),band=bandInt,id=factor(id)) %>% group_by(id) %>% mutate(nd=n_distinct(lowBound))
# unique(dtest[dtest$nd==4,]$sbjCode) # 7 in wrong condition
dtest <- dtest %>% group_by(id,lowBound) %>% filter(nBand>=5 & nd==6)
# for any id that has at least 1 nBand >=5, remove all rows with that id. 
dtest <- dtest %>% group_by(id) %>% filter(!id %in% unique(dtest$id[dtest$nBand<5]))

dtestAgg <- dtest %>% group_by(id,condit,catOrder,feedbackType,vb,band,lowBound,highBound,input) %>% mutate(vxCapped=ifelse(vx>1600,1600,vx)) %>%
  summarise(vxMean=mean(vx),devMean=mean(dist),vxMed=median(vx),devMed=median(dist),
            vxMeanCap=mean(vxCapped),.groups = "keep")

ds <- d %>% filter(expMode %in% c("train","train-Nf","test-Nf","test-train-nf")) %>% 
filter(!id %in% unique(dtest$id[dtest$nBand<5])) %>% 
select(id,condit,catOrder,feedbackType,expMode,trial,gt.train,vb,band,bandInt,lowBound,highBound,input,vx,dist,vxb) 
head(ds,6)

dst <- ds %>% filter(expMode=="train",catOrder=="orig")
head(dst)
colnames(dst)






vTrainTrial <- dst %>% filter(condit=="Varied",gt.train<=84) %>% group_by(gt.train,vb) %>% summarise(sdVx=sd(vx),vx=mean(vx),sdDist=sd(dist),dist=mean(dist)) %>% 
  group_by(vb) %>% mutate(gt.trainBin=cut(gt.train,breaks=5,labels=c(1:5)))

binTrainTrial <- dst %>% filter(gt.train<=83) %>% group_by(gt.train,vb,condit) %>% summarise(sdVx=sd(vx),vx=mean(vx),sdDist=sd(dist),dist=mean(dist)) %>% 
  group_by(vb) %>% mutate(gt.trainBin=cut(gt.train,breaks=6,labels=c(1:6)))


tMax=84
bandVec <- rep(c(800,1000,1200),each=tMax/3)
bandVec <- bandVec[sample(1:length(bandVec),tMax,replace=FALSE)]

trainTrials <- dst %>% filter(gt.train<=tMax) %>% group_by(condit,gt.train,vb,bandInt,input) %>% summarise(vx=mean(vx)) 

tv <- trainTrials %>% filter(condit=="Varied") %>% group_by(gt.train) %>% mutate(bandInt2=bandVec[gt.train]) %>% filter(bandInt==bandInt2) %>% select(-bandInt2) %>% 
  rbind(.,trainTrials %>% filter(condit=="Constant"))

head(dst)
dst[1:2,]

```

### Empirical Learning Patterns - Group Level

```{r fig.width=11, fig.height=12}
#| eval: true
#| warning: false
# grouping by condit, display count of each trial in dst
# dst %>% group_by(condit,trial) %>% summarise(n=n())
# dst %>% group_by(condit,gt.train) %>% summarise(n=n())

# display histogram with count on x axis
#dst %>% ggplot(aes(gt.train)) + geom_histogram() + facet_wrap(~condit)

ggplot(dst %>% filter(gt.train<=84), aes(x = gt.train, y = vx,color=vb)) +
  geom_point(aes(color = vb), stat = "summary", fun = mean) + 
  stat_summary(aes(color = vb), geom = "line", fun = mean) +
  stat_summary(geom="errorbar",fun.data=mean_se,width=.4,alpha=.7)+facet_wrap(~condit)

ggplot(binTrainTrial, aes(x = gt.trainBin, y = vx,color=vb)) +
  geom_point(aes(color = vb), stat = "summary", fun = mean) + 
  stat_summary(aes(color = vb), geom = "line", fun = mean) +
  stat_summary(geom="errorbar",fun.data=mean_se,width=.4,alpha=.7)+facet_wrap(~condit)


# ggplot(vTrainTrial, aes(x = gt.train, y = vx,color=vb)) +
#   geom_point(aes(color = vb), stat = "summary", fun = mean) + 
#   stat_summary(aes(color = vb), geom = "line", fun = mean) +
#   stat_summary(geom="errorbar",fun.data=mean_se,width=.4,alpha=.7)
# 
# # plot vx and sdVx over gt.train, separate facets for vx and sdVx
# vTrainTrial %>% pivot_longer(cols=c(vx,sdVx),names_to="vxType",values_to="vxVal") %>% 
#   ggplot(aes(gt.train,vxVal,color=vb)) + geom_line() + 
#   stat_summary(aes(color = vb), geom = "line", fun = mean) +
#   stat_summary(geom="errorbar",fun.data=mean_se,width=.4,alpha=.7)+
#   facet_wrap(~vxType, scale="free_y")


  
```

```{r}
input.activation<-function(x.target, c){
  return(exp(-1*c*(x.target-inputNodes)^2))
}

output.activation<-function(x.target, weights, c){
  return(weights%*%input.activation(x.target, c))
}

mean.prediction<-function(x.target, weights, association.parameter){
  probability<-output.activation(x.target, weights, c)/sum(output.activation(x.target, weights, c))
  return(outputNodes%*%probability) # integer prediction
}
# function to generate exam predictions
exam.prediction<-function(x.target, weights, c,trainVec){
  #trainVec = sort(unique(x.learning))
  nearestTrain = trainVec[which.min(abs(trainVec-x.target))]
  aresp = mean.prediction(nearestTrain, weights, c)
  xUnder = ifelse(min(trainVec) == nearestTrain, nearestTrain, trainVec[which(trainVec == nearestTrain) - 1])
  xOver = ifelse(max(trainVec) == nearestTrain, nearestTrain, trainVec[which(trainVec == nearestTrain) + 1])
  mUnder = mean.prediction(xUnder, weights, c)
  mOver = mean.prediction(xOver, weights, c)
  exam.output = round(aresp + ((mOver - mUnder) / (xOver - xUnder)) * (x.target - nearestTrain), 3)
  exam.output
}
  
update.weights<-function(x.new, y.new, weights, c, lr){
  y.feedback.activation<-exp(-1*c*(y.new-outputNodes)^2)
  x.feedback.activation<-output.activation(x.new, weights, c)
  return(weights+lr*(y.feedback.activation-x.feedback.activation)%*%t(input.activation(x.new, c)))
}

train.alm<-function(dat, c=0.05, lr=0.5, weights){
   alm.train<-rep(NA,nrow(dat))  
  for (i in 1:nrow(dat)){
    weights <- update.weights(dat$input[i], dat$vx[i], weights, c, lr)
    resp = mean.prediction(dat$input[i], weights, c)
    alm.train[i]=resp
    weights[weights<0]=0
  }
  alm.train
}

trainTest.alm<-function(dat, c=0.05, lr=0.5, weights,testVec){
   alm.train<-rep(NA,nrow(dat))  
  for (i in 1:nrow(dat)){
    weights <- update.weights(dat$input[i], dat$vx[i], weights, c, lr)
    resp = mean.prediction(dat$input[i], weights, c)
    alm.train[i]=resp
    weights[weights<0]=0
  }
  almPred <- sapply(testVec,mean.prediction,weights,c)
  examPred <- sapply(testVec,exam.prediction,weights,c,trainVec=c(1,sort(unique(dat$input))))
  list(almPred=almPred,examPred=examPred)
}

wrap_alm <- function(parms,dat, weights){
    c=parms[1]; lr=parms[2]
   pred=train.alm(dat, c=c, lr=lr, weights=weights)
   sqrt(mean((dat$vx -pred)^2))
}

wrap_optim <- function(dat,wm){
  bounds_lower <- c(.0000001, .00001)
  bounds_upper <- c(5, 5)

 optim(c(.1, .2),
   fn = wrap_alm,
   dat = dat, weights = wm,
   method = "L-BFGS-B",
   lower = bounds_lower,
   upper = bounds_upper,
   control = list(maxit = 1e4, pgtol = 0, factr = 0)
 )
}

```

```{r}

inputNodes = seq(1,7,1)  # 
outputNodes = seq(50,1600,50)
wm=matrix(.00001,nrow=length(outputNodes),ncol=length(inputNodes))


fitVaried <- tv %>% filter(condit=="Varied") %>% wrap_optim(.,wm)
fitConstant <- tv %>% filter(condit=="Constant") %>% wrap_optim(.,wm)

# call train.alm with the optimized parameters
vPred <- tv %>% filter(condit=="Varied") %>% cbind(., pred=train.alm(., c=fitVaried$par[1], lr=fitVaried$par[2], weights=wm))

cPred <- tv %>% filter(condit=="Constant") %>% cbind(., pred=train.alm(., c=fitConstant$par[1], lr=fitConstant$par[2], weights=wm))

# plot the results, showing gt.train on x axis, and separate lines with vx and pred
vPred %>% ggplot(aes(gt.train,vx)) + geom_line() + 
geom_line(aes(y=pred),color="red") + facet_wrap(~bandInt, scale="free_y")+ggtitle("Varied training and predictions")

cPred %>% ggplot(aes(gt.train,vx)) + geom_line() +
geom_line(aes(y=pred),color="red") + facet_wrap(~bandInt, scale="free_y")+ggtitle("Constant training and predictions")

# pred <- train.alm(dat, c = .05, lr = .2, wm)
# sqrt(mean((dat$vx - pred)^2))

# pred <- train.alm(dat, c = .113, lr = .048, wm)
# sqrt(mean((dat$vx - pred)^2))

# c=.05; lr=.5; y.new=dat$vx[1]; x.new=dat$input[1]; weights=wm; i=1
#dat=tv %>% filter(condit=="Varied")

```

```{r}

testVec=seq(2,7)
varTestPred <- tv %>% filter(condit=="Varied") %>%
  trainTest.alm(.,.113,.0488,wm,testVec)
  #trainTest.alm(.,fitVaried$par[1],fitVaried$par[2],wm,testVec)
varTestPred


testVec=seq(2,7)
conTestPred <- tv %>% filter(condit=="Constant") %>%
  trainTest.alm(.,.0008,.144,wm,testVec)
  #trainTest.alm(.,fitVaried$par[1],fitVaried$par[2],wm,testVec)
conTestPred





```

```{r}



```
