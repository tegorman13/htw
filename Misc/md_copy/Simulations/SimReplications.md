---
title: General Simulations
date: last-modified
categories: [Simulation, ALM, EXAM, R]
code-fold: true
code-tools: true
execute: 
  warning: false
--- 

### Functions
```{r}
pacman::p_load(tidyverse,reshape2)

input.activation<-function(x.target, association.parameter){
  return(exp(-1*association.parameter*(x.target-x.plotting)^2))
}

output.activation<-function(x.target, weights, association.parameter){
  return(weights%*%input.activation(x.target, association.parameter))
}

mean.prediction<-function(x.target, weights, association.parameter){
  probability<-output.activation(x.target, weights, association.parameter)/sum(output.activation(x.target, weights, association.parameter))
  return(y.plotting%*%probability) # integer prediction
}
# function to generate exam predictions
exam.prediction<-function(x.target, weights, association.parameter){
  trainVec = sort(unique(x.learning))
  nearestTrain = trainVec[which.min(abs(trainVec-x.target))]
  aresp = mean.prediction(nearestTrain, weights, association.parameter)
  xUnder = ifelse(min(trainVec) == nearestTrain, nearestTrain, trainVec[which(trainVec == nearestTrain) - 1])
  xOver = ifelse(max(trainVec) == nearestTrain, nearestTrain, trainVec[which(trainVec == nearestTrain) + 1])
  mUnder = mean.prediction(xUnder, weights, association.parameter)
  mOver = mean.prediction(xOver, weights, association.parameter)
  exam.output = round(aresp + ((mOver - mUnder) / (xOver - xUnder)) * (x.target - nearestTrain), 3)
  exam.output
}
  
update.weights<-function(x.new, y.new, weights, association.parameter, update.parameter){
  y.feedback.activation<-exp(-1*association.parameter*(y.new-y.plotting)^2)
  x.feedback.activation<-output.activation(x.new, weights, association.parameter)
  return(weights+update.parameter*(y.feedback.activation-x.feedback.activation)%*%t(input.activation(x.new, association.parameter)))
}

learn.alm<-function(y.learning, association.parameter=0.05, update.parameter=0.5){
  weights<-matrix(rep(0.00, length(y.plotting)*length(x.plotting)), nrow=length(y.plotting), ncol=length(x.plotting))
  for (i in 1:length(y.learning)){
    weights<-update.weights(x.learning[i], y.learning[i], weights, association.parameter, update.parameter)
    resp=mean.prediction(x.learning[i],weights,association.parameter)
    weights[weights<0]=0
  }
  alm.predictions<-sapply(x.plotting, mean.prediction, weights=weights, association.parameter=association.parameter)
  exam.predictions <- sapply(x.plotting, exam.prediction, weights=weights, association.parameter=association.parameter)
  return(list(alm.predictions=alm.predictions, exam.predictions=exam.predictions))
}

```



### No noise, 1 training rep

Red dots are training points - gray lines are individual simulations, black line is average of simulations

```{r fig.height=9,fig.width=10}
# | eval: false

trainRep=1

x.plotting<<-seq(0,100, .5)
y.plotting<<-seq(0, 210, by=2)
f.plotting<-as.numeric(x.plotting*2.2+30)
x.learning<-rep(x.plotting[20*c(4:7)+1])
f.learning<-rep(f.plotting[20*c(4:7)+1])

parmVec <- expand.grid(assoc=c(.1,0.5),update=c(0.2,1),noise=c(0),trainRep=c(1))
#parmVec <- expand.grid(assoc=c(.01),update=c(0.5),noise=c(30),trainRep=c(1,2,3,4))

parmVec$sim <- 1:nrow(parmVec)
nSim=nrow(parmVec)

nRep=5
output <- list()
for (i in 1:nrow(parmVec)){
  x.learning<-rep(x.plotting[20*c(4:7)+1],times=parmVec$trainRep[i])
  f.learning<-rep(f.plotting[20*c(4:7)+1],times=parmVec$trainRep[i])
  #noise.learning <- rnorm(n_distinct(f.learning),sd=parmVec$noise[i])
  output[[i]] <- replicate(nRep, list(learn.alm(f.learning+rep(rnorm(n_distinct(f.learning),sd=parmVec$noise[i]),times=parmVec$trainRep[i]), 
                                         association.parameter=parmVec$assoc[i], update.parameter=parmVec$update[i])))
}


# convert list of dataframes to a list of lists, each list is a simulation, each element is a dataframe
output1 <- lapply(output, function(x) lapply(x, as.data.frame)) # 10 dfs x 9 lists
output2 <- lapply(output1, function(x) Reduce(rbind,x))# 1 df x 9 lists
output3 <- lapply(output2, function(x) mutate(x, x=rep(x.plotting,nRep),y=rep(f.plotting,nRep),
                                              repN=rep(seq(1,nRep),each=length(x.plotting))))
o4 <- Reduce(rbind,output3) %>% 
  mutate(sim=rep(seq(1,nrow(parmVec)),each=nRep*length(x.plotting))) %>%
  left_join(.,parmVec,by="sim") %>%
  mutate(pvec=paste0("c=",assoc,"_lr=",update,"_noise=",noise,"_nrep=",trainRep),pv=factor(pvec),rn=factor(repN)) 

oMeans <- o4 %>% group_by(pv,x,y) %>% 
  summarise(alm.predictions=mean(alm.predictions),exam.predictions=mean(exam.predictions),.groups="keep")

o4 %>% ggplot(aes(x=x,y=alm.predictions,color=rn))+geom_line(alpha=.7)+
   scale_color_manual(values=rep("grey",nRep))+
  theme(legend.position="none")+
  geom_point(data=data.frame(x=x.learning,y=f.learning),aes(x=x,y=y),color="red")+
  geom_line(data=o4,aes(x=x,y=y),color="black",alpha=.5, linetype=2)+
  geom_line(data=oMeans,aes(x=x,y=alm.predictions),color="black")+
  facet_wrap(~pv, scales="free_y")+ggtitle("ALM predictions")

o4 %>% ggplot(aes(x=x,y=exam.predictions,color=rn))+ geom_line()+ #geom_line(color="grey",alpha=.4)+
  scale_color_manual(values=rep("grey",nRep))+
  theme(legend.position="none")+
  geom_point(data=data.frame(x=x.learning,y=f.learning),aes(x=x,y=y),color="red")+
  geom_line(data=o4,aes(x=x,y=y),color="black",alpha=.5,linetype=2)+
  geom_line(data=oMeans,aes(x=x,y=exam.predictions),color="black")+
  facet_wrap(~pv, scales="free_y")+ggtitle("EXAM predictions")
```





### High noise, 1 training rep

Red dots are training points - gray lines are individual simulations, black line is average of simulations
```{r fig.height=9,fig.width=10}
#| eval: false

trainRep=1

x.plotting<<-seq(0,100, .5)
y.plotting<<-seq(0, 210, by=2)
f.plotting<-as.numeric(x.plotting*2.2+30)
x.learning<-rep(x.plotting[20*c(4:7)+1])
f.learning<-rep(f.plotting[20*c(4:7)+1])


parmVec <- expand.grid(assoc=c(.1,0.5),update=c(0.2,1),noise=c(30),trainRep=c(1))
#parmVec <- expand.grid(assoc=c(.01),update=c(0.5),noise=c(30),trainRep=c(1,2,3,4))

parmVec$sim <- 1:nrow(parmVec)
nSim=nrow(parmVec)

nRep=10
output <- list()
for (i in 1:nrow(parmVec)){
  x.learning<-rep(x.plotting[20*c(4:7)+1],times=parmVec$trainRep[i])
  f.learning<-rep(f.plotting[20*c(4:7)+1],times=parmVec$trainRep[i])
  #noise.learning <- rnorm(n_distinct(f.learning),sd=parmVec$noise[i])
  output[[i]] <- replicate(nRep, list(learn.alm(f.learning+rep(rnorm(n_distinct(f.learning),sd=parmVec$noise[i]),times=parmVec$trainRep[i]), 
                                         association.parameter=parmVec$assoc[i], update.parameter=parmVec$update[i])))
}

# 
# nRep=3
# output <- list()
# for (i in 1:nrow(parmVec)){
#   output[[i]] <- replicate(nRep, list(learn.alm(f.learning + rnorm(length(f.learning), sd=10), 
#                                          association.parameter=parmVec$assoc[i], update.parameter=parmVec$update[i])))
# }



#output[[i]] <- replicate(nRep, list(learn.alm(f.learning + rnorm(length(f.learning), sd=10)

# convert list of dataframes to a list of lists, each list is a simulation, each element is a dataframe
output1 <- lapply(output, function(x) lapply(x, as.data.frame)) # 10 dfs x 9 lists
output2 <- lapply(output1, function(x) Reduce(rbind,x))# 1 df x 9 lists
output3 <- lapply(output2, function(x) mutate(x, x=rep(x.plotting,nRep),y=rep(f.plotting,nRep),
                                              repN=rep(seq(1,nRep),each=length(x.plotting))))
o4 <- Reduce(rbind,output3) %>% 
  mutate(sim=rep(seq(1,nrow(parmVec)),each=nRep*length(x.plotting))) %>%
  left_join(.,parmVec,by="sim") %>%
  mutate(pvec=paste0("c=",assoc,"_lr=",update,"_noise=",noise,"_nrep=",trainRep),pv=factor(pvec),rn=factor(repN)) 

oMeans <- o4 %>% group_by(pv,x,y) %>% 
  summarise(alm.predictions=mean(alm.predictions),exam.predictions=mean(exam.predictions),.groups="keep")

o4 %>% ggplot(aes(x=x,y=alm.predictions,color=rn))+geom_line(alpha=.7)+
   scale_color_manual(values=rep("grey",nRep))+
  theme(legend.position="none")+
  geom_point(data=data.frame(x=x.learning,y=f.learning),aes(x=x,y=y),color="red")+
  geom_line(data=o4,aes(x=x,y=y),color="black",alpha=.5, linetype=2)+
  geom_line(data=oMeans,aes(x=x,y=alm.predictions),color="black")+
  facet_wrap(~pv, scales="free_y")+ggtitle("ALM predictions")+ylim(0,300)

o4 %>% ggplot(aes(x=x,y=exam.predictions,color=rn))+ geom_line()+ #geom_line(color="grey",alpha=.4)+
  scale_color_manual(values=rep("grey",nRep))+
  theme(legend.position="none")+
  geom_point(data=data.frame(x=x.learning,y=f.learning),aes(x=x,y=y),color="red")+
  geom_line(data=o4,aes(x=x,y=y),color="black",alpha=.5,linetype=2)+
  geom_line(data=oMeans,aes(x=x,y=exam.predictions),color="black")+
  facet_wrap(~pv, scales="free_y")+ggtitle("EXAM predictions")+ylim(0,300)
```




### High noise, 60 training rep

Red dots are training points - gray lines are individual simulations, black line is average of simulations

```{r fig.height=9,fig.width=10}
#| eval: false

trainRep=1

x.plotting<<-seq(0,100, .5)
y.plotting<<-seq(0, 210, by=2)
f.plotting<-as.numeric(x.plotting*2.2+30)
x.learning<-rep(x.plotting[20*c(4:7)+1])
f.learning<-rep(f.plotting[20*c(4:7)+1])


parmVec <- expand.grid(assoc=c(.1,0.5),update=c(0.2,1),noise=c(30),trainRep=c(60))
#parmVec <- expand.grid(assoc=c(.01),update=c(0.5),noise=c(30),trainRep=c(1,2,3,4))

parmVec$sim <- 1:nrow(parmVec)
nSim=nrow(parmVec)

nRep=10
output <- list()
for (i in 1:nrow(parmVec)){
  x.learning<-rep(x.plotting[20*c(4:7)+1],times=parmVec$trainRep[i])
  f.learning<-rep(f.plotting[20*c(4:7)+1],times=parmVec$trainRep[i])
  #noise.learning <- rnorm(n_distinct(f.learning),sd=parmVec$noise[i])
  output[[i]] <- replicate(nRep, list(learn.alm(f.learning+rep(rnorm(n_distinct(f.learning),sd=parmVec$noise[i]),times=parmVec$trainRep[i]), 
                                         association.parameter=parmVec$assoc[i], update.parameter=parmVec$update[i])))
}


# convert list of dataframes to a list of lists, each list is a simulation, each element is a dataframe
output1 <- lapply(output, function(x) lapply(x, as.data.frame)) # 10 dfs x 9 lists
output2 <- lapply(output1, function(x) Reduce(rbind,x))# 1 df x 9 lists
output3 <- lapply(output2, function(x) mutate(x, x=rep(x.plotting,nRep),y=rep(f.plotting,nRep),
                                              repN=rep(seq(1,nRep),each=length(x.plotting))))
o4 <- Reduce(rbind,output3) %>% 
  mutate(sim=rep(seq(1,nrow(parmVec)),each=nRep*length(x.plotting))) %>%
  left_join(.,parmVec,by="sim") %>%
  mutate(pvec=paste0("c=",assoc,"_lr=",update,"_noise=",noise,"_nrep=",trainRep),pv=factor(pvec),rn=factor(repN)) 

oMeans <- o4 %>% group_by(pv,x,y) %>% 
  summarise(alm.predictions=mean(alm.predictions),exam.predictions=mean(exam.predictions),.groups="keep")

o4 %>% ggplot(aes(x=x,y=alm.predictions,color=rn))+geom_line(alpha=.7)+
   scale_color_manual(values=rep("grey",nRep))+
  theme(legend.position="none")+
  geom_point(data=data.frame(x=x.learning,y=f.learning),aes(x=x,y=y),color="red")+
  geom_line(data=o4,aes(x=x,y=y),color="black",alpha=.5, linetype=2)+
  geom_line(data=oMeans,aes(x=x,y=alm.predictions),color="black")+
  facet_wrap(~pv, scales="free_y")+ggtitle("ALM predictions")+ylim(0,300)

o4 %>% ggplot(aes(x=x,y=exam.predictions,color=rn))+ geom_line()+ #geom_line(color="grey",alpha=.4)+
  scale_color_manual(values=rep("grey",nRep))+
  theme(legend.position="none")+
  geom_point(data=data.frame(x=x.learning,y=f.learning),aes(x=x,y=y),color="red")+
  geom_line(data=o4,aes(x=x,y=y),color="black",alpha=.5,linetype=2)+
  geom_line(data=oMeans,aes(x=x,y=exam.predictions),color="black")+
  facet_wrap(~pv, scales="free_y")+ggtitle("EXAM predictions")+ylim(0,300)
```




### [Shiny App](ALM_Shiny.html)






```{r}
#| eval: false

x.plotting<<-seq(0,90, .5)
y.plotting<<-seq(0, 210, by=2)
f.plotting<-as.numeric(x.plotting * 2.2 + 30)
x.learning<-x.plotting[10*c(3:9)+1]
f.learning<-f.plotting[10*c(3:9)+1]

# Single Simulation
# get alm and exam predictions for full range of x.plotting
output<-learn.alm(f.learning)
alm.predictions<-output$alm.predictions
exam.predictions<-output$exam.predictions

# plot the results
plot(x.plotting, f.plotting, type="l", col="blue", lwd=.5, xlab="x", ylab="f(x)")
points(x.learning, f.learning, col="red", pch=19)
lines(x.plotting, alm.predictions, col="green", lwd=2)
lines(x.plotting, exam.predictions, col="purple", lwd=2)
legend("topright", legend=c("f(x)", "training data", "ALM", "Exam"), col=c("blue", "red", "green", "purple"), lty=1.5, cex=0.8)

#function to plot in greyscale
plot.grey<-function(predictions){
  lines(x.plotting, predictions, col="grey")
}

```


```{r}
#| eval: false

# Average of 100 simulations:
# get alm and exam predictions for full range of x.plotting, averaged over 100 simulations
nSim<-10
output <- replicate(nSim, list(learn.alm(f.learning + rnorm(length(f.learning), sd=0.1), 
                                         association.parameter=0.05, update.parameter=0.5)))

#alm.predictions<-do.call(rbind, lapply(output, function(x) x$alm.predictions))
alm.predictions <- Reduce(rbind,output %>% map("alm.predictions"))
exam.predictions <- Reduce(rbind,output %>% map("exam.predictions"))

alm.predictions.avg<-apply(alm.predictions, 2, mean)
exam.predictions.avg<-apply(exam.predictions, 2, mean)
dfAvg<-data.frame(x=x.plotting, f=f.plotting, alm=alm.predictions.avg, exam=exam.predictions.avg)
dfAvg<-reshape2::melt(dfAvg, id.vars="x")
dfAvg$model<-factor(dfAvg$variable, levels=c("f", "alm", "exam"))
ggplot(dfAvg, aes(x=x, y=value, color=model)) + geom_line() + geom_point(data=data.frame(x=x.learning, f=f.learning), aes(x=x, y=f), color="red", size=2) + theme_bw() + theme(legend.position="topright")


alm.predictions<-as.data.frame(alm.predictions) %>% mutate(sim=seq(1:nSim))
alm.predictions<-pivot_longer(alm.predictions, cols=1:ncol(alm.predictions)-1, 
                              names_to=c("sim"), values_to="alm",names_repair = "unique") 
colnames(alm.predictions)=c("sim","x","pred")
alm.predictions <- alm.predictions %>% mutate(stim = as.numeric(gsub("V", "", x)),model="alm",x=x.plotting[stim])


exam.predictions<-as.data.frame(exam.predictions) %>% mutate(sim=seq(1:nSim))
exam.predictions<-pivot_longer(exam.predictions, cols=1:ncol(exam.predictions)-1, 
                               names_to=c("sim"), values_to="exam",names_repair = "unique")
colnames(exam.predictions)=c("sim","x","pred")
exam.predictions <- exam.predictions %>% mutate(stim = as.numeric(gsub("V", "", x)),model="exam",x=x.plotting[stim])


df<- rbind(alm.predictions,exam.predictions)

ggplot(df, aes(x=x, y=pred, color=sim)) + geom_line(alpha=.2) + facet_wrap(~model) + theme_bw() + 
  geom_point(data=data.frame(x=x.learning, f=f.learning), aes(x=x, y=f), color="red", size=2)+
  geom_line(data=data.frame(x=x.plotting, f=f.plotting),aes(x=x,y=f),color="black")+
  geom_line(data=dfAvg %>% filter(model!="f"),aes(x=x,y=value),color="green")


```

