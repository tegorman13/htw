


library(tidyverse)
input.activation<-function(x.target, c){return(exp((-1*c)*(x.target-inputNodes)^2))}
output.activation<-function(x.target, weights, c){return(weights%*%input.activation(x.target, c))}
mean.prediction<-function(x.target, weights, c){
  probability<-output.activation(x.target, weights, c)/sum(output.activation(x.target, weights, c))

  return(t(outputNodes)%*%probability) # integer prediction
}
update.weights<-function(x.new, y.new, weights, c, lr){
  y.feedback.activation<-exp(-1*c*(y.new-outputNodes)^2)
  x.feedback.activation<-output.activation(x.new, weights, c)
  return(weights+lr*(y.feedback.activation-x.feedback.activation)%*%t(input.activation(x.new, c)))
}
train.alm<-function(dat, c=0.05, lr=0.5, weights){
  alm.train<-rep(NA,nrow(dat))
  for (i in 1:nrow(dat)){
    weights <- update.weights(dat$x[i], dat$y[i], weights, c, lr)
    resp = mean.prediction(dat$x[i], weights, c)
    alm.train[i]=resp
    weights[weights<0]=0
  }
  list(almTrain = alm.train, weights = weights)
}
alm.response <- function(input=1,c) {
  input.activation <- exp(-c*(input.layer - input)^2)
  input.activation <<- input.activation/sum(input.activation)
  #print(length(input.activation)); print(dim(weight.mat))
  output.activation <<- weight.mat %*% input.activation
  output.probability <<- output.activation/sum(output.activation)
  mean.response <<- sum(output.layer * output.probability)
  mean.response
}
alm.trial <- function(input, corResp,c,lr){
  alm.response(input,c)
  alm.update(corResp,c,lr)
  # print(paste0("input=",input,"; corResp=",corResp,"; mean.response=",mean.response))
  mean.response
}


generate.data <- function(x, type = 'linear', noise = NA) {
  if (type == 'linear') {
    y <- round(2.2*x + 30,0)
  }
  else if (type == 'exponential') {
    y <- round(200*(1-exp(-x/25)),0)
  }
  else if (type == 'quadratic') {
    y <- round(210 - ((x-50)^2)/12,0)
  }
  # if noise is specified, add noise to the y values
  if(!is.na(noise)) {
    y <- y + round(rnorm(length(y), 0, noise),2)
  }
  data.frame(x, y,type)
}
envTypes <- c('linear', 'exponential', 'quadratic')
lowDensityTrainBlock <- c(30.5, 36.0, 41.0, 46.5, 53.5, 59.0, 64.0, 69.5)
medDensityTrainBlock <- c(
  30.0, 31.5, 33.0, 34.5, 36.5, 38.5, 41.0, 43.5, 46.0,
  48.5, 51.5, 54.0, 56.5, 59.0, 61.5, 63.5, 65.5, 67.0, 68.5, 70.0
)
#  low density has 25 training blocks, medium has 10 blocks
# generate training data, for each combination of environment type and density.

lowTrain <- map_dfr(envTypes, ~ generate.data(rep(lowDensityTrainBlock,25), type = .x)) %>% group_by(type) %>% mutate(block = rep(1:25, each = 8),trial=seq(1,200))


medTrain <- map_dfr(envTypes, ~ generate.data(rep(medDensityTrainBlock,10), type = .x)) %>% 
  group_by(type) %>% mutate(block = rep(1:10, each = 20),trial=seq(1,200))


# Initialize weights
inputNodes <<- matrix(seq(0,100,.5) ) # half step units for inputs, from 0 to 100
outputNodes <<- matrix(seq(0,250,1)) # single step units for outputs, from 0 to 250
weight.mat <- matrix(0.0000000,nrow=length(outputNodes),ncol=length(inputNodes )) 
# Train the model on the low density training data
lowTrainResults <- train.alm(lowTrain, c = 0.05, lr = 0.5, weights = weight.mat)
# Train the model on the medium density training data
medTrainResults <- train.alm(medTrain, c = 0.05, lr = 0.5, weights = weight.mat)
# Plot the results
plot(lowTrainResults$almTrain, type = 'l', main = 'Learning Performance on Low Density Training Data')
plot(medTrainResults$almTrain, type = 'l', main = 'Learning Performance on Medium Density Training Data')




alm.sim <- function(dat, c, lr,testRange=seq(0,100,.5)){
  input.layer <<- matrix(seq(0,100,.5) ) # half step units for inputs, from 0 to 100
  output.layer <<- matrix(seq(0,250,1)) # single step units for outputs, from 0 to 250
  weight.mat <<- matrix(0.0000000,nrow=length(output.layer),ncol=length(input.layer )) # weights initialized to 0 (as in Delosh 1997)
  xt<<-dat$x
  st <- map2_dbl(dat$x, dat$y, ~alm.trial(.x,.y,c,lr))
  dat <- dat %>% mutate(almResp = st)
  
  return(list(d=dat,wm=weight.mat,c=c,lr=lr)) # final weightmat is probs incorrect for all but last
}





input.activation<-function(x.target, c){
  return(exp((-1*c)*(x.target-inputNodes)^2))
}
output.activation<-function(x.target, weights, c){
  return(weights%*%input.activation(x.target, c))
}
mean.prediction<-function(x.target, weights, c){
  probability<-output.activation(x.target, weights, c)/sum(output.activation(x.target, weights, c))
  return(outputNodes%*%probability) # integer prediction
}
update.weights<-function(x.new, y.new, weights, c, lr, noise_sd = NULL){
  y.feedback.activation<-exp(-1*c*(y.new-outputNodes)^2)
  x.feedback.activation<-output.activation(x.new, weights, c)
  return(weights+lr*(y.feedback.activation-x.feedback.activation)%*%t(input.activation(x.new, c)))
}

train.alm<-function(dat, c=0.05, lr=0.5, weights){
  alm.train<-rep(NA,nrow(dat))  
  for (i in 1:nrow(dat)){
    weights <- update.weights(dat$x[i], dat$y[i], weights, c, lr)
    resp = mean.prediction(dat$x[i], weights, c)
    alm.train[i]=resp
    weights[weights<0]=0
  }
  alm.train
}


# Modify the sim_data function to accept the dataset as an argument
sim_data <- function(dat, c=0.5, lr=0.2, inNodes=7, outNodes=32, trainVec=c(5,6,7)) {
  inputNodes <<- seq(1,7,length.out=inNodes*1)  
  outputNodes <<- seq(50,1600,length.out=outNodes*1) 
  wm=matrix(.0000,nrow=length(outputNodes),ncol=length(inputNodes))
  tt<-train.alm(dat, c, lr, wm)
}

ld <- generate.data(rep(lowDensityTrainBlock,25))

lowTrain <- map_dfr(envTypes, ~ generate.data(rep(lowDensityTrainBlock,25), type = .x)) %>% 
  group_by(type) %>% mutate(block = rep(1:25, each = 8),trial=seq(1,200))


map(lowTrain,~sim_data(dat=.x))


tibble(crossing(
  c = c(.5,5),lr = c(.05,1),noise = c(0),
  inNodes=c(7),outNodes=c(32),
  trainVec=list(list(5,6,7)),trainRep = c(9),
  lossFun = list("MAE"),
  simNum = 1:1,
  update_func = list("update.weights"),update_func_name = c("uW"),
  noise_sd = c(0)
)) %>%   mutate(id=seq(1,nrow(.)),td = pmap(list(trainVec,trainRep,noise),~gen_train(trainVec=.x,trainRep=..2,noise=..3) )) %>% 
  ungroup() %>%
  mutate(d = pmap(list(td, c, lr, update_func,noise_sd,inNodes,outNodes), 
                  ~sim_data(dat = .x, c = ..2, lr = ..3, update_func = ..4, noise_sd = ..5,inNodes=..6,outNodes=..7)),
         almTrainDat = map(d, "almTrain"),weights=map(d,"weights"))%>%
  unnest(c(almTrainDat, td)) %>% select(-d) %>% mutate(input=as.factor(input)) %T>%
  {pf(.) } %>% trainTab
