
pacman::p_load(tidyverse,data.table)
options(dplyr.summarise.inform=FALSE)

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

# Add a memory parameter (0 <= memory_coeff <= 1) to the function arguments
update.weightsMemory<-function(x.new, y.new, weights, c, lr, memory_coeff=.1){
  y.feedback.activation<-exp(-1*c*(y.new-outputNodes)^2)
  x.feedback.activation<-output.activation(x.new, weights, c)
  
  # Modify the weights update equation to include the memory component
  delta_weights <- lr * (y.feedback.activation - x.feedback.activation) %*% t(input.activation(x.new, c))
  new_weights <- (1 - memory_coeff) * delta_weights + memory_coeff * weights
  return(new_weights)
  
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


trainTest.alm<-function(dat, c=0.05, lr=0.5, weights,testVec, update_func, noise_sd){
  update_func=get(update_func)
  alm.train<-rep(NA,nrow(dat))  
  for (i in 1:nrow(dat)){
    weights <- update.weights(dat$input[i], dat$vx[i], weights, c, lr)
    resp = mean.prediction(dat$input[i], weights, c)
    alm.train[i]=resp
    weights[weights<0]=0
  }
  almPred <- sapply(testVec,mean.prediction,weights,c)
  examPred <- sapply(testVec,exam.prediction,weights,c,trainVec=c(1,sort(unique(dat$input))))
  list(almTrain=alm.train,almPred=almPred,examPred=examPred)
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

# sim_data <- function(c=0.5, lr=0.2,noise=0,inNodes=7,outNodes=32,
#                      trainVec=c(5,6,7),trainRep=10,testVec=seq(2,7)) {
#   
#   inputNodes <<- seq(1,7,length.out=inNodes*1)  
#   outputNodes <<- seq(50,1600,length.out=outNodes*1) 
#   wm=matrix(.0000,nrow=length(outputNodes),ncol=length(inputNodes))
#   dat<-gen_train(trainVec,trainRep,noise)
#   #trainDat <- train.alm(dat,c,lr,wm)
#   tt<-trainTest.alm(dat,c,lr,wm,testVec)
# }

# Modify the sim_data function to accept the dataset as an argument
sim_data <- function(dat, c=0.5, lr=0.2, inNodes=7, outNodes=32, trainVec=c(5,6,7)) {
  inputNodes <<- seq(1,7,length.out=inNodes*1)  
  outputNodes <<- seq(50,1600,length.out=outNodes*1) 
  wm=matrix(.0000,nrow=length(outputNodes),ncol=length(inputNodes))
  tt<-trainTest.alm(dat, c, lr, wm, trainVec)
}

gen_train <- function(trainVec=c(5,6,7),trainRep=10,noise=0){
  bandVec=c(0,100,350,600,800,1000,1200)
  if(class(trainVec)=="list"){trainVec=unlist(trainVec)}
  ts <- rep(seq(1,length(trainVec)),trainRep)
  # print(trainVec)
  # print(length(ts)); print(length(trainRep))
  noiseVec=rnorm(length(ts),mean=0)*noise
  if(noise==0) {noiseVec=noiseVec*0}
  tibble(trial=seq(1,length(ts)),input=trainVec[ts],vx=bandVec[trainVec[ts]]+noiseVec)
}



wrap_alm <- function(par,dat, weights,lossFun){
  c=par[1]; lr=par[2]
  pred=train.alm(dat, c=c, lr=lr, weights=weights)
  #sqrt(mean((dat$vx -pred)^2))
  lossFun(dat$vx,pred)
}

wrap_optim <- function(dat,lossFun=RMSE){
  if(class(lossFun)=="character"){lossFun=get(lossFun)}
  inputNodes = seq(1,7,1)  # 
  outputNodes = seq(50,1600,50)
  wm=matrix(.00001,nrow=length(outputNodes),ncol=length(inputNodes))
  testVec=seq(2,7)
  
  bounds_lower <- c(.0000001, .00001)
  bounds_upper <- c(10, 10)
  parmsLab <- c("c","lr")
  
  fit=optim(par=c(.1, .2),
            fn = wrap_alm,
            dat = dat, weights = wm,lossFun=lossFun,
            method = "L-BFGS-B",
            lower = bounds_lower,
            upper = bounds_upper,
            control = list(maxit = 1e5, pgtol = 0, factr = 0)
  )
  
  l=reduce(list(list(fit),fit$par,fit$value),append)
  names(l)=c("Fit",parmsLab,"Value")
  return(l)
}

# implement wrap_optim using grid search
wrap_grid <- function(dat,lossFun=RMSE){
  if(class(lossFun)=="character"){lossFun=get(lossFun)}  
  inputNodes = seq(1,7,1)  # 
  outputNodes = seq(50,1600,50)
  wm=matrix(.0000,nrow=length(outputNodes),ncol=length(inputNodes))
  testVec=seq(2,7)
  # define grid boundaries
  cRange <- seq(.000001, 5, length.out = 30)
  lrRange <- seq(.05, 5, length.out = 20)
  # create grid
  grid <- expand.grid(c = cRange, lr = lrRange)
  grid$Value <- NA
  # loop through grid
  for (i in 1:nrow(grid)) {
    grid$Value[i] <- wrap_alm(par=c(grid[i,c("c") ],grid[i,c("lr") ]), dat=dat, weights=wm,lossFun=lossFun)
  }
  
  # find best fit
  bestFit <- grid[which.min(grid$Value), ]
  bestFit$Value <- min(grid$Value)
  
  # return best fit, and best c and lr
  return(list(Fit = bestFit, c = bestFit$c, lr = bestFit$lr, Value = bestFit$Value))
  
}


RMSE <- function(x,y){
  # print("rmseTrial")
  sqrt(mean((x-y)^2))
}

RMSE.blocked <- function(x,y,blocks=6){
  #print("rmseBlocked")
  data.table(x=x,y=y,t=seq(1,length(x))) %>% 
    .[, `:=`(fitBins = cut(t, breaks = ..blocks, labels = c(1:..blocks)))] %>%
    .[, .(predMean = mean(x), obsMean = mean(y)), keyby = .(fitBins)] %>%
    .[, RMSE(predMean,obsMean)] %>% as.numeric()
}

MAE <- function(x, y) {
  mean(abs(x - y))
}

MAPE <- function(x, y) {
  mean(abs((x - y) / y)) * 100
}

MedAE <- function(x, y) {
  median(abs(x - y))
}

HuberLoss <- function(x, y, delta = 1) {
  error <- x - y
  abs_error <- abs(error)
  loss <- ifelse(abs_error <= delta, 0.5 * error^2, delta * (abs_error - 0.5 * delta))
  mean(loss)
}


sigmoid <- function(x) {1/(1+exp(-x))}
