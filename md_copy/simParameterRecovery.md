---
title: Parameter Recovery Simulations
date: last-modified
categories: [Simulation, ALM, EXAM, R]
page-layout: full
code-fold: true
code-tools: true
execute: 
  warning: false
  eval: false 
--- 

# Parameter Recovery Simulations

Assessing the parameter recovery of ALM and EXAM on synthetic data. 

```{r}
#| code-fold: show
pacman::p_load(tidyverse,data.table)
options(dplyr.summarise.inform=FALSE)

input.activation<-function(x.target, c) { return(exp((-1*c)*(x.target-inputNodes)^2))}

output.activation<-function(x.target, weights, c){
  return(weights%*%input.activation(x.target, c))
}
mean.prediction<-function(x.target, weights, c){
  probability<-output.activation(x.target, weights, c)/sum(output.activation(x.target, weights, c))
  return(outputNodes%*%probability) # integer prediction
}

update.weights<-function(x.new, y.new, weights, c, lr){t
  yt.feedback.activation<-exp(-1*c*(y.new-outputNodes)^2)
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

# Modify the sim_data function to accept the dataset as an argument
sim_data <- function(dat, c=0.5, lr=0.2, inNodes=7, outNodes=32, trainVec=c(5,6,7)) {
  inputNodes <<- seq(1,7,length.out=inNodes*1)  
  outputNodes <<- seq(50,1600,length.out=outNodes*1) 
  wm=matrix(.0000,nrow=length(outputNodes),ncol=length(inputNodes))
  tt<-trainTest.alm(dat, c, lr, wm, trainVec)
}

gen_train <- function(trainVec=c(5,6,7),trainRep=10,noise=0){
   bandVec=c(0,100,350,600,800,1000,1200)
   ts <- rep(seq(1,length(trainVec)),trainRep)
   noiseVec=rnorm(length(ts),mean=0)*noise
   if(noise==0) {noiseVec=noiseVec*0}
   tibble(input=trainVec[ts],vx=bandVec[trainVec[ts]]+noiseVec)
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
  list(almTrain=alm.train,almPred=almPred,examPred=examPred)
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
```


### Loss Functions
```{r}

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


```









```{r}

k= parmVec %>% group_by(simNum,c,lr) %>% mutate(td = list(gen_train(trainRep = first(trainRep), noise = first(noise))),
                                                 o=map(td,head,1),vx1=map(o,"vx"))

parmVec <- tibble(crossing(c = c(0.1), lr = c(0.1,0.4,1), noise = c(10), trainRep = c(20), lossFun = list("RMSE", "MAE","MAPE", "MedAE", "HuberLoss"), simNum = 1:10))

sdp <- parmVec %>% group_by(simNum,c,lr) %>% mutate(td = list(gen_train(trainRep = first(trainRep), noise = first(noise)))) %>% ungroup() %>%
  mutate(d = pmap(list(td, c, lr), ~sim_data(dat = .x, c = ..2, lr = ..3)),
         almTrainDat = map(d, "almTrain"),
         almTestDat = map(d, "almPred"),
         examTestDat = map(d, "examPred"),
         fitO = map2(td, lossFun, ~wrap_optim(.x, .y)),
         cFitO = map_dbl(fitO, "c"),
         lrFitO = map_dbl(fitO, "lr"),
         optimValO = map_dbl(fitO, "Value")) 

sdpResults <- sdp %>% 
  mutate(lossFun=as.character(lossFun),cDiff=cFitO-c,lrDiff=lrFitO-lr) %>%
  relocate(simNum,lossFun, c,lr,cFitO,lrFitO,optimValO,cDiff,lrDiff) %>% 
  arrange(lossFun,lr,c)

averaged_sdp <- sdpResults %>%
  group_by(lossFun, c, lr) %>%
  summarise(
    avg_cFitO = mean(cFitO),
    avg_lrFitO = mean(lrFitO),
    avg_optimValO = mean(optimValO),
    avg_cDiff = mean(cDiff),
    avg_lrDiff = mean(lrDiff),
    .groups = "drop"
  )


averaged_sdp <- sdpResults %>%
  group_by(lossFun, c, lr) %>%
  summarise(
    avg_cFitO = mean(cFitO),
    var_cFitO = var(cFitO),
    avg_lrFitO = mean(lrFitO),
    var_lrFitO = var(lrFitO),
    avg_optimValO = mean(optimValO),
    var_optimValO = var(optimValO),
    avg_cDiff = mean(cDiff),
    var_cDiff = var(cDiff),
    avg_lrDiff = mean(lrDiff),
    var_lrDiff = var(lrDiff),
    .groups = "drop"
  )

# averaged_sdp <- sdp %>%
#   group_by(c, lr, noise, trainRep, lossFun) %>%
#   summarise(across(starts_with("cFit") | starts_with("lrFit") | starts_with("optimVal"), list(mean = mean), .names = "mean_{.col}")) %>% 
#   mutate(lossFun=as.character(lossFun),
#          diff_cFitO = abs(c - mean_cFitO),
#          diff_lrFitO = abs(lr - mean_lrFitO)) %>%
#   relocate(noise, trainRep, lossFun, c, diff_cFitO, lr, diff_lrFitO) %>%
#   dplyr::arrange(c,lr)


#k=parmVec %>% group_by(simNum) %>% mutate(td = list(gen_train(trainRep = first(trainRep), noise = first(noise))))
#sdp <- parmVec %>% mutate(d = pmap(list(c, lr, noise, trainRep), ~sim_data(c = ..1, lr = ..2, noise = ..3, trainRep = ..4)),

```










```{r}


```


```{r}
library(plotly)

g2= grid %>% filter(Value<160) %>% arrange(Value)
#plot_ly() %>% add_trace(data=g2,x=grid$c,y=grid$lr,z=grid$Value,type="mesh3d")
 
plot_ly(data=g2,x=~c,y=~lr,z=~Value,type = 'mesh3d')
plot_ly(g2,type = 'surface')





fig <- plot_ly(x = grid2$lr, y = grid2$c, z = grid2$Value) %>% add_surface()

p <- ggplot(g2, aes(c, lr, z= Value)) +
  stat_contour(geom="polygon",aes(fill=stat(level))) +
  scale_fill_distiller(palette = "Spectral", direction = 1)
ggplotly(p)


p <- ggplot(grid, aes(c, lr, z= Value,colour=stat(level))) +
  geom_contour() 
ggplotly(p)


plot_ly(g2, x = ~c, y = ~lr, z = ~Value, type = 'scatter3d', mode = 'lines+markers',
        opacity = 7, 
        line = list(width = 6, colorscale = 'Viridis', reverscale = FALSE)
        )


#install.packages("echarts4r")
library(echarts4r)
g2 |> 
  e_charts(c) |> 
  e_surface(lr, Value, wireframe = list(show = FALSE)) |> 
  e_visual_map(Value)

```



 c    lr noise trainRep lossFun   mean_cFitO mean_cFitG mean_lrFitO mean_lrFitG mean_optimValO mean_optimValG
  <dbl> <dbl> <dbl>    <dbl> <chr>          <dbl>      <dbl>       <dbl>       <dbl>          <dbl>          <dbl>
1   0.1   0.4   500       20 RMSE           4.94           5       1.00         1.09       3.30e- 5        0.115  
2   0.1   0.4   500       20 MAE            0.100          5       0.268        1.09       2.56e+ 1        0.0921 
3   0.1   0.4   500       20 MAPE           0.101          5       0.268        1.09       2.31e+ 0        0.00954
4   0.1   0.4   500       20 MedAE          0.288          5       0.429        1.09       3.65e- 1        0.124  
5   0.1   0.4   500       20 HuberLoss      8.09           5       1.00         1.09       1.98e-17        0.00662



```{r}



```

