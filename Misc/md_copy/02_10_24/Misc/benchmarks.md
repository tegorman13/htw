---
title: Benchmarking 
date: last-modified
categories: [Simulation, ALM, R]
code-fold: show
code-tools: true
execute: 
  warning: false
  eval: false
editor: 
  markdown: 
    wrap: 72
---






```{r}
source(here::here("Functions", "packages.R"))
test <- readRDS(here("data/e1_08-21-23.rds")) |>  filter(expMode2 == "Test") |>
  select(id,condit,bandInt,vb,vx,dist,sdist,bandType)

e1_vxBMM <- brm(vx ~ condit * bandInt + (1 + bandInt|id),
                        data=test,file=paste0(here::here("data/model_cache", "e1_testVxBand_RF_5k")))
new_data_grid=map_dfr(1, ~data.frame(unique(test[,c("id","condit","bandInt")])))


t1 <- system.time({
tidy_pred <- test |> add_predicted_draws(e1_vxBMM) |> mutate(.residual = vx - .prediction); 
predict_per_row <- tidy_pred |> group_by(.row) |> mean_hdi(.prediction,.residual)

})

m1 <- system.time({
  predict_vx <- test |> cbind(predict(e1_vxBMM, test)) |> mutate(resid=vx-Estimate)
}
)

cat(paste0("tidybayes preds = ",round(t1["elapsed"],2)," \n ",
           "predict.brmsfit = ",round(m1["elapsed"],2)))


predict_vx |> group_by(id,condit,bandInt) |> summarise(mean(vx),mean(Estimate),mean(resid))

```
tidybayes preds = 169.94 
predict.brmsfit = 10.27




https://cran.r-project.org/web/packages/brms/vignettes/brms_threading.html

https://discourse.mc-stan.org/t/using-the-apple-m1-gpus-question-from-a-noob/23089/39?page=2

```{r}
#| warning: false

pacman::p_load(tidyverse,tidybayes,brms,broom,broom.mixed,lme4,here,knitr,gt,gghalves,patchwork,ggdist,microbenchmark)
e1 <- readRDS(here("data/e1_08-21-23.rds"))
test <- e1 |> filter(expMode2 == "Test")  

options(mc.cores = 4, brms.backend = "cmdstanr")

n_cores=4
bayes_seed <- 1234
n_iter=20000
n_threads=2


gprior<- c( prior(normal(800, 100), class = Intercept),
    prior(normal(400, 10), class = sigma)
  )


##############

tf3 <- system.time(
fit_int_norm <- brm(vx ~ 1 + condit, 
  data = test,
  family = gaussian(),
  iter=n_iter,
  prior = gprior,
  silent=2,
  cores=n_cores,
  threads = threading(n_threads),
  seed=bayes_seed
)
)
tf3

tf4 <- system.time(
fit_int_norm <- brm(vx ~ 1 + condit +(1|bandInt) + (1 + bandInt|id),
  data = test |> filter(id %in% 1:15),
  family = gaussian(),
  iter=n_iter,
  prior = gprior,
  silent=2,
  cores=n_cores,
  threads = threading(n_threads),
  seed=bayes_seed
) )
tf4


#########


tf <- system.time(
fit_int_norm <- brm(vx ~ 1 + condit, 
  data = test,
  family = gaussian(),
  iter=n_iter,
  prior = gprior,
  silent=2,
  cores=n_cores,
  seed=bayes_seed
)
)
tf

tf2 <- system.time(
fit_int_norm <- brm(vx ~ 1 + condit +(1|bandInt) + (1 + bandInt|id),
  data = test |> filter(id %in% 1:15),
  family = gaussian(),
  iter=n_iter,
  prior = gprior,
  silent=2,
  cores=n_cores,
  seed=bayes_seed
)
)
tf2


##########


cat(paste0("int only gaussian = ",round(tf["elapsed"],2)," \n ",
          # "foreach parallel time=",round(tfP["elapsed"],2)," \n ",
           "irt version = ",round(tf2["elapsed"],2), " \n ",
           "int only threading = ",round(tf3["elapsed"],2)," \n ",
          "irt  threading = ",round(tf4["elapsed"],2)," \n "
            
          ))

```
### default brms settings (2K iterations)
int only gaussian = 6.85 
irt version = 12.15

int only gaussian = 7.51 
irt version = 12

int only gaussian = 6.67 
irt version = 11.82 
int only threading = 7.49 
irt  threading = 11.33 


### default brms settings (4K iterations) - 2 threads
 int only gaussian = 8.1 
irt version = 14.69 
int only threading = 6.99 
irt  threading = 15.24 


### default brms settings (4K iterations) - 4 threads
brms int only gaussian = 8.11 
irt version = 15.05 
int only threading = 7.28 
irt  threading = 16.95 

### default brms settings (4K iterations) - 4 threads - 2 cores
int only gaussian = 11.27 
irt version = 21.72 
int only threading = 8.39 
irt  threading = 22.67 

### default brms settings (4K iterations) - 2 threads - 2 cores
int only gaussian = 11.27 
irt version = 21.81 
int only threading = 9.38 
irt  threading = 20.63 


### default brms settings (4K iterations) - 8 threads - 4 cores
int only gaussian = 8.27 
irt version = 15.4 
int only threading = 7.69 
irt  threading = 21.28 
 
int only gaussian = 8.37 
irt version = 15.18 
int only threading = 8.02 
irt  threading = 21.47  


### default brms settings (4K iterations) - 8 threads - 8 cores

int only gaussian = 8.11 
irt version = 14.88 
int only threading = 7.92 
irt  threading = 20.72 




### default brms settings (10K iterations) - 8 threads - 8 cores

int only gaussian = 11.38 
irt version = 23.64 
int only threading = 11.58 
irt  threading = 37.1 

int only gaussian = 11.53 
irt version = 23.76 
int only threading = 12 
irt  threading = 37.45 

### default brms settings (10K iterations) - 8 threads - 4 cores

int only gaussian = 11.48 
irt version = 23.29 
int only threading = 10.91 
irt  threading = 36.54


### default brms settings (10K iterations) - 4 threads - 4 cores

int only gaussian = 11.99 
irt version = 24.25 
int only threading = 10.13 
irt  threading = 29.13 

### default brms settings (10K iterations) - 4 threads - 1 cores

int only gaussian = 29.35 
irt version = 66.61 
int only threading = 18.97 
irt  threading = 71.21 


### default brms settings (10K iterations) - 2 threads - 2 cores

int only gaussian = 17.98 
irt version = 38.31 
int only threading = 13.87 
irt  threading = 37.78 


### default brms settings (10K iterations) - 1 threads - 4 cores

int only gaussian = 11.14 
irt version = 22.89 
int only threading = 12.25 
irt  threading = 27.03 

int only gaussian = 11.35 
 irt version = 22.99 
 int only threading = 12.28 
 irt  threading = 27.38 


### default brms settings (10K iterations) - 2 threads - 4 cores

int only gaussian = 11.22 
irt version = 23.07 
int only threading = 9 
irt  threading = 22.34 

int only gaussian = 11.31 
 irt version = 23.14 
 int only threading = 9.49 
 irt  threading = 22.21 


### default brms settings (10K iterations) - 3 threads - 4 cores

int only gaussian = 11.15 
irt version = 22.88 
int only threading = 9.4 
irt  threading = 26.75 

int only gaussian = 11 
 irt version = 22.88 
 int only threading = 9.5 
 irt  threading = 25.7 



### (20K iterations) - 1 threads - 4 cores

int only gaussian = 18.24 
 irt version = 40.37 
 int only threading = 20.32 
 irt  threading = 48.94 

int only gaussian = 17.8 
 irt version = 39.26 
 int only threading = 20.54 
 irt  threading = 48.14 
 
 
 int only gaussian = 18.12 
 irt version = 40.31  
 int only threading = 20.6 
 irt  threading = 40.02 - when I commented out the threading specification
  


### (20K iterations) - 2 threads - 4 cores

int only gaussian = 17.82 
irt version = 39.34 
int only threading = 14.22 
irt  threading = 37.73 
 
 
### (20K iterations) - adding condit factor -  2 threads - 4 cores
int only gaussian = 7.78 
 irt version = 45.66 
 int only threading = 8.57 
 irt  threading = 45.97 
 
 
### (20K iterations) - adding condit factor -  3 threads - 4 cores

int only gaussian = 8.58 
irt version = 47.07 
int only threading = 10.71 
irt  threading = 57.87 


### (20K iterations) - condit factor and bandInt RF slope -  2 threads - 4 cores

int only gaussian = 7.85 
 irt version = 80.71 
 int only threading = 8.69 
 irt  threading = 88.3 




 
```{r}

pacman::p_load(tidyverse,foreach,doParallel,future,furrr,here)
purrr::walk(c(here("Functions/alm_functions.R","Functions/Display_Functions.R")),source)

n_cores <- parallel::detectCores()

param_grid <- tibble(crossing(
  c = seq(.5,5,.25),
  lr = seq(0.01, 1,.1),
  noise_sd = c(0,.0001,0.001, 0.01),
  inNodes = c(5, 7,14,28),
  outNodes = c(16, 32,64)
))
nrow(param_grid)

gen_train <- function(trainVec=c(5,6,7),trainRep=3,noise=0){
  bandVec=c(0,100,350,600,800,1000,1200)
  if(class(trainVec)=="list"){trainVec=unlist(trainVec)}
  ts <- rep(seq(1,length(trainVec)),trainRep)
  noiseVec=rnorm(length(ts),mean=0)*noise
  if(noise==0) {noiseVec=noiseVec*0}
  tibble(trial=seq(1,length(ts)),input=trainVec[ts],vx=bandVec[trainVec[ts]]+noiseVec)
}
fit_alm <- function(data, c, lr, noise_sd, inNodes, outNodes) {
  mse_list <- replicate(5, {
    train_data <- data[, c("trial", "input", "cor")] %>% rename("vx" = cor)
    sim_result <- sim_train(
      dat = train_data,
      c = c,
      lr = lr,
      inNodes = inNodes,
      outNodes = outNodes,
      noise_sd = noise_sd
    )
    train_data$almTrain <- sim_result$almTrain
    mse <- mean((data$vx - train_data$almTrain)^2)
    mse
  })
  avg_mse <- mean(mse_list)
  return(avg_mse)
}

gt <- gen_train(trainVec=c(5,6,7),trainRep=8) %>% mutate(cor=vx,err=(800-0)*exp(-.1*seq(1,n()))+0,vx=cor-err)


furrr::furrr_options(seed = TRUE)
plan(multisession, workers = n_cores-1)
tff <- system.time({

param_grid <- param_grid %>% mutate(performance = future_map_dbl(seq_len(nrow(.)), function(idx) {
    fit_alm(gt, c = c[idx], lr = lr[idx], noise_sd = noise_sd[idx], inNodes = inNodes[idx], outNodes = outNodes[idx])
  },.options = furrr_options(seed = T)))
  best_paramsF <- param_grid %>%
    arrange((performance)) 
  bestF <- head(best_paramsF,1)
})


# cluster <- parallel::makeCluster(n_cores-1)                 
# doParallel::registerDoParallel(cluster)
# tfP <- system.time({
#   param_grid <- param_grid %>%
#     mutate(performance = foreach(idx = seq_len(nrow(.)), .combine = c) %dopar% {
#       fit_alm(gt, c = c[idx], lr = lr[idx], noise_sd = noise_sd[idx], inNodes = inNodes[idx], outNodes = outNodes[idx])
#     })
#   best_paramsP <- param_grid %>%
#     arrange((performance)) 
#   bestP <- head(best_paramsP,1)
#   })
# stopImplicitCluster()
#   tfP


tI <- system.time({
gt <- gen_train(trainVec=c(5,6,7),trainRep=8) %>% mutate(cor=vx,err=(600-0)*exp(-.1*seq(1,n()))+0,vx=cor-err)

param_grid <- param_grid %>%
  mutate(performance = map_dbl(seq_len(nrow(.)), ~ {
    fit_alm(gt, c = c[.x], lr = lr[.x], noise_sd = noise_sd[.x], inNodes = inNodes[.x], outNodes = outNodes[.x])
  }))
best_paramsI <- param_grid %>%
  arrange((performance)) 
bestI <- head(best_paramsI,1)
})

tI




cat(paste0("furr time=",round(tff["elapsed"],2)," \n ",
          # "foreach parallel time=",round(tfP["elapsed"],2)," \n ",
           " Standard Time=",round(tI["elapsed"],2)))

```

4.2 -6.2 times faster with furr on m1

| Runs | Machine   | Cores | Method   | Time (seconds) |
|------|-----------|-------|----------|----------------|
| 9120 | M1        | 9/10  | Furr     | 9.56           |
| 9120 | M1        | 9/10  | Standard | 58.8           |
| 912  | M1        | 8/10  | Standard | 5.7            |
| 912  | M1        | 8/10  | Parallel | 1.2            |
| 912  | M1        | 9/10  | Furr     | 1.3            |
| 912  | 2015 iMac | 2/4   | Standard | 26.2           |
| 912  | 2015 iMac | 2/4   | Parallel | 18.5           |
| 912  | 2015 iMac | 3/4   | Furr     | 12.88          |
| 912  | GTX       | 3/4   | Standard | 25.32          |
| 912  | GTX       | 3/4   | Furr     | 13.1           |

------------------------------------------------------------------------

Furrr and standard version both work equally fast, tested with up to 120
simulation repetitions

```{r}
library(furrr)
furrr::furrr_options(seed = TRUE)
plan(multisession, workers = parallel::detectCores())

parmVec <- tibble(crossing(c = c(0.1,.5), lr = c(0.4), noise = c(500), trainRep = c(20), lossFun = list("RMSE", "RMSE.blocked"), simNum = 1:30))
tf <- system.time(
sdpf <- parmVec %>% 
  mutate(d = future_pmap(list(c, lr, noise, trainRep), ~sim_data(c = ..1, lr = ..2, noise = ..3, trainRep = ..4),
                         .options = furrr_options(seed = T)),
                           almTrainDat = map(d, "almTrain"),
                          almTestDat = map(d, "almPred"),
                          examTestDat = map(d, "examPred"),
                          td = map(trainRep, ~gen_train(trainRep = .)),
                          fitO = map2(td, lossFun, ~wrap_optim(.x, .y)),
                          fitG = map2(td, lossFun, ~wrap_grid(.x, .y)),
                          cFitO = map_dbl(fitO, "c"),
                          lrFitO = map_dbl(fitO, "lr"),
                          optimValO = map_dbl(fitO, "Value"),
                          cFitG = map_dbl(fitG, "c"),
                          lrFitG = map_dbl(fitG, "lr"),
                          optimValG = map_dbl(fitG, "Value"))
)
tf

ts <- system.time(
sdp <- parmVec %>% mutate(d = pmap(list(c, lr, noise, trainRep), ~sim_data(c = ..1, lr = ..2, noise = ..3, trainRep = ..4)),
                          almTrainDat = map(d, "almTrain"),
                          almTestDat = map(d, "almPred"),
                          examTestDat = map(d, "examPred"),
                          td = map(trainRep, ~gen_train(trainRep = .)),
                          fitO = map2(td, lossFun, ~wrap_optim(.x, .y)),
                          fitG = map2(td, lossFun, ~wrap_grid(.x, .y)),
                          cFitO = map_dbl(fitO, "c"),
                          lrFitO = map_dbl(fitO, "lr"),
                          optimValO = map_dbl(fitO, "Value"),
                          cFitG = map_dbl(fitG, "c"),
                          lrFitG = map_dbl(fitG, "lr"),
                          optimValG = map_dbl(fitG, "Value"))
)
ts


```

M1 Max times: tf user system elapsed 221.326 4.517 226.233

ts user system elapsed 221.140 4.288 225.330

Mac Pro times: tf user system elapsed 1125.102 76.859 1474.612

ts user system elapsed 1132.716 76.260 1493.154

# Benchmarking ALM Model Fit Functions

```{r}
#| code-fold: true
pacman::p_load(tidyverse,data.table,microbenchmark())


d <- readRDS(here('dPrune-01-19-23.rds'))

dtest <- d %>% filter(expMode %in% c("test-Nf","test-train-nf")) %>% group_by(id,lowBound) %>% 
  mutate(nBand=n(),band=bandInt,id=factor(id)) %>% group_by(id) %>% mutate(nd=n_distinct(lowBound))
# unique(dtest[dtest$nd==4,]$sbjCode) # 7 in wrong condition
dtest <- dtest %>% group_by(id,lowBound) %>% filter(nBand>=5 & nd==6)
# for any id that has at least 1 nBand >=5, remove all rows with that id. 
dtest <- dtest %>% group_by(id) %>% filter(!id %in% unique(dtest$id[dtest$nBand<5]))

dtestAgg <- dtest %>% group_by(id,condit,catOrder,feedbackType,vb,band,lowBound,highBound,input) %>% mutate(vxCapped=ifelse(vx>1600,1600,vx)) %>%
  summarise(vxMean=mean(vx),devMean=mean(dist),vxMed=median(vx),devMed=median(dist),
            vxMeanCap=mean(vxCapped),.groups = "keep")

# select first row for each id in d, then create histogram for nTrain
#  d  %>% group_by(id) %>% slice(1) %>% ggplot(aes(nTrain)) + geom_histogram() + facet_wrap(~condit)
  
ds <- d %>% filter(expMode %in% c("train","train-Nf","test-Nf","test-train-nf")) %>% 
filter(!id %in% unique(dtest$id[dtest$nBand<5])) %>% 
select(id,condit,catOrder,feedbackType,expMode,trial,gt.train,vb,band,bandInt,lowBound,highBound,input,vx,dist,vxb) 

dst <- ds %>% filter(expMode=="train",catOrder=="orig")

vTrainTrial <- dst %>% filter(condit=="Varied",gt.train<=84) %>% group_by(gt.train,vb) %>% summarise(sdVx=sd(vx),vx=mean(vx),sdDist=sd(dist),dist=mean(dist)) %>% 
  group_by(vb) %>% mutate(gt.trainBin=cut(gt.train,breaks=5,labels=c(1:5)))

binTrainTrial <- dst %>% filter(gt.train<=83) %>% group_by(gt.train,vb,condit) %>% summarise(sdVx=sd(vx),vx=mean(vx),sdDist=sd(dist),dist=mean(dist)) %>% 
  group_by(vb) %>% mutate(gt.trainBin=cut(gt.train,breaks=6,labels=c(1:6)))


tMax=84
bandVec <- rep(c(800,1000,1200),each=tMax/3)
bandVec <- bandVec[sample(1:length(bandVec),tMax,replace=FALSE)]

trainTrials <- dst %>% filter(gt.train<=tMax) %>% group_by(condit,gt.train,vb,bandInt,input) %>% summarise(vx=mean(vx)) 

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

wrap_alm <- function(parms,dat, weights,lossFun){
    c=parms[1]; lr=parms[2]
   pred=train.alm(dat, c=c, lr=lr, weights=weights)
   #sqrt(mean((dat$vx -pred)^2))
   lossFun(dat$vx,pred)
}

wrap_optim <- function(dat,wm,lossFun){
  bounds_lower <- c(.0000001, .00001)
  bounds_upper <- c(5, 5)

 optim(c(.1, .2),
   fn = wrap_alm,
   dat = dat, weights = wm,lossFun=lossFun,
   method = "L-BFGS-B",
   lower = bounds_lower,
   upper = bounds_upper,
   control = list(maxit = 1e4, pgtol = 0, factr = 0)
 )
}

RMSE <- function(x,y){
  sqrt(mean((x-y)^2))
}

## First average observed and predicted data into blocks, then compute RMSE
RMSE.tb <- function(x,y,blocks=6){
  data.frame(x,y) %>% mutate(t=row_number(),fitBins=cut(t,breaks=blocks,labels=c(1:blocks))) %>%
    group_by(fitBins) %>% 
    summarise(predMean=mean(x),obsMean=mean(y)) %>% 
    summarise(RMSE(predMean,obsMean)) %>% as.numeric()
}

## Recode RMSE.tb using data.table functions rather than dplyr
RMSE.tb2 <- function(x,y,blocks=6){
  data.table(x=x,y=y,t=seq(1,length(x))) %>% 
    .[, `:=`(fitBins = cut(t, breaks = ..blocks, labels = c(1:..blocks)))] %>%
    .[, .(predMean = mean(x), obsMean = mean(y)), keyby = .(fitBins)] %>%
    .[, RMSE(predMean,obsMean)] %>% as.numeric()
}

```

### dplyr RMSE vs data.table RMSE

```{r}
dpVsDt=microbenchmark(
  dplyrMethod={
  fitVaried <- tv %>% filter(condit=="Varied") %>% wrap_optim(.,wm,lossFun=RMSE.tb);
  fitConstant <- tv %>% filter(condit=="Constant") %>% wrap_optim(.,wm,lossFun=RMSE.tb)},
  dtMethod={
  fitVaried2 <- tv %>% filter(condit=="Varied") %>% wrap_optim(.,wm,lossFun=RMSE.tb2);
  fitConstant2 <- tv %>% filter(condit=="Constant") %>% wrap_optim(.,wm,lossFun=RMSE.tb2)},
  times=5
)
knitr::kable(summary(dpVsDt),format="markdown")
```

| expr        |       min |        lq |      mean |    median |        uq |       max | neval | cld |
|:-------|-------:|-------:|-------:|-------:|-------:|-------:|-------:|:-------|
| dplyrMethod | 11.282291 | 12.996857 | 13.558623 | 13.970367 | 14.658515 | 14.885086 |     5 | a   |
| dtMethod    |  4.933999 |  5.306232 |  5.681235 |  5.515397 |  6.131193 |  6.519355 |     5 | b   |

The data.table version seems to be consistently more than 2x faster

### Separate fits, vs. nesting method vs. split method

```{r}
nestSplit<-microbenchmark(
separate={fitVaried <- tv %>% filter(condit=="Varied") %>% wrap_optim(.,wm,lossFun=RMSE.tb2);
fitConstant <- tv %>% filter(condit=="Constant") %>% wrap_optim(.,wm,lossFun=RMSE.tb2) },
nestGroups = tv %>% group_by(condit) %>% nest() %>% mutate(fit=map(data,~wrap_optim(.,wm,RMSE.tb2))),
splitGroups = tv %>% split(.$condit) %>% map(~wrap_optim(.,wm,RMSE.tb2)),
times=5
)
knitr::kable(summary(nestSplit),format="markdown") # 03/02/23 - Mac Pro
```

| expr        |      min |       lq |     mean |   median |       uq |      max | neval | cld |
|:-------|-------:|-------:|-------:|-------:|-------:|-------:|-------:|:-------|
| separate    | 6.394459 | 6.483235 | 6.709882 | 6.508353 | 7.059420 | 7.103941 |     5 | a   |
| nestGroups  | 6.253850 | 6.457956 | 6.789962 | 6.756942 | 7.044965 | 7.436094 |     5 | a   |
| splitGroups | 6.117184 | 6.455866 | 6.605814 | 6.461572 | 6.640798 | 7.353648 |     5 | a   |

Not much of an effect for nesting method.

### Computing RMSE over Raw trials vs. RMSE of blocked training performance

```{r}
trialBlock<-microbenchmark(
  trialFit = tv %>% split(.$condit) %>% map(~wrap_optim(.,wm,RMSE)),
  blockFit = tv %>% split(.$condit) %>% map(~wrap_optim(.,wm,RMSE.tb2)),
  times=5
)
knitr::kable(summary(trialBlock),format="markdown")
```

| expr     |      min |       lq |     mean |   median |       uq |      max | neval | cld |
|:-------|-------:|-------:|-------:|-------:|-------:|-------:|-------:|:-------|
| trialFit | 2.665184 | 2.674186 | 2.806476 | 2.803491 | 2.843806 | 3.045712 |     5 | a   |
| blockFit | 5.169167 | 5.264218 | 5.350769 | 5.285151 | 5.491828 | 5.543483 |     5 | b   |

The models are fit about 2x faster when computing RMSE over trials. This
shouldn't be surprising, since fitting to blocked data requires several
additional computations (e.g. grouping, computing means)
