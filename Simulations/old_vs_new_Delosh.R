
library(dplyr); library(purrr); library(tidyr)
library(patchwork)

new_version_simulation <- function() {

# Constants
INPUT_LAYER_DEFAULT <- seq(0, 100, 0.5)
OUTPUT_LAYER_DEFAULT <- seq(0, 250, 1)

alm.response <- function(input = 1, c, input.layer, output.layer,weight.mat) {
  input.activation <- exp(-c * (input.layer - input)^2) / sum(exp(-c * (input.layer - input)^2))
  output.activation <- weight.mat %*% input.activation
  output.probability <- output.activation / sum(output.activation)
  mean.response <- sum(output.layer * output.probability)
  list(mean.response = mean.response, input.activation = input.activation, output.activation = output.activation)
}

alm.update <- function(corResp, c, lr, output.layer, input.activation, output.activation, weight.mat) {
  fz <- exp(-c * (output.layer - corResp)^2)
  teacherSignal <- (fz - output.activation) * lr
  wChange <- teacherSignal %*% t(input.activation)
  weight.mat <- weight.mat + wChange
  weight.mat[weight.mat < 0] = 0
  return(weight.mat)
}

alm.trial <- function(input, corResp, c, lr, input.layer, output.layer, weight.mat) {
  alm_resp <- alm.response(input, c, input.layer,output.layer, weight.mat)
  updated_weight.mat <- alm.update(corResp, c, lr, output.layer, alm_resp$input.activation, alm_resp$output.activation, weight.mat)
  return(list(mean.response = alm_resp$mean.response, weight.mat = updated_weight.mat))
}

alm.sim <- function(dat, c, lr, input.layer = INPUT_LAYER_DEFAULT, output.layer = OUTPUT_LAYER_DEFAULT) {
  weight.mat <- matrix(0.00, nrow = length(output.layer), ncol = length(input.layer))
  xt <- dat$x
  n <- nrow(dat)
  st <- numeric(n) # Initialize the vector to store mean responses
  for(i in 1:n) {
    trial <- alm.trial(dat$x[i], dat$y[i], c, lr, input.layer, output.layer, weight.mat)
    weight.mat <- trial$weight.mat
    st[i] <- trial$mean.response
  }
  
  dat <- dat %>% mutate(almResp = st)
  return(list(d = dat, wm = weight.mat, c = c, lr = lr))
}

exam.response <- function(input, c, trainVec, input.layer = INPUT_LAYER_DEFAULT,output.layer = OUTPUT_LAYER_DEFAULT, weight.mat) {
  nearestTrain <- trainVec[which.min(abs(input - trainVec))]
  aresp <- alm.response(nearestTrain, c, input.layer = input.layer,output.layer = OUTPUT_LAYER_DEFAULT,weight.mat)$mean.response
  
  xUnder <- ifelse(min(trainVec) == nearestTrain, nearestTrain, trainVec[which(trainVec == nearestTrain) - 1])
  xOver <- ifelse(max(trainVec) == nearestTrain, nearestTrain, trainVec[which(trainVec == nearestTrain) + 1])
  
  mUnder <- alm.response(xUnder, c, input.layer = input.layer, output.layer, weight.mat)$mean.response
  mOver <- alm.response(xOver, c, input.layer = input.layer,output.layer, weight.mat)$mean.response
  
  exam.output <- round(aresp + ((mOver - mUnder) / (xOver - xUnder)) * (input - nearestTrain), 3)
  exam.output
}



simOrganize <- function(simOut) {
  dat <- simOut$d
  weight.mat <- simOut$wm
  c <- simOut$c
  lr <- simOut$lr
  trainX <- unique(dat$x)
  
  almResp <- generate.data(seq(0,100,.5), type = first(dat$type)) %>% rowwise() %>% 
    mutate(model = "ALM", resp = alm.response(x, c, input.layer = INPUT_LAYER_DEFAULT,output.layer = OUTPUT_LAYER_DEFAULT, weight.mat = weight.mat)$mean.response)
  
  examResp <- generate.data(seq(0,100,.5), type = first(dat$type)) %>% rowwise() %>% 
    mutate(model = "EXAM", resp = exam.response(x, c, trainVec = trainX, input.layer = INPUT_LAYER_DEFAULT,output.layer = OUTPUT_LAYER_DEFAULT, weight.mat))
  
  organized_data <- bind_rows(almResp, examResp) %>% 
    mutate(type = first(dat$type),
           error = abs(resp - y),
           c = c,
           lr = lr,
           type = factor(type, levels = c("linear", "exponential", "quadratic")),
           test_region = ifelse(x %in% trainX, "train", 
                                ifelse(x > min(trainX) & x < max(trainX), "interpolate", "extrapolate")))
  organized_data
}

#### Run Simulation

envTypes <- c("linear", "exponential", "quadratic")
noise=0

lowDensityTrainBlock <- c(30.5, 36.0, 41.0, 46.5, 53.5, 59.0, 64.0, 69.5)
medDensityTrainBlock <- seq(30.0, 70.0, length.out=20)
highDensityTrainBlock <- seq(30.0, 70.0, length.out=50)

lowTrain <- map_dfr(envTypes, ~ generate.data(rep(lowDensityTrainBlock,25), type = .x, noise)) %>% group_by(type) %>% mutate(block = rep(1:25, each = 8),trial=seq(1,200))
medTrain <- map_dfr(envTypes, ~ generate.data(rep(medDensityTrainBlock,10), type = .x, noise)) %>% group_by(type) %>% mutate(block = rep(1:10, each = 20),trial=seq(1,200))
highTrain <- map_dfr(envTypes, ~ generate.data(rep(highDensityTrainBlock,4), type = .x, noise)) %>% group_by(type) %>% mutate(block = rep(1:4, each = 50),trial=seq(1,200))

lowSim <- map(envTypes, ~ alm.sim(lowTrain %>% filter(type == .x), c = .2, lr = .2))
medSim <- map(envTypes, ~ alm.sim(medTrain %>% filter(type == .x), c = .2, lr = .2))
highSim <- map(envTypes, ~ alm.sim(highTrain %>% filter(type == .x), c = .2, lr = .2))

lowSimTest <- map_dfr(lowSim,simOrganize) %>% mutate(density = "low")
medSimTest <- map_dfr(medSim,simOrganize) %>% mutate(density = "med")
highSimTest <- map_dfr(highSim,simOrganize) %>% mutate(density = "high")

simTestAll <- rbind(lowSimTest,medSimTest,highSimTest) %>% group_by(type,density,model) %>%
  mutate(type=factor(type,levels=c("linear","exponential","quadratic")),
         density=factor(density,levels=c("low","med","high"))) %>%
  dplyr::relocate(density,type,test_region)

print(sim_results(lowSim,medSim,highSim, simTestAll))
list(test=simTestAll,train=list(low=lowSim,med=medSim,high=highSim))

}



old_version_simulation <- function() {

alm.response <- function(input=1,c) {
  # input.activation <- exp(-c*(input.layer - input)^2)
  # input.activation <<- input.activation/sum(input.activation)
  input.activation <<- exp(-c * (input.layer - input)^2) / sum(exp(-c * (input.layer - input)^2))
  
  output.activation <<- weight.mat %*% input.activation
  output.probability <<- output.activation/sum(output.activation)
  mean.response <<- sum(output.layer * output.probability)
  mean.response
}

alm.update <- function(corResp,c,lr){
  fz <- exp(-c*(output.layer - corResp)^2)
  teacherSignal <- (fz - output.activation)*lr
  wChange <- teacherSignal %*% t(input.activation)
  weight.mat <<- weight.mat + (wChange)
  weight.mat[weight.mat<0]=0 
  weight.mat <<- weight.mat
}

alm.trial <- function(input, corResp,c,lr){
  alm.response(input,c)
  alm.update(corResp,c,lr)
  mean.response
}

exam.response <- function(input,c){
  trainVec = sort(unique(xt))
  nearestTrain <- trainVec[which.min(abs(input - trainVec))]
  aresp <- alm.response(nearestTrain,c)
  xUnder = ifelse(min(trainVec) == nearestTrain, nearestTrain, trainVec[which(trainVec == nearestTrain) - 1])
  xOver = ifelse(max(trainVec) == nearestTrain, nearestTrain, trainVec[which(trainVec == nearestTrain) + 1])
  mUnder <- alm.response(xUnder,c)
  mOver <- alm.response(xOver,c)
  exam.output = round(aresp + ((mOver - mUnder) / (xOver - xUnder)) * (input - nearestTrain), 3)
  exam.output
}


alm.sim <- function(dat, c, lr,testRange=seq(0,100,.5)){
  
  input.layer <<- matrix(seq(0,100,.5) ) # half step units for inputs, from 0 to 100
  output.layer <<- matrix(seq(0,250,1)) # single step units for outputs, from 0 to 250
  weight.mat <<- matrix(0.00,nrow=length(output.layer),ncol=length(input.layer )) # weights initialized to 0 (as in Delosh 1997)
  xt<<-dat$x
  # run training
  st <- map2_dbl(dat$x, dat$y, ~alm.trial(.x,.y,c,lr))
  # append training data to the data frame
  dat <- dat %>% mutate(almResp = st)
  
  return(list(d=dat,wm=weight.mat,c=c,lr=lr)) # final weightmat is probs incorrect for all but last
}


simOrganize <- function(simOut){
  dat <- simOut$d
  weight.mat <<- simOut$wm
  c <- simOut$c
  lr <- simOut$lr
  
  trainX <- unique(dat$x)
  xt <<- trainX
  
  almResp <- generate.data(seq(0,100,.5), type = first(dat$type)) %>% rowwise() %>% 
    mutate(model="ALM",resp = alm.response(x,c))
  
  examResp <- generate.data(seq(0,100,.5), type = first(dat$type)) %>% rowwise() %>% 
    mutate(model="EXAM",resp = exam.response(x,c))
  
  bind_rows(almResp,examResp) %>% 
    mutate(type=first(dat$type),
           error = abs(resp - y),
           c=c,lr=lr,
           type=factor(type,levels=c("linear","exponential","quadratic"))) %>%
    mutate(test_region = ifelse(x %in% trainX, "train", ifelse(x > min(trainX) & x < max(trainX), "interpolate", "extrapolate")))
    
}

#######


#### Run Simulation

envTypes <- c("linear", "exponential", "quadratic")
noise=0

lowDensityTrainBlock <- c(30.5, 36.0, 41.0, 46.5, 53.5, 59.0, 64.0, 69.5)
medDensityTrainBlock <- seq(30.0, 70.0, length.out=20)
highDensityTrainBlock <- seq(30.0, 70.0, length.out=50)

lowTrain <- map_dfr(envTypes, ~ generate.data(rep(lowDensityTrainBlock,25), type = .x, noise)) %>% group_by(type) %>% mutate(block = rep(1:25, each = 8),trial=seq(1,200))
medTrain <- map_dfr(envTypes, ~ generate.data(rep(medDensityTrainBlock,10), type = .x, noise)) %>% group_by(type) %>% mutate(block = rep(1:10, each = 20),trial=seq(1,200))
highTrain <- map_dfr(envTypes, ~ generate.data(rep(highDensityTrainBlock,4), type = .x, noise)) %>% group_by(type) %>% mutate(block = rep(1:4, each = 50),trial=seq(1,200))

lowSim <- map(envTypes, ~ alm.sim(lowTrain %>% filter(type == .x), c = .2, lr = .2))
medSim <- map(envTypes, ~ alm.sim(medTrain %>% filter(type == .x), c = .2, lr = .2))
highSim <- map(envTypes, ~ alm.sim(highTrain %>% filter(type == .x), c = .2, lr = .2))

lowSimTest <- map_dfr(lowSim,simOrganize) %>% mutate(density = "low")
medSimTest <- map_dfr(medSim,simOrganize) %>% mutate(density = "med")
highSimTest <- map_dfr(highSim,simOrganize) %>% mutate(density = "high")

simTestAll <- rbind(lowSimTest,medSimTest,highSimTest) %>% group_by(type,density,model) %>%
  mutate(type=factor(type,levels=c("linear","exponential","quadratic")),
         density=factor(density,levels=c("low","med","high"))) %>%
  dplyr::relocate(density,type,test_region)

print(sim_results(lowSim,medSim,highSim, simTestAll))
list(test=simTestAll,train=list(low=lowSim,med=medSim,high=highSim))

# simTestAll |> ungroup() |> filter(model=="EXAM", test_region=="extrapolate") %>% 
#   summarise(mean_error=round(mean(error),3), p=paste("old: ", mean_error)) |> as.character() |> print()
# 
# 
# simTestAll %>% ggplot(aes(x=block,y=dev,color=type)) + stat_summary(geom="line",fun=mean,alpha=.4)+
#   stat_summary(geom="point",fun=mean,alpha=.4)+
#   stat_summary(geom="errorbar",fun.data=mean_cl_normal,alpha=.4)+facet_wrap(~density, scales="free_x")

}



# library(microbenchmark)
# 
# results <- microbenchmark(
#   new_version_simulation(),
#   old_version_simulation(),
#   times = 5L
# )
# 
# print(results)



# Unit: seconds
# expr      min       lq     mean   median       uq      max neval cld
# new_version_simulation() 1.163419 1.168211 1.237266 1.197497 1.313758 1.401563    10  a 
# old_version_simulation() 1.311915 1.328216 1.440362 1.478216 1.525288 1.566956    10   b
