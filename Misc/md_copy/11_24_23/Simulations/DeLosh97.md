library(tidyverse)

source(here::here("Functions","deLosh_data.R"))

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
  trainX <- round(unique(dat$x)/.5)*.5
  
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

generateSimData <- function(density, envTypes, noise) {
  reps <- 200 / length(trainingBlocks[[density]])
  map_dfr(envTypes, ~ 
            generate.data(rep(trainingBlocks[[density]], reps), type = .x, noise)) |>
    group_by(type) |>
    mutate(block = rep(1:reps, each = length(trainingBlocks[[density]])),
           trial=seq(1,200))
}

simulateAll <- function(density,envTypes, noise, c = .2, lr = .2) {
  trainMat <- generateSimData(density, envTypes, noise)
  trainData <- map(envTypes, ~ alm.sim(trainMat %>% filter(type == .x), c = c, lr = lr))
  assign(paste(density),list(train=trainData, test=map_dfr(trainData, simOrganize) %>% mutate(density = density)))
}

#### Run Simulation
noise=0
envTypes <- c("linear", "exponential", "quadratic")
densities <- c("low", "med", "high")

results <- map(densities, ~ simulateAll(.x, envTypes, noise)) |>
  set_names(densities) 

trainAll <- results %>%
  map_df(~ map_df(.x$train, pluck, "d"), .id = "density") |>
  mutate(stage=as.numeric(cut(trial,breaks=20,labels=seq(1,20))),
         dev=sqrt((y-almResp)^2),
         density=factor(density,levels=c("low","med","high")),
         type=factor(type,levels=c("linear","exponential","quadratic"))) |>
  dplyr::relocate(density,type,stage)

simTestAll <- results |> map("test") |> bind_rows() |>
  group_by(type,density,model) %>%
  mutate(type=factor(type,levels=c("linear","exponential","quadratic")),
         density=factor(density,levels=c("low","med","high"))) %>%
  dplyr::relocate(density,type,test_region)



trainAll %>% ggplot(aes(x=block,y=dev,color=type)) + stat_summary(geom="line",fun=mean,alpha=.4)+
  stat_summary(geom="point",fun=mean,alpha=.4)+
  stat_summary(geom="errorbar",fun.data=mean_cl_normal,alpha=.4)+facet_wrap(~density, scales="free_x")



simTestAll %>% ggplot(aes(x=x,y=y)) + 
  geom_point(aes(x=x,y=resp,shape=model,color=model),alpha=.7,size=1) + 
  geom_line(aes(x=x,y=y),alpha=.4)+ 
  geom_point(data=simTestAll %>% filter(test_region=="train"),aes(x=x,y=y),color="black",size=1,alpha=1) +
  facet_grid(density~type) + 
  theme_bw() + theme(legend.position="bottom")


# t <- map(results, pluck, "train") 
# str(t,max.level=2)
# head(t$low[[1]]$d)

# sim_results <- function(lowSim,medSim,highSim, simTestAll)
# {
#   simAll <- rbind(bind_rows(lowSim %>% map("d")) %>% mutate(density = "low"), 
#                   bind_rows(medSim %>% map("d")) %>% mutate(density = "med"), 
#                   bind_rows(highSim %>% map("d")) %>% mutate(density = "high"))
#   
#   simAll <- simAll %>% mutate(stage=as.numeric(cut(trial,breaks=20,labels=seq(1,20))),
#                               dev=sqrt((y-almResp)^2),
#                               #reorder density factor levels
#                               density=factor(density,levels=c("low","med","high")),
#                               type=factor(type,levels=c("linear","exponential","quadratic"))) %>%
#     dplyr::relocate(density,type,stage)
#   
#   tp <- simAll %>% ggplot(aes(x=block,y=dev,color=type)) + stat_summary(geom="line",fun=mean,alpha=.4)+
#     stat_summary(geom="point",fun=mean,alpha=.4)+
#     stat_summary(geom="errorbar",fun.data=mean_cl_normal,alpha=.4)+facet_wrap(~density, scales="free_x")
#   
#   testp <- simTestAll %>% ggplot(aes(x=x,y=y)) + 
#     geom_point(aes(x=x,y=resp,shape=model,color=model),alpha=.7,size=1) + 
#     geom_line(aes(x=x,y=y),alpha=.4)+ 
#     geom_point(data=simTestAll %>% filter(test_region=="train"),aes(x=x,y=y),color="black",size=1,alpha=1) +
#     facet_grid(density~type) + 
#     theme_bw() + theme(legend.position="bottom")
#   
#   simTestAll %>% group_by(test_region,model,density) %>%
#     summarise(mean_error=mean(error)) %>%
#     pivot_wider(names_from = test_region,values_from = mean_error) |> pander::pandoc.table()
#   
#   tp + testp
# }

# lowDensityTrainBlock <- c(30.5, 36.0, 41.0, 46.5, 53.5, 59.0, 64.0, 69.5)
# medDensityTrainBlock <- seq(30.0, 70.0, length.out=20)
# highDensityTrainBlock <- seq(30.0, 70.0, length.out=50)



# lowTrain <- map_dfr(envTypes, ~ generate.data(rep(lowDensityTrainBlock,25), type = .x, noise)) %>% group_by(type) %>% mutate(block = rep(1:25, each = 8),trial=seq(1,200))
# medTrain <- map_dfr(envTypes, ~ generate.data(rep(medDensityTrainBlock,10), type = .x, noise)) %>% group_by(type) %>% mutate(block = rep(1:10, each = 20),trial=seq(1,200))
# highTrain <- map_dfr(envTypes, ~ generate.data(rep(highDensityTrainBlock,4), type = .x, noise)) %>% group_by(type) %>% mutate(block = rep(1:4, each = 50),trial=seq(1,200))
# 
# lowSim <- map(envTypes, ~ alm.sim(lowTrain %>% filter(type == .x), c = .2, lr = .2))
# medSim <- map(envTypes, ~ alm.sim(medTrain %>% filter(type == .x), c = .2, lr = .2))
# highSim <- map(envTypes, ~ alm.sim(highTrain %>% filter(type == .x), c = .2, lr = .2))
# 
# lowSimTest <- map_dfr(lowSim,simOrganize) %>% mutate(density = "low")
# medSimTest <- map_dfr(medSim,simOrganize) %>% mutate(density = "med")
# highSimTest <- map_dfr(highSim,simOrganize) %>% mutate(density = "high")
# 
# simTestAll <- rbind(lowSimTest,medSimTest,highSimTest) %>% group_by(type,density,model) %>%
#   mutate(type=factor(type,levels=c("linear","exponential","quadratic")),
#          density=factor(density,levels=c("low","med","high"))) %>%
#   dplyr::relocate(density,type,test_region)

# print(sim_results(lowSim,medSim,highSim, simTestAll))
# DeLosh_Results=list(test=simTestAll,train=list(low=lowSim,med=medSim,high=highSim))

# trainAll <- rbind(bind_rows(lowResults$train %>% map("d")) %>% mutate(density = "low"), 
#                   bind_rows(medResults$train %>% map("d")) %>% mutate(density = "med"), 
#                   bind_rows(highResults$train %>% map("d")) %>% mutate(density = "high")) |>
#   mutate(stage=as.numeric(cut(trial,breaks=20,labels=seq(1,20))),
#          dev=sqrt((y-almResp)^2),
#          density=factor(density,levels=c("low","med","high")),
#          type=factor(type,levels=c("linear","exponential","quadratic"))) |>
#   dplyr::relocate(density,type,stage)
