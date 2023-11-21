#
#
#
#
#
#
#
#
#
#
pacman::p_load(tidyverse,data.table,patchwork,glue,knitr,kableExtra,here)
options(dplyr.summarise.inform=FALSE)
purrr::walk(here(c("Functions/alm_functions.R","Functions/Display_Functions.R")),source)


update.weights.with_noise <- function(x.new, y.new, weights, c, lr, noise_sd){
  y.feedback.activation <- exp(-1 * c * (y.new - outputNodes)^2)
  x.feedback.activation <- output.activation(x.new, weights, c)
  weight_updates <- lr * (y.feedback.activation - x.feedback.activation) %*% t(input.activation(x.new, c))
  noise <- matrix(rnorm(nrow(weight_updates) * ncol(weight_updates), sd = noise_sd), 
                  nrow = nrow(weight_updates), ncol = ncol(weight_updates))
  updated_weights <- weights + weight_updates + noise
  return(updated_weights)
}


update.weights<-function(x.new, y.new, weights, c, lr, noise_sd = NULL){
  y.feedback.activation<-exp(-1*c*(y.new-outputNodes)^2)
  x.feedback.activation<-output.activation(x.new, weights, c)
  return(weights+lr*(y.feedback.activation-x.feedback.activation)%*%t(input.activation(x.new, c)))
}

sim_data <- function(dat, c=0.5, lr=0.2, inNodes=7, outNodes=32, trainVec=c(5,6,7),noise_sd=0,update_func="update.weights" ) {
  inputNodes <<- seq(1,7,length.out=inNodes*1)  
  outputNodes <<- seq(50,1600,length.out=outNodes*1) 
  #print(length(outputNodes))
  wm=matrix(.0000,nrow=length(outputNodes),ncol=length(inputNodes))
  # wm=matrix(rnorm(length(outputNodes)*length(inputNodes),.1,5),nrow=length(outputNodes),ncol=length(inputNodes))
  tt<-trainTest.alm(dat, c, lr, wm, trainVec, update_func, noise_sd)
}   

trainTest.alm <- function(dat, c=0.05, lr=0.5, weights, testVec, update_func, noise_sd) {
   alm.train <- rep(NA, nrow(dat))  
   update_func=get(update_func)
   decay_factor = 0.79
  for (i in 1:nrow(dat)) {
    #lr = lr * decay_factor^i
    resp = mean.prediction(dat$input[i], weights, c)
    weights <- update_func(dat$input[i], dat$vx[i], weights, c, lr, noise_sd)
    alm.train[i] = resp
    weights[weights<0] = 0
  }
  almPred <- sapply(testVec, mean.prediction, weights, c)
  examPred <- sapply(testVec, exam.prediction, weights, c, trainVec=c(1,sort(unique(dat$input))))
  list(almTrain=alm.train, almPred=almPred, examPred=examPred,weights=weights)
}

#
#
#
#
#
#| warning: false
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

#
#
#
#
# Define parameters
params <- tibble(crossing(
  c = c(.5,5),
  lr = c(.05,1),
  noise = c(0),
  inNodes = c(7),
  outNodes = c(32),
  trainVec = list(list(5,6,7)),
  trainRep = c(9),
  lossFun = list("MAE"),
  simNum = 1:1,
  update_func = list("update.weights"),
  update_func_name = c("uW"),
  noise_sd = c(0)
))

# Generate training data
params <- params %>% 
  mutate(
    id = seq(1, nrow(.)),
    td = pmap(list(trainVec, trainRep, noise), ~gen_train(trainVec = .x, trainRep = ..2, noise = ..3))
  )

# Run simulations
params <- params %>% 
  mutate(
    d = pmap(list(td, c, lr, update_func, noise_sd, inNodes, outNodes), 
              ~sim_data(dat = .x, c = ..2, lr = ..3, update_func = ..4, noise_sd = ..5, inNodes = ..6, outNodes = ..7)),
    almTrainDat = map(d, "almTrain"),
    weights = map(d, "weights")
  )

# Unnest and select relevant columns
params <- params %>% 
  unnest(c(almTrainDat, td)) %>% 
  select(-d) %>% 
  mutate(input = as.factor(input))

# Apply pf function and trainTab
result <- params %T>% 
  {pf(.) } %>% 
  trainTab

#
#
#
#
# Define a function to fit the model and return the parameters
fit_model <- function(data, initial_c, initial_lr) {
  # Fit the model here and extract the parameters
  # This is a placeholder and should be replaced with your actual model fitting code
  fit_c = initial_c #+ rnorm(1, 0, 0.1)
  fit_lr = initial_lr #+ rnorm(1, 0, 0.1)
  
  tibble(
    gen_c = initial_c,
    gen_lr = initial_lr,
    fit_c = fit_c,
    fit_lr = fit_lr,
    error_c = fit_c - initial_c,
    error_lr = fit_lr - initial_lr
  )
}

params <- tibble(crossing(
  c = seq(.01,2,length.out=10),
  lr = seq(.01,2,length.out=10)
))
# Run the simulations
results <- params %>%
  mutate(simulation = map2(c, lr, ~fit_model(td, .x, .y))) %>%
  unnest(simulation)


fit_model <- function(data, initial_c, initial_lr) {
  # Simulate data from the ALM model
  sim_data <- sim_data(dat = data, c = initial_c, lr = initial_lr, 
                       update_func = "update.weights", noise_sd = 0, 
                       inNodes = 7, outNodes = 32)
  
  # Extract the fitted parameters
  fit_c = sim_data$c
  fit_lr = sim_data$lr
  
  tibble(
    gen_c = initial_c,
    gen_lr = initial_lr,
    fit_c = fit_c,
    fit_lr = fit_lr,
    error_c = fit_c - initial_c,
    error_lr = fit_lr - initial_lr
  )
}

# Run the simulations
results <- params %>%
  mutate(simulation = map2(c, lr, ~fit_model(td, .x, .y))) %>%
  unnest(simulation)

# Print the results
print(results)

# Print the results
print(results %>% select(c,lr,gen_c,gen_lr,fit_c,fit_lr,error_c,error_lr))

ggplot(results, aes(x = gen_c, y = fit_c)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  labs(x = "Generating c", y = "Fitted c", title = "Parameter Recovery for c")

# Plot for lr parameter
ggplot(results, aes(x = gen_lr, y = fit_lr)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  labs(x = "Generating lr", y = "Fitted lr", title = "Parameter Recovery for lr")

#
#
#
#
#
#
#
#
#| eval: false
pacman::p_load(tidyverse,data.table,knitr,kableExtra,glue)
input.activation<-function(x.target, c){return(exp((-1*c)*(x.target-inputNodes)^2))}
output.activation<-function(x.target, weights, c){return(weights%*%input.activation(x.target, c))}
mean.prediction<-function(x.target, weights, c){
  probability<-output.activation(x.target, weights, c)/sum(output.activation(x.target, weights, c))
  return(outputNodes%*%probability) # integer prediction
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
  list(almTrain = alm.train, weights = weights)
  }
sim_train <- function(dat, c=0.5, lr=0.2, inNodes=7, outNodes=32, trainVec=c(5,6,7),noise_sd=0,update_func="update.weights" ) {
  inputNodes <<- seq(1,7,length.out=inNodes*1)  
  outputNodes <<- seq(50,1600,length.out=outNodes*1) 
  wm=matrix(.0000,nrow=length(outputNodes),ncol=length(inputNodes))
  tt<-train.alm(dat, c, lr, wm)
} 
gen_train <- function(trainVec = c(5, 6, 7), trainRep = 3, noise = 0) {
    bandVec <- c(0, 100, 350, 600, 800, 1000, 1200)
    if (class(trainVec) == "list") {trainVec <- unlist(trainVec)}
    ts <- rep(seq(1, length(trainVec)), trainRep)
    noiseVec <- rnorm(length(ts), mean = 0) * noise
    if (noise == 0) {noiseVec <- noiseVec * 0}
    tibble(trial = seq(1, length(ts)), input = trainVec[ts], vx = bandVec[trainVec[ts]] + noiseVec)
}

tibble(crossing(
    c = c(.5, 5), lr = c(.05, 1), noise = c(0),
    inNodes = c(7), outNodes = c(32),
    trainVec = list(list(5, 6, 7)), trainRep = c(9),
    lossFun = list("MAE"),
    simNum = 1:1,
)) %>%
    mutate(id = seq(1, nrow(.)), td = pmap(list(trainVec, trainRep, noise), ~ gen_train(trainVec = .x, trainRep = ..2, noise = ..3))) %>%
    ungroup() %>%
    mutate(
        d = pmap(
            list(td, c, lr,inNodes, outNodes),
            ~ sim_train(dat = .x, c = ..2, lr = ..3, inNodes = ..4, outNodes = ..5)
        ),
        almTrainDat = map(d, "almTrain"), weights = map(d, "weights")
    ) %>%
    unnest(c(almTrainDat, td)) %>%
    select(-d) %>%
    mutate(input = as.factor(input)) %>%
    trainTab()

tibble(crossing(
    c = c(.5, 5), lr = c(.05, 1), noise = c(0),
    inNodes = c(7), outNodes = c(32),
    trainVec = list(list(5, 6, 7)), trainRep = c(9),
    lossFun = list("MAE"),
    simNum = 1:1,
)) %>%
    mutate(id = seq(1, nrow(.)), td = pmap(list(trainVec, trainRep, noise), ~ gen_train(trainVec = .x, trainRep = ..2, noise = ..3))) %>%
    ungroup() %>%
    mutate(
        d = pmap(
            list(td, c, lr,inNodes, outNodes),
            ~ sim_train(dat = .x, c = ..2, lr = ..3, inNodes = ..4, outNodes = ..5)
        ),
        almTrainDat = map(d, "almTrain"), weights = map(d, "weights")
    ) %>%
    unnest(c(almTrainDat, td)) %>%
    select(-d) %>%
    mutate(input = as.factor(input)) %>%
    trainTab()
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
gt <- gen_train(trainVec=c(5,6,7),trainRep=18) %>% mutate(cor=vx,err=(400-0)*exp(-.1*seq(1,n()))+0,vx=cor-err)

bias <- 1000; 
gt <- gen_train(trainVec=c(5,6,7),trainRep=228,noise=0) %>% mutate(
  cor = vx,
  err = (bias - 0) * exp(-.005 * seq(1, n())) + 0,
  en = map2_dbl(err,cor, ~rnorm(n = 1, mean = .y, sd = .x/2)),
  enAvg = map2_dbl(err,cor, ~mean(rnorm(n = 1, mean = .y, sd = .x))),
  weight = (seq(1, n()) - 1) / (n() - 1),
  vx = (weight*en)+bias*(1-weight),
  vx=en
)
gt %>% ggplot(aes(x = trial, y = vx, color = as.factor(input))) +
  geom_line() + ylim(c(-10,1600))


k=wrap_optim(gt,lossFun = "MAE")
s=sim_data(dat=mutate(gt,vx=cor),c=k %>% pluck("c"),lr= k %>% pluck("lr"))
ggp <- gt %>% mutate(pred=s %>% pluck("almTrain"),c=k %>% pluck("c"),lr= k %>% pluck("lr"),input=as.factor(input)) %>%  
  ggplot(aes(x = trial, y = pred, color = input)) +
  geom_line() + ylim(c(0,1600))
ggo <-  gt %>% ggplot(aes(x = trial, y = vx, color = as.factor(input))) +
  geom_line() + ylim(c(-400,1600))

ggo+ggp
#
#
#
#
#
#
#
#| eval: false
#k = t[1,]
#image(matrix(unlist(k$weights),nrow=7))
# mutate(md=map(weights,~matrix(unlist(.),nrow=7)))

wms <- t %>% filter(trial==1) %>% 
  mutate(k=map(weights,~ as.data.frame(matrix(unlist(.),nrow=7)) %>% 
                 rownames_to_column("Input") %>%
                 pivot_longer(-c(Input), names_to = "Output",
                              names_prefix = "V", values_to = "weight") %>%  mutate(Output= fct_relevel(Output,unique(.$Output)))))

wms %>% unnest(k) %>% ggplot(.,aes(x=Input, y=Output, fill=weight)) + 
  geom_raster() + 
  scale_fill_viridis_c()+facet_wrap(~id+c)



md=matrix(unlist(k$weights),nrow=7)
kd=as.data.frame(matrix(unlist(k$weights),nrow=7))%>%
  rownames_to_column("Input") %>%
  pivot_longer(-c(Input), names_to = "Output",names_prefix = "V", values_to = "weight")

kd %>% 
  mutate(Output= fct_relevel(Output,unique(kd$Output))) %>%
  ggplot(aes(x=Input, y=Output, fill=weight)) + 
  geom_raster() + 
  scale_fill_viridis_c()

#
#
#
#
#
#
pv <- t <- parmVec <- tibble(crossing(
  c = c(0.00003),
  lr = c(0.051),
  noise = c(0),
  inNodes=c(7),
  outNodes=c(32),
  trainVec=list(list(1,5,6,7),list(5,6,7)),
  trainRep = c(4),
  lossFun = list("MAE"),
  simNum = 1:1,
  update_func = list("update.weights"),
  update_func_name = c("update.weights"),
  noise_sd = c(0)
)) %>% mutate(id=seq(1,nrow(.)))

#
#
#
#

parmVec <- tibble(crossing(
  c = c(0.1),
  lr = c(0.1, 0.4, 1),
  noise = c(10),
  trainRep = c(20),
  lossFun = list("MAE"),
  simNum = 1:10,
  update_func = list("update.weights", "update.weights.with_noise"),
  update_func_name = c("update.weights", "update.weights.with_noise"),
  noise_sd = c(0.1, 0.5)
))



t <- parmVec %>%
  group_by(simNum, c, lr, update_func,noise_sd) %>%
  mutate(td = list(gen_train(trainRep = first(trainRep), noise = first(noise)))) %>%
  ungroup() %>%
  mutate(d = pmap(list(td, c, lr, update_func,noise_sd), 
                  ~sim_data(dat = .x, c = ..2, lr = ..3, update_func = ..4, noise_sd = ..5)),
         almTrainDat = map(d, "almTrain"))%>%
  unnest(almTrainDat, td) %>%
  select(-d) %>% 
  mutate(trial = rep(seq(1, nrow(.)/10), 10),input=as.factor(input))

t %>% group_by(lr, update_func_name, simNum, trial, input) %>%
  summarize(almTrainDat = mean(almTrainDat), .groups = "drop") %>%
  ggplot(aes(x = trial, y = almTrainDat, color = input)) +
  geom_line() + ylim(c(0,1300))+
  facet_grid(lr ~ update_func_name)
#
#
#
#
#
#
#| eval: false 
#| 
parmVec <- tibble(crossing(c = c(0.1), lr = c(0.1,0.4,1), 
                           noise = c(10), trainRep = c(20), lossFun = list("MAE"), simNum = 1))

t <- parmVec %>% group_by(simNum,c,lr) %>% mutate(td = list(gen_train(trainRep = first(trainRep), noise = first(noise)))) %>% ungroup() %>%
  mutate(d = pmap(list(td, c, lr), ~sim_data(dat = .x, c = ..2, lr = ..3)),
         almTrainDat = map(d, "almTrain"))

# extract and plot each almTrainDat 
almTrainDat <- parmVec %>% group_by(simNum,c,lr) %>% mutate(td = list(gen_train(trainRep = first(trainRep), noise = first(noise)))) %>% ungroup() %>%
  mutate(d = pmap(list(td, c, lr), ~sim_data(dat = .x, c = ..2, lr = ..3)),
         almTrainDat = map(d, "almTrain")) %>% unnest(almTrainDat) %>% select(-d)

# For each unique value of lr, plot the learning curve, showing almTrainDat as a function of trial number, color by input, and facet by lr. Convert input to factor first. 

almTrainDat %>% group_by(lr) %>% mutate(trial = seq(1,n()),input= as.factor(input)) %>% ggplot(aes(x=trial,y=almTrainDat,color=input)) + geom_line() + facet_wrap(~lr)



parmVec <- tibble(crossing(c = c(0.1), lr = c(.01,.05,0.1), noise = c(5), trainRep = c(20), lossFun = list("MAE"), simNum = 1:30))

# The rest of the code remains the same
almTrainDat <- parmVec %>% group_by(simNum,c,lr) %>% mutate(td = list(gen_train(trainRep = first(trainRep), noise = first(noise)))) %>% ungroup() %>%
  mutate(d = pmap(list(td, c, lr), ~sim_data(dat = .x, c = ..2, lr = ..3)),
         almTrainDat = map(d, "almTrain")) %>% unnest(almTrainDat,td) %>% select(-d) %>%
  mutate(trial = rep(seq(1, nrow(.)/30), 30),input=as.factor(input))


mean_sd_almTrainDat <- almTrainDat %>%
  group_by(lr,input, trial) %>%
  summarise(avg_almTrainDat = mean(almTrainDat), sd_almTrainDat = sd(almTrainDat), .groups = 'drop')

# Check the results
head(mean_sd_almTrainDat)


ggplot(mean_sd_almTrainDat, aes(x = trial, y = avg_almTrainDat,color=input)) +
   geom_line() +
  geom_errorbar(aes(ymin = avg_almTrainDat - sd_almTrainDat, ymax = avg_almTrainDat + sd_almTrainDat), width = 0.2) +
  facet_wrap(~lr) +
  labs(title = "ALM Train Data: Mean and Variance Over Trials",
       x = "Trial",
       y = "Average ALM Train Data") + ylim(c(0,1300))+
  theme_minimal()



#
#
#
#
#
#
#
#

#
#
#
#

#
#
#
#plot(seq(1,90), (200-.50)*exp(-.1*seq(1,90))+.50)
gt <- gen_train(trainVec=c(5,6,7),trainRep=18) %>% mutate(cor=vx,err=(400-0)*exp(-.1*seq(1,n()))+0,vx=cor-err)
gt <- gen_train(trainVec=c(5,6,7),trainRep=28,noise=10) %>% mutate(cor=vx,err=(100-0)*exp(-.03*seq(1,n()))+0,en=map_dbl(err,~rnorm(n=1,mean=0,sd=.x)),vx1=cor-err,vx2=cor-en)

gt <- gen_train(trainVec=c(5,6,7),trainRep=28,noise=0) %>% group_by(input) %>% mutate(cor=vx,err=(200--10)*exp(-.01*seq(1,n()))+(-10),en=map(err,~rnorm(n=1,mean=0,sd=.x)),vx=cor-err)

gt <- gen_train(trainVec=c(5,6,7),trainRep=28,noise=10) %>% mutate(cor=vx,
            err=ifelse(seq(1, n()) <= n()/2, (700-0)*exp(-.01*seq(1,n()))+0, (700-0)*exp(-.06*(seq(1,n())-n()/2))+0),
            en=map(err,~rnorm(n=1,mean=0,sd=.x)),
            vx=cor-err)

gt <- gen_train(trainVec=c(5,6,7),trainRep=28,noise=10) %>%
  mutate(
    cor = vx,
    err = (700 - 0) * exp(-0.03 * seq(1, n()) * (input / max(input))) + 0,
    en = map(err, ~rnorm(n = 1, mean = 0, sd = .x)),
    vx = cor - err
  )
gt <- gen_train(trainVec=c(5,6,7),trainRep=28,noise=10) %>% mutate(
  cor = vx,
  err = (700 - 1) * exp(-.03 * seq(1, n())) + 1,
  en = map(err, ~rnorm(n = 1, mean = 0, sd = .x)),
  weight = (seq(1, n()) - 1) / (n() - 1),
  vx = cor * weight - err + abs(min(cor * weight - err)) + 1
)

#
#
#
#
#
#
