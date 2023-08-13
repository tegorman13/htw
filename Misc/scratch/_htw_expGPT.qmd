Consider this study description, then propose a statistical analysis approach that can be conducted with R. The main interest of the study is on the difference between the varied and constant training conditions, but the effect of catOrder and feedbacktype will also need to be considered. \# Introduction In project 1, we applied model-based techniques to quantify and control for the similarity between training and testing experience, which in turn enabled us to account for the difference between varied and constant training via an extended version of a similarity based generalization model. In project 2, we will go a step further, implementing a full process model capable of both 1) producing novel responses and 2) modeling behavior in both the learning and testing stages of the experiment. Project 2 also places a greater emphasis on extrapolation performance following training - as varied training has often been purported to be particularly beneficial in such situations. Extrapolation has long been a focus of the literature on function learning [@brehmerHypothesesRelationsScaled1974; @carrollFunctionalLearningLearning1963]. Central questions of the function learning literature have included the relative difficulties of learning various functional forms (e.g. linear vs.bilinear vs. quadratic), and the relative effectiveness of rule-based vs. association-based exemplar models vs. various hybrid models [@bottNonmonotonicExtrapolationFunction2004; @deloshExtrapolationSineQua1997; @jonesActiveFunctionLearning2018; @kalishPopulationLinearExperts2004; @mcdanielConceptualBasisFunction2005; @mcdanielPredictingTransferPerformance2009]. However the issue of training variation has received surprisingly little attention in this area. \# Methods \## Participants Data was collected from 647 participants (after exclusions). The results shown below consider data from subjects in our initial experiment, which consisted of 196 participants (106 constant, 90 varied). The follow-up experiments entailed minor manipulations: 1) reversing the velocity bands that were trained on vs. novel during testing; 2) providing ordinal rather than numerical feedback during training (e.g. correct, too low, too high). The data from these subsequent experiments are largely consistently with our initial results shown below. \## Task We developed a novel visuomotor extrapolation task, termed the ["Hit The Wall" (HTW]{.underline}) task, wherein participants learned to launch a projectile such that it hit a rectangle at the far end of the screen with an appropriate amount of force (specified by the velocity band values). Although the projectile had both x and y velocity components, only the x-dimension was relevant for the task. On each trial, we record the vx and vy produced by the participant. We calculate the distance between the vx value and the closest edge of the current velocity band. Note that lower distance values correspond to better performance. \## Design (expMode in dataframe) \### training stage 1) train - 90 training trials split evenly divided between velocity bands. Varied training with 3 velocity bands and Constant training with 1 band. 2) train-Nf - interleaved no feedback testing that happens throughout training \### Testing Stage 1) No-feedback testing from 3 novel extrapolation bands. 15 trials each. (test-Nf) 2) No-feedbacd testing from the 3 bands used during the training phase (2 of which were novel for the constant group). 9 trials each. (test-train-Nf) 3) Feedback testing for each of the 3 extrapolation bands. 10 trials each. (test-feedback) The is also a between subjects order manipulation (catOrder). In the orig condition, constant subjects train on band 800-1000, and varied subjects train on bands 800-1000, 1000-1200, and 1200-1400. In the rev condition, constant subjects train on band 600-800, and varied subjects train on bands 600-800, 350-550, and 100-300. There is also the feedbackType between subjects manipulation. Subjects in the continuous condition receive condition feedback, and subjects in the ordinal condition receive ordinal feedback. d \<- readRDS('dPrune-01-19-23.rds') colnames(d) \[1\] "id" "sbjCode""condit""trial""nGoodTrial""goodThrow""gt.bandStage""expMode" "gt.stage"\[10\] "gt.train" "trainStage" "expStage""band" "input" "vb" "bandInt" "lowBound" "highBound" "feedback" "fb""runTotal" "mode" "stage" "catOrder" "dist" "vx" "vxb" "vxi" "vy" "nTrain" "nTestNf""nInt" "nTestF" "nTotal" "lastTrain" "lastTrial""trainVec" "feedbackType" "fullCond"\
dtest \<- d %\>% filter(expMode %in% c("test-Nf","test-train-nf")) %\>% group_by(id,lowBound) %\>% mutate(nBand=n(),band=bandInt,id=factor(id)) %\>% group_by(id) %\>% mutate(nd=n_distinct(lowBound)) dtest \<- dtest %\>% group_by(id,lowBound) %\>% filter(nBand\>=5 & nd==6) dtest \<- dtest %\>% group_by(id) %\>% filter(!id %in% unique(dtest$id[dtest$nBand\<5\])) dtestAgg \<- dtest %\>% group_by(id,condit,catOrder,feedbackType,vb,band,lowBound,highBound,input) %\>% mutate(vxCapped=ifelse(vx\>1600,1600,vx)) %\>% summarise(vxMean=mean(vx),devMean=mean(dist),vxMed=median(vx),devMed=median(dist), vxMeanCap=mean(vxCapped),.groups = "keep") ds \<- d %\>% filter(expMode %in% c("train","train-Nf","test-Nf","test-train-nf")) %\>% filter(!id %in% unique(dtest$id[dtest$nBand\<5\])) %\>% select(id,condit,catOrder,feedbackType,expMode,trial,gt.train,vb,band,bandInt,lowBound,highBound,input,vx,dist,vxb) head(ds) \# A tibble: 6 × 16 id condit catOrder feedbackType expMode trial gt.train vb band bandInt lowBound highBound input vx dist vxb <dbl> <fct> <fct> <fct> <fct> <dbl> <int> <fct> <fct> <dbl> <fct> <dbl> <dbl> <dbl> <dbl> <dbl> 1 1 Varied orig continuous train 2 1 1000-1200 5 1000 1000 1200 6 2614. 1414. 1600 2 1 Varied orig continuous train 3 2 1200-1400 6 1200 1200 1400 7 782. 418. 800 3 1 Varied orig continuous train 4 3 800-1000 4 800 800 1000 5 1214. 214. 1200

# Simulation GPT Prompt

Assess whether the code below correctly implements the ALM model \## ALM (Associative Learning Model) The basis of the model is an associative neural network that uses a delta learning rule to associate inputs to outputs. Responses to stimuli within the training domain are obtained by simple activation of the network. The underlying associative learning model is a simple two-layer connectionist network that updates weight strengths using the delta learning rule. The $M$ nodes along the input of the model represent the possible values along some input domain. The $N$ output nodes represent the possible responses along the output range. A real stimulus/response will be mapped to the nodes via some psychophysical function $\psi(x)$. The equations presented will ignore this transformation to simplify things. When an input stimulus $X$ is presented, it will cause a Gaussian activation of input nodes according to the function $a_i(X)=e^{-\gamma \cdot\left(X-X_i\right)^2}$, where $\gamma$ is a scaling parameter describing the steepness of the gradient. The input activation is then normalized so that the area under the curve is equal to 1 by dividing all values by the sum of all input activation values. The activation of each output node is calculated by summing the products of the input nodes and the weights that connect them to a particular output node. This is shown in the formula $O_j(X)=\sum_{i=1}^M w_{j i} \cdot a_i(X)$, where $w_{j i}$ designates the strength of association between input node $X_i$ and output node $Y_j$. At this point we have a distribution that should describe the responses in a stochastic simulation. In order to obtain deterministic predictions, the expected value of the output activation is calculated. First the probability of choosing a particular output given a specific input is calculated according to $P\left[Y_j \mid X\right]=\frac{O_j(X)}{\sum_{k=1}^{\mathrm{L}} O_k(X)}$. The mean output given stimulus $X$ is the weighted average $m(X)=\sum_{j=1}^L Y_j \cdot P\left[Y_j \mid X\right]$. In addition to creating associations, the model allows for learning. This step adjusts weight connections to improve accuracy. The first step is obtaining an error signal. When participants perform the task, they are shown a feedback signal $Z$ that demonstrates what the proper response should have been. This signal activates the output nodes using another Gaussian similarity function $f_j(Z)=e^{-\gamma \cdot\left(Z-Y_j\right)^2}$. This feedback signal is used to update the weights according to the delta learning rule: $w_{j i}(t+1)=w_{j i}(t)+\alpha \cdot\left\{f_j[Z(t)]-O_j[X(t)]\right\} \cdot a_i[X(t)]$

```{r fig.width=12,fig.height=10}
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
    c = c(.1,.5, 2), lr = c(.05, 1), noise = c(0),
    inNodes = c(7), outNodes = c(32),
    trainVec = list(list(5, 6, 7)), trainRep = c(8),
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
    mutate(input = as.factor(input)) %T>%
  {pf(.) } %>% trainTab


```

### HTW design + ABC prompt

## Participants

Data was collected from 647 participants (after exclusions). The results shown below consider data from subjects in our initial experiment, which consisted of 196 participants (106 constant, 90 varied). The follow-up experiments entailed minor manipulations: 1) reversing the velocity bands that were trained on vs. novel during testing; 2) providing ordinal rather than numerical feedback during training (e.g. correct, too low, too high). The data from these subsequent experiments are largely consistently with our initial results shown below. \## Task We developed a novel visuomotor extrapolation task, termed the ["Hit The Wall" (HTW]{.underline}) task, wherein participants learned to launch a projectile such that it hit a rectangle at the far end of the screen with an appropriate amount of force (specified by the velocity band values). Although the projectile had both x and y velocity components, only the x-dimension was relevant for the task. On each trial, we record the vx and vy produced by the participant. We calculate the distance between the vx value and the closest edge of the current velocity band. Note that lower distance values correspond to better performance. \## Design (expMode in dataframe) \### training stage 1) train - 90 training trials split evenly divided between velocity bands. Varied training with 3 velocity bands and Constant training with 1 band. 2) train-Nf - interleaved no feedback testing that happens throughout training \### Testing Stage 1) No-feedback testing from 3 novel extrapolation bands. 15 trials each. (test-Nf) 2) No-feedbacd testing from the 3 bands used during the training phase (2 of which were novel for the constant group). 9 trials each. (test-train-Nf) 3) Feedback testing for each of the 3 extrapolation bands. 10 trials each. (test-feedback) The is also a between subjects order manipulation (catOrder). In the orig condition, constant subjects train on band 800-1000, and varied subjects train on bands 800-1000, 1000-1200, and 1200-1400. In the rev condition, constant subjects train on band 600-800, and varied subjects train on bands 600-800, 350-550, and 100-300. There is also the feedbackType between subjects manipulation. Subjects in the continuous condition receive condition feedback, and subjects in the ordinal condition receive ordinal feedback.

```{r}
# load and view data
pacman::p_load(tidyverse,data.table,lme4)
d <- readRDS("dPrune-01-19-23.rds")
levels(d$condit)
# Prepare the data for analysis
ds <- d %>% filter(expMode %in% c("train","train-Nf","test-Nf","test-train-nf")) %>% 
filter(!id %in% unique(dtest$id[dtest$nBand<5])) %>% 
select(id,condit,catOrder,feedbackType,expMode,trial,gt.train,vb,band,bandInt,lowBound,highBound,input,vx,dist,vxb) 
dst <- ds %>% filter(expMode=="train")
dst <- dst %>%
  group_by(id, vb) %>%
  mutate(trial_band = row_number())
head(dst)
```

| id  | condit | catOrder | feedbackType | expMode | trial | gt.train | vb        | band | bandInt | lowBound | highBound | input | vx   | dist | vxb | trial_band |
|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|
| 1   | Varied | orig     | continuous   | train   | 2     | 1        | 1000-1200 | 5    | 1000    | 1200     | 6         | 2614  | 1414 | 1600 | 1   |            |
| 2   | Varied | orig     | continuous   | train   | 3     | 2        | 1200-1400 | 6    | 1200    | 1400     | 7         | 782   | 418  | 800  | 1   |            |

```{r}
input.activation <- function(x.target, c, noise_sd) {
  noise <- rnorm(length(inputNodes), mean = 0, sd = noise_sd)
  return(exp((-1 * c) * (x.target - inputNodes + noise)^2))
}
output.activation <- function(x.target, weights, c, noise_sd) {
  noise <- rnorm(length(outputNodes), mean = 0, sd = noise_sd)
  return(weights %*% input.activation(x.target, c, noise_sd) + noise)
}
mean.prediction <- function(x.target, weights, c, noise_sd) {
  probability <- output.activation(x.target, weights, c, noise_sd) / sum(output.activation(x.target, weights, c, noise_sd))
  return(outputNodes %*% probability) # integer prediction
}
update.weights <- function(x.new, y.new, weights, c, lr, noise_sd) {
  y.feedback.activation <- exp(-1 * c * (y.new - outputNodes)^2)
  x.feedback.activation <- output.activation(x.new, weights, c, noise_sd)
  return(weights + lr * (y.feedback.activation - x.feedback.activation) %*% t(input.activation(x.new, c, noise_sd)))
}
train.alm <- function(dat, c = 0.05, lr = 0.5, noise_sd = 0, weights) {
  alm.train <- rep(NA, nrow(dat))
  for (i in 1:nrow(dat)) {
    weights <- update.weights(dat$input[i], dat$vx[i], weights, c, lr, noise_sd)
    resp = mean.prediction(dat$input[i], weights, c, noise_sd)
    alm.train[i] = resp
    weights[weights < 0] = 0
  }
  list(almTrain = alm.train, weights = weights)
}
sim_train <- function(dat, c = 0.5, lr = 0.2, inNodes = 7, outNodes = 32, noise_sd = 0) {
  inputNodes <<- seq(1, 7, length.out = inNodes * 1)
  outputNodes <<- seq(50, 1600, length.out = outNodes * 1)
  wm = matrix(.0000, nrow = length(outputNodes), ncol = length(inputNodes))
  tt <- train.alm(dat, c, lr, noise_sd, wm)
}
MAE <- function(x, y) {mean(abs(x - y))}
RMSE <- function(x,y){sqrt(mean((x-y)^2))}

library(abc)
library(mvtnorm)
fit_alm_sim <- function(data, c, lr, noise_sd, inNodes, outNodes) {
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
  return(train_data$almTrain)
}
generate_prior <- function(n) {
  prior_samples <- tibble(
    c = runif(n, 0.001, 5),
    lr = runif(n, 0.001, 3),
    noise_sd = runif(n, 0, 0.0001),
    inNodes = sample(c( 3,7, 14,21,28,35), size = n, replace = TRUE),
    outNodes = sample(c(8,16, 32,48,64,80,96), size = n, replace = TRUE)
  )
  return(prior_samples)
}
n_prior_samples <- 10000
prior_samples <- generate_prior(n_prior_samples)

gt <- gen_train(trainVec=c(5,6,7),trainRep=36) %>% mutate(cor=vx,err=(1500-0)*exp(-.05*seq(1,n()))+0,vx=cor-err)

simulated_data <- future_map_dfc(seq_len(nrow(prior_samples)), function(idx) {
  params <- prior_samples[idx, ]
  fit_alm_sim(gt, params$c, params$lr, params$noise_sd, params$inNodes, params$outNodes)
},.options = furrr_options(seed = T))

gt_obs <- gt$vx
tolerance <- 0.1 * sd(gt_obs)*500
abc_result <- abc(
  target = gt_obs,
  param = prior_samples,
  sumstat = do.call(rbind, simulated_data),
  tol = .1,
  method = "rejection",
  names=colnames(gt_obs)
)

  


```