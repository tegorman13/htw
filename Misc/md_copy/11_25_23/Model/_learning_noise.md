---
title: Noisy Learning
date: last-modified
date-format: "MM/D/YYYY, h:mm A"
categories: [Simulation, ALM, R, Bayesian, Model-Structure]
page-layout: full
execute: 
  warning: false
  eval: false 
---



```{r}
pacman::p_load(tidyverse,data.table,knitr,kableExtra,glue,future,furrr,here)
purrr::walk(here(c("Functions/alm_functions.R","Functions/Display_Functions.R","Functions/Noisy_Functions.R")),source)
```







It appears that the model may be too simple or not flexible enough to capture the slow, gradual learning observed in the human data. Here are some suggestions to improve the model's ability to better fit the human learning data:

Introduce a momentum term: Adding a momentum term to the weight update rule can help the model exhibit more gradual learning behavior. The momentum term is a fraction of the previous weight update that is added to the current weight update, which can help to smooth out the learning process.

Vary learning rate over time: Instead of using a fixed learning rate, you can try to gradually decrease the learning rate over time. This approach, called learning rate annealing, can help the model learn more slowly at the beginning of training and become more fine-tuned as training progresses.

Incorporate noise: You can add noise to the input activation, output activation, or weight updates to introduce some stochasticity into the learning process. This can lead to more gradual learning, as the model will not be able to rely solely on deterministic updates.

Hierarchical or recurrent structure: Instead of using a simple two-layer connectionist network, you can try a more complex network architecture, such as a hierarchical or recurrent neural network. This can potentially help capture the slow, gradual learning behavior observed in human data.

Regularization: Add regularization techniques like L1 or L2 regularization to the learning process. This can help prevent overfitting and create a more constrained model, which may lead to slower and more gradual learning.

Task-specific features: Incorporate task-specific features, constraints, or biases into the model that may better represent human cognitive processes. This may help the model to better capture the observed learning patterns in human data.

Model comparison: Compare the performance of your current model with alternative models, such as reinforcement learning or Bayesian models, which may exhibit different learning dynamics.

By trying out these suggestions, you can potentially improve your model's ability to fit the slow, gradual learning observed in human data from the velocity production task. Experiment with different combinations of these techniques to find the best approach for your specific problem.





## Noisy layers version

```{r}

```


```{r fig.width=12,fig.height=9}

tibble(crossing(
    c = c(.5, 2), lr = c(0.005, 0.05, 0.1, 0.5), noise_sd = c(0.001, 0.005, 0.01, 0.05),
    noise = c(0),
    inNodes = c(7), outNodes = c(32),
    trainVec = list(list(5, 6, 7)), trainRep = c(16),
    lossFun = list("MAE"),
    simNum = 1:1
)) %>%
    mutate(id = seq(1, nrow(.)), td = pmap(list(trainVec, trainRep, noise), ~ gen_train(trainVec = .x, trainRep = ..2, noise = ..3))) %>%
    ungroup() %>%
    mutate(
        d = pmap(
            list(td, c, lr, inNodes, outNodes, noise_sd),
            ~ sim_train(dat = .x, c = ..2, lr = ..3, inNodes = ..4, outNodes = ..5, noise_sd = ..6)
        ),
        almTrainDat = map(d, "almTrain"), weights = map(d, "weights")
    ) %>%
    unnest(c(almTrainDat, td)) %>%
    select(-d) %>%
    mutate(input = as.factor(input)) %T>%
  {pf(.) } %>% trainTab %>% {. ->> tt}

```




## Grid Search
```{r}
n_cores <- parallel::detectCores()
param_grid <- tibble(crossing(
  c = seq(.5,5,.5),
  lr = seq(0.01, 1,.1),
  noise_sd = c(0,.0001,0.001, 0.01),
  inNodes = c(5, 7,14,28),
  outNodes = c(16, 32,64)
))
nrow(param_grid)

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

param_grid <- param_grid %>% mutate(performance = future_map_dbl(seq_len(nrow(.)), function(idx) {
    fit_alm(gt, c = c[idx], lr = lr[idx], noise_sd = noise_sd[idx], inNodes = inNodes[idx], outNodes = outNodes[idx])
  },
  .options = furrr_options(seed = T)))

best_params <- param_grid %>%
  arrange((performance)) 
best <- head(best_params,1)

s=sim_train(dat=mutate(gt,vx=cor), c = best$c,
  lr = best$lr,inNodes = best$inNodes,outNodes = best$outNodes,
  noise_sd = best$noise_sd
)

ggp <- gt %>% mutate(pred=s %>% pluck("almTrain"),input=as.factor(input))
ggp %>% ggplot(aes(x = trial, y = pred, color = input)) +
  geom_line() + ylim(c(0,1600))
ggp  %>% ggplot(aes(x = trial, y = vx, color = input)) +
  geom_line() + ylim(c(0,1600))

```





```r
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

#gt <- gen_train(trainVec=c(5,6,7),trainRep=8) %>% mutate(cor=vx,err=(800-0)*exp(-.1*seq(1,n()))+0,vx=cor-err)
gt <- gen_train(trainVec=c(5,6,7),trainRep=36) %>% mutate(cor=vx,err=(1500-0)*exp(-.05*seq(1,n()))+0,vx=cor-err)



simulated_data <- future_map_dfc(seq_len(nrow(prior_samples)), function(idx) {
  params <- prior_samples[idx, ]
  fit_alm_sim(gt, params$c, params$lr, params$noise_sd, params$inNodes, params$outNodes)
},.options = furrr_options(seed = T))

gt_obs <- gt$vx
tolerance <- 0.1 * sd(gt_obs)*500

length(gt_obs)
str(simulated_data)

abc_result <- abc(
  target = gt_obs,
  param = prior_samples,
  sumstat = do.call(rbind, simulated_data),
  tol = .1,
  method = "rejection",
  names=colnames(gt_obs)
)
str(abc_result)





posterior_samples <- abc_result$unadj.values
colnames(posterior_samples) <- c("c", "lr", "noise_sd", "inNodes", "outNodes")
posterior_samples_long <- tidyr::pivot_longer(as.data.frame(posterior_samples), everything())

ggplot(posterior_samples_long, aes(x=value)) +
  geom_density() +
  facet_wrap(~name, scales="free") +
  theme_minimal() +
  labs(x="Value", y="Density", title="Posterior Density Plots")


summary_statistics <- data.frame(
  mean = apply(posterior_samples, 2, mean),
  median = apply(posterior_samples, 2, median),
  q025 = apply(posterior_samples, 2, function(x) quantile(x, 0.025)),
  q975 = apply(posterior_samples, 2, function(x) quantile(x, 0.975))
)
summary_statistics=rownames_to_column(summary_statistics,var="parameter")
summary_statistics
head(summary_statistics)


s = sim_train(dat=mutate(gt,vx=cor),c=summary_statistics$mean[1],lr=summary_statistics$mean[2],inNodes=summary_statistics$mean[4],outNodes=summary_statistics$mean[5],noise_sd=summary_statistics$mean[3])

ggp <- gt %>% mutate(pred=s %>% pluck("almTrain"),input=as.factor(input))
ggp %>% ggplot(aes(x = trial, y = pred, color = input)) +
  geom_line() + ylim(c(0,1600))
ggp  %>% ggplot(aes(x = trial, y = vx, color = input)) +
  geom_line() + ylim(c(0,1600))


simulated_data_accepted <- abc_result$ss
target_data <- as.vector(gt_obs)

# Plot target data and accepted simulations
ggplot() +
  geom_density(data = as.data.frame(simulated_data_accepted), aes(x = V1), alpha = 0.3, color = "blue") + geom_vline(aes(xintercept = target_data), color = "red", linetype = "dashed") +
  labs(x = "Summary statistics", y = "Density", title = "Diagnostic Plot: Marginal Posterior Distribution of Summary Statistics") +
  theme_minimal()





posterior_predictive_simulations <- future_map_dfc(seq_len(nrow(posterior_samples)), function(idx) {
  params <- prior_samples[idx, ]
  fit_alm_sim(gt, params$c, params$lr, params$noise_sd, params$inNodes, params$outNodes)
},.options = furrr_options(seed = T))



# Calculate summary statistics for the posterior predictive simulations
posterior_predictive_summary <- apply(posterior_predictive_simulations, 2, function(x) c(mean = mean(x), sd = sd(x)))

# Compare with the observed data summary statistics
observed_data_summary <- as.data.frame(t(c(mean = mean(gt_obs), sd = sd(gt_obs))))
colnames(observed_data_summary) <- colnames(posterior_predictive_summary)

comparison <- rbind(posterior_predictive_summary, observed_data_summary)
rownames(comparison)[3] <- "Observed"
comparison




library(GGally)
ggpairs(as.data.frame(posterior_samples))



```








#ABC Binned
```{r}
nTrain=90
gt_con <- gen_train(trainVec=c(5,6,7),trainRep=nTrain/3) %>% mutate(cor=vx,err=(500-0)*exp(-.05*seq(1,n()))+0,vx=cor-err)
gt_con <- gen_train(trainVec=c(5),trainRep=nTrain) %>% mutate(cor=vx,err=(500-0)*exp(-.09*seq(1,n()))+0,vx=cor-err)

gt <- gt_con

bin_size <- 8
binned_data <- gt %>%
  mutate(bin=cut(trial,breaks=bin_size,labels=c(1:bin_size))) %>%
  group_by(bin,cor,input) %>%
  summarize(vx_mean = mean(vx), .groups = "drop")

fit_alm_sim_binned <- function(data, c, lr, noise_sd, inNodes, outNodes) {
  train_sims <- replicate(5, {
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
    
   binned_data <- train_data %>%
  mutate(bin=cut(trial,breaks=bin_size,labels=c(1:bin_size))) %>%
  group_by(bin,vx,input) %>%
  summarize(vx_mean = mean(almTrain), .groups = "drop")
  binned_data$vx_mean
  })
  train_avg <- rowMeans(train_sims )
  return( train_avg)
}

n_prior_samples <- 500
prior_samples <- generate_prior(n_prior_samples)

simulated_data <- future_map_dfc(seq_len(nrow(prior_samples)), function(idx) {
  params <- prior_samples[idx, ]
  fit_alm_sim_binned(gt, params$c, params$lr, params$noise_sd, params$inNodes, params$outNodes)
},.options = furrr_options(seed = T))

gt_obs <- binned_data$vx_mean

abc_result <- abc(
  target = gt_obs,param = prior_samples,
  sumstat = do.call(rbind, simulated_data),
  tol = .1, method = "rejection",
  names=colnames(gt_obs)
)

posterior_samples <- abc_result$unadj.values
colnames(posterior_samples) <- c("c", "lr", "noise_sd", "inNodes", "outNodes")
posterior_samples_long <- tidyr::pivot_longer(as.data.frame(posterior_samples), everything())

ggplot(posterior_samples_long, aes(x=value)) +
  geom_density() +
  facet_wrap(~name, scales="free") +
  theme_minimal() +
  labs(x="Value", y="Density", title="Posterior Density Plots")

summary_statistics <- data.frame(
  mean = apply(posterior_samples, 2, mean),
  median = apply(posterior_samples, 2, median),
  q025 = apply(posterior_samples, 2, function(x) quantile(x, 0.025)),
  q975 = apply(posterior_samples, 2, function(x) quantile(x, 0.975))
)
( summary_statistics=rownames_to_column(summary_statistics,var="parameter") )

s = sim_train(dat=mutate(gt,vx=cor),c=summary_statistics$mean[1],lr=summary_statistics$mean[2],inNodes=summary_statistics$mean[4],outNodes=summary_statistics$mean[5],noise_sd=summary_statistics$mean[3])

ggp <- gt %>% mutate(pred=s %>% pluck("almTrain"),input=as.factor(input)) %>%
  mutate(bin=cut(trial,breaks=bin_size,labels=c(1:bin_size)),input=as.factor(input),binN=as.numeric(bin)) %>%
  group_by(bin,cor,input,binN) %>% summarize(vx_mean = mean(vx), pred_mean=mean(pred), .groups = "drop")
ggp %>% ggplot(aes(x = binN, y = pred_mean,color=input)) +
  geom_line() + ylim(c(-200,1600))
ggp  %>% ggplot(aes(x = binN, y = vx_mean, color = input)) +
  geom_line() + ylim(c(-200,1600))

```




s=sim_data(dat=mutate(gt,vx=cor),c=best %>% pluck("c"),best= k %>% pluck("lr"))


ggp <- gt %>% mutate(pred=s %>% pluck("almTrain"),c=k %>% pluck("c"),lr= k %>% pluck("lr"),input=as.factor(input)) %>%
  ggplot(aes(x = trial, y = pred, color = input)) +
  geom_line() + ylim(c(0,1600))
ggo <-  gt %>% ggplot(aes(x = trial, y = vx, color = as.factor(input))) +
  geom_line() + ylim(c(-400,1600))

ggo+ggp

```


### Optimize for single decay curve
generate data that follows an exponetial decay function of error over trials, inspect
ability of model to fit that data. 
```{r}
gt <- gen_train(trainVec=c(5,6,7),trainRep=36) %>% mutate(cor=vx,err=(1500-0)*exp(-.05*seq(1,n()))+0,vx=cor-err)
gt_con <- gen_train(trainVec=c(5),trainRep=nTrain) %>% mutate(cor=vx,err=(500-0)*exp(-.09*seq(1,n()))+0,vx=cor-err)

gt_con %>% ggplot(aes(x = trial, y = vx, color = as.factor(input))) +
  geom_line() + ylim(c(-10,1600))
head(gt,10)
# bias <- 1000; 
# gt <- gen_train(trainVec=c(5,6,7),trainRep=228,noise=0) %>% mutate(
#   cor = vx,
#   err = (bias - 0) * exp(-.005 * seq(1, n())) + 0,
#   en = map2_dbl(err,cor, ~rnorm(n = 1, mean = .y, sd = .x/2)),
#   enAvg = map2_dbl(err,cor, ~mean(rnorm(n = 1, mean = .y, sd = .x))),
#   weight = (seq(1, n()) - 1) / (n() - 1),
#   vx = (weight*en)+bias*(1-weight),
#   vx=en
# )
gt %>% ggplot(aes(x = trial, y = vx, color = as.factor(input))) +
  geom_line() + ylim(c(-10,1600))


wrap_optim <- function(dat,lossFun=RMSE){
  
    wrap_alm <- function(par,dat, weights,lossFun){
      c=par[1]; lr=par[2]; noise_sd=par[3]
      pred=train.alm(dat, c=c, lr=lr, noise_sd=noise_sd, weights=weights)
      lossFun(dat$vx,pred)
    }
  if(class(lossFun)=="character"){lossFun=get(lossFun)}
  inputNodes = seq(1,7,1)  # 
  outputNodes = seq(50,1600,50)
  wm=matrix(.00001,nrow=length(outputNodes),ncol=length(inputNodes))
  testVec=seq(2,7)
  
  bounds_lower <- c(.0000001, .00001,.000000001)
  bounds_upper <- c(10, 10,10)
  parmsLab <- c("c","lr","noise_sd")
  
  fit=optim(par=c(.1, .2,.05),
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


k=wrap_optim(gt,lossFun = "MAE")
s=sim_data(dat=mutate(gt,vx=cor),c=k %>% pluck("c"),lr= k %>% pluck("lr"))
ggp <- gt %>% mutate(pred=s %>% pluck("almTrain"),c=k %>% pluck("c"),lr= k %>% pluck("lr"),input=as.factor(input)) %>%  
  ggplot(aes(x = trial, y = pred, color = input)) +
  geom_line() + ylim(c(0,1600))
ggo <-  gt %>% ggplot(aes(x = trial, y = vx, color = as.factor(input))) +
  geom_line() + ylim(c(-400,1600))

ggo+ggp
```


