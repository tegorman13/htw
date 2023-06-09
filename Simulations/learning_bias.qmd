---
title: learning bias
date: last-modified
categories: [Simulation, ALM, R]
---



```{r}
pacman::p_load(tidyverse, data.table, patchwork, glue, knitr, kableExtra,here)
purrr::walk(here(c("Functions/alm_functions.R", "Functions/Display_Functions.R")), source)
```

```{r}


mean.prediction <- function(x.target, weights, bias, c, noise_sd) {
    probability <- output.activation(x.target, weights, c, noise_sd) / sum(output.activation(x.target, weights, c, noise_sd))
    return(outputNodes %*% probability + bias) # integer prediction with bias
}
input.activation <- function(x.target, c, noise_sd) {
  noise <- rnorm(length(inputNodes), mean = 0, sd = noise_sd)
  return(exp((-1 * c) * (x.target - inputNodes + noise)^2))
}
output.activation <- function(x.target, weights, c, noise_sd) {
    noise <- rnorm(length(outputNodes), mean = 0, sd = noise_sd)
    return(weights %*% input.activation(x.target, c, noise_sd) + noise)
}

update.weights <- function(x.new, y.new, weights, bias, c, lr, noise_sd) {
  y.feedback.activation <- exp(-1 * c * (y.new - outputNodes)^2)
  x.feedback.activation <- output.activation(x.new, weights, c, noise_sd)
  # Update weights
  new_weights <- weights + lr * (y.feedback.activation - x.feedback.activation) %*% t(input.activation(x.new, c, noise_sd))
  # Update bias
  new_bias <- bias + blr * (y.new - mean.prediction(x.new, weights, bias, c, noise_sd))
  return(list(weights = new_weights, bias = new_bias))
}

train.alm <- function(dat, c = 0.05, lr = 0.5, noise_sd = 0, weights, bias = 0) {
    alm.train <- rep(NA, nrow(dat))
    for (i in 1:nrow(dat)) {
        result <- update.weights(dat$input[i], dat$vx[i], weights, bias, c, lr, noise_sd)
        weights <- result$weights
        bias <- result$bias
        resp <- mean.prediction(dat$input[i], weights, bias, c, noise_sd)
        alm.train[i] <- resp
        weights[weights < 0] <- 0
    }
    list(almTrain = alm.train, weights = weights, bias = bias)
}

sim_train <- function(dat, c = 0.5, lr = 0.2, inNodes = 7, outNodes = 32, noise_sd = 0,bias=0) {
    inputNodes <<- seq(1, 7, length.out = inNodes * 1)
    outputNodes <<- seq(50, 1600, length.out = outNodes * 1)
    wm <- matrix(.000001, nrow = length(outputNodes), ncol = length(inputNodes))
    tt <- train.alm(dat, c, lr, noise_sd, wm, bias)
    return(tt)
}

gen_train <- function(trainVec = c(5, 6, 7), trainRep = 3, noise = 0) {
    bandVec <- c(0, 100, 350, 600, 800, 1000, 1200)
    if (class(trainVec) == "list") {
        trainVec <- unlist(trainVec)
    }
    ts <- rep(seq(1, length(trainVec)), trainRep)
    noiseVec <- rnorm(length(ts), mean = 0) * noise
    if (noise == 0) {
        noiseVec <- noiseVec * 0
    }
    tibble(trial = seq(1, length(ts)), input = trainVec[ts], vx = bandVec[trainVec[ts]] + noiseVec)
}


```


```{r fig.width=13, fig.height=12}
#| eval: FALSE

run_simulations <- function(param_grid, data) {
    simulations <- list()

    for (i in 1:nrow(param_grid)) {
        c_value <- param_grid$c[i]
        lr_value <- param_grid$lr[i]
        noise_sd_value <- param_grid$noise_sd[i]
        inNodes_value <- param_grid$inNodes[i]
        outNodes_value <- param_grid$outNodes[i]
        bias_value <- param_grid$bias[i]

        sim_result <- sim_train(
            dat = data,
            c = c_value,
            lr = lr_value,
            inNodes = inNodes_value,
            outNodes = outNodes_value,
            noise_sd = noise_sd_value,
            bias = bias_value
        )

        simulations[[i]] <- cbind(data,alm_train=sim_result$almTrain,biasEnd=sim_result$bias)
    }

    return(simulations)
}

plot_learning_curves <- function(simulations,param_grid) {
    curves <- data.frame(simulations[[1]], Simulation=1,c=param_grid$c[1],lr=param_grid$lr[1],noise_sd=param_grid$noise_sd[1],inNodes=param_grid$inNodes[1],outNodes=param_grid$outNodes[1],biasStart=param_grid$bias[1])

    for (i in 2:length(simulations)) {
        curve_data <- data.frame(simulations[[i]],Simulation = i,
        c=param_grid$c[i],lr=param_grid$lr[i],noise_sd=param_grid$noise_sd[i],inNodes=param_grid$inNodes[i],outNodes=param_grid$outNodes[i],biasStart=param_grid$bias[i])
        curves <- rbind(curves, curve_data)
    }
    curves <- curves %>% mutate(simLab=paste0("c=",c," lr=",lr," n_sd=",noise_sd,
                                
                                              " biasEnd=",round(biasEnd,1)," biasStart=",biasStart),
                                input=as.factor(input))
    

    p <- ggplot(curves, aes(x = trial, y = alm_train, color=input)) +
        geom_line(alpha = 0.8) + facet_wrap(~simLab,ncol=3)+
        theme_minimal() + ylim(c(0,1600))
        labs(title = "Learning Curves for Simulated Learners", x = "Trial", y = "Velocity")

    return(p)
}


param_grid <- tibble(crossing(
    c = c(1.5, 2, 3),
    lr = c(.01,.5),
    noise_sd = c(.0010),
    inNodes = c( 7),
    outNodes = c(32),
    bias=c(-800,100,500)
))

blr<<-0
simulations <- run_simulations(param_grid[1:35,], gen_train(trainRep=10))
learning_curve_plot <- plot_learning_curves(simulations, param_grid)
print(learning_curve_plot)
blr<<-.1
simulations <- run_simulations(param_grid[1:35,], gen_train(trainRep=10))
learning_curve_plot <- plot_learning_curves(simulations, param_grid)
print(learning_curve_plot)
blr<<-.01
simulations <- run_simulations(param_grid[1:35,], gen_train(trainRep=10))
learning_curve_plot <- plot_learning_curves(simulations, param_grid)
print(learning_curve_plot)


```

