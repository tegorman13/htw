---
title: Node Manipulations
date: last-modified
categories: [Simulation, ALM, R, Model-Structure]
page-layout: full
code-fold: true
code-tools: true
execute: 
  warning: false
  eval: false 
--- 



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


```


### Loss Functions
```{r}
RMSE <- function(x,y){
 # print("rmseTrial")
  sqrt(mean((x-y)^2))
}
MAE <- function(x, y) {
  mean(abs(x - y))
}
```




```{r}
# Test varying input and output nodes, c, lr parameters and noise
test_nodes_exhaustive <- function(dat, c_values = seq(0.1, 1, 0.1), lr_values = seq(0.1, 1, 0.1),
                                  input_nodes = seq(3, 10), output_nodes = seq(16, 64, 8),
                                  noise_values = seq(0, 1, 0.1), n_simulations = 10) {
  results <- list()

  for (c in c_values) {
    for (lr in lr_values) {
      for (inNodes in input_nodes) {
        for (outNodes in output_nodes) {
          for (noise in noise_values) {
            rmse_values <- c()

            for (i in 1:n_simulations) {
              dat_noise <- gen_train(noise = noise)
              result <- sim_data(dat_noise, c, lr, inNodes, outNodes)
              rmse <- RMSE(dat_noise$vx, result$almPred)
              rmse_values <- c(rmse_values, rmse)
            }

            results[[paste0("c", c, "_lr", lr, "_in", inNodes, "_out", outNodes, "_noise", noise)]] <-
              list(c = c, lr = lr, inNodes = inNodes, outNodes = outNodes, noise = noise,
                   mean_rmse = mean(rmse_values), var_rmse = var(rmse_values))
          }
        }
      }
    }
  }

  results
}

# Visualize the results with mean and variance
visualize_results_exhaustive <- function(results) {
  df <- do.call(rbind, lapply(results, function(x) data.frame(t(unlist(x)))))
  ggplot(df, aes(x = inNodes, y = outNodes, fill = mean_rmse)) +
    geom_tile() +
    facet_grid(c ~ lr) +
    scale_fill_gradient(low = "green", high = "red") +
    labs(title = "Mean RMSE for varying input and output nodes, c, lr, and noise",
         x = "Number of Input Nodes",
         y = "Number of Output Nodes",
         fill = "Mean RMSE") +
    theme_minimal()
}

# Generate test data
dat <- gen_train()

# Test and visualize results
results_exhaustive <- test_nodes_exhaustive(dat)
visualize_results_exhaustive(results_exhaustive)



```

```{r}
# Load necessary libraries
library(shiny)
library(tidyverse)

# Copy the required functions and data from your previous code here

# Shiny UI
ui <- fluidPage(
    titlePanel("ALM Model Performance"),
    sidebarLayout(
        sidebarPanel(
            sliderInput("c_value", "C value:", min = 0.1, max = 1, value = 0.5, step = 0.1),
            sliderInput("lr_value", "Learning rate:", min = 0.1, max = 1, value = 0.5, step = 0.1),
            sliderInput("inNodes", "Number of Input Nodes:", min = 3, max = 10, value = 7, step = 1),
            sliderInput("outNodes", "Number of Output Nodes:", min = 16, max = 64, value = 32, step = 8),
            sliderInput("noise", "Noise:", min = 0, max = 1, value = 0, step = 0.1),
            actionButton("run_model", "Run Model")
        ),
        mainPanel(
            plotOutput("model_performance")
        )
    )
)

# Shiny server
server <- function(input, output) {
    observeEvent(input$run_model, {
        dat_noise <- gen_train(noise = input$noise)
        result <- sim_data(dat_noise, input$c_value, input$lr_value, input$inNodes, input$outNodes)
        rmse <- RMSE(dat_noise$vx, result$almPred)

        combined_data <- dat_noise %>%
            mutate(
                almPred = result$almPred,
                examPred = result$examPred
            )

        output$model_performance <- renderPlot({
            ggplot(combined_data, aes(x = input, y = vx)) +
                geom_point(aes(y = almPred, color = "ALM Prediction")) +
                geom_line(aes(y = examPred, color = "EXAM Prediction")) +
                labs(
                    title = paste0("ALM Model Performance (RMSE: ", round(rmse, 4), ")"),
                    x = "Input",
                    y = "Prediction"
                ) +
                scale_color_manual(values = c("ALM Prediction" = "blue", "EXAM Prediction" = "red")) +
                theme_minimal()
        })
    })
}

# Run the application
shinyApp(ui = ui, server = server)



```





```{r, fig.width=11}
# Test varying input and output nodes more exhaustively
test_nodes_exhaustive <- function(dat, c = 0.5, lr = 0.2, input_nodes = seq(3, 15,3), output_nodes = seq(16, 128, 8)) {
  results <- list()

  for (inNodes in input_nodes) {
    for (outNodes in output_nodes) {
      result <- sim_data(dat, c, lr, inNodes, outNodes)
      rmse <- RMSE(dat$vx, result$almPred)
      results[[paste0("in", inNodes, "_out", outNodes)]] <- list(inNodes = inNodes, outNodes = outNodes, rmse = rmse)
    }
  }

  results
}

# Visualize the results using bar plots
visualize_bar_plots <- function(results) {
  df <- do.call(rbind, lapply(results, function(x) data.frame(t(unlist(x)))))
  
  ggplot(df, aes(x = interaction(inNodes, outNodes, sep = "_"), y = rmse)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    labs(title = "RMSE for varying input and output nodes",
         x = "Combinations of Input and Output Nodes",
         y = "RMSE") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8))
}

# Generate test data
dat <- gen_train()

# Test and visualize results using more exhaustive parameters
results_exhaustive <- test_nodes_exhaustive(dat)
visualize_bar_plots(results_exhaustive)



```










```{r}

# Generate synthetic data for testing
gen_train <- function(trainVec = c(5, 6, 7), trainRep = 10, noise = 0) {
    bandVec <- c(0, 100, 350, 600, 800, 1000, 1200)
    ts <- rep(seq(1, length(trainVec)), trainRep)
    noiseVec <- rnorm(length(ts), mean = 0) * noise
    if (noise == 0) {
        noiseVec <- noiseVec * 0
    }
    tibble(input = trainVec[ts], vx = bandVec[trainVec[ts]] + noiseVec)
}

# Test varying input and output nodes
test_nodes <- function(dat, c = 0.5, lr = 0.2, input_nodes = seq(3, 10), output_nodes = seq(16, 64, 8)) {
    results <- list()

    for (inNodes in input_nodes) {
        for (outNodes in output_nodes) {
            result <- sim_data(dat, c, lr, inNodes, outNodes)
            rmse <- RMSE(dat$vx, result$almPred)
            results[[paste0("in", inNodes, "_out", outNodes)]] <- list(inNodes = inNodes, outNodes = outNodes, rmse = rmse)
        }
    }

    results
}

# Visualize the results
visualize_results <- function(results) {
    df <- do.call(rbind, lapply(results, function(x) data.frame(t(unlist(x)))))
    ggplot(df, aes(x = inNodes, y = outNodes, fill = rmse)) +
        geom_tile() +
        scale_fill_gradient(low = "green", high = "red") +
        labs(
            title = "RMSE for varying input and output nodes",
            x = "Number of Input Nodes",
            y = "Number of Output Nodes",
            fill = "RMSE"
        ) +
        theme_minimal()
}

# Generate test data
dat <- gen_train()

# Test and visualize results
results <- test_nodes(dat)
visualize_results(results)



```