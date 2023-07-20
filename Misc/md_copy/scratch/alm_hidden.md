---
title: alm hidden
date: last-modified
---





gen_train <- function(trainVec=c(5,6,7),trainRep=3,noise=0){
  bandVec=c(0,100,350,600,800,1000,1200)
  if(class(trainVec)=="list"){trainVec=unlist(trainVec)}
  ts <- rep(seq(1,length(trainVec)),trainRep)

  noiseVec=rnorm(length(ts),mean=0)*noise
  if(noise==0) {noiseVec=noiseVec*0}
  tibble(trial=seq(1,length(ts)),input=trainVec[ts],vx=bandVec[trainVec[ts]]+noiseVec)
}




```{r}

pacman::p_load(tidyverse)
# Modify the sim_data function to accept the dataset as an argument
sim_data <- function(dat, c=0.5, lr=0.2, inNodes=7, outNodes=32, trainVec=c(5,6,7)) {
  inputNodes <<- seq(1,inNodes,length.out=inNodes*1)  
  outputNodes <<- seq(50,1600,length.out=outNodes*1) 
  wm=matrix(.0000,nrow=length(outputNodes),ncol=length(inputNodes))
  tt<-trainTest.alm(dat, c, lr, wm, trainVec)
}

# Define a function to calculate the ideal output for a given input value
ideal_output <- function(x) {
  outputNodes[which.min(abs(outputNodes - (x * 50)))]
}

input_activation <- function(input_nodes,x, c) {
    exp(-c * (x - input_nodes)^2)
  }
  # Sigmoid activation function for the hidden layer
sigmoid_activation <- function(x) {
    1 / (1 + exp(-x))
}
  # Output layer activation function
  output_activation <- function(hidden_activations, weights) {
    hidden_activations %*% weights
  }
  

# Modify the gen_train function to return a data frame with input, output, and trial columns
gen_train <- function(trainVec=c(5,6,7),trainRep=3,noise=0){
  bandVec=c(0,100,350,600,800,1000,1200)
  if(class(trainVec)=="list"){trainVec=unlist(trainVec)}
  ts <- rep(seq(1,length(trainVec)),trainRep)
  # print(trainVec)
  # print(length(ts)); print(length(trainRep))
  noiseVec=rnorm(length(ts),mean=0)*noise
  if(noise==0) {noiseVec=noiseVec*0}
  df <- tibble(input = trainVec[ts], trial = seq_along(ts),vx=bandVec[trainVec[ts]]+noiseVec)
  return(df)
}

# New ALM with Hidden Layer
New_ALM <- function(input_hidden_weights, hidden_output_weights, bias_hidden_weights, bias_output_weights, inNodes,outNodes,hidNodes, dat, param) {
  # Input layer and hidden layer positions
  input_nodes <- seq(1, inNodes)
  hidden_nodes <- seq(1, hidNodes)
  output_nodes <- seq(1, outNodes)
  # Gaussian activation function for the input layer
   c=param["c"]
   lrA = param[2]
   lrW = param[3]
  
  # Training the model
  n = nrow(dat)
  predictions = rep(0, n)
  for(i in 1:n) {
    x_target = dat$input[i]
    
    # Calculate input activations
    input_activations = input_activation(input_nodes,x_target, c)
    # Calculate hidden layer activations
    #hidden_activations =sigmoid_activation(input_hidden_weights %*% input_activations + (bias_hidden_weights))
    hidden_activations =input_hidden_weights %*% input_activations #+ (bias_hidden_weights)

    # Calculate output activations
    out_activations = (hidden_output_weights %*% hidden_activations) #+ (bias_output_weights)
    # Generate a response by taking the weighted average of the output units
    responseH = sum(hidden_activations * output_nodes) / sum(hidden_activations)
    response = sum(out_activations * output_nodes) / sum(out_activations)
    predictions[i] = response
    error = dat$vx[i] - response
    
    y.feedback.activation<-exp(-1*c*(dat$vx[i]-output_nodes)^2)
    act_error <- y.feedback.activation - out_activations
    
    hidden_output_weights <- hidden_output_weights + lrW*act_error %*% t(hidden_activations)
    input_hidden_weights = input_hidden_weights +lrA*hidden_activations %*% t(input_activations)
    
    # Update bias weights
    bias_hidden_weights = (bias_hidden_weights) + (lrA * error * hidden_activations * (1 - hidden_activations))
    bias_output_weights = (bias_output_weights) + lrW * error
    }

return(predictions)
}

init_weights <- function(nrow, ncol, mu, sigma) {
  random_vector = mu + sigma * runif(nrow * ncol, 0, 1)
  temp_weights <- matrix(random_vector, nrow = nrow, ncol = ncol)
  return(temp_weights)
}
#Generate training data
set.seed(123)
train_data <- gen_train(trainVec = c(2, 4, 6), trainRep = 6, noise = 0)

#Initialize weights for input to hidden layer, hidden to output layer, and bias terms
input_hidden_weights <- init_weights(5, 7, 0.0025, 0.001)
hidden_output_weights <- init_weights(32, 5, 0.0025, 0.001)
bias_hidden_weights <- runif(5,0,.01)
bias_output_weights <- runif(32,0,.01)

#Parameter settings
c = 0.1
phi = 6
lrA = .4
lrW = 0.4
param = c(c=sc, lrA=lrA, lrW=lrW, phi)
train_data <- gen_train(trainVec = c(2, 4, 6), trainRep = 6, noise = 0)

#Run the new ALM model
result <- New_ALM(input_hidden_weights, hidden_output_weights, bias_hidden_weights, bias_output_weights, 7,32,7, train_data, param)
result


#Print the results
cat("Input Ideal Predicted\n")
for (i in 1:nrow(train_data)) {
cat(train_data[i, "input"], " ", train_data[i, "output"], " ", result[i], "\n")
}










```











```{r}

pacman::p_load(tidyverse,data.table)
options(dplyr.summarise.inform=FALSE)

input.activation<-function(x.target, c){
  return(exp((-1*c)*(x.target-inputNodes)^2))
}

# Add hidden layer activation function
hidden.activation <- function(x.target, input_weights, c){
  return(input_weights %*% input.activation(x.target, c))
}

# Update output activation function to include hidden layer
output.activation <- function(x.target, input_weights, hidden_weights, c){
  hidden_activations <- hidden.activation(x.target, input_weights, c)
  return(hidden_weights %*% hidden_activations)
}

# Modify mean.prediction function to include hidden layer
mean.prediction <- function(x.target, input_weights, hidden_weights, c){
  probability <- output.activation(x.target, input_weights, hidden_weights, c) / sum(output.activation(x.target, input_weights, hidden_weights, c))
  return(outputNodes %*% probability)
}

# Update weights function to include hidden layer
update.weights <- function(x.new, y.new, input_weights, hidden_weights, c, lr, noise_sd = NULL){
  y.feedback.activation <- exp(-1 * c * (y.new - outputNodes)^2)
  x.feedback.activation <- output.activation(x.new, input_weights, hidden_weights, c)
  
  hidden_activations <- hidden.activation(x.new, input_weights, c)
  
  delta_hidden_weights <- lr * (y.feedback.activation - x.feedback.activation) %*% t(hidden_activations)
  new_hidden_weights <- hidden_weights + delta_hidden_weights
  
  #delta_input_weights <- lr * t(hidden_activations * (y.feedback.activation - x.feedback.activation)) %*% t(input.activation(x.new, c))
  delta_input_weights <- lr * outer(hidden_activations, (y.feedback.activation - x.feedback.activation)) %*% t(input.activation(x.new, c))

  new_input_weights <- input_weights + delta_input_weights
  
  return(list(input_weights = new_input_weights, hidden_weights = new_hidden_weights))
}

# Modify train.alm function to include hidden layer
train.alm <- function(dat, c=0.05, lr=0.5, input_weights, hidden_weights){
  alm.train <- rep(NA, nrow(dat))
  
  for (i in 1:nrow(dat)){
    weights <- update.weights(dat$input[i], dat$vx[i], input_weights, hidden_weights, c, lr)
    input_weights <- weights$input_weights
    hidden_weights <- weights$hidden_weights
    
    resp = mean.prediction(dat$input[i], input_weights, hidden_weights, c)
    alm.train[i] = resp
    
    input_weights[input_weights < 0] = 0
    hidden_weights[hidden_weights < 0] = 0
  }
  
  alm.train
}

# Modify sim_data function to include hidden layer
sim_data <- function(dat, c=0.5, lr=0.2, inNodes=7, hiddenNodes=5, outNodes=32, trainVec=c(5,6,7)) {
  inputNodes <<- seq(1, 7, length.out=inNodes*1)
  outputNodes <<- seq(50, 1600, length.out=outNodes*1)
  
  input_weights <- matrix(0, nrow=hiddenNodes, ncol=inNodes)
  hidden_weights <- matrix(0, nrow=length(outputNodes), ncol=hiddenNodes)
  
  tt <- train.alm(dat, c, lr, input_weights, hidden_weights)
}

trainTest.alm<-function(dat, c=0.05, lr=0.5, weights){
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
```



```{r}
# Generate data from a population of ALM learning agents
generate_population_data <- function(dat, c_values, lr_values, inNodes=7, hiddenNodes=5, outNodes=32, trainVec=c(5,6,7)) {
  population_data <- list()
  
  for (c_val in c_values) {
    for (lr_val in lr_values) {
      cat("Running simulation for c =", c_val, "and lr =", lr_val, "\n")
      
      # Run simulation for the current combination of c and lr values
      sim_result <- sim_data(dat, c=c_val, lr=lr_val, inNodes=inNodes, hiddenNodes=hiddenNodes, outNodes=outNodes, trainVec=trainVec)
      
      # Store the result in the population_data list
      param_combo <- paste("c", c_val, "lr", lr_val, sep="_")
      population_data[[param_combo]] <- sim_result
    }
  }
  
  return(population_data)
}

# Example usage of generate_population_data
c_values <- seq(0.1, 1, by=0.1)
lr_values <- seq(0.1, 1, by=0.1)

# Assume 'dat' is the dataset for the ALM model
population_data <- generate_population_data(dat=gen_train(), c_values, lr_values)


c=.05; lr=.5
i=1
y.new=dat$vx[i]
x.new=dat$input[i]
```




```{r}

#Here's a new version of ALM with 7 input nodes, 32 output nodes, and a hidden layer with 10 nodes:
input.activation <- function(x.target, association.parameter) {
  return(exp(-1 * association.parameter * (x.target - x.plotting)^2))
}

hidden.activation <- function(input.activations, weights, association.parameter) {
  return(exp(-1 * association.parameter * (weights %*% input.activations)))
}

output.activation <- function(hidden.activations, weights, association.parameter) {
  return(weights %*% hidden.activations)
}

mean.prediction <- function(x.target, weights, association.parameter) {
  hidden.activations <- hidden.activation(input.activation(x.target, association.parameter), weights, association.parameter)
  probability <- output.activation(hidden.activations, weights, association.parameter) / sum(output.activation(hidden.activations, weights, association.parameter))
  return(y.plotting %*% probability) # integer prediction
}

update.weights <- function(x.new, y.new, weights, association.parameter, update.parameter) {
  hidden.activations <- hidden.activation(input.activation(x.new, association.parameter), weights, association.parameter)
  y.feedback.activation <- exp(-1 * association.parameter * (y.new - y.plotting)^2)
  x.feedback.activation <- output.activation(hidden.activations, weights, association.parameter)
  delta.weights <- update.parameter * (y.feedback.activation - x.feedback.activation) %*% t(hidden.activations %*% input.activation(x.new, association.parameter))
  return(weights + delta.weights)
}

learn.alm <- function(y.learning, association.parameter = 0.05, update.parameter = 0.5) {
  weights1 <- matrix(rnorm(7 * 10, mean = 0, sd = 0.1), nrow = 7, ncol = 10)
  weights2 <- matrix(rnorm(10 * 32, mean = 0, sd = 0.1), nrow = 32, ncol = 10)
  for (i in 1:length(y.learning)) {
    hidden.activations <- hidden.activation(input.activation(x.learning[i], association.parameter), weights1, association.parameter)
    weights2 <- update.weights(x.learning[i], y.learning[i], weights2, association.parameter, update.parameter, hidden.activations)
    resp = mean.prediction(x.learning[i], weights2, association.parameter)
    weights2[weights2 < 0] = 0
  }
  alm.predictions <- sapply(x.plotting, mean.prediction, weights = weights2, association.parameter = association.parameter)
  exam.predictions <- sapply(x.plotting, exam.prediction, weights = weights2, association.parameter = association.parameter)
  return(list(alm.predictions = alm.predictions, exam.predictions = exam.predictions))
}



```




```{r}
# Define input and output variables
x.learning <- seq(0, 100, length.out = 100)
y.learning <- sin(x.learning) + rnorm(length(x.learning), mean = 0, sd = 0.1)
x.plotting <- seq(0, 100, length.out = 10)
y.plotting <- sin(x.plotting)


# Train ALM model
results <- learn.alm(y.learning)

# Plot predictions
plot(x.learning, y.learning, type = "l", col = "black")
lines(x.learning, results$alm.predictions, col = "red")
lines(x.learning, results$exam.predictions, col = "blue")
legend("topright", legend = c("ALM", "Exemplar"), col = c("red", "blue"), lty = 1)
```

