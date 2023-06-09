---
title: RNN
---




```{r}
pacman::p_load(tidyverse, glue, knitr, patchwork,kableExtra)
purrr::walk(c("Functions/alm_functions.R","Functions/Display_Functions.R"),source)

gen_train <- function(trainVec = c(5, 6, 7), trainRep = 3, noise = 0) {
    bandVec <- c(0, 100, 350, 600, 800, 1000, 1200)
    if (class(trainVec) == "list") {trainVec <- unlist(trainVec)}
    ts <- rep(seq(1, length(trainVec)), trainRep)
    noiseVec <- rnorm(length(ts), mean = 0) * noise
    if (noise == 0) {noiseVec <- noiseVec * 0}
    tibble(trial = seq(1, length(ts)), input = trainVec[ts], vx = bandVec[trainVec[ts]] + noiseVec)
}

input.activation<-function(x.target, c){return(exp((-1*c)*(x.target-inputNodes)^2))}

# Generate network layers and matrices
generate_network <- function(input_size, hidden_size, output_size) {
  input_weights <- matrix(runif(input_size * hidden_size, min = 0, max = .001), nrow = hidden_size, ncol = input_size)
  hidden_weights <- matrix(runif(hidden_size * hidden_size, min = 0, max = .001), nrow = hidden_size, ncol = hidden_size)
  output_weights <- matrix(runif(hidden_size * output_size, min = 0, max = .001), nrow = output_size, ncol = hidden_size)
  
  list(input_weights = input_weights, hidden_weights = hidden_weights, output_weights = output_weights)
}
hidden_rbf_activation <- function(input_act, hidden_state, gamma=.5) {
  distances <- matrix(0, nrow = nrow(hidden_state), ncol = 1)
  # for (i in 1:nrow(hidden_weights)) {
  #   distances[i] <- sum((input_act - hidden_weights[i,]) ^ 2)
  # }
   for (i in 1:nrow(hidden_state)) {
    distances[i] <- sum((input_act - hidden_state[i]) ^ 2)
  }
  exp(-gamma * distances)
}
rnn_alm_train <- function(dat, c = 0.05, lr = 0.5, network, hidden_size = 5,feedback_c = 0.5,gamma=.5) {
  input_weights <- network$input_weights
  hidden_weights <- network$hidden_weights
  output_weights <- network$output_weights
  
  hidden_state <- matrix(0, nrow = hidden_size, ncol = 1)
  #hidden_state <- matrix(seq(min(inputNodes),max(inputNodes),length.out=hidden_size))
  hidden_nodes <- matrix(seq(min(inputNodes),max(inputNodes),length.out=hidden_size))

  alm_train <- rep(NA, nrow(dat))
  
  # could have input connect to output & hidden, and hidden -> output
  # could transform hidden into radial basis function
  # just a feedforward version with the hidden layer
  for (i in 1:nrow(dat)) {
    # Input and output activations
    raw_input_act <- input.activation(dat$input[i], c)
    input_act <- raw_input_act %*% t(input_weights)
   # hidden_state <- sigmoid(hidden_weights %*% hidden_state + t(input_act))
   hidden_state= (hidden_state*.3) + hidden_rbf_activation(input_act, hidden_nodes, gamma)
   # hidden_state <- (hidden_weights %*% hidden_state + t(input_act))
    output_act <- t(t(hidden_state*1) %*% t(output_weights))
    
    # Mean prediction
    prob <- output_act/(sum(output_act)+.001)
    resp= outputNodes %*% prob
    alm_train[i] <- resp
    
    
  # Calculate output layer error
output_error <- (exp(-1 * c * (dat$vx[i] - outputNodes)^2) - output_act)

# Compute the derivatives of the RBF activation functions
hidden_rbf_derivative <- -2 * gamma * (hidden_state - hidden_nodes) * hidden_state
input_rbf_derivative <- -2 * gamma * ( input.activation(dat$input[i],c) - dat$input[i]) * input.activation(dat$input[i],c)

input.activation(dat$input[i],c) %*% t(input_weights)

# Backpropagate the error
hidden_error <- (output_weights %*% output_error) * hidden_rbf_derivative
input_error <- (hidden_weights %*% hidden_error) * input_rbf_derivative

# Update weights
output_weights <- output_weights + lr * output_error %*% t(hidden_state)
hidden_weights <- hidden_weights + lr * hidden_error %*% t(input_act)
input_weights <- input_weights + lr * input_error %*% t(input.activation(dat$input[i], c))

    
    # error <- (as.numeric(dat$vx[i] - resp))
    # feedback.activation<-exp(-1*c*(dat$vx[i]-outputNodes)^2)
    # 
    # output_weights <- output_weights + lr * (feedback.activation - output_act) %*% t(hidden_state)
    # hidden_weights <- hidden_weights + lr * error * (hidden_state * (1 - hidden_state)) %*% (input_act)
    # input_weights <- input_weights + lr * error * (hidden_state * (1 - hidden_state)) %*% (input.activation(dat$input[i], c))
    # output_weights[output_weights < 0] <- 0
    # input_weights[input_weights < 0] <- 0
  }
  
  list(almTrain = alm_train, network = list(input_weights = input_weights, hidden_weights = hidden_weights, output_weights = output_weights))
}


# Modify sim_train function to use the RNN version of ALM
sim_train_rnn <- function(dat, c = 0.5, lr = 0.2, inNodes = 7, outNodes = 32, hidden_size = 5, feedback_c = 0.05) {
  inputNodes <<- seq(1, 7, length.out = inNodes)
  outputNodes <<- seq(50, 1600, length.out = outNodes)
  network <- generate_network(length(inputNodes), hidden_size, length(outputNodes))
  tt <- rnn_alm_train(dat, c, lr, network, hidden_size,feedback_c)
}

dat=gen_train(trainRep=10)
tt <- rnn_alm_train(dat, c=.5, lr=.2, network, hidden_size)
dat %>% mutate(almTrain=tt$almTrain)
df %>% knitr::kable()
```


mean.prediction<-function(x.target, weights, c){
  probability<-output.activation(x.target, weights, c)/sum(output.activation(x.target, weights, c))
  return(outputNodes%*%probability) # integer prediction
}


```{r, fig.width=12,fig.height=10}
# c=1; lr=.5; noise_sd=.001; inNodes=7; outNodes=32; i=1; hidden_size=5;  trainVec = c(5, 6, 7); dat=gen_train(); gamma=.5; feedback_c=.5


# Activation functions
sigmoid <- function(x) {
  return(1 / (1 + exp(-x)))
}

# Replace sim_train with sim_train_rnn in the pipeline
tibble(crossing(
  c = c(.005,.05,1), lr = c(.0005,1,2), noise = c(0),
  inNodes = c(7), outNodes = c(32),
  hidden_size = c(1,3,7),
  trainVec = list(list(1,3,7)), trainRep = c(1000),
  lossFun = list("RMSE"),
  feedback_c=c(.01,.5,2),
  simNum = 1:1,
)) %>%
  mutate(id = seq(1, nrow(.)), td = pmap(list(trainVec, trainRep, noise), ~ gen_train(trainVec = .x, trainRep = ..2, noise = ..3))) %>%
  ungroup() %>%
  mutate(
    d = pmap(
      list(td, c, lr, inNodes, outNodes, hidden_size,feedback_c),
      ~ sim_train_rnn(dat = .x, c = ..2, lr = ..3, inNodes = ..4, outNodes = ..5, hidden_size = ..6,feedback_c = ..7)
    ),
    almTrainDat = map(d, "almTrain"),weights = map(d, "weights")
) %>%
  unnest(c(almTrainDat, td)) %>%
  select(-d) %>%
  mutate(input = as.factor(input)) %T>%
  {pfLc(.,groupVars=c("id","input","c","lr","feedback_c","hidden_size")) } %>% trainTab(nBlocks =12,groupVars=c("id","c","lr","feedback_c","hidden_size")) %>% {. ->> tt}


    #trainTab(nBlocks =12,groupVars=c("id","c","lr","hidden_size"))

```






To implement a Recurrent Neural Network (RNN) version of the ALM model, we can modify the existing ALM model by incorporating hidden states to capture temporal dependencies in the data. In this RNN-ALM version, we will use the same input and output nodes but include a hidden layer with recurrent connections. Here is an outline of the modifications and additional functions needed to implement the RNN-ALM:

Define the RNN architecture: We will modify the existing architecture by adding a hidden layer with recurrent connections, updating the weight matrices to include connections from input to hidden layer, hidden to output layer, and hidden layer to itself.

Update the output activation function: The output activation function will be updated to include the hidden layer.

Modify the weight update function: The weight update function should be adapted to include updates for the new weight matrices connecting the input, hidden, and output layers.

Adapt the training and simulation functions: The training and simulation functions will need to be modified to accommodate the new RNN-ALM architecture and weight update function.

Here is an RNN-ALM implementation based on the provided ALM code:
```{r}
#| eval: false


# Add additional libraries
pacman::p_load(Rcpp)
# Register the Rcpp function to be used in R
sourceCpp("Functions/compute_gradients_bptt.cpp")


# Set the number of hidden units
H <- 10

# 1. Define the RNN architecture
init_weights <- function(M, H, N) {
  list(
    w_xh = matrix(runif(M * H, min = -1, max = 1), nrow = H, ncol = M),
    w_hh = matrix(runif(H * H, min = -1, max = 1), nrow = H, ncol = H),
    w_ho = matrix(runif(H * N, min = -1, max = 1), nrow = N, ncol = H)
  )
}

# 2. Update the output activation function
output_activation_rnn <- function(x_target, weights, c) {
  a_x <- input.activation(x_target, c)
  a_h <- tanh(weights$w_xh %*% a_x + weights$w_hh %*% h_prev)
  h_prev <<- a_h
  return(weights$w_ho %*% a_h)
}

# 3. Modify the weight update function
update_weights_rnn <- function(x_new, y_new, weights, c, lr) {
  y_feedback_activation <- exp(-1 * c * (y_new - outputNodes)^2)
  x_feedback_activation <- output_activation_rnn(x_new, weights, c)
  
  # Calculate gradients using backpropagation through time (BPTT)
  # For simplicity, we use Rcpp package for efficient computation
  # Add Rcpp code to calculate gradients for w_xh, w_hh, and w_ho
  # (code not provided here)
  
  # Update the weights using the gradients
  weights$w_xh <- weights$w_xh + lr * dw_xh
  weights$w_hh <- weights$w_hh + lr * dw_hh
  weights$w_ho <- weights$w_ho + lr * dw_ho
  
  return(weights)
}

# 4. Adapt the training and simulation functions
train_alm_rnn <- function(dat, c = 0.05, lr = 0.5, weights) {
  alm_train <- rep(NA, nrow(dat))
  h_prev <<- rep(0, H)
  
  for (i in 1:nrow(dat)) {
    weights <- update_weights_rnn(dat$input[i], dat$vx[i], weights, c, lr)
    resp <- mean.prediction
    (dat$input[i], weights, c)
alm_train[i] = resp
weights$w_xh[weights$w_xh < 0] = 0
weights$w_hh[weights$w_hh < 0] = 0
weights$w_ho[weights$w_ho < 0] = 0
}

list(almTrain = alm_train, weights = weights)
}

sim_train_rnn <- function(dat, c = 0.5, lr = 0.2, inNodes = 7, outNodes = 32, hiddenNodes = 10, trainVec = c(5, 6, 7), noise_sd = 0, update_func = "update.weights_rnn") {
inputNodes <<- seq(1, 7, length.out = inNodes * 1)
outputNodes <<- seq(50, 1600, length.out = outNodes * 1)
#Initialize the weights
weights <- init_weights(length(inputNodes), hiddenNodes, length(outputNodes))

#Train the RNN-ALM model
tt <- train_alm_rnn(dat, c, lr, weights)
}

#Modify the main script to use the RNN-ALM functions

tibble(crossing(
c = c(.5, 5), lr = c(.05, 1), noise = c(0),
inNodes = c(7), outNodes = c(32), hiddenNodes = c(10),
trainVec = list(list(5, 6, 7)), trainRep = c(9),
lossFun = list("MAE"),
simNum = 1:1,
)) %>%
mutate(id = seq(1, nrow(.)), td = pmap(list(trainVec, trainRep, noise), ~ gen_train(trainVec = .x, trainRep = ..2, noise = ..3))) %>%
ungroup() %>%
mutate(
d = pmap(
list(td, c, lr, inNodes, outNodes, hiddenNodes),
~ sim_train_rnn(dat = .x, c = ..2, lr = ..3, inNodes = ..4, outNodes = ..5, hiddenNodes = ..6)
),
almTrainDat = map(d, "almTrain"), weights = map(d, "weights")
) %>%
unnest(c(almTrainDat, td)) %>%
select(-d) %>%
mutate(input = as.factor(input)) %>%
trainTab()








```


This RNN-ALM implementation is based on the provided ALM code and adds a hidden layer with recurrent connections to capture temporal dependencies in the data. The code modifies the architecture, output activation function, weight update function, and adapts the training and simulation functions to accommodate the new RNN-ALM model. Please note that the backpropagation through time (BPTT) gradients computation using Rcpp is not provided here. You can implement this part using existing Rcpp examples and tutorials or find an existing R package that handles RNN training.
