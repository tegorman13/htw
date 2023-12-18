INPUT_LAYER_DEFAULT <- seq(0, 100, 0.5)
OUTPUT_LAYER_DEFAULT <- seq(0, 250, 1)




adjust_layer <- function(input.layer, k){
  if(k == 1) return(input.layer)
  new.layer <- c()
  for(i in 1:(length(input.layer) - 1)){
    new.layer <- c(new.layer, seq(input.layer[i], input.layer[i+1], length.out = k + 1))
  }
  new.layer <- unique(new.layer)
  return(new.layer)
}






exam.response3 <- function(input, c, trainData, input.layer = INPUT_LAYER_DEFAULT, output.layer = OUTPUT_LAYER_DEFAULT, weight.mat) {
  # Compute ALM response for the current input
  almResp <- alm.response(input, c, input.layer, output.layer, weight.mat).response
  
  # Find associated y-values above and below the ALM response
  associatedYsAbove <- trainData[trainData > almResp]
  associatedYsBelow <- trainData[trainData < almResp]
  
  # Compute slopes for the associated y-values above and below the ALM response
  slopeAbove <- if (length(associatedYsAbove) > 1) sd(associatedYsAbove) / 100 else 0
  slopeBelow <- if (length(associatedYsBelow) > 1) sd(associatedYsBelow) / 100 else 0
  
  # Compute the overall slope as the average of the slopes above and below
  slope <- (slopeAbove + slopeBelow) / 2
  
  # Compute exam output based on estimated slope
  exam.output <- round(almResp + slope * (input - almResp), 3)
  
  exam.output
}



#trainData = dsAvg |> filter(condit=="Constant",expMode2=="Train")
exam.response2 <- function(input, c, trainData, input.layer = INPUT_LAYER_DEFAULT, output.layer = OUTPUT_LAYER_DEFAULT, weight.mat) {
  # Find nearest training x-value
  nearestTrainX <- trainData$x[which.min(abs(input - trainData$x))]

  # For the nearest training x-value, fetch all experienced y-values
  associatedYs <- trainData$y[trainData$x == nearestTrainX]

  # If we have more than one associated y-value, we can estimate a slope
  if (length(associatedYs) > 1) {
    # Calculate mean and standard deviation of associated y-values
    meanY <- mean(associatedYs)
    sdY <- sd(associatedYs)

    # Slope is approximated as the standard deviation of y-values divided by some constant (this can be adjusted/tuned)
    slope <- sdY / 100
  } else {
    # If only one y-value, default slope to zero
    slope <- 0
  }

  # Compute ALM response for nearestTrainX
  aresp <- alm.response(nearestTrainX, c, input.layer, output.layer, weight.mat)$mean.response

  # Compute exam output based on estimated slope
  exam.output <- round(aresp + slope * (input - nearestTrainX), 3)
  
  exam.output
}


# Define the sigmoid function
sigmoid <- function(x) {
  1 / (1 + exp(-x))
}

# Initialize the weights
weights <- matrix(runif(3), nrow = 2)

# Define the neural network
neural_network <- function(x, y, weights, learning_rate = 0.1) {
  # Calculate the activations of the hidden layer nodes
  hidden_activations <- sigmoid(weights %*% c(x, y))
  
  # Calculate the predicted criterion value
  y_pred <- sum(hidden_activations)
  
  # Calculate the error
  error <- y - y_pred
  
  # Update the weights
  weights <- weights + learning_rate * error * c(x, y)
  
  # Calculate the slope estimate
  delta <- (y - y_pred) / (x - x_prev)
  
  # Update the previous predictor value
  x_prev <- x
  
  return(list(weights = weights, y_pred = y_pred, delta = delta))
}

# Use the neural network to make predictions
predictions <- sapply(1:length(X), function(i) neural_network(X[i], Y[i], weights))









library(tidyverse)

# Define the fully connectionist architecture function
connectionist_model <- function(input, output_layer, weight_mat, gamma, c, alpha, X_train, Y_train) {
  # Input Activation Function
  input_activation <- function(X) {
    exp(-gamma * (X_train - X)^2)
  }
  
  # Output Activation Function
  output_activation <- function(input_act) {
    weight_mat %*% input_act
  }
  
  # Compute input activations and normalize
  act <- input_activation(input)
  normalized_act <- act / sum(act)
  
  # Compute output activations
  out_act <- output_activation(normalized_act)
  output_prob <- out_act / sum(out_act)
  mean_response <- sum(Y_train * output_prob)
  
  # Compute the slope within the network
  slope <- (mean_response[which.max(X_train > input)] - mean_response[which.min(X_train < input)]) /
           (X_train[which.max(X_train > input)] - X_train[which.min(X_train < input)])
  
  # Update the weights
  feedback_signal <- exp(-c * (output_layer - mean_response)^2)
  weight_updates <- (feedback_signal - out_act) * alpha
  weight_mat <<- weight_mat + weight_updates
  
  # Clip the weights to avoid negative values
  weight_mat[weight_mat < 0] <- 0
  
  # Output Calculation
  expected_output <- mean_response + slope * (input - X_train[which.min(abs(input - X_train))])
  
  return(expected_output)
}

# Initialize parameters and layers
gamma <- 1  # for example
c <- 1      # for example
alpha <- 0.1 # learning rate
X_train <- seq(1:10)
Y_train <- 2* X_train
input_layer <- seq(from = min(X_train), to = max(X_train), length.out=length(X_train)) # Adjust as necessary
output_layer <- Y_train
weight_mat <- matrix(runif(length(input_layer) * length(output_layer)), nrow = length(output_layer))

# Simulate the connectionist model for a new input
new_input <- 1.5 # Example input
connectionist_output <- connectionist_model(new_input, output_layer, weight_mat, gamma, c, alpha, X_train, Y_train)

# connectionist_output now holds the predicted value from the integrated model
print(connectionist_output)










library(tidyverse)

# Define the slope estimation layer within the network
slope_estimation_layer <- function(weight.mat, input.layer, output.layer) {
  # Calculate the slopes between each adjacent pair of inputs and outputs
  slopes <- numeric(length(output.layer) - 1)
  for (i in 1:(length(output.layer) - 1)) {
    deltaY <- output.layer[i + 1] - output.layer[i]
    deltaX <- input.layer[i + 1] - input.layer[i]
    slopes[i] <- deltaY / deltaX
  }
  # Pad the slopes to align with the number of outputs
  slopes <- c(slopes[1], slopes, slopes[length(slopes)])
  return(slopes)
}

# Modify the alm.trial function to include the slope estimation layer
alm_trial_with_slope <- function(input, corResp, c, lr, input.layer, output.layer, weight.mat) {
  alm_resp <- alm.response(input, c, input.layer, output.layer, weight.mat)
  updated_weight.mat <- alm.update(corResp, c, lr, output.layer, alm_resp$input.activation, alm_resp$output.activation, weight.mat)
  
  # Compute the slopes using the updated weight matrix
  slopes <- slope_estimation_layer(updated_weight.mat, input.layer, output.layer)
  
  return(list(mean.response = alm_resp$mean.response, weight.mat = updated_weight.mat, slopes = slopes))
}

# Define the full connectionist model
full_connectionist_model <- function(dat, c, lr, input.layer = INPUT_LAYER_DEFAULT, output.layer = OUTPUT_LAYER_DEFAULT) {
  weight.mat <- matrix(0.00, nrow = length(output.layer), ncol = length(input.layer))
  xt <- dat$x
  n <- nrow(dat)
  st <- numeric(n) # Initialize the vector to store mean responses
  slopes_list <- list() # Initialize the list to store slope estimates for each trial
  
  for(i in 1:n) {
    trial <- alm_trial_with_slope(dat$x[i], dat$y[i], c, lr, input.layer, output.layer, weight.mat)
    weight.mat <- trial$weight.mat
    st[i] <- trial$mean.response
    slopes_list[[i]] <- trial$slopes # Store the slopes for each trial
  }
  dat <- cbind(dat, almResp = st)
  return(list(d = dat, wm = weight.mat, c = c, lr = lr, slopes = slopes_list))
}

# Example usage with a dataset 'dat'
# Assume INPUT_LAYER_DEFAULT and OUTPUT_LAYER_DEFAULT are defined
# c and lr are hyperparameters for the ALM model
# dat is a dataframe with columns x and y representing the training data
result <- full_connectionist_model(dat, c = 0.5, lr = 0.1)



# Define the sigmoid function for monotonically tuned neurons
sigmoid <- function(x) {
  return(1 / (1 + exp(-x)))
}

# This function computes the network response based on the input
guigon.response <- function(input = 1, steepness, input.layer, output.layer, weight.mat) {
  input.activation <- sigmoid((input.layer - input) / steepness)
  output.activation <- weight.mat %*% input.activation
  output.probability <- output.activation / sum(output.activation)
  mean.response <- sum(output.layer * output.probability)
  
  list(mean.response = mean.response, input.activation = input.activation, output.activation = output.activation)
}

# This function updates the weight matrix using least-square error learning
guigon.update <- function(corResp, steepness, lr, output.layer, input.activation, output.activation, weight.mat) {
  desired_output <- rep(0, length(output.layer))
  desired_output[which.min(abs(output.layer - corResp))] <- 1 # One-hot encoding of desired response
  error <- desired_output - output.activation
  wChange <- lr * (error %*% t(input.activation))
  weight.mat <- weight.mat + wChange
  
  return(weight.mat)
}

# This function performs a single trial of the Guigon model
guigon.trial <- function(input, corResp, steepness, lr, input.layer, output.layer, weight.mat) {
  resp <- guigon.response(input, steepness, input.layer, output.layer, weight.mat)
  updated_weight.mat <- guigon.update(corResp, steepness, lr, output.layer, resp$input.activation, resp$output.activation, weight.mat)
  
  return(list(mean.response = resp$mean.response, weight.mat = updated_weight.mat))
}

# This function simulates the Guigon model across multiple trials
guigon.sim <- function(dat, steepness, lr, input.layer, output.layer) {
  weight.mat <- matrix(0.1, nrow = length(output.layer), ncol = length(input.layer)) # Initialize weights
  xt <- dat$x
  n <- nrow(dat)
  st <- numeric(n) # Initialize the vector to store mean responses
  
  for(i in 1:n) {
    trial <- guigon.trial(dat$x[i], dat$y[i], steepness, lr, input.layer, output.layer, weight.mat)
    weight.mat <- trial$weight.mat
    st[i] <- trial$mean.response
  }
  dat <- cbind(dat, guigonResp = st)
  
  return(list(d = dat, wm = weight.mat, steepness = steepness, lr = lr))
}

gs <- guigon.sim(i1,.4,.1,input.layer,output.layer )

