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
   dat <- cbind(dat,almResp=st)
  return(list(d = dat, wm = weight.mat, c = c, lr = lr))
}

exam.response <- function(input, c, trainVec, input.layer = INPUT_LAYER_DEFAULT,output.layer = OUTPUT_LAYER_DEFAULT, weight.mat) {
  nearestTrain <- trainVec[which.min(abs(input - trainVec))]
  aresp <- alm.response(nearestTrain, c, input.layer = input.layer,output.layer,weight.mat)$mean.response
  
  xUnder <- ifelse(min(trainVec) == nearestTrain, nearestTrain, trainVec[which(trainVec == nearestTrain) - 1])
  xOver <- ifelse(max(trainVec) == nearestTrain, nearestTrain, trainVec[which(trainVec == nearestTrain) + 1])
  
  mUnder <- alm.response(xUnder, c, input.layer = input.layer, output.layer, weight.mat)$mean.response
  mOver <- alm.response(xOver, c, input.layer = input.layer,output.layer, weight.mat)$mean.response
  
  exam.output <- round(aresp + ((mOver - mUnder) / (xOver - xUnder)) * (input - nearestTrain), 3)
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