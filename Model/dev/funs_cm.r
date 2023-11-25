alm.responseOnly <- function(input = 1, c, input.layer, output.layer,weight.mat, trainVec=NULL) {
  input.activation <- exp(-c * (input.layer - input)^2) / (sum(exp(-c * (input.layer - input)^2)) +.01)
  output.activation <- (weight.mat %*% input.activation) 
  output.probability <- output.activation / sum(output.activation) +.01
 mean.response <- sum(output.layer * output.probability)
}


alm.response <- function(input = 1, c, input.layer, output.layer,weight.mat, trainVec=NULL) {
  input.activation <- exp(-c * (input.layer - input)^2) / (sum(exp(-c * (input.layer - input)^2)) +.0001)
  output.activation <- (weight.mat %*% input.activation) 
  output.probability <- output.activation / (sum(output.activation) +.0001)
  mean.response <- sum(output.layer * output.probability)
  list(mean.response = mean.response, input.activation = input.activation, output.activation = output.activation)
}

alm.update <- function(corResp, c, lr, output.layer, input.activation, output.activation, weight.mat) {
  fz <- exp(-c * (output.layer - corResp)^2) + .001
  teacherSignal <- (fz - output.activation) * lr
  wChange <- teacherSignal %*% t(input.activation)
  weight.mat <- weight.mat + wChange
  return(weight.mat)
}

alm.trial <- function(input, corResp, c, lr, input.layer, output.layer, weight.mat) {
  alm_resp <- alm.response(input, c, input.layer,output.layer, weight.mat)
  updated_weight.mat <- alm.update(corResp, c, lr, output.layer, alm_resp$input.activation, alm_resp$output.activation, weight.mat)
  return(list(mean.response = alm_resp$mean.response, weight.mat = updated_weight.mat))
}

alm.sim <- function(dat, c, lr, 
  input.layer = INPUT_LAYER_DEFAULT, output.layer = OUTPUT_LAYER_DEFAULT, 
  weight.mat=NULL) {

  # initialize weight.mat only if it is not provided
  if(is.null(weight.mat)){
  weight.mat <- matrix(0.000001, nrow = length(output.layer), ncol = length(input.layer))
  }

  xt <- dat$x
  n <- nrow(dat)
  st <- numeric(n) 
  for(i in 1:n) {
    trial <- alm.trial(dat$x[i], dat$y[i], c, lr, input.layer, output.layer, weight.mat)
    weight.mat <- trial$weight.mat
    st[i] <- trial$mean.response
  }
  dat <- cbind(dat,almResp=st)
  return(list(d = dat, wm = weight.mat, c = c, lr = lr))
}



exam.response <- function(input, c, input.layer = INPUT_LAYER_DEFAULT,output.layer = OUTPUT_LAYER_DEFAULT, weight.mat, trainVec) {
  nearestTrain <- trainVec[which.min(abs(input - trainVec))]
  aresp <- alm.response(nearestTrain, c, input.layer = input.layer,output.layer,weight.mat)$mean.response
  
  xUnder <- ifelse(min(trainVec) == nearestTrain, nearestTrain, trainVec[which(trainVec == nearestTrain) - 1])
  xOver <- ifelse(max(trainVec) == nearestTrain, nearestTrain, trainVec[which(trainVec == nearestTrain) + 1])
  
  mUnder <- alm.response(xUnder, c, input.layer = input.layer, output.layer, weight.mat)$mean.response
  mOver <- alm.response(xOver, c, input.layer = input.layer,output.layer, weight.mat)$mean.response
  
  exam.output <- round(aresp + ((mOver - mUnder) / (xOver - xUnder)) * (input - nearestTrain), 3)
  exam.output
}

# Describe EXAM's response process: 
# 1. Calculate the nearest training item to the test item
# 2. Calculate the mean ALM response for the nearest training item
# 3. Identify the nearest training item to the left and right of the nearest training item
# 4. Calculate the mean ALM response for the nearest training item to the left and right
# 5. Calculate the slope of the line between the two nearest training items
# 6. Compute the mean response for the test item

# Alternative EXAM response process for when training only contains one item:
# 1. Calculate the mean ALM response for the nearest training item
# 2. Use output values from training that were above the mean ALM response to infer additional input value. 
# 3. Use the ALM response for the nearest item, and for the inferred item, to compute a slope. 
alt_exam <- function(input, c, input.layer = input_layer, output.layer = output_layer, weight.mat, trainVecX, trainVecY) {
    if (length(input) > 1) {
        stop("Input must be scalar")
    }
    if (length(trainVecX) == 0) {
        stop("No training data")
    }

    nearestTrain <- trainVecX[which.min(abs(input - trainVecX))]
    aresp <- alm.response(nearestTrain, c, input.layer, output.layer, weight.mat)$mean.response
    
    if (nearestTrain>input) {
        upper_output <- trainVecY[which(trainVecY > aresp)]
        if (length(upper_output) == 0) {
            return(alm.response(input, c, input.layer, output.layer, weight.mat)$mean.response)
        }
        upper_output <- upper_output[which.max(abs(upper_output - nearestTrain))]
        infer_input <- nearestTrain + (upper_output - aresp)
        slope <- (alm.response(infer_input, c, input.layer, output.layer, weight.mat)$mean.response - aresp) / (infer_input - nearestTrain)
        return(round(aresp + slope * (input - nearestTrain), 3))
    }

    lower_output <- trainVecY[which(trainVecY < aresp)]
    if (length(lower_output) == 0) {
        return(alm.response(input, c, input.layer, output.layer, weight.mat)$mean.response)
    }
    lower_output <- lower_output[which.min(abs(lower_output - nearestTrain))]
    infer_input <- nearestTrain + (aresp - lower_output)
    slope <- (alm.response(infer_input, c, input.layer, output.layer, weight.mat)$mean.response - aresp) / (infer_input - nearestTrain)
    return(round(aresp + slope * (input - nearestTrain), 3))
}