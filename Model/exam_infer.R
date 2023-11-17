
library(dplyr)
library(purrr)
library(tidyr)



alm.response <- function(input = 1, c, input.layer, output.layer,weight.mat, trainVec=NULL) {
  input.activation <- exp(-c * (input.layer - input)^2) / (sum(exp(-c * (input.layer - input)^2)) +.01)
  output.activation <- (weight.mat %*% input.activation) 
  output.probability <- output.activation / sum(output.activation) +.01
  mean.response <- sum(output.layer * output.probability)
  list(mean.response = mean.response, input.activation = input.activation, output.activation = output.activation)
}

alm.update <- function(corResp, c, lr, output.layer, input.activation, output.activation, weight.mat) {
  fz <- exp(-c * (output.layer - corResp)^2)
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

alm.sim <- function(dat, c, lr, input.layer = INPUT_LAYER_DEFAULT, output.layer = OUTPUT_LAYER_DEFAULT) {
  weight.mat <- matrix(0.000001, nrow = length(output.layer), ncol = length(input.layer))
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


exam_infer.response <- function(input, c, input_layer, output_layer, weight.mat, trainVec, trainOutputs) {
  # Get the ALM response to the nearest training item
  nearestTrain <- trainVec[which.min(abs(input - trainVec))]
  aresp <- alm.response(nearestTrain, c, input_layer, output_layer, weight.mat = weight.mat)$mean.response
  
  # Retrieve training outputs above and below the learned response to the training item
  sortedTrainOutputs <- sort(trainOutputs)
  yRetrieveLow <- sortedTrainOutputs[which.min(abs(aresp - sortedTrainOutputs))]
  yRetrieveHigh <- sortedTrainOutputs[which.max(abs(aresp - sortedTrainOutputs))]
  
  xInferLow <- (yRetrieveLow - aresp) / 2 + input
  xInferHigh <- (yRetrieveHigh - aresp) / 2 + input
  
  slope <- (yRetrieveHigh - yRetrieveLow) / (xInferHigh - xInferLow)
  exam.output <- aresp + slope * (input-nearestTrain)
  return(exam.output)
}


exam_infer.response <- function(input, c, input_layer, output_layer, weight.mat, trainVec, trainOutputs) {
  # Get the ALM response to the nearest training item
  nearestTrain <- trainVec[which.min(abs(input - trainVec))]
  aresp <- alm.response(nearestTrain, c, input_layer, output_layer, weight.mat = weight.mat)$mean.response
  
  # Retrieve training outputs above and below the learned response to the training item
  sortedTrainOutputs <- sort(trainOutputs)
  
  # Calculate similarity and dissimilarity scores
  similarity_scores <- exp(-c * (sortedTrainOutputs - aresp)^2)
  dissimilarity_scores <- 1 - similarity_scores
  
  
  # Rescale similarity and dissimilarity scores
  similarityScores <- (similarityScores - min(similarityScores)) / (max(similarityScores) - min(similarityScores))
  dissimilarityScores <- (dissimilarityScores - min(dissimilarityScores)) / (max(dissimilarityScores) - min(dissimilarityScores))
  
  balanceWeight=.9
  combinedScores <- balanceWeight * dissimilarityScores + (1 - balanceWeight) * similarityScores
  combinedScores
  
  outputsBelow <- sortedTrainOutputs[sortedTrainOutputs < aresp]
  outputsAbove <- sortedTrainOutputs[sortedTrainOutputs > aresp]
  
  if(length(outputsBelow) > 0) {
    scoresBelow <- combinedScores[sortedTrainOutputs < aresp]
    yRetrieveLow <- outputsBelow[which.max(scoresBelow)]
  } else {
    yRetrieveLow <- min(sortedTrainOutputs)  # default if no outputs are below
  }
  
  if(length(outputsAbove) > 0) {
    scoresAbove <- combinedScores[sortedTrainOutputs > aresp]
    yRetrieveHigh <- outputsAbove[which.min(scoresAbove)]  # Corrected to which.min
  } else {
    yRetrieveHigh <- max(sortedTrainOutputs)  # default if no outputs are above
  }
  
  similarityLow <- exp(-c * (yRetrieveLow - aresp)^2)
  similarityHigh <- exp(-c * (yRetrieveHigh - aresp)^2)
  # Infer X value for yRetrieveLow and yRetrieveHigh
  
  xInferLow <- input - (input - nearestTrain) * similarityLow
  xInferHigh <- input + (nearestTrain - input) * similarityHigh
  
  
  xInferLow <- input - (aresp - yRetrieveLow) * (input - nearestTrain) / (aresp - min(sortedTrainOutputs))
  xInferHigh <- input + (yRetrieveHigh - aresp) * (input - nearestTrain) / (max(sortedTrainOutputs) - aresp)
  
  
  # Use InferLow, InferHigh, yRetrieveLow, and yRetrieveHigh, and aresp to compute slope
  slope <- (yRetrieveHigh - yRetrieveLow) / (xInferHigh - xInferLow)
  
  # Calculate the exam output using the inferred slope
  exam.output <- aresp + slope * (input - nearestTrain)
  
  return(exam.output)
}





input_layer <- c(100, 350, 600, 800, 1000, 1200)
output_layer <- c(100, 350, 600, 800, 1000, 1200)

ds <- readRDS(here::here("data/e1_md_11-06-23.rds")) |> mutate(sbj=id) |> relocate(sbj,.after=id)


split_data <- split(ds, ds$id)
split_data <- ds %>% split(.$id) |> head(1)
c = .0008; lr=2.0; 
#split_data[[1]] |> filter(expMode2=="Train")

train_results <- alm.sim(split_data[[1]] |> filter(expMode2=="Train"), 
                         c, lr, input_layer, output_layer)

trainVec <- unique(train_results$d$x) |> sort()
trainOutputs <- unique(train_results$d$y) |> sort()
# generate predictions from exam_infer.response
map_dbl(c(100, 350, 600, 800, 1000, 1200), ~ exam_infer.response(.x, c, input_layer, output_layer, train_results$wm,  trainVec=trainVec, trainOutputs))

test_prediction <- map_dbl(input_layer, ~ exam_infer.response(.x, c, input_layer, output_layer, train_results$wm,  trainVec=trainVec, trainOutputs))
map_dbl(input_layer, ~ alm.responseOnly(.x, c, input_layer, output_layer, train_results$wm))
map_dbl(input_layer, ~ exam.response(.x, c, input_layer, output_layer, train_results$wm,  trainVec=trainVec))




split_data <- split(ds, ds$id)
split_data <- split_data[c(3)]
c = .0008; lr=2.0; 
train_results <- alm.sim(split_data[[1]] |> filter(expMode2=="Train"), 
                         c, lr, input_layer, output_layer)

trainVec <- c(0,unique(train_results$d$x)) |> sort()
trainOutputs <- unique(train_results$d$y) |> sort()
map_dbl(350, ~ exam_infer.response(.x, c, input_layer, output_layer, train_results$wm,  trainVec=trainVec, trainOutputs))
map_dbl(input_layer, ~ alm.responseOnly(.x, c, input_layer, output_layer, train_results$wm))
map_dbl(input_layer, ~ exam.response(.x, c, input_layer, output_layer, train_results$wm,  trainVec=trainVec))

