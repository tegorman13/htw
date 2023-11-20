
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
 fdissimilarity_scores <- 1 - similarity_scores
  
 f
  # Rescale similarity and dissimilarity scores
  similarityScores <- (similarityScores - min(similarityScores)) / (max(similarityScores) - min(similarityScores))
  dissimilarityScores <- (dissimilarityScores -fmin(dissimilarityScores)) / (max(dissimilarityScores) - min(dissimilarityScores))
  
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
    nearestTrain <- trainVecX[which.min(abs(input - trainVecX))]
    aresp <- alm.response(nearestTrain, c, input.layer, output.layer, weight.mat)$mean.response
    # case when only single training item
    if (length(trainVecX) >0) {
        # use output values from training (trainVecY) that were above the mean ALM response to infer additional input value
        upper_output <- trainVecY[which(trainVecY > aresp)]
        # randomly select one of the values in upper_output, only 1 value should be in upper_output
        #upper_output <- upper_output[sample(1:length(upper_output),size=1)]
        # randomly select the upper_output value furthest away from the nearest training item
        upper_output <- upper_output[which.max(abs(upper_output - nearestTrain))]



        # infer input value based ut <- upper_output[sample(1:length(upper_output))]on distance between alm.response and upper_output
        infer_input <- nearestTrain + (upper_output - aresp)
        # use the ALM response for the nearest item, and for the inferred item, to compute a slope
        slope <- (alm.response(infer_input, c, input.layer, output.layer, weight.mat)$mean.response - aresp) / (infer_input - nearestTrain)
        exam_output <- round(aresp + slope * (input - nearestTrain), 3)
    }

}


## Toy scenario 
# input space ranges from X=1:7. True Function is Y=X. 
# Trained from X=3:5
# Assume perfect learning - initalize weight matrix diagonals to 1 for trained pairs. 

input_layer <- c(1,2,3,4,5,6,7)
output_layer <- c(1,2,3,4,5,6,7)

# generate weight matrix reflective of learning from X=3:5
# value of 1 for trained pairs, 0 else where
weight_mat <- matrix(0, nrow=7, ncol=7)
# only portion of diag that was trained set to 1, e.g. (3,3), (4,4), (5,5). 
weight_mat[3,3] <- 1
weight_mat[4,4] <- 1
weight_mat[5,5] <- 1

# generate training data
#train_vec <- c(3,4,5)
train_vec <- c(3)
train_outputs <- train_vec

n=5

td <- tibble(x=rep(train_vec, each=n), y=rep(train_outputs, each=n)) |> 
  mutate(trial=sample(seq(1,n*length(train_vec)), replace=FALSE)) |> 
  arrange(trial)

train_dat <- alm.sim(dat=td,
   c=.5, lr=.01, input_layer, output_layer, weight_mat)
train_dat$d

# test alt_exam on extrapolation items x=1 & x=2. 
map_dbl(c(1,2,4,5), ~ alt_exam(.x, c=1, input_layer, output_layer, train_dat$wm, train_vec, train_dat$d$almResp))

















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

