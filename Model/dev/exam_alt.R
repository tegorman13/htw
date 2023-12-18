



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

alt_exam <- function(input, c, input.layer = input_layer, output.layer = output_layer, weight.mat, trainVecX, trainVecY) {
    if (length(input) > 1) {
        stop("Input must be scalar")
    }
    if (length(trainVecY) == 0) {
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
