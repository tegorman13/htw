


adjust_layer <- function(input.layer, k){
  if(k == 1) return(input.layer)
  new.layer <- c()
  for(i in 1:(length(input.layer) - 1)){
    new.layer <- c(new.layer, seq(input.layer[i], input.layer[i+1], length.out = k + 1))
  }
  new.layer <- unique(new.layer)
  return(new.layer)
}


### Loss Functions

#  mse <- mean((predictions - test_data$y)^2)

RMSE <- function(x,y){
 # print("rmseTrial")
  sqrt(mean((x-y)^2, na.rm=TRUE))
}

RMSE.blocked <- function(x,y,blocks=6){
  #print("rmseBlocked")
  data.table(x=x,y=y,t=seq(1,length(x))) %>% 
    .[, `:=`(fitBins = cut(t, breaks = ..blocks, labels = c(1:..blocks)))] %>%
    .[, .(predMean = mean(x), obsMean = mean(y)), keyby = .(fitBins)] %>%
    .[, RMSE(predMean,obsMean)] %>% as.numeric()
}

MAE <- function(x, y) {
  mean(abs(x - y))
}

MAPE <- function(x, y) {
  mean(abs((x - y) / y)) * 100
}

MedAE <- function(x, y) {
  median(abs(x - y))
}

HuberLoss <- function(x, y, delta = 1) {
  error <- x - y
  abs_error <- abs(error)
  loss <- ifelse(abs_error <= delta, 0.5 * error^2, delta * (abs_error - 0.5 * delta))
  mean(loss)
}


