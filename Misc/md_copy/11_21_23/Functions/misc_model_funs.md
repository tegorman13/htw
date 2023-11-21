

weight_plot <- function(weight.mat){
  tibbleFromMat <-
    weight.mat %>%
    tibble::as_tibble() %>%
    tibble::rownames_to_column("Var1") %>%
    tidyr::pivot_longer(-Var1, names_to = "Var2", values_to = "value") %>%
    dplyr::mutate(
      Output_Layer = factor(Var1, levels = 1:dim(weight.mat)[1]),
      Input_Layer = factor(gsub("V", "", Var2), levels = 1:dim(weight.mat)[2])
    )
  ggplot(tibbleFromMat, aes(Output_Layer, Input_Layer)) +
    geom_tile(aes(fill = value)) +
    geom_text(aes(label = round(value, 3))) +
    scale_fill_gradient(low = "white", high = "red")
}





round_tibble <- function(tbl, rn) {
  tbl %>% 
    mutate(across(where(is.numeric), ~round(., rn)))
}

strip_list_notation <- function(str) {
  # Remove 'list(' at the beginning
  str <- gsub("^list\\(", "", str)
  # Remove ')' at the end
  str <- gsub("\\)$", "", str)
  return(str)
}


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

nll <- function(obsv,pred,sigma)
{
 -sum(dnorm(obsv,mean=pred,sd=sigma,log=TRUE))
}


nll2 <- function(obsv,pred,sigma)
{
  nll= -sum(dnorm(obsv,mean=pred,sd=sigma,log=TRUE)) 
  #print(nll)
  if (is.nan(nll)) {
    nll <- 1e4 # Large penalty
  }
  return(nll)
}



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


