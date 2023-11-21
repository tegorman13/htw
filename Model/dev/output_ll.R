pacman::p_load(dplyr,tidyr,purrr,tibble)

input_layer <- c(100, 350, 600, 800, 1000, 1200)
output_layer <- c(100, 350, 600, 800, 1000, 1200)

ds <- readRDS(here::here("data/e1_md_11-06-23.rds")) |> mutate(sbj=id) |> relocate(sbj,.after=id) |>
  mutate(y=ifelse(y>1600,1600,y))
split_data <- split(ds, ds$id)
d1 <- split_data[[1]]

alm.response <- function(input = 1, c, input.layer, output.layer,weight.mat, trainVec=NULL) {
  input.activation <- exp(-c * (input.layer - input)^2) / (sum(exp(-c * (input.layer - input)^2)) +.0001)
  output.activation <- (weight.mat %*% input.activation) 
  output.probability <- output.activation / (sum(output.activation) +.0001)
  mean.response <- sum(output.layer * output.probability)
  list(mean.response = mean.response, input.activation = input.activation, output.activation = output.activation,op=output.probability)
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

alm.update <- function(corResp, c, lr, output.layer, input.activation, output.activation, weight.mat) {
  fz <- exp(-c * (output.layer - corResp)^2)
  teacherSignal <- (fz - output.activation) * lr
  weight.mat + ( teacherSignal %*% t(input.activation))
  #pmax(weight.mat,0)
  
}





alm.sim <- function(train, c, lr, input.layer = input_layer, output.layer = output_layer) {
  weight.mat <- matrix(0.000001, nrow = length(output.layer), ncol = length(input.layer))
  st <- numeric(nrow(train)) 
  for(i in 1:nrow(train)) {
    alm_resp <- alm.response(train$x[i], c, input_layer, output_layer, weight.mat)
    weight.mat <- alm.update(train$y[i], c, lr, output_layer, alm_resp$input.activation, alm_resp$output.activation, weight.mat)
    st[i] <- alm_resp$mean.response
  }
  return(list(d = cbind(train,almResp=st,dev=train$y-st) , wm = weight.mat, c = c, lr = lr))
}




alm_nll <- function(par, data,add_exam_nll=FALSE,test=NULL,delta=NULL) {
  c <- par[1]  # Sensitivity parameter
  lr <- par[2] # Learning rate parameter
  weight.mat <- matrix(0.000001, nrow = length(output_layer), ncol = length(input_layer))
  total_nll <- 0
  print("start training")
  for (i in 1:nrow(data)) {
    input <- data$x[i]
    corResp <- data$y[i]
    alm_resp <- alm.response(input, c, input_layer, output_layer, weight.mat)
    mean_response = sum(alm_resp$op*output_layer)
    variance_response = sum(alm_resp$op * (output_layer - mean_response)^2)
    if(variance_response<0) {print(c); print(lr)}
    observed_prob = dnorm(corResp, mean_response, sqrt(variance_response)) +.001
    print(paste0(corResp," ; ", mean_response," ; ",variance_response," ; ", observed_prob))
    total_nll <- total_nll + -log(observed_prob + 1e-5)
    weight.mat <- alm.update(corResp, c, lr, output_layer, alm_resp$input.activation, alm_resp$output.activation, weight.mat)
  }
  
  if(add_exam_nll){
    total_nll = total_nll + 10*exam_nll(test,c=c,input_layer,output_layer,weight.mat,trainVec=sort(unique(train$x)),delta)
  }
  return(total_nll)
}


exam_nll <- function(data, c, input.layer = INPUT_LAYER_DEFAULT, output.layer = OUTPUT_LAYER_DEFAULT, weight.mat, trainVec,delta=.000001 ) {
  total_log_likelihood <- 0
  
  for (i in 1:nrow(data)) {
    input <- data$x[i]
    corResp <- data$y[i]
    nearestTrain <- trainVec[which.min(abs(input - trainVec))]
    alm_out <- alm.response(nearestTrain, c, input.layer, output.layer, weight.mat)
    aresp <- alm_out$mean.response
    xUnder <- ifelse(min(trainVec) == nearestTrain, nearestTrain, trainVec[which(trainVec == nearestTrain) - 1])
    xOver <- ifelse(max(trainVec) == nearestTrain, nearestTrain, trainVec[which(trainVec == nearestTrain) + 1])
    mUnder <- alm.response(xUnder, c, input.layer, output.layer, weight.mat)$mean.response
    mOver <- alm.response(xOver, c, input.layer, output.layer, weight.mat)$mean.response
    exam.output <- round(aresp + ((mOver - mUnder) / (xOver - xUnder)) * (input - nearestTrain), 3)
    # Compute P[Y | X_i] for each X_i
    response_prob <- exp(-delta * (corResp - exam.output)^2)
    # Sum over M to get P[y_t | X(t)]
    P_yt_given_Xt <- sum(alm_out$input.activation * response_prob)  # aresp is P[X_i | X(t)]
    
    # Compute log likelihood for this trial and accumulate
    total_log_likelihood <- total_log_likelihood + -log(P_yt_given_Xt + 1e-5)
  }
  print(total_log_likelihood)
  return(total_log_likelihood)
}

weight.mat <- matrix(runif(nrow * ncol, min = 0.1, max = 1), nrow = nrow, ncol = ncol)


predict_full <- function(dat, c, lr, input.layer = input_layer, output.layer = output_layer) {
  weight.mat <- matrix(0.000001, nrow = length(output.layer), ncol = length(input.layer))
  train <- dat$train
  test <- dat$test
  st <- numeric(nrow(train)) 
  for(i in 1:nrow(train)) {
    alm_resp <- alm.response(train$x[i], c, input_layer, output_layer, weight.mat)
    weight.mat <- alm.update(train$y[i], c, lr, output_layer, alm_resp$input.activation, alm_resp$output.activation, weight.mat)
    st[i] <- alm_resp$mean.response
  }
  
  
  train = train |> mutate(ALM=st,dev=y-ALM,c=c,lr=lr) |> 
    rename(Observed=y) |> 
    pivot_longer(c("Observed","ALM"),names_to = "Resp",values_to = "vx") 
  
  trainVec=sort(unique(train$x))
  test_extrap <- map_dbl(test$x, ~exam.response(.x,c,input.layer,output.layer,weight.mat,trainVec))
  test_alm <- map_dbl(test$x, ~alm.response(.x,c,input.layer,output.layer,weight.mat)$mean.response)
  test <- test |> mutate(c=c,lr=lr,ALM=test_alm,EXAM=test_extrap) |> rename(Observed=y) |> 
    pivot_longer(c("Observed","ALM","EXAM"),names_to = "Resp",values_to = "vx") 
  
  return(tibble::lst(train,test,test_extrap,test_alm, wm = weight.mat, c = c, lr = lr))
}

fit_ex_testOnly <- function(par, dat=dat_params)
{
  c <- par[1]  # Sensitivity parameter
  lr <- par[2] # Learning rate parameter
  delta=par[3]
  train <- dat$train
  test <- dat$test
  
  train_alm <- alm.sim(train,c,lr,input_layer,output_layer)
  exam_nll(test,c=c,input_layer,output_layer,train_alm$wm,trainVec=sort(unique(train$x)),delta)
}

fit_ex_testTrain <- function(par, dat=dat_params)
{
  c <- par[1]  # Sensitivity parameter
  lr <- par[2] # Learning rate parameter
  delta=par[3]
  train <- dat$train
  test <- dat$test
  alm_nll(c(c,lr), train,TRUE,test,delta)
}



init_params <- c(.08, .5,.0001) 
lower_bounds <- c(1e-6, 1e-2,1e-12)   
upper_bounds <- c(7, 10,9)


dat_params <- tibble::lst(train = d1 |> filter(expMode2=="Train",tr>20), test = d1 |> filter(expMode2=="Test"))

ex_testOnly <- optim(par = init_params, fn = fit_ex_testOnly, dat=dat_params,
              method = 'L-BFGS-B', lower = lower_bounds, upper = upper_bounds,control=list(trace=2))


pf <- predict_full(dat_params,ex_testOnly$par[1], ex_testOnly$par[2])


pf$test |> 
  ggplot(aes(x=x,y=vx,fill=Resp))+stat_summary(geom="bar",fun="mean",position=position_dodge()) + 
  stat_summary(geom="errorbar",fun.data="mean_se",position=position_dodge())


pf$train |> 
  ggplot(aes(x=tr,y=vx,col=Resp))+stat_summary(geom="line",fun="mean") + 
  stat_summary(geom="errorbar",fun.data="mean_se") +facet_wrap(~x)

weight_plot(pf$wm)


ex_testTrain <- optim(par = init_params, fn = fit_ex_testTrain, dat=dat_params,
                     method = 'L-BFGS-B', lower = lower_bounds, upper = upper_bounds,control=list(trace=2))



pf <- predict_full(dat_params,ex_testTrain$par[1], ex_testTrain$par[2])

pf$test |> 
  ggplot(aes(x=x,y=vx,fill=Resp))+stat_summary(geom="bar",fun="mean",position=position_dodge()) + 
  stat_summary(geom="errorbar",fun.data="mean_se",position=position_dodge())


pf$train |> 
  ggplot(aes(x=tr,y=vx,col=Resp))+stat_summary(geom="line",fun="mean") + 
  stat_summary(geom="errorbar",fun.data="mean_se") +facet_wrap(~x)


pf$train |> 
  ggplot(aes(x=x,y=vx,col=Resp))+stat_summary(geom="bar",fun="mean") + 
  stat_summary(geom="errorbar",fun.data="mean_se") +facet_wrap(~x)

weight_plot(pf$wm)








pf <- predict_full(dat_params,.0006, .1)

pf$test |> 
  ggplot(aes(x=x,y=vx,fill=Resp))+stat_summary(geom="bar",fun="mean",position=position_dodge()) + 
  stat_summary(geom="errorbar",fun.data="mean_se",position=position_dodge())














fit1 <- optim(par = init_params, fn = alm_neg_log_likelihood, data = ,
             method = 'L-BFGS-B', lower = lower_bounds, upper = upper_bounds,control=list(trace=2))

sim_best <- alm.sim(d1,fit1$par[1],fit1$par[2])





sim_best$d |> ggplot(aes(x=tr,y=almResp))+geom_line()+geom_line(aes(y=y),col="red")+facet_wrap(~x)
sim_best$d |> ggplot(aes(x=almResp))+geom_histogram()+facet_wrap(~x)

sim_best$d |> pivot_longer(c(y,almResp),names_to="Resp",values_to = "vx") |> 
  ggplot(aes(x=x,y=vx,fill=Resp))+stat_summary(geom="bar",fun="mean",position=position_dodge()) + 
  stat_summary(geom="errorbar",fun.data="mean_se",position=position_dodge())



weight_plot(sim_best$wm)







trainVec=sort(unique(train$x))
testVec <- c(100,350,600,800,1000,1200) [!(c(100,350,600,800,1000,1200) %in% train$x)]

test_extrap <- map_dbl(testVec, ~exam.response(.x,c,input.layer,output.layer,weight.mat,trainVec))
test_train <- map_dbl(trainVec, ~alm.response(.x,c,input.layer,output.layer,weight.mat)$mean.response)




weight_plot(weight.mat)



fit2 <- optim(par = init_params, fn = alm_neg_log_likelihood, data = d1,
              method = 'Nelder-Mead',control=list(trace=3))

fit3 <- optim(par = init_params, fn = alm_neg_log_likelihood, data = d1,
              method = 'SANN',control=list(trace=3))


# Calculate the probability of the observed response
#observed_prob <- alm_resp$op[which.min(abs(corResp - output_layer))]
# sum(alm_resp$op*output_layer)
# 1/sum(alm_resp$op^2)


init_params <- c(.008, .5) 
lower_bounds <- c(1e-8, 1e-4)   
upper_bounds <- c(5, 5)

exam_neg_log_likelihood <- function(par, data, trainVec, delta) {
  c <- par[1]  # Sensitivity parameter
  weight.mat <- par[2:(length(input_layer) * length(output_layer) + 1)]  # Weight matrix parameters
  weight.mat <- matrix(weight.mat, nrow = length(output_layer), ncol = length(input_layer))
  
  total_neg_log_likelihood <- 0
  
  for (i in 1:nrow(data)) {
    input <- data$x[i]
    corResp <- data$y[i]
    
    # Calculate the probability for each training value
    prob_sum <- 0
    for (j in 1:length(trainVec)) {
      input_activation <- exp(-c * (trainVec[j] - input)^2) / sum(exp(-c * (trainVec - input)^2))
      exam_response <- exam.response(trainVec[j], c, input.layer, output.layer, weight.mat, trainVec)
      response_prob <- exp(-delta * (corResp - exam_response)^2)
      prob_sum <- prob_sum + input_activation * response_prob
    }
    
    # Add to the total negative log likelihood
    total_neg_log_likelihood <- total_neg_log_likelihood - log(prob_sum + 1e-5)
  }
  
  return(total_neg_log_likelihood)
}
