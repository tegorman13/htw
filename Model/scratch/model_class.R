# Load required packages
pacman::p_load(dplyr, tidyr, purrr, tibble)

# Define input and output layers
input_layer <- c(100, 350, 600, 800, 1000, 1200)
output_layer <- c(100, 350, 600, 800, 1000, 1200)



# Function to calculate exam response for a given input
exam.response <- function(input, c, input.layer = input_layer, output.layer = output_layer, weight.mat, trainVec) {
    nearestTrain <- trainVec[which.min(abs(input - trainVec))]
    aresp <- alm.response(nearestTrain, c, input.layer, output.layer, weight.mat)$mean.response
    
    xUnder <- ifelse(min(trainVec) == nearestTrain, nearestTrain, trainVec[which(trainVec == nearestTrain) - 1])
    xOver <- ifelse(max(trainVec) == nearestTrain, nearestTrain, trainVec[which(trainVec == nearestTrain) + 1])
    
    mUnder <- alm.response(xUnder, c, input.layer, output.layer, weight.mat)$mean.response
    mOver <- alm.response(xOver, c, input.layer, output.layer, weight.mat)$mean.response
    
    exam.output <- round(aresp + ((mOver - mUnder) / (xOver - xUnder)) * (input - nearestTrain), 3)
    exam.output
}




# Function to calculate activation and response for a given input
alm.response <- function(input, c, input.layer, output.layer, weight.mat, trainVec = NULL) {
    input.activation <- exp(-c * (input.layer - input)^2) / (sum(exp(-c * (input.layer - input)^2)) + .0001)
    output.activation <- (weight.mat %*% input.activation) 
    output.probability <- output.activation / (sum(output.activation) + .0001)
    mean.response <- sum(output.layer * output.probability)
    list(mean.response = mean.response, input.activation = input.activation, output.activation = output.activation, op = output.probability)
}

# Function to update weight matrix based on correct response
alm.update <- function(corResp, c, lr, output.layer, input.activation, output.activation, weight.mat) {
    fz <- exp(-c * (output.layer - corResp)^2)
    teacherSignal <- (fz - output.activation) * lr
    weight.mat + (teacherSignal %*% t(input.activation))
}

# Function to simulate ALM learning
alm.sim <- function(train, c, lr, input.layer = input_layer, output.layer = output_layer) {
    weight.mat <- matrix(0.000001, nrow = length(output.layer), ncol = length(input.layer))
    st <- numeric(nrow(train)) 
    for(i in 1:nrow(train)) {
        alm_resp <- alm.response(train$x[i], c, input.layer, output.layer, weight.mat)
        weight.mat <- alm.update(train$y[i], c, lr, output.layer, alm_resp$input.activation, alm_resp$output.activation, weight.mat)
        st[i] <- alm_resp$mean.response
    }
    return(list(d = cbind(train, almResp = st, dev = train$y - st), wm = weight.mat, c = c, lr = lr))
}


# Function to calculate negative log likelihood for ALM learning
alm_nll <- function(par, data, add_exam_nll = FALSE, test = NULL, delta = NULL) {
    c <- par[1]  # Sensitivity parameter
    lr <- par[2] # Learning rate parameter
    weight.mat <- matrix(0.000001, nrow = length(output_layer), ncol = length(input_layer))
    total_nll <- 0
    for (i in 1:nrow(data)) {
        input <- data$x[i]
        corResp <- data$y[i]
        alm_resp <- alm.response(input, c, input_layer, output_layer, weight.mat)
        mean_response <- sum(alm_resp$op * output_layer)
        variance_response <- sum(alm_resp$op * (output_layer - mean_response)^2)
        observed_prob <- dnorm(corResp, mean_response, sqrt(variance_response)) + .001
        total_nll <- total_nll + -log(observed_prob + 1e-5)
        weight.mat <- alm.update(corResp, c, lr, output_layer, alm_resp$input.activation, alm_resp$output.activation, weight.mat)
    }
    
    if (add_exam_nll) {
        total_nll <- total_nll + 10 * exam_nll(test, c = c, input_layer, output_layer, weight.mat, trainVec = sort(unique(train$x)), delta)
    }
    return(total_nll)
}

# Function to calculate negative log likelihood for exam response
exam_nll <- function(data, c, input.layer = input_layer, output.layer = output_layer, weight.mat, trainVec, delta = .000001) {
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
        response_prob <- exp(-delta * (corResp - exam.output)^2)
        P_yt_given_Xt <- sum(alm_out$input.activation * response_prob)
        total_log_likelihood <- total_log_likelihood + -log(P_yt_given_Xt + 1e-5)
    }
    return(total_log_likelihood)
}


# Function to fit the model and calculate negative log likelihood for test data only
fit_ex_testOnly <- function(par, dat = dat_params) {
    c <- par[1]  # Sensitivity parameter
    lr <- par[2] # Learning rate parameter
    delta <- par[3]
    train <- dat$train
    test <- dat$test
    
    train_alm <- alm.sim(train, c, lr, input_layer, output_layer)
    exam_nll(test, c = c, input_layer, output_layer, train_alm$wm, trainVec = sort(unique(train$x)), delta)
}

# Function to fit the model and calculate negative log likelihood for both train and test data
fit_ex_testTrain <- function(par, dat = dat_params) {
    c <- par[1]  # Sensitivity parameter
    lr <- par[2] # Learning rate parameter
    delta <- par[3]
    train <- dat$train
    test <- dat$test
    alm_nll(c(c, lr), train, TRUE, test, delta)
}


# Define initial parameters, lower bounds, and upper bounds
init_params <- c(.08, .5, .0001) 
lower_bounds <- c(1e-6, 1e-2, 1e-12)   
upper_bounds <- c(7, 10, 9)

# Define data parameters
dat_params <- tibble::lst(train = d1 |> filter(expMode2 == "Train", tr > 20), test = d1 |> filter(expMode2 == "Test"))

# Fit the model and calculate negative log likelihood for test data only
ex_testOnly <- optim(par = init_params, fn = fit_ex_testOnly, dat = dat_params,
                                         method = 'L-BFGS-B', lower = lower_bounds, upper = upper_bounds, control = list(trace = 2))


# Function to predict responses for train and test data
predict_full <- function(dat, c, lr, input.layer = input_layer, output.layer = output_layer) {
    weight.mat <- matrix(0.000001, nrow = length(output.layer), ncol = length(input.layer))
    train <- dat$train
    test <- dat$test
    st <- numeric(nrow(train)) 
    for(i in 1:nrow(train)) {
        alm_resp <- alm.response(train$x[i], c, input.layer, output.layer, weight.mat)
        weight.mat <- alm.update(train$y[i], c, lr, output.layer, alm_resp$input.activation, alm_resp$output.activation, weight.mat)
        st[i] <- alm_resp$mean.response
    }
    
    train <- train |> mutate(ALM = st, dev = y - ALM, c = c, lr = lr) |> 
        rename(Observed = y) |> 
        pivot_longer(c("Observed", "ALM"), names_to = "Resp", values_to = "vx") 
    
    trainVec <- sort(unique(train$x))
    test_extrap <- map_dbl(test$x, ~exam.response(.x, c, input.layer, output.layer, weight.mat, trainVec))
    test_alm <- map_dbl(test$x, ~alm.response(.x, c, input.layer, output.layer, weight.mat)$mean.response)
    test <- test |> mutate(c = c, lr = lr, ALM = test_alm, EXAM = test_extrap) |> 
        rename(Observed = y) |> 
        pivot_longer(c("Observed", "ALM", "EXAM"), names_to = "Resp", values_to = "vx") 
    
    return(tibble::lst(train, test, test_extrap, test_alm, wm = weight.mat, c = c, lr = lr))
}


# Predict responses for train and test data
pf <- predict_full(dat_params, ex_testOnly$par[1], ex_testOnly$par[2])

# Plot the predicted responses for test data
pf$test |> 
    ggplot(aes(x = x, y = vx, fill = Resp)) + 
    stat_summary(geom = "bar", fun = "mean", position = position_dodge()) + 
    stat_summary(geom = "errorbar", fun.data = "mean_se", position = position_dodge())

# Plot the predicted responses for train data
pf$train |> 
    ggplot(aes(x = tr, y = vx, col = Resp)) + 
    stat_summary(geom = "line", fun = "mean") + 
    stat_summary(geom = "errorbar", fun.data = "mean_se") +
    facet_wrap(~x)


