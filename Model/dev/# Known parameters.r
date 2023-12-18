

library(tibble)

# Known parameters
c_true <- 0.5
lr_true <- 0.5

input_layer <- seq(1,5,1)
output_layer <- input_layer
weight.mat <- matrix(0.000001, nrow = length(output_layer), ncol = length(input_layer))
diag(weight.mat) <- 1

# Simulate data
set.seed(123) # for reproducibility
# generate train data (20 observations from linear identity function)
train <-  tibble(x = rep(seq(1, 5, 1),3), y = x)
# 5 reps of each each training item


sim_data <- alm.sim(train, c = c_true, lr = lr_true)

c_values <- seq(0.1, 2, length.out = 5)


alm.sim <- function(train, c, lr, input.layer = input_layer, output.layer = output_layer) {
    weight.mat <- matrix(0.000001, nrow = length(output.layer), ncol = length(input.layer))
    diag(weight.mat) <- 1
    st <- numeric(nrow(train)) 
    for(i in 1:nrow(train)) {
        alm_resp <- alm.response(train$x[i], c, input.layer, output.layer, weight.mat)
        weight.mat <- alm.update(train$y[i], c, lr, output.layer, alm_resp$input.activation, alm_resp$output.activation, weight.mat)
        st[i] <- alm_resp$mean.response
    }
    return(list(d = cbind(train, almResp = st, dev = train$y - st), wm = weight.mat, c = c, lr = lr))
}

alm_rmsd <- function(par, data) {
    c <- par[1]  # Sensitivity parameter
    lr <- par[2] # Learning rate parameter
    weight.mat <- matrix(0.000001, nrow = length(output_layer), ncol = length(input_layer))
    diag(weight.mat) <- 1
    rmsd <- 0
    for (i in 1:nrow(data)) {
        input <- data$x[i]
        corResp <- data$y[i]
        alm_resp <- alm.response(input, c, input_layer, output_layer, weight.mat)
        mean_response <- sum(alm_resp$op * output_layer)
        rmsd <- rmsd + (corResp - mean_response)^2
        weight.mat <- alm.update(corResp, c, lr, output_layer, alm_resp$input.activation, alm_resp$output.activation, weight.mat)
    }
    return(sqrt(rmsd/nrow(data)))
}



# Fit the model for each c value
for (c_val in c_values) {
    # simulate model
    sim_data <- alm.sim(train, c = c_val, lr = lr_true)
    # Fit the model
    fit <- optim(c(c_val, lr_true), alm_rmsd, data = sim_data$d, method = "Nelder-Mead")

    # Estimated parameters
    c_est <- fit$par[1]
    lr_est <- fit$par[2]

    # Compare estimated parameters to true parameters
    cat("True c: ", c_val, "\n")
    cat("Estimated c: ", c_est, "\n")
    cat("True lr: ", lr_true, "\n")
    cat("Estimated lr: ", lr_est, "\n")
    cat("\n")
}



input.activation <- exp(-c * (input_layer - input)^2) / (sum(exp(-c * (input_layer - input)^2)) + .0001)
output.activation <- (weight.mat %*% input.activation) 

output.probability <- output.activation / (sum(output.activation) + .0001)
mean.response <- sum(output_layer * output.probability)
