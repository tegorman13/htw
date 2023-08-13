
input.activation <- function(x.target, c, noise_sd) {
  noise <- rnorm(length(inputNodes), mean = 0, sd = noise_sd)
  return(exp((-1 * c) * (x.target - inputNodes + noise)^2))
}

output.activation <- function(x.target, weights, c, noise_sd) {
  noise <- rnorm(length(outputNodes), mean = 0, sd = noise_sd)
  return(weights %*% input.activation(x.target, c, noise_sd) + noise)
}

mean.prediction <- function(x.target, weights, c, noise_sd) {
  probability <- output.activation(x.target, weights, c, noise_sd) / sum(output.activation(x.target, weights, c, noise_sd))
  return(outputNodes %*% probability) # integer prediction
}

update.weights <- function(x.new, y.new, weights, c, lr, noise_sd) {
  y.feedback.activation <- exp(-1 * c * (y.new - outputNodes)^2)
  x.feedback.activation <- output.activation(x.new, weights, c, noise_sd)
  return(weights + lr * (y.feedback.activation - x.feedback.activation) %*% t(input.activation(x.new, c, noise_sd)))
}

train.alm <- function(dat, c = 0.05, lr = 0.5, noise_sd = 0, weights) {
  alm.train <- rep(NA, nrow(dat))
  for (i in 1:nrow(dat)) {
    weights <- update.weights(dat$input[i], dat$vx[i], weights, c, lr, noise_sd)
    resp = mean.prediction(dat$input[i], weights, c, noise_sd)
    alm.train[i] = resp
    weights[weights < 0] = 0
  }
  list(almTrain = alm.train, weights = weights)
}

sim_train <- function(dat, c = 0.5, lr = 0.2, inNodes = 7, outNodes = 32, noise_sd = 0) {
  inputNodes <<- seq(1, 7, length.out = inNodes * 1)
  outputNodes <<- seq(50, 1600, length.out = outNodes * 1)
  wm = matrix(.000000001, nrow = length(outputNodes), ncol = length(inputNodes))
  tt <- train.alm(dat, c, lr, noise_sd, wm)
}


# loss functions
MAE <- function(x, y) {mean(abs(x - y))}
RMSE <- function(x,y){sqrt(mean((x-y)^2))}