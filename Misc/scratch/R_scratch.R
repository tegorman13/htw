



k="c(5, 6, 7)"
k
k2="5, 6, 7"

k2

# insert "." between each number in k2 (and remove commas and spaces)
gsub(" ", ".", gsub(",", "", k2))


library(RSNNS)

basePath <- ("./")

#inputs <- as.matrix(seq(0,1,0.01))
#outputs <- as.matrix(sin(x))

set.seed(2)

inputs <- as.matrix(seq(0,10,0.1))
outputs <- as.matrix(sin(inputs) + runif(inputs*0.2))
#outputs <- as.matrix(sin(inputs))
outputs <- normalizeData(outputs, "0_1")

numHiddenUnits <- 40

snnsObject <- SnnsRObjectFactory()

snnsObject$setLearnFunc('RadialBasisLearning')
snnsObject$setUpdateFunc('Topological_Order')
snnsObject$setUnitDefaults(0,0,1,0,1,'Act_RBF_Gaussian','Out_Identity')
snnsObject$createNet(c(ncol(inputs),numHiddenUnits,ncol(outputs)), TRUE)

snnsObject$setTTypeUnitsActFunc("UNIT_INPUT", "Act_Identity")
snnsObject$setTTypeUnitsActFunc("UNIT_HIDDEN", "Act_RBF_Gaussian")
snnsObject$setTTypeUnitsActFunc("UNIT_OUTPUT", "Act_IdentityPlusBias")
#snnsObject$setTTypeUnitsActFunc("UNIT_OUTPUT", "Act_Logistic")

patset <- snnsObject$createPatSet(inputs, outputs)
snnsObject$setCurrPatSet(patset$set_no)

snnsObject$shufflePatterns(TRUE)
snnsObject$DefTrainSubPat()

snnsObject$initializeNet(c(0,0,0,0,0), "RBF_Weights_Kohonen")

snnsObject$initializeNet(c(0, 1, 0, 0.01, 0.01), "RBF_Weights")

parameters <- c(1e-8, 0, 1e-8, 0.1, 0.8)
maxit <- 1000

error <- vector()
for(i in 1:maxit) {
  res <- snnsObject$learnAllPatterns(parameters)
  error[i] <- res[[2]]
}

par(mfrow=c(2,1))

plot(error, type="l")

predictions <- snnsObject$predictCurrPatSet("output", c(0))

plot(inputs, outputs)
lines(inputs, predictions, col="green")
# 
# snnsObject$saveNet(paste(basePath,"rbf_sinSnnsR.net",sep=""),"rbf_sinSnnsR")
# snnsObject$saveNewPatterns(paste(basePath,"rbf_sinSnnsR.pat",sep=""), patset$set_no);









#https://mpopov.com/blog/2021/02/28/animation-of-optimization-in-torch/
library(torch)
library(gganimate)
library(tidyverse)


f <- function(x) (6 * x - 2) ^ 2 * sin(12 * x - 4)

adam_iters <- (function(n_iters, learn_rate) {
  
  x <- torch_zeros(1, requires_grad = TRUE)
  
  f <- function(x) (6 * x - 2)^2 * torch_sin(12 * x - 4)
  
  optimizer <- optim_adam(x, lr = learn_rate)
  
  iters <- tibble(
    iter = 1:n_iters,
    
    x = replicate(n_iters, {
      # Evaluate at current value of x:
      y <- f(x)
      # Zero out the gradients before the backward pass:
      optimizer$zero_grad()
      # Compute gradient on evaluation tensor:
      y$backward()
      # Update value of x:
      optimizer$step()
      # Remember updated value of x:
      as.numeric(x)
    })
  )
  
  # Add starting value and return:
  bind_rows(tibble(iter = 0, x = 0), iters)
  
})(n_iters = 50, learn_rate = 0.25)


ggplot(adam_iters) +
  geom_function(fun = f, size = 1, n = 100) +
  geom_point(aes(x = x, y = f(x)), size = 5)

anim <- ggplot(adam_iters) +
  geom_function(fun = f, size = 1, n = 100) +
  geom_point(aes(x = x, y = f(x)), size = 5) +
  transition_manual(iter)

anim <- anim +
  scale_y_continuous(name = NULL, breaks = NULL, minor_breaks = NULL) +
  scale_x_continuous(name = NULL, breaks = NULL, minor_breaks = NULL)

anim



#https://torch.mlverse.org/docs/articles/examples/basic-nn-module.html
library(torch)

rbf <- nn_module(
  clasname = "rbf",
  initialize = function(in_features, out_features) {
    self$centers <- nn_parameter(torch_randn(out_features, in_features))
    self$betas <- nn_parameter(torch_randn(out_features))
    self$weights <- nn_parameter(torch_randn(out_features))
  },
  forward = function(x) {
    # calculate the distance between x and each center
    distances <- torch_norm(x - self$centers, p = 2, dim = 2)
    # apply the radial basis function to the distances
    activations <- torch_exp(-self$betas * distances)
    # multiply the activations by the weights and sum them
    torch_sum(self$weights * activations, dim = 1, keepdim = FALSE)
  }
)

model <- rbf(1, 40)
model <- nn_sequential(
  rbf(1, 64),
  nn_linear(64, 1)
)

model <- nn_sequential(
  nn_linear(1, 64),
  rbf(64,64),
  nn_linear(64, 1)
)

model$parameters
model$weights

x <- torch_randn(1, 1)
y_pred <- model(x)
y_pred




#https://jonnylaw.rocks/posts/2021-02-02-neural-networks-in-r/
n <- 100
x <- seq(-5, 5, length.out=n)
y <- 4 * sin(x)
y_obs <- 4 * sin(x) + rnorm(n, sd = 1)
non_linear <- tibble(x, y_obs) %>% sample_n(20)
observed <- tibble(x, y) %>% sample_n(50)
loss <- function(y, pred) {
  (y - pred)$pow(2)$sum()
}


qplot(x, y, geom = "line") +
  geom_point(data = non_linear, aes(x = x, y = y_obs, colour = "Observed")) +
  theme(legend.position = "none")




loss <- function(y, pred) {
  (y - pred)$pow(2)$sum()
}

# train_torch <- function(x, y, model, loss, epochs, learning_rate) {
#   x_in <- x
#   optimiser <- optim_adam(model$parameters, lr = learning_rate)
#   for (i in seq_len(epochs)) {
#     pred <- model(x_in)
#     losses <- loss(y, pred)
#     model$zero_grad()
#     losses$backward()
#     if (i %% 10 == 0)
#       cat("Epoch: ", i, "   Loss: ", losses$item(), "\n")
#     optimiser$step()
#   }
#   model
# }

#trained_model <- train_torch(scale(observed$x), scale(observed$y), model, loss, 200)


trained_model <-
  train_torch(torch_tensor(as.matrix(non_linear$x)),
              torch_tensor(as.matrix(non_linear$y_obs)),
              model,
              nnf_mse_loss,
              100, 
              0.1)



model <- rbf(1, 20)
# Convert the data to tensors
x_obs <- torch_tensor(non_linear$x, requires_grad = FALSE)
y_obs <- torch_tensor(non_linear$y_obs, requires_grad = FALSE)
x_test <- torch_tensor(observed$x, requires_grad = FALSE)
y_test <- torch_tensor(observed$y, requires_grad = FALSE)

# Train the model
optimizer <- optim_adam(model$parameters, lr = 0.01)
for (i in 1:1000) {
  optimizer$zero_grad()
  output <- model(x_obs)
  loss_val <- loss(y_obs, output)
  loss_val$backward()
  optimizer$step()
}
output_test <- model(x_test)
test_loss <- loss(y_test, output_test)$item()
cat("Test loss: ", test_loss, "\n")



df <- tibble(non_linear$x, prediction = as_array(model(torch_tensor(x)$unsqueeze(-1))))
ggplot() +
  geom_point(data = non_linear, aes(x = x, y = y_obs), color = "blue") +
  geom_line(data = tibble(x = x, y = y), aes(x = x, y = y), color = "black") +
  geom_line(data = df, aes(x = x, y = y), color = "red")



predictions <- tibble(x, prediction = as_array(trained_model(torch_tensor(x)$unsqueeze(-1))))
qplot(x, y, geom = "line", colour = "truth") +
  geom_line(data = predictions, aes(x, prediction, colour = "fitted")) +
  geom_point(data = non_linear, aes(x, y_obs, colour = "observed")) +
  labs(colour = "") +
  theme(legend.position = c(0.9, 0.9))




library(dplyr)

fun_df <- tibble(cond = c("cond1", "cond2"), 
                 fun = c(function(x, p) (1 + exp(-p[2] * (x - p[1])))^(-1), 
                         function(x, p) (1 + exp(-p[2] * (x - p[3])))^(-1)))






