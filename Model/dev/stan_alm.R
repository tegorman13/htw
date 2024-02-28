
library(rstan)

stan_code_alm <- "
data {
  int<lower=0> N; // number of data points
  int<lower=0> M; // number of input nodes
  int<lower=0> L; // number of output nodes
  matrix[N, M] X; // input data
  int<lower=1, upper=L> Y[N]; // output data, assuming categorical with labels 1 through L
  real<lower=0> c; // constant for Gaussian function
  real<lower=0, upper=1> alpha; // learning rate
}

parameters {
  matrix[L, M] w; // weights from input nodes to output nodes
}

model {
  // Priors
  for (l in 1:L) {
    for (m in 1:M) {
      w[l, m] ~ normal(0, 1); // Assuming a simple Gaussian prior for each weight
    }
  }

  for (n in 1:N) {
    vector[L] logits;
    for (l in 1:L) {
      // Compute the weighted sum for each output node
      real weighted_sum = 0;
      for (m in 1:M) {
        weighted_sum += w[l, m] * X[n, m];
      }
      logits[l] = weighted_sum;
    }
    
    Y[n] ~ categorical_logit(logits); // Assuming Y is categorical
  }
    

}
"


# Generate data
set.seed(123)
N <- 100 # number of data points
M <- 2 # number of input nodes
L <- 3 # number of output nodes
X <- matrix(runif(N * M, -1, 1), nrow = N) # input data
Y <- sample(1:L, N, replace = TRUE) # output data
c <- 1 # constant for Gaussian function
alpha <- 0.1 # learning rate

# Prepare data for Stan
stan_data <- list(N = N, M = M, L = L, X = X, Y = Y, c = c, alpha = alpha)

# Compile and fit Stan model

model_alm <- stan_model(model_code = stan_code_alm)
fit_alm <- sampling(model_alm, data = stan_data)

# Generate test data
X_test <- matrix(runif(N * M, -1, 1), nrow = N) # test input data

# Make predictions on test data
pred_alm <- extract(fit_alm, permuted = TRUE)$m
Y_test <- apply(pred_alm, 1, which.max) # predicted output for test data

# Print accuracy on test data
cat("Accuracy on test data:", mean(Y_test == Y), "\n")
