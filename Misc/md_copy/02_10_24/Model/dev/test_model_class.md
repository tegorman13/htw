
library(testthat)
library(tibble)

# Set up some test parameters
test_inp_layer <- c(1)
test_out_layer <- c(2)
test_c_range <- seq(from = 0.1, to = 1, by = 0.1)
test_lr_range <- seq(from = 0.2, to = 1, by = 0.2)

tests <- list()

# First test: check that the function returns a list
test1 <- function() {
  result <- generate_simulated_data(n_agents = 1, n_train_trials = 1, trainVec=c(1), 
                                    input_layer=test_inp_layer, output_layer=test_out_layer)
  testthat::expect_is(result, "list")
}
tests[["Test #1: Should return a list"]] <- test1

# Second test: check the length of output list equals to 'n_agents'
test2 <- function() {
  n_agents_trial <- 10
  result <- generate_simulated_data(n_agents = n_agents_trial, n_train_trials = 1, trainVec=c(1), 
                                    input_layer=test_inp_layer, output_layer=test_out_layer)
  testthat::expect_equal(length(result), n_agents_trial)
}
tests[["Test #2: Returns a list of correct length"]] <- test2

# Third test: check that every item in the list is a data frame
test3 <- function() {
  n_agents_trial <- 3
  result <- generate_simulated_data(n_agents = n_agents_trial, n_train_trials = 1, trainVec=c(1), 
                                    input_layer=test_inp_layer, output_layer=test_out_layer)
  lapply(result, function(x) testthat::expect_s3_class(x, "data.frame"))
}
tests[["Test #3: Each list item is a data frame"]] <- test3

# Run the tests
for (name in names(tests)) {
  cat("\n", name, "\n")
  tests[[name]]()
}












library(testthat)

# Test case 1: alm_nll function
test_that("alm_nll calculates negative log likelihood correctly", {
  # Create sample data
  data <- data.frame(x = c(1, 2, 3), y = c(0.5, 0.6, 0.7))
  input_layer <- c(0.1, 0.2, 0.3)
  output_layer <- c(0.4, 0.5, 0.6)
  par <- c(0.8, 0.9)
  
  # Call alm_nll function
  nll <- alm_nll(par, data, add_exam_nll = FALSE)
  
  # Check if the calculated negative log likelihood is correct
  expect_equal(nll, 6.907755)
})

# Test case 2: exam_nll function
test_that("exam_nll calculates negative log likelihood correctly", {
  # Create sample data
  data <- data.frame(x = c(1, 2, 3), y = c(0.5, 0.6, 0.7))
  input_layer <- c(0.1, 0.2, 0.3)
  output_layer <- c(0.4, 0.5, 0.6)
  weight_mat <- matrix(0.000001, nrow = length(output_layer), ncol = length(input_layer))
  trainVec <- c(0.9, 1, 1.1)
  c <- 0.8
  
  # Call exam_nll function
  nll <- exam_nll(data, c, input.layer = input_layer, output.layer = output_layer, weight.mat = weight_mat, trainVec = trainVec)
  
  # Check if the calculated negative log likelihood is correct
  expect_equal(nll, 6.907755)
})

# Run the tests
testthat::test_report()


# Define a simple linear regression model
simple_model <- function(x, par) {
  par[1] * x + par[2]
}
# Generate some sample data
x <- 1:3
y <- simple_model(x, c(0.8, 0.9))
# Calculate the negative log likelihood
nll <- -sum(dnorm(y, mean = simple_model(x, c(0.8, 0.9)), sd = 1, log = TRUE))
# Print the calculated nll
print(nll)
