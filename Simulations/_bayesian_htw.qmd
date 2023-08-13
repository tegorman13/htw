---
title: Bayesian Simulation
date-modified: last-modified
categories: [Simulation, ALM, R]
page-layout: full
--- 




To implement the Bayesian model, we first need to define the prior beliefs and likelihood function. The prior beliefs could be represented by a Gaussian distribution with a certain mean and standard deviation. The likelihood function would represent the probability of observing the training data given the current beliefs. We can update the beliefs using Bayes' rule, which combines the prior beliefs with the likelihood function to produce a posterior distribution.

Here's an implementation of the Bayesian model in R:


```{r}
pacman::p_load(tidyverse)

library(tidyverse)

# Define the prior beliefs
prior_beliefs <- function(mu, sigma) {
  list(mu = mu, sigma = sigma)
}

# Define the likelihood function
likelihood <- function(data, mu, sigma) {
  dnorm(data, mean = mu, sd = sigma)
}

# Update beliefs using Bayes' rule
update_beliefs <- function(prior, data) {
  likelihoods <- likelihood(data, prior$mu, prior$sigma)
  posterior_mu <- (sum(likelihoods * prior$mu) + prior$mu) / (sum(likelihoods) + 1)
  posterior_sigma <- sqrt((sum(likelihoods * (prior$sigma ^ 2)) + prior$sigma ^ 2) / (sum(likelihoods) + 1))
  
  prior_beliefs(posterior_mu, posterior_sigma)
}

# Generate training data
generate_training_data <- function(n, constant_training) {
  if (constant_training) {
    velocities <- rnorm(n, mean = 50, sd = 10)
  } else {
    velocities <- runif(n, min = 0, max = 100)
  }
  
  velocities
}

# Simulate the model for an artificial population
simulate_population <- function(n_subjects, n_trials, constant_training) {
  subjects <- vector("list", n_subjects)
  
  for (i in 1:n_subjects) {
    # Generate training data
    training_data <- generate_training_data(n_trials, constant_training)
    
    # Initialize prior beliefs
    beliefs <- prior_beliefs(mu = 50, sigma = 30)
    
    # Update beliefs based on the training data
    for (data_point in training_data) {
      beliefs <- update_beliefs(beliefs, data_point)
    }
    
    subjects[[i]] <- beliefs
  }
  
  subjects
}

# Simulate populations for constant and varied learners
n_subjects <- 100
n_trials <- 100
constant_learners <- simulate_population(n_subjects, n_trials, constant_training = TRUE)
varied_learners <- simulate_population(n_subjects, n_trials, constant_training = FALSE)

# Plot learning curves and generalization patterns
plot_beliefs <- function(population, title) {
  df <- tibble(
    subject = rep(1:length(population), each = 100),
    mu = rep(sapply(population, function(x) x$mu), each = 100),
    sigma = rep(sapply(population, function(x) x$sigma), each = 100),
    velocity = unlist(lapply(population, function(x) rnorm(100, mean = x$mu, sd = x$sigma)))
  )
  
  ggplot(df, aes(x = subject, y = velocity, group = subject)) +
    geom_line(alpha = 0.5) +
    geom_point(aes(color = mu), size = 3, alpha = 0.5) +
    scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = mean(df$mu)) +
    labs(title = title, x = "Subject", y = "Generated Velocity") +
    theme_minimal()
}

plot_beliefs(constant_learners, "Constant Learners")

```


We can modify the simulation to match the design of your experiment. I've implemented this by creating separate functions for each stage of the experiment and combining them within a single function for the entire simulation.

```{r}

library(tidyverse)

# Training stage
training_stage <- function(n_subjects, n_trials, constant_training, orig_condition) {
  subjects <- vector("list", n_subjects)
  
  for (i in 1:n_subjects) {
    if (orig_condition) {
      if (constant_training) {
        training_data <- rep(900, n_trials)
      } else {
        training_data <- c(rep(900, n_trials/3), rep(1100, n_trials/3), rep(1300, n_trials/3))
      }
    } else {
      if (constant_training) {
        training_data <- rep(700, n_trials)
      } else {
        training_data <- c(rep(700, n_trials/3), rep(450, n_trials/3), rep(200, n_trials/3))
      }
    }
    subjects[[i]] <- training_data
  }
  
  subjects
}

# Testing stage
testing_stage <- function(training_subjects, n_trials, feedback_type) {
  subjects <- vector("list", length(training_subjects))
  test_bands <- c(seq(100, 1400, by = 250))
  
  for (i in 1:length(training_subjects)) {
    training_data <- training_subjects[[i]]
    test_data <- vector("list", length(test_bands))
    
    for (j in 1:length(test_bands)) {
      if (feedback_type == "continuous") {
        test_data[[j]] <- rnorm(n_trials, mean = test_bands[j], sd = 100)
      } else { # ordinal
        test_data[[j]] <- round(rnorm(n_trials, mean = test_bands[j], sd = 100), -2)
      }
    }
    
    subjects[[i]] <- list(training = training_data, testing = test_data)
  }
  
  subjects
}

# Entire simulation
simulate_experiment <- function(n_subjects, n_training_trials, n_testing_trials, constant_training, orig_condition, feedback_type) {
  training_subjects <- training_stage(n_subjects, n_training_trials, constant_training, orig_condition)
  testing_subjects <- testing_stage(training_subjects, n_testing_trials, feedback_type)
  
  testing_subjects
}

# Simulate populations for varied and constant learners
n_subjects <- 100
n_training_trials <- 90
n_testing_trials <- 15

constant_orig_continuous <- simulate_experiment(n_subjects, n_training_trials, n_testing_trials, constant_training = TRUE, orig_condition = TRUE, feedback_type = "continuous")
constant_orig_ordinal <- simulate_experiment(n_subjects, n_training_trials, n_testing_trials, constant_training = TRUE, orig_condition = TRUE, feedback_type = "ordinal")
constant_rev_continuous <- simulate_experiment(n_subjects, n_training_trials, n_testing_trials, constant_training = TRUE, orig_condition = FALSE, feedback_type = "continuous")
constant_rev_ordinal <- simulate_experiment(n_subjects, n_training_trials, n_testing_trials, constant_training = TRUE, orig_condition = FALSE, feedback_type = "ordinal")

varied_orig_continuous <- simulate_experiment(n_subjects, n_training_trials, n_testing_trials, constant_training = FALSE, orig_condition = TRUE, feedback_type = "continuous")
varied_orig_ordinal <- simulate_experiment(n_subjects, n_training_trials, n_testing_trials, constant_training = FALSE, orig_condition = TRUE, feedback_type = "ordinal")
varied_rev_continuous <- simulate_experiment(n_subjects, n_training_trials, n_testing_trials, constant_training = FALSE, orig_condition = FALSE, feedback_type = "continuous")
varied_rev_ordinal <- simulate_experiment(n_subjects, n_training_trials, n_testing_trials, constant_training = FALSE, orig_condition = FALSE, feedback_type = "ordinal")

# Combine the simulation results
combine_results <- function(simulation, training_type, condition, feedback_type) {
  tibble(
    subject = rep(1:length(simulation), each = 15 * 6),
    training_type = rep(training_type, each = 15 * 6 * length(simulation)),
    condition = rep(condition, each = 15 * 6 * length(simulation)),
    feedback_type = rep(feedback_type, each = 15 * 6 * length(simulation)),
    test_band = rep(rep(seq(100, 1400, by = 250), each = 15), length(simulation)),
    generated_velocity = unlist(lapply(simulation, function(x) unlist(x$testing)))
  )
}

results <- bind_rows(
  combine_results(constant_orig_continuous, "constant", "orig", "continuous"),
  combine_results(constant_orig_ordinal, "constant", "orig", "ordinal"),
  combine_results(constant_rev_continuous, "constant", "rev", "continuous"),
  combine_results(constant_rev_ordinal, "constant", "rev", "ordinal"),
  combine_results(varied_orig_continuous, "varied", "orig", "continuous"),
  combine_results(varied_orig_ordinal, "varied", "orig", "ordinal"),
  combine_results(varied_rev_continuous, "varied", "rev", "continuous"),
  combine_results(varied_rev_ordinal, "varied", "rev", "ordinal")
)

# Visualize the results with ggplot
ggplot(results, aes(x = test_band, y = generated_velocity)) +
  geom_point(aes(color = training_type), alpha = 0.5, position = position_jitter(width = 10, height = 0)) +
  facet_grid(condition ~ feedback_type) +
  labs(
    title = "Simulation Results",
    x = "Test Band",
    y = "Generated Velocity",
    color = "Training Type"
  ) +
  theme_minimal()


  results_summary <- results %>%
  group_by(training_type, condition, feedback_type, test_band) %>%
  summarise(
    mean_velocity = mean(generated_velocity),
    sd_velocity = sd(generated_velocity),
    n = n()
  ) %>%
  arrange(training_type, condition, feedback_type, test_band)

# Display the summary table
results_summary

```

Based on the first 20 rows of the results_summary table, we can interpret the simulation results as follows:

For participants in the constant training type and orig condition with continuous feedback, their mean generated velocities are close to the center of each test band (e.g., 102 for the 100-300 band, 353 for the 350-550 band, etc.), with standard deviations around 97 to 103.
For participants in the constant training type and orig condition with ordinal feedback, their mean generated velocities are also close to the center of each test band, but the standard deviations are slightly larger (around 103 to 108), indicating more variability in their responses.
For participants in the constant training type and rev condition with continuous feedback, their mean generated velocities are again close to the center of each test band, and the standard deviations are similar to the orig condition with continuous feedback (around 98 to 101).
For participants in the constant training type and rev condition with ordinal feedback, their mean generated velocities are generally close to the center of each test band, but there is a bit more variability in the 100-300 band with a mean of 97.7. The standard deviations are around 102 to 107, which is slightly larger than the rev condition with continuous feedback.
These results suggest that participants in the constant training group generally perform well in estimating the center of the test bands, regardless of the condition and feedback type. The ordinal feedback appears to introduce more variability in the responses compared to continuous feedback. The rev condition doesn't seem to have a large impact on the performance of the constant training group.



To introduce variability in the noisiness of simulated agents and adjust the noisiness during the training stage between constant and varied training, we can modify the training stage function to include a noise_factor parameter. This parameter will control the amount of noise in the generated velocities for each participant.

```{r}
library(tidyverse)

# Training stage with noise factor
training_stage <- function(n_subjects, n_trials, constant_training, orig_condition, noise_factor) {
  subjects <- vector("list", n_subjects)
  
  for (i in 1:n_subjects) {
    if (orig_condition) {
      if (constant_training) {
        training_data <- rep(900, n_trials) + rnorm(n_trials, mean = 0, sd = 10 * noise_factor)
      } else {
        training_data <- c(rep(900, n_trials/3), rep(1100, n_trials/3), rep(1300, n_trials/3)) +
          rnorm(n_trials, mean = 0, sd = 10 * noise_factor)
      }
    } else {
      if (constant_training) {
        training_data <- rep(700, n_trials) + rnorm(n_trials, mean = 0, sd = 10 * noise_factor)
      } else {
        training_data <- c(rep(700, n_trials/3), rep(450, n_trials/3), rep(200, n_trials/3)) +
          rnorm(n_trials, mean = 0, sd = 10 * noise_factor)
      }
    }
    subjects[[i]] <- training_data
  }
  
  subjects
}

# Entire simulation with noise factor
simulate_experiment <- function(n_subjects, n_training_trials, n_testing_trials, constant_training, orig_condition, feedback_type, noise_factor) {
  training_subjects <- training_stage(n_subjects, n_training_trials, constant_training, orig_condition, noise_factor)
  testing_subjects <- testing_stage(training_subjects, n_testing_trials, feedback_type)
  
  testing_subjects
}

# Simulate populations for varied and constant learners with different noise factors
n_subjects <- 100
n_training_trials <- 90
n_testing_trials <- 15

constant_noise_factor <- 1
varied_noise_factor <- 2

constant_orig_continuous <- simulate_experiment(n_subjects, n_training_trials, n_testing_trials, constant_training = TRUE, orig_condition = TRUE, feedback_type = "continuous", noise_factor = constant_noise_factor)
constant_orig_ordinal <- simulate_experiment(n_subjects, n_training_trials, n_testing_trials, constant_training = TRUE, orig_condition = TRUE, feedback_type = "ordinal", noise_factor = constant_noise_factor)
constant_rev_continuous <- simulate_experiment(n_subjects, n_training_trials, n_testing_trials, constant_training = TRUE, orig_condition = FALSE, feedback_type = "continuous", noise_factor = constant_noise_factor)
constant_rev_ordinal <- simulate_experiment(n_subjects, n_training_trials, n_testing_trials, constant_training = TRUE, orig_condition = FALSE, feedback_type = "ordinal", noise_factor = constant_noise_factor)

varied_orig_continuous <- simulate_experiment(n_subjects, n_training_trials, n_testing_trials, constant_training = FALSE, orig_condition = TRUE, feedback_type = "continuous", noise_factor = varied_noise_factor)
varied_orig_ordinal <- simulate_experiment(n_subjects, n_training_trials, n_testing_trials, constant_training = FALSE, orig_condition = TRUE, feedback_type = "ordinal", noise_factor = varied_noise_factor)
varied_rev_continuous <- simulate_experiment(n_subjects, n_training_trials, n_testing_trials, constant_training = FALSE, orig_condition = FALSE, feedback_type = "continuous", noise_factor = varied_noise_factor



```



```{r}


library(shiny)
library(ggplot2)
library(tidyverse)

# Include the modified training_stage, testing_stage, and simulate_experiment functions here

# Shiny UI
ui <- fluidPage(
  titlePanel("Simulation Experiment"),

  sidebarLayout(
    sidebarPanel(
      sliderInput("n_subjects", "Number of Subjects", min = 10, max = 200, value = 100, step = 10),
      sliderInput("n_training_trials", "Number of Training Trials", min = 10, max = 200, value = 90, step = 10),
      sliderInput("n_testing_trials", "Number of Testing Trials", min = 5, max = 50, value = 15, step = 5),
      sliderInput("constant_noise_factor", "Constant Noise Factor", min = 1, max = 5, value = 1, step = 0.5),
      sliderInput("varied_noise_factor", "Varied Noise Factor", min = 1, max = 5, value = 2, step = 0.5)
    ),

    mainPanel(
      tableOutput("results_table"),
      plotOutput("learning_curve_plot"),
      plotOutput("generalization_plot")
    )
  )
)

# Shiny Server
server <- function(input, output) {
  results <- reactive({
    constant_orig_continuous <- simulate_experiment(input$n_subjects, input$n_training_trials, input$n_testing_trials, constant_training = TRUE, orig_condition = TRUE, feedback_type = "continuous", noise_factor = input$constant_noise_factor)
    
    # Process the results for the table
    results_table <- data.frame(
      Group = c("Constant Orig Continuous"),
      Mean_Training_Performance = mean(unlist(lapply(constant_orig_continuous, function(x) x$training))),
      Mean_Testing_Performance = mean(unlist(lapply(constant_orig_continuous, function(x) unlist(x$testing))))
    )
    
    results_table
  })

  output$results_table <- renderTable({
    req(results())
  })

  output$learning_curve_plot <- renderPlot({
    req(results())
    
    # Process the results for the learning curve plot
    learning_curve_data <- tibble(
      Subject = rep(1:length(results()), each = input$n_training_trials),
      Group = rep("Constant Orig Continuous", input$n_subjects * input$n_training_trials),
      Trial = rep(1:input$n_training_trials, input$n_subjects),
      Velocity = unlist(lapply(results(), function(x) x$training))
    )

    ggplot(learning_curve_data, aes(x = Trial, y = Velocity, color = Group)) +
      geom_line(alpha = 0.5) +
      labs(title = "Learning Curve", x = "Trial", y = "Velocity") +
      theme_minimal()
  })

  output$generalization_plot <- renderPlot({
    req(results())
    
    # Process the results for the generalization plot
    generalization_data <- tibble(
      Subject = rep(1:length(results()), each = input$n_testing_trials * length(c(seq(100, 1400, by = 250)))),
      Group = rep("Constant Orig Continuous", input$n_subjects * input$n_testing_trials * length(c(seq(100, 1400, by = 250)))),
      Band = rep(c(seq(100, 


```










```{r}
install.packages("keras")
install.packages("tensorflow")

library(keras)
library(tensorflow)

# Create a simple neural network model
create_nn_model <- function(n_hidden_nodes) {
  model <- keras_model_sequential() %>%
    layer_dense(units = n_hidden_nodes, activation = "relu", input_shape = c(1)) %>%
    layer_dense(units = 1)
  
  model %>%
    compile(
      loss = "mean_squared_error",
      optimizer = optimizer_adam(lr = 0.01),
      metrics = c("mean_absolute_error")
    )
  
  model
}

# Modify the training stage to use the neural network model
training_stage_nn <- function(n_subjects, n_trials, constant_training, orig_condition, n_hidden_nodes) {
  subjects <- vector("list", n_subjects)
  
  for (i in 1:n_subjects) {
    if (orig_condition) {
      if (constant_training) {
        training_data <- rep(900, n_trials)
      } else {
        training_data <- c(rep(900, n_trials/3), rep(1100, n_trials/3), rep(1300, n_trials/3))
      }
    } else {
      if (constant_training) {
        training_data <- rep(700, n_trials)
      } else {
        training_data <- c(rep(700, n_trials/3), rep(450, n_trials/3), rep(200, n_trials/3))
      }
    }
    
    # Train the neural network model
    model <- create_nn_model(n_hidden_nodes)
    history <- model %>%
      fit(as.matrix(training_data[-length(training_data)]), as.matrix(training_data[-1]), epochs = 100, verbose = 0)
    
    subjects[[i]] <- list(model = model, history = history)
  }
  
  subjects
}

# Modify the testing stage to use the neural network model
testing_stage_nn <- function(training_subjects, n_trials) {
  subjects <- vector("list", length(training_subjects))
  test_bands <- c(seq(100, 1400, by = 250))
  
  for (i in 1:length(training_subjects)) {
    model <- training_subjects[[i]]$model
    test_data <- vector("list", length(test_bands))
    
    for (j in 1:length(test_bands)) {
      test_input <- rep(test_bands[j], n_trials)
      test_output <- predict(model, as.matrix(test_input))
      test_data[[j]] <- test_output
    }
    
    subjects[[i]] <- list(training = training_subjects[[i]]$history, testing = test_data)
  }
  
  subjects
}

# Modify the entire simulation to use the neural network model
simulate_experiment_nn <- function(n_subjects, n_training_trials, n_testing_trials, constant_training, orig_condition, n_hidden_nodes) {
  training_subjects <- training_stage_nn(n_subjects, n_training_trials, constant_training, orig_condition, n_hidden_nodes)
  testing_subjects <- testing_stage_nn(training_subjects, n_testing_trials)
  
  testing_subjects
}

# Simulate populations for varied and constant learners with different numbers of hidden



```























### Turbo
To implement the computational cognitive Bayesian model in R, we first need to specify the model structure and the parameters. Here is an outline of the model:

Model structure:
The learner's belief about the correct velocity to produce is modeled by a Gaussian distribution.
The mean of the distribution represents the learner's best estimate of the correct velocity.
The standard deviation of the distribution represents the learner's uncertainty about the correct velocity.
Parameters:
The mean of the prior distribution represents the learner's initial belief about the correct velocity before any training.
The standard deviation of the prior distribution represents the learner's uncertainty about the correct velocity before any training.
The learning rate determines how much the learner updates their belief about the correct velocity after each trial.
The noise parameter represents the amount of noise in the learner's motor system.
With this outline in mind, here is the R code to simulate the model for an artificial population of varied and constant learners:

```{r}

library(tidyverse)

# Set up the model parameters
prior_mean <- 1000   # Mean of the prior distribution
prior_sd <- 200      # Standard deviation of the prior distribution
learning_rate <- 0.1 # Learning rate
noise <- 50          # Noise parameter

# Define a function to simulate the model for one participant
simulate_participant <- function(training_data) {
  # Initialize the posterior distribution with the prior distribution
  posterior_mean <- prior_mean
  posterior_sd <- prior_sd
  
  # Simulate the participant's performance during training
  for (i in 1:nrow(training_data)) {
    # Compute the likelihood of the observed data given the current belief
    likelihood <- dnorm(training_data$vx[i], mean = posterior_mean, sd = sqrt(posterior_sd^2 + noise^2))
    
    # Update the posterior distribution using Bayes' rule
    posterior_mean <- (1 - learning_rate) * posterior_mean + learning_rate * training_data$vx[i]
    posterior_sd <- sqrt((1 - learning_rate) * posterior_sd^2 + learning_rate * (training_data$vx[i] - posterior_mean)^2)
  }
  
  # Return the final posterior mean and condition label for the participant
  return(list(mean = posterior_mean, condit = unique(training_data$condit)))
}

# Set up the training data
n_trials <- 90
n_participants <- 196
training_data <- data.frame(id = rep(1:n_participants, each = n_trials),
                             condit = rep(c("Constant", "Varied"), each = n_trials / 2),
                             vb = rep(c("800-1000", "1000-1200", "1200-1400"), times = n_trials / 3))
training_data$vx <- ifelse(training_data$condit == "Constant", 1000, as.numeric(gsub("-.*", "", training_data$vb)))

# Simulate the model for each participant
participants <- training_data %>%
  group_by(id) %>%
  nest() %>%
  mutate(participant = map(data, simulate_participant)) %>%
  unnest(participant) %>%
  left_join(training_data[, c("id", "condit")], by = "id")

# Compute the learning curves
learning_curves <- participants %>%
  group_by(condit) %>%
  summarize(mean_posterior_mean = mean(mean))

# Plot the learning curves
ggplot(training_data, aes(x = condit, y = vx)) +
  ylim(0, 1400) +
  geom_point(aes(y = mean_posterior_mean), data = participants, color = "red") +
  geom_point(color = "blue") +
  facet_wrap(~ vb) +
  scale_x_discrete(labels = c("Constant" = "Constant Velocity", "Varied" = "Varied Velocity")) +
  labs(x = "Velocity Condition", y = "Velocity") +
  ggtitle("Learning Curves") +
  theme(plot.title = element_text(hjust = 0.5))

# Set up the testing data
n_trials <- 15
testing_data <- data.frame(id = rep(1:n_participants, each = 2 * n_trials),
                            condit = rep(c("Constant", "Varied"), each = n_trials),
                            vb = rep(c("800-1000", "1000-1200", "1200-1400"), each = n_trials))
testing_data$vx <- ifelse(testing_data$condit == "Constant", 1000, as.numeric(gsub("-.*", "", testing_data$vb))))

# Simulate the model for each participant during testing
testing_participants <- testing_data %>% 
  group_by(id) %>% 
  nest() %


```


```{r}
library(data.table)

prior_mean <- 1000 # Mean of the prior distribution
prior_sd <- 200 # Standard deviation of the prior distribution
learning_rate <- 0.1 # Learning rate
noise <- 50 # Noise parameter

#Define a function to simulate the model for one participant
simulate_participant <- function(training_data) {
  
  # Initialize the posterior distribution with the prior distribution
  posterior_mean <- prior_mean
  posterior_sd <- prior_sd
  
  # Create an empty vector to store the posterior means for each trial
  posterior_means <- numeric(nrow(training_data))
  
  # Define a function to simulate the model for one participant
  for (i in 1:nrow(training_data)) {
    # Compute the likelihood of the observed data given the current belief
    likelihood <- dnorm(training_data$vx[i], mean = posterior_mean, sd = sqrt(posterior_sd^2 + noise^2))
    
    # Update the posterior distribution using Bayes' rule
    posterior_mean <- (1 - learning_rate) * posterior_mean + learning_rate * training_data$vx[i]
    posterior_sd <- sqrt((1 - learning_rate) * posterior_sd^2 + learning_rate * (training_data$vx[i] - posterior_mean)^2)
    
    # Store the posterior mean for this trial
    posterior_means[i] <- posterior_mean
  }
  
  # Return the final posterior mean and condition label for the participant
  list(mean = posterior_mean, condit = unique(training_data$condit))
}
# Set up the training data
n_trials <- 90
n_participants <- 196
training_data <- data.table(id = rep(1:n_participants, each = n_trials),
                             condit = rep(c("Constant", "Varied"), each = n_trials / 2),
                             vb = rep(c("800-1000", "1000-1200", "1200-1400"), times = n_trials / 3))
training_data[, vx := ifelse(condit == "Constant", 1000, as.numeric(gsub("-.*", "", vb)))]

# Simulate the model for each participant
participants <- training_data[, .(simulate_participant(.SD)), by = id]$V1

head(participants)

# Extract the posterior means and condition labels for each participant
posterior_means <- sapply(participants, function(x) x[[1]])
condits <- sapply(participants, function(x) x[[2]])

# Compute the learning curves
learning_curves <- tapply(posterior_means, condits, mean)

# Plot the learning curves
ggplot(training_data, aes(x = vb, y = vx, color = condit)) +
  geom_point() +
  stat_summary(aes(y = posterior_means), fun = mean, geom = "line", size = 1.2, linetype = "dashed") +
  ylim(0, 1400) +
  scale_color_manual(values = c("red", "blue")) +
  labs(title = "Learning Curves", x = "VB", y = "VX", color = "Condition") +
  theme_minimal()



```