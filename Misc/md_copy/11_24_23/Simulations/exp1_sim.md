


Design The experiment employed a 2 (Training Condition: varied vs. constant).

Procedure Upon arrival at the laboratory, participants were provided with a description of the experiment and signed informed consent forms. They were then seated in front of a computer equipped with a mouse and were given instructions on how to perform the “Hit The Wall” (HTW) visuomotor extrapolation task.

The HTW task involved launching projectiles to hit a target displayed on the computer screen. Participants completed a total of 90 trials during the training stage. In the varied training condition, participants encountered three velocity bands (800-1000, 1000-1200, and 1200-1400). In contrast, participants in the constant training condition encountered only one velocity band (800-1000).

During the training stage, participants in both conditions also completed “no feedback” trials, where they received no information about their performance. These trials were randomly interleaved with the regular training trials.

Following the training stage, participants proceeded to the testing stage, which consisted of three phases. In the first phase, participants completed “no-feedback” testing from three novel extrapolation bands (100-300, 350-550, and 600-800), with each band consisting of 15 trials.

In the second phase of testing, participants completed “no-feedback” testing from the three velocity bands used during the training stage (800-1000, 1000-1200, and 1200-1400). In the constant training condition, two of these bands were novel, while in the varied training condition, all three bands were encountered during training.

The third and final phase of testing involved “feedback” testing for each of the three extrapolation bands (100-300, 350-550, and 600-800), with each band consisting of 10 trials. Participants received feedback on their performance during this phase.

Throughout the experiment, participants’ performance was measured by calculating the distance between the produced x-velocity of the projectiles and the closest edge of the current velocity band. Lower distances indicated better performance.

```{r}
library(tidyverse)

# simulate design of experiment
set.seed(123)
n_participants <- 20
num_constant_subjects <- 10
num_varied_subjects <- 10
n_trials <- 90
n_bands_varied <- 3
n_bands_constant <- 1
# Set the training positions
training_pos_constant <- 800
training_pos_varied <- c(800, 1000, 1200)

test_positions <- c(200, 400, 600, 800, 1000, 1200)
# Set the noise level (standard deviation of the noise)
noise_level <- 0.1
# Set the initial performance for each subject
initial_performance <- 0.5
# Set the learning rate (how much performance improves with each trial)
learning_rate <- 0.01

# Set the number of trials for each phase
n_training_trials <- 90
n_testing_trials <- 15
n_feedback_trials <- 10

# Set the number of phases
n_phases <- 3


# Initialize vectors to store the performance of each subject
performance_constant <- rep(initial_performance, num_constant_subjects)
performance_varied <- rep(initial_performance, num_varied_subjects)





# Simulate the training phase for the constant condition
for (i in 1:num_constant_subjects) {
  for (j in 1:n_training_trials) {
    # Simulate the participant's response
    response <- rnorm(1, training_pos_constant, noise_level)
    # Calculate the distance between the response and the target
    distance <- abs(response - training_pos_constant)
    # Update the performance
    performance_constant[i] <- performance_constant[i] + learning_rate * (1 - distance)
  }
}




```

