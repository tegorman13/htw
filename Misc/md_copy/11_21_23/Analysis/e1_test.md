

#source(here::here("Functions", "packages.R"))
library(pacman)
pacman::p_load(tidyverse,tidybayes,brms, lme4, bayesplot,bayestestR,parameters,marginaleffects,
                emmeans, equatiomatic, here, pacman,  broom,
               broom.mixed,lme4,emmeans,here,knitr,kableExtra,gt,
                wesanderson,glue, ggdist,ggforce,patchwork,gghalves,ggh4x,
                pander,
                install = TRUE,
                update = FALSE
               )

test <- readRDS(here("data/e1_08-21-23.rds")) |> filter(expMode2 == "Test") 
options(brms.backend="cmdstanr",mc.cores=4)
e1Sbjs <- test |> group_by(id,condit) |> summarise(n=n())
testAvg <- test %>% group_by(id, condit, vb, bandInt,bandType,tOrder) %>%
  summarise(nHits=sum(dist==0),vx=mean(vx),dist=mean(dist),sdist=mean(sdist),n=n(),Percent_Hit=nHits/n)


result <- test_summary_table(test, "dist", mfun = list(mean = mean, median = median, sd = sd))
result$constant
result$varied





library(mvtnorm)
n <- 100 # number of trials 
k <- 6 # number of groups 
m <- 5 # number of lines 
mu <- seq(0.5, 0.9, length.out = m) # mean reward for each line 
sigma <- 0.1 # standard deviation of reward for each line 
theta <- 0.1 # decay rate of reward for each line 
labels <- c("human", "GP-UCB", "lesion", "lesion", "random") # labels for each line

#Generate data
data <- list() 

for (i in 1:k) {

##Generate trial numbers for each group
x <- seq(from = i * 5, to = i * 5 + 4, length.out = n)

#Generate rewards for each line
y <- matrix(NA, nrow = n, ncol = m) 
for (j in 1:m) { 
  # Generate Gaussian noise with correlation structure 
  noise <- rmvnorm(n, mean = rep(0, n), sigma = sigma^2 * exp(-abs(outer(x, x, "-")) * theta))
  # Generate reward values with decay function
y[, j] <- mu[j] * exp(-theta * x) + noise[, j]

}
data[[i]] <- data.frame(x = rep(x, m), y = c(y), id = rep(labels, each = n)) }

ggplot(data = do.call(rbind, data), aes(x = x, y = y, color = id)) + 
geom_line() + 
facet_wrap(~ x %/% n + 1) + 
labs(x = "Trial", y = "Normalized mean reward ± s.e.m.", title = "Synthetic data") + 
theme_bw() + 
theme(legend.position = "top")




library(latex2exp)


eqns <- c( 
  "a(x) = \\sum_i w_i x_i", 
  "o(x) = \\sum_i w_i x_i", 
  "p(x) = \\frac{\\exp(o(x))}{\\sum_j \\exp(o(j))}", 
  "m(x) = \\sum_j p(j) r(j)", 
  "f(x) = r(x) - m(x)", 
  "\\Delta w_i = \\alpha f(x) x_i", 
  "e(x) = \\sum_j p(j) a(j)", 
  "i(x) = a(x) - e(x)" 
)

labels <- c( "Input activation", "Output activation", "Output probability", "Mean output", "Feedback activation", "Weight adjustment", "Extrapolation", "Interpolation" )

#Create matrix with labels and equations
mat <- matrix(c(labels, eqns), ncol = 2, byrow = TRUE)

#Convert matrix to data frame
df <- as.data.frame(mat)

#Set column names
colnames(df) <- c(" ", "EXAM")

#Convert equations to expressions
df$EXAM <- TeX(df$EXAM)

print(df, row.names = FALSE, right = FALSE) 



table1 <- tibble(
  Exemplar = c("T1", "T2", "T3", "M1", "M2", "M3", "M4", "N1", "N2", "N3"),
  Symptoms = c("1", "1", "1", "0", "0", "0", "1", "1", "1", "0"),
  Observed = c(".88", ".89", ".83", ".79", ".78", ".82", ".76", ".88", ".82", ".86"),
  ALCOVE = c(".76", ".78", ".82", ".76", ".75", ".76", ".76", ".82", ".86", ".82")
)



# | Exemplar | Symptoms | Observed | ALCOVE | Config cue |
# |---|---|---|---|---|
# | T1 | 1 | .88 | .82 | .76 |
# | T2 | 0 | .89 | .78 | .76 |
# | T3 | 1 | .73 | .72 | .75 |
# | M1 | 0 | .79 | .82 | .76 |
# | M2 | 1 | .83 | .82 | .76 |
# | M3 | 1 | .77 | .82 | .76 |
# | M4 | 0 | .82 | .82 | .76 |
# | N1 | 1 | .88 | .82 | .76 |
# | N2 | 0 | .88 | .82 | .76 |
# | N3 | 1 | .88 | .82 | .76 |



library(tidyverse)
library(ggplot2)

# Define the number of epochs
epochs <- 100

# Define the learning curve patterns for each category type
curve_I <- 1 - exp(-0.05 * seq(1, epochs))
curve_II <- 1 - exp(-0.1 * seq(1, epochs))
curve_III <- 1 - exp(-0.15 * seq(1, epochs))
curve_IV <- 1 - exp(-0.2 * seq(1, epochs))
curve_V <- 1 - exp(-0.25 * seq(1, epochs))
curve_VI <- 1 - exp(-0.3 * seq(1, epochs))

# Combine the learning curve patterns into a data frame
learning_curves <- data.frame(
  epoch = rep(seq(1, epochs), times = 6),
  category_type = rep(c("I", "II", "III", "IV", "V", "VI"), each = epochs),
  probability_correct = c(curve_I, curve_II, curve_III, curve_IV, curve_V, curve_VI)
)

# Plot the learning curves
ggplot(data = learning_curves, aes(x = epoch, y = probability_correct)) +
  geom_line(aes(color = category_type)) +
  scale_color_manual(values = c("red", "blue", "green", "purple", "orange", "black")) +
  labs(x = "Epoch", y = "Pr(correct)", title = "Learning Curves for Different Category Types")
