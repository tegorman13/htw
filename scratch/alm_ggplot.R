

library(tidyverse)
library(patchwork)
library(ggforce)

theme_set(theme_grey() +
            theme_void() +
            theme(plot.margin = margin(0, 5.5, 0, 5.5)))
p1 <-
  tibble(x = seq(from = .01, to = .99, by = .01),
         d = (dbeta(x, 2, 2)) / max(dbeta(x, 2, 2))) %>% 
  ggplot(aes(x = x, y = d)) +
  geom_area(fill = "skyblue", size = 0) +
  annotate(geom = "text",
           x = .5, y = .2,
           label = "beta",
           size = 7) +
  annotate(geom = "text",
           x = .5, y = .6,
           label = "italic(A[omega])*', '*italic(B[omega])", 
           size = 7, family = "Times", parse = TRUE) +
  scale_x_continuous(expand = c(0, 0)) +
  theme(axis.line.x = element_line(size = 0.5))

# a gamma density
p2 <-
  tibble(x = seq(from = 0, to = 5, by = .01),
         d = (dgamma(x, 1.75, .85) / max(dgamma(x, 1.75, .85)))) %>% 
  ggplot(aes(x = x, y = d)) +
  geom_area(fill = "skyblue", size = 0) +
  annotate(geom = "text",
           x = 2.5, y = .2,
           label = "gamma",
           size = 7) +
  annotate(geom = "text",
           x = 2.5, y = .6,
           label = "list(italic(S)[kappa], italic(R)[kappa])",
           size = 7, family = "Times", parse = TRUE) +
  scale_x_continuous(expand = c(0, 0)) +
  theme(axis.line.x = element_line(size = 0.5))

p3 <-
  tibble(x = c(.5, .475, .26, .08, .06,
               .5, .55, .85, 1.15, 1.175,
               1.5, 1.4, 1, .25, .2,
               1.5, 1.49, 1.445, 1.4, 1.39),
         y = c(1, .7, .6, .5, .2,
               1, .7, .6, .5, .2,
               1, .7, .6, .5, .2,
               1, .75, .6, .45, .2),
         line = rep(letters[2:1], each = 5) %>% rep(., times = 2),
         plot = rep(1:2, each = 10)) %>% 
  
  ggplot(aes(x = x, y = y, group = interaction(plot, line))) +
  geom_bspline(aes(color = line),
               size = 2/3, show.legend = F) + 
  annotate(geom = "text",
           x = 0, y = .1,
           label = "omega(kappa-2)+1*', '*(1-omega)(kappa-2)+1",
           size = 7, parse = T, family = "Times", hjust = 0) +
  annotate(geom = "text",
           x = c(1/3, 1.15), y = .7,
           label = "'~'",
           size = 10, parse = T, family = "Times") +
  scale_color_manual(values = c("grey75", "black")) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 2)) +
  ylim(0, 1)

p4 <-
  tibble(x = seq(from = .01, to = .99, by = .01),
         d = (dbeta(x, 2, 2)) / max(dbeta(x, 2, 2))) %>% 
  ggplot(aes(x = x, y = d)) +
  geom_area(fill = "skyblue", size = 0) +
  annotate(geom = "text",
           x = .5, y = .2,
           label = "beta",
           size = 7) +
  scale_x_continuous(expand = c(0, 0)) +
  theme(axis.line.x = element_line(size = 0.5))

# an annotated arrow
p5 <-
  tibble(x     = c(.375, .625),
         y     = c(1/3, 1/3),
         label = c("'~'", "italic(s)")) %>% 
  
  ggplot(aes(x = x, y = y, label = label)) +
  geom_text(size = c(10, 7), parse = T, family = "Times") +
  geom_segment(x = 0.5, xend = 0.5,
               y = 1, yend = 0,
               arrow = my_arrow) +
  xlim(0, 1)

p6 <-
  tibble(x = 0:1,
         d = (dbinom(x, size = 1, prob = .6)) / max(dbinom(x, size = 1, prob = .6))) %>% 
  
  ggplot(aes(x = x, y = d)) +
  geom_col(fill = "skyblue", width = .4) +
  annotate(geom = "text",
           x = .5, y = .2,
           label = "Bernoulli",
           size = 7) +
  annotate(geom = "text",
           x = .5, y = .94,
           label = "theta", 
           size = 7, family = "Times", parse = T) +
  xlim(-.75, 1.75) +
  theme(axis.line.x = element_line(size = 0.5))

p7 <-
  tibble(x     = c(.35, .65),
         y     = c(1/3, 1/3),
         label = c("'~'", "italic(i)*'|'*italic(s)")) %>% 
  
  ggplot(aes(x = x, y = y, label = label)) +
  geom_text(size = c(10, 7), parse = T, family = "Times") +
  geom_segment(x = .5, xend = .5,
               y = 1, yend = 0,
               arrow = my_arrow) +
  xlim(0, 1)

p8 <-
  tibble(x     = .5,
         y     = .5,
         label = "italic(y[i])['|'][italic(s)]") %>% 
  
  ggplot(aes(x = x, y = y, label = label)) +
  geom_text(size = 7, parse = T, family = "Times") +
  xlim(0, 1)

layout <- c(
  area(t = 1, b = 2, l = 1, r = 1),
  area(t = 1, b = 2, l = 2, r = 2),
  area(t = 4, b = 5, l = 1, r = 1),
  area(t = 3, b = 4, l = 1, r = 2),
  area(t = 6, b = 6, l = 1, r = 1),
  area(t = 7, b = 8, l = 1, r = 1),
  area(t = 9, b = 9, l = 1, r = 1),
  area(t = 10, b = 10, l = 1, r = 1)
)

(p1 + p2 + p4 + p3 + p5 + p6 + p7 + p8) + 
  plot_layout(design = layout) &
  ylim(0, 1)










library(tidyverse)
library(ggforce)
library(latex2exp)

theme_set(theme_grey() +
            theme_void() +
            theme(plot.margin = margin(0, 5.5, 0, 5.5)))

# Parameters
input_x <- c(1, 2)
output_x <- c(1.5, 2.5, 3.5)
input_y <- 3
output_y <- 1
activation_levels <- c(0.2, 0.7, 0.5)

# Input Layer
input_layer <- tibble(x = input_x, y = rep(input_y, length(input_x)))

# Output Layer
output_layer <- tibble(x = output_x, y = rep(output_y, length(output_x)), activation = activation_levels)

# Gaussian Activation
gaussian_activation <- tibble(
  x = seq(from = 0, to = 3, by = 0.01),
  y = exp(-2 * (x - 1)^2) + 2.5,
  x2 = x,
  y2 = exp(-2 * (x - 2)^2) + 2.5
)

# Plot
p <- ggplot() +
  # Input nodes
  geom_point(data = input_layer, aes(x = x, y = y), size = 4, color = 'green') +
  annotate("text", x = input_x, y = rep(input_y, length(input_x)) + 0.3, label = c("Input #1", "Input #2")) +
  
  # Gaussian Activations
  geom_line(data = gaussian_activation, aes(x = x, y = y), color = 'blue') +
  geom_line(data = gaussian_activation, aes(x = x2, y = y2), color = 'blue') +
  
  # Connections
  geom_segment(data = expand.grid(input_x, output_x), aes(x = Var1, xend = Var2, y = input_y, yend = output_y), arrow = arrow(type = 'closed', length = unit(0.2, 'inches'))) +
  
  # Output nodes
  geom_point(data = output_layer, aes(x = x, y = y), size = 4, color = 'red') +
  geom_bar(data = output_layer, aes(x = x, y = activation), stat = 'identity', position = 'dodge', fill = 'red', alpha = 0.3, width = 0.3) +
  annotate("text", x = output_x, y = rep(output_y, length(output_x)) - 0.3, label = c("Output #1", "Output #2", "Output #3")) +
  
  # Equation annotations
  annotate("text", x = 0.5, y = input_y + 1, label = TeX(r"( $\gamma^2 = \alpha^2 + \beta^2$ )")) +
  annotate("text", x = 1.5, y = output_y - 1, 
           label = TeX(r"($S = \{z \in \bf{C}\, |\, |z|<1 \} \, \textrm{and} \, S_2=\partial{S}$)")) +
  
  # Input stimulus and output response
  annotate("text", x = mean(input_x), y = input_y + 1.3, label = "Input Stimulus") +
  annotate("text", x = mean(output_x), y = output_y - 1.3, label = "Output Response") +
  
  # Coordinate limits and axis labels
  coord_cartesian(xlim = c(0, 4), ylim = c(-1, 5)) +
  labs(x = "", y = "") +
  theme_void()

# Show the plot
print(p)








theme_set(theme_grey() +
            theme_void() +
            theme(plot.margin = margin(0, 5.5, 0, 5.5)))

# Parameters
input_x <- c(1, 2)
output_x <- c(1.5, 2.5, 3.5)
input_y <- 3
output_y <- 1
activation_levels <- c(0.2, 0.7, 0.5)

# Input Layer
input_layer <- tibble(x = input_x, y = rep(input_y, length(input_x)))

# Output Layer
output_layer <- tibble(x = output_x, y = rep(output_y, length(output_x)), activation = activation_levels)

# Gaussian Activation
gaussian_activation <- tibble(
  x = seq(from = 0, to = 3, by = 0.01),
  y = exp(-8 * (x - 1)^2)+ 2.5,
  x2 = x,
  y2 = exp(-8 * (x - 2)^2)+2.5
)

# Random connection weights
connection_weights <- runif(length(input_x) * length(output_x), 0.1, 1)

# Plot
p <- ggplot() +
  # Input nodes
  geom_point(data = input_layer, aes(x = x, y = y), size = 4, color = 'green') +
  annotate("text", x = input_x, y = rep(input_y, length(input_x)) + 0.3, label = c("Input #1", "Input #2")) +
  # Gaussian Activations
  geom_line(data = gaussian_activation, aes(x = x, y = y), color = 'blue') +
  geom_line(data = gaussian_activation, aes(x = x2, y = y2), color = 'blue') +
  # Connections
  geom_segment(data = expand.grid(input_x, output_x), aes(x = Var1, xend = Var2, y = input_y, yend = output_y, size = connection_weights), arrow = arrow(type = 'closed', length = unit(0.2, 'inches'))) +
  # Output nodes
  geom_point(data = output_layer, aes(x = x, y = y), size = 4, color = 'red') +
  geom_bar(data = output_layer, aes(x = x, y = activation), stat = 'identity', position = 'dodge', fill = 'red', alpha = 0.3, width = 0.3) +
  annotate("text", x = output_x, y = rep(output_y, length(output_x)) - 0.3, label = c("Output #1", "Output #2", "Output #3")) +
  # Equation annotations
  annotate("text", x = 0.2, y = input_y + 1, label = TeX(r"( $a_i(X)=e^{-\gamma \cdot\left(X-X_i\right)^2$})"), parse = TRUE) +
  annotate("text", x = 0.5, y = output_y - 1, label = TeX(r"( O_j(X)=\sum_{i=1}^M w_{j i} \cdot a_i(X))"), parse = TRUE) +
  annotate("text", x = 0.5, y = output_y - 1.4, 
           label = TeX("$\\m(X)=\\sum_{j=1}^L Y_j \\cdot P\\left[Y_j \\mid X\\right]$",output='character'), parse = TRUE) + 
  
  # Input stimulus and output response
  annotate("text", x = mean(input_x), y = input_y + 1.3, label = "Input Stimulus") +
  annotate("text", x = mean(output_x), y = output_y - 1.3, label = "Output Response") +
  # Coordinate limits and axis labels
  coord_cartesian(xlim = c(0, 4), ylim = c(-1, 5)) +
  labs(x = "", y = "") + theme_void() +
  scale_size_continuous(range = c(1, 3), guide = FALSE)  # Line sizes based on weights

# Show the plot
print(p)







#Random connection weights
connection_weights_matrix <- matrix(runif(length(input_x) * length(output_x), 0.1, 1), nrow=length(output_x))
# Calculating the activations
activation_levels <- rowSums(connection_weights_matrix)
# Average activation
average_activation <- mean(activation_levels)

# Input Layer
input_layer <- tibble(x = input_x, y = rep(input_y, length(input_x)))

# Output Layer
output_layer <- tibble(x = output_x, y = rep(output_y, length(output_x)), activation = activation_levels)

# Plot
p <- ggplot() +
  # Input nodes
  geom_point(data = input_layer, aes(x = x, y = y), size = 4, color = 'green') +
  annotate("text", x = input_x, y = rep(input_y, length(input_x)) - 0.2, label = c("Input #1", "Input #2")) +
  
  # Connections
  geom_segment(data = expand.grid(input_x, output_x), aes(x = Var1, xend = Var2, y = input_y, yend = output_y, size = as.vector(connection_weights_matrix)), arrow = arrow(type = 'closed', length = unit(0.2, 'inches'))) +
  
  # Output nodes
  geom_point(data = output_layer, aes(x = x, y = y), size = 4, color = 'red') +
  geom_bar(data = output_layer, aes(x = x, y = activation), stat = 'identity', position = 'dodge', fill = 'red', alpha = 0.3, width = 0.3) +
  annotate("text", x = output_x, y = rep(output_y, length(output_x)) + 0.2, label = c("Output #1", "Output #2", "Output #3")) +
  
  # Response
  geom_point(aes(x = mean(output_x), y = output_y - 0.7), size = average_activation * 5, color = 'purple') +
  annotate("text", x = mean(output_x), y = output_y - 1, label = "Response (average of activations)", size = 3) +
  
  # Coordinate limits and axis labels
  coord_cartesian(xlim = c(0, 4), ylim = c(-1, 3)) +
  labs(x = "", y = "") +
  theme_void() +
  scale_size_continuous(range = c(1, 3), guide = FALSE)  # Line sizes based on weights

# Show the plot
print(p)