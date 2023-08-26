
source(here::here("Functions", "packages.R"))
test <- readRDS(here("data/e1_08-21-23.rds")) |> 
  filter(expMode2 == "Test") 



biModel <- glmer(result == "Hit" ~ condit*bandInt + (1 | id), 
                 data = test, 
                 family = binomial)
summary(biModel)
library(multcomp)
summary(glht(biModel))

biModel2 <- glmer(result == "Hit" ~ condit*bandInt + (bandInt | id), data = test, family = binomial)
summary(biModel2)

biModel3 <- glmer(result == "Hit" ~ condit*vb + (vb | id), data = test, family = binomial)
summary(biModel3)
plot(ranef(biModel3))
# Refitting the model using different optimizer
biModel3_refit <- glmer(result == "Hit" ~ condit*vb + (vb | id), 
                        data = test, family = binomial,
                        control=glmerControl(optimizer="bobyqa"))
summary(biModel3_refit)

agg_data <- test %>%
  group_by(vb, condit) %>%
  summarise(hits = sum(result == "Hit"), trials = n()) %>%
  mutate(hit_rate = hits / trials)
print(agg_data)




ggplot(agg_data, aes(x = vb, y = hit_rate, color = condit)) +
  geom_point() + geom_smooth(method = "lm", se = FALSE)

residuals <- residuals(biModel3, type = "pearson")
hist(residuals, breaks = 20)


# Extract fixed effect estimates
fix_eff <- summary(biModel3)$coefficients

# Create a data frame for plotting
plot_data <- data.frame(
  vb = rep(levels(test$vb)[2:6], each = 2),
  condit = rep(c("Constant", "Varied"), times = length(levels(test$vb))-1),
  estimate = c(fix_eff[3:7, 1], fix_eff[3:7, 1] + fix_eff[8:12, 1])
)
ggplot(plot_data, aes(x = vb, y = estimate, color = condit)) +
  geom_point() +
  geom_line(aes(group = condit)) +
  labs(title = "Fixed Effect Estimates Across `vb` Categories",
       x = "vb Category",
       y = "Estimate")





resid_vals <- resid(biModel3)
fitted_vals <- fitted(biModel3)
ggplot() + geom_point(aes(x=fitted_vals, y=resid_vals)) +
  ggtitle('Residuals vs Fitted Values') +
  xlab('Fitted Values') + ylab('Residuals')


ranef_plot <- dotplot(ranef(biModel3, condVar=TRUE))
print(ranef_plot)



biModel_simple <- glmer(result == "Hit" ~ condit + vb + (1 | id), data = test, family = binomial)
anova(biModel3, biModel_simple)



conf_matrix_breakdown <- test %>%
  group_by(condit, vb, result) %>%
  summarise(count = n())




ggplot(conf_matrix_breakdown, aes(x = vb, y = result)) +
  geom_tile(aes(fill = count), color = "white") +
  facet_grid(~condit) +
  scale_fill_gradient(low = "white", high = "red") +
  labs(title = "Confusion Matrix Breakdown by 'condit' and 'vb'",
       x = "vb",
       y = "Result",
       fill = "Count")




table(test$result,test$vb)
table(test$vxCat,test$vb)


conf_matrix_constant <- test %>% 
  filter(condit == "Constant") %>% 
  select(vxCat, vb) %>% 
  count(vxCat, vb) %>%
  spread(vb, n, fill = 0)

# Creating a confusion matrix for the 'Varied' condition
conf_matrix_varied <- test %>% 
  filter(condit == "Varied") %>% 
  select(vxCat, vb) %>% 
  count(vxCat, vb) %>% 
  spread(vb, n, fill = 0)


common_cols <- intersect(names(conf_matrix_constant), names(conf_matrix_varied))
conf_matrix_constant <- conf_matrix_constant %>% select(all_of(common_cols))
conf_matrix_varied <- conf_matrix_varied %>% select(all_of(common_cols))
# Compute the difference matrix
conf_matrix_difference <- data.table::setDT(conf_matrix_varied) - data.table::setDT(conf_matrix_constant)
# Display the difference matrix
print(conf_matrix_difference)


create_normalized_conf_matrix <- function(data, condition) {
  conf_matrix <- data %>% 
    filter(condit == condition) %>% 
    select(vxCat, vb) %>% 
    count(vxCat, vb) %>%
    spread(vb, n, fill = 0)
  
  total_trials <- sum(conf_matrix[,-1])
  conf_matrix[,-1] <- conf_matrix[,-1] / total_trials
  return(conf_matrix)
}

# Creating and normalizing confusion matrices for each condition
conf_matrix_constant <- create_normalized_conf_matrix(test, "Constant")
conf_matrix_varied <- create_normalized_conf_matrix(test, "Varied")










compute_metrics <- function(mat){
  # Accuracy
  accuracy <- sum(diag(mat)) / sum(mat)
  # Precision
  precision <- diag(mat) / rowSums(mat)
  # Recall
  recall <- diag(mat) / colSums(mat)
  # F1 score
  f1_score <- 2 * (precision * recall) / (precision + recall)
  return(list(accuracy = accuracy, precision = precision, recall = recall, f1_score = f1_score))
}

conf_matrix2 <- matrix(c(61, 514, 138, 593, 109, 295, 208, 137, 96, 121,
                         5, 200, 105, 564, 132, 521, 339, 197, 108, 123,
                         0, 60, 55, 396, 165, 558, 429, 280, 186, 157,
                         0, 1, 3, 55, 22, 172, 215, 179, 124, 129,
                         0, 5, 6, 29, 6, 88, 166, 186, 166, 222,
                         0, 7, 3, 25, 10, 63, 133, 172, 156, 296), nrow = 10)

rownames(conf_matrix2) <- c("0", "100", "btw1", "350", "btw2", "600", "800", "1000", "1200", ">=1401")
colnames(conf_matrix2) <- c("100-300", "350-550", "600-800", "800-1000", "1000-1200", "1200-1400")


calculate_metrics <- function(conf_matrix) {
  metrics_list <- list()
  for (col in 1:ncol(conf_matrix)) {
    TP <- conf_matrix[, col]
    FN <- sum(conf_matrix[, col]) - TP
    FP <- sum(conf_matrix[, ]) - TP
    TN <- sum(conf_matrix) - TP - FN - FP
    
    Accuracy <- (TP + TN) / (TP + TN + FP + FN)
    Precision <- TP / (TP + FP)
    Recall <- TP / (TP + FN)
    F1 <- 2 * ((Precision * Recall) / (Precision + Recall))
    
    metrics <- data.frame(Accuracy, Precision, Recall, F1)
    metrics_list[[col]] <- metrics
  }
  
  names(metrics_list) <- colnames(conf_matrix)
  return(metrics_list)
}


metrics1 <- calculate_metrics(conf_matrix1)
metrics2 <- calculate_metrics(conf_matrix2)


metrics1
metrics2

library(caret)
pred <- predict(biModel, test, type="response")
confusionMatrix(table(round(pred), test$result), positive = "Hit")


test$binary_result <- as.numeric(test$result == "Hit")
roc_obj <- roc(test$binary_result, test$vx) 
roc_obj  




library(pROC)
# Get predicted probabilities
pred_probs <- predict(biModel3, type = "response")
# Compute ROC curve
roc_obj <- roc(test$result == "Hit", pred_probs)
# Calculate AUC
auc(roc_obj)


