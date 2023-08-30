
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



# takes a long time to fit
mn_model1 <- brm(vxCat ~ bandInt + condit + (1 | id), family = categorical, data=test)
# 
# mn_model1
# Family: categorical 
# Links: mu100 = logit; mubtw1 = logit; mu350 = logit; mubtw2 = logit; mu600 = logit; mu800 = logit; mu1000 = logit; mu1200 = logit; mu1401 = logit 
# Formula: vxCat ~ bandInt + condit + (1 | id) 
# Data: test (Number of observations: 9491) 
# Draws: 4 chains, each with iter = 2000; warmup = 1000; thin = 1;
# total post-warmup draws = 4000
# 
# Group-Level Effects: 
#   ~id (Number of levels: 156) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(mu100_Intercept)      2.21      0.19     1.86     2.61 1.00      933     1899
# sd(mubtw1_Intercept)     1.37      0.16     1.09     1.69 1.01      901     2047
# sd(mu350_Intercept)      0.92      0.09     0.76     1.11 1.01      761     1644
# sd(mubtw2_Intercept)     0.54      0.10     0.34     0.73 1.01      756     1330
# sd(mu600_Intercept)      0.32      0.07     0.17     0.45 1.01      494      570
# sd(mu800_Intercept)      0.58      0.07     0.44     0.73 1.01      581     1105
# sd(mu1000_Intercept)     0.94      0.09     0.78     1.12 1.00      718     1622
# sd(mu1200_Intercept)     1.33      0.11     1.13     1.56 1.00      897     1532
# sd(mu1401_Intercept)     1.78      0.14     1.52     2.07 1.00     1069     1803
# 
# Population-Level Effects: 
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# mu100_Intercept         1.20      0.36     0.50     1.90 1.01     1149     2173
# mubtw1_Intercept       -0.04      0.30    -0.63     0.56 1.00     2093     2799
# mu350_Intercept         1.61      0.25     1.12     2.11 1.00     1595     2391
# mubtw2_Intercept       -0.28      0.26    -0.78     0.22 1.00     1947     2347
# mu600_Intercept         0.80      0.24     0.33     1.26 1.00     1764     2368
# mu800_Intercept        -0.16      0.25    -0.65     0.32 1.00     1684     2426
# mu1000_Intercept       -1.39      0.27    -1.93    -0.86 1.00     1788     2465
# mu1200_Intercept       -2.60      0.30    -3.21    -2.00 1.00     1789     2373
# mu1401_Intercept       -3.48      0.35    -4.17    -2.80 1.00     1348     2274
# mu100_bandInt           0.00      0.00    -0.00     0.00 1.00     1027     1309
# mu100_conditVaried     -0.06      0.45    -0.94     0.83 1.01      970     1652
# mubtw1_bandInt          0.01      0.00     0.00     0.01 1.00     1004     1334
# mubtw1_conditVaried     0.32      0.36    -0.37     1.03 1.00     1495     2050
# mu350_bandInt           0.01      0.00     0.00     0.01 1.00     1008     1330
# mu350_conditVaried      0.30      0.29    -0.26     0.87 1.00     1340     2150
# mubtw2_bandInt          0.01      0.00     0.01     0.01 1.00     1017     1306
# mubtw2_conditVaried     0.44      0.28    -0.11     0.99 1.00     1505     2360
# mu600_bandInt           0.01      0.00     0.01     0.01 1.00     1008     1347
# mu600_conditVaried      0.31      0.26    -0.19     0.81 1.00     1395     2017
# mu800_bandInt           0.01      0.00     0.01     0.01 1.00      997     1315
# mu800_conditVaried      0.51      0.27    -0.02     1.04 1.00     1411     1796
# mu1000_bandInt          0.01      0.00     0.01     0.01 1.00     1010     1391
# mu1000_conditVaried     0.91      0.30     0.36     1.50 1.00     1344     1990
# mu1200_bandInt          0.01      0.00     0.01     0.01 1.00     1015     1371
# mu1200_conditVaried     1.19      0.34     0.55     1.84 1.00     1391     2314
# mu1401_bandInt          0.01      0.00     0.01     0.02 1.00     1013     1402
# mu1401_conditVaried     1.24      0.38     0.49     2.00 1.00      996     1678
# 
# Draws were sampled using sample(hmc). For each parameter, Bulk_ESS
# and Tail_ESS are effective sample size measures, and Rhat is the potential
# scale reduction factor on split chains (at convergence, Rhat = 1).

mn_model1 |> linpred_draws(new_data_grid) |> 
  summarize(across(.linpred, lst(mean, sd, median), .names = "{.fn}")) |> print(n=54)
# A tibble: 54 × 8
# Groups:   id, condit, bandInt, .row [6]
# id    condit bandInt  .row .category   mean    sd median
# <fct> <fct>    <dbl> <int> <fct>      <dbl> <dbl>  <dbl>
#   1 1     Varied     100     4 1         -0.938 1.46  -0.770
# 2 1     Varied     100     4 2          0.271 0.829  0.321
# 3 1     Varied     100     4 3          2.20  0.454  2.20 
# 4 1     Varied     100     4 4          0.713 0.467  0.710
# 5 1     Varied     100     4 5          2.15  0.326  2.14 
# 6 1     Varied     100     4 6          1.84  0.382  1.83 
# 7 1     Varied     100     4 7          0.364 0.476  0.365
# 8 1     Varied     100     4 8         -1.27  0.684 -1.23 
# 9 1     Varied     100     4 9         -0.763 0.523 -0.750
# 10 1     Varied     350     5 1         -0.415 1.49  -0.259
# 11 1     Varied     350     5 2          1.55  0.867  1.60 
# 12 1     Varied     350     5 3          3.77  0.517  3.77 
# 13 1     Varied     350     5 4          2.66  0.528  2.66 
# 14 1     Varied     350     5 5          4.31  0.407  4.30 
# 15 1     Varied     350     5 6          4.33  0.455  4.32 
# 16 1     Varied     350     5 7          3.12  0.538  3.12 
# 17 1     Varied     350     5 8          1.69  0.722  1.72 
# 18 1     Varied     350     5 9          2.45  0.576  2.45 
# 19 1     Varied     600     6 1          0.108 1.57   0.256
# 20 1     Varied     600     6 2          2.83  0.995  2.82 
# 21 1     Varied     600     6 3          5.34  0.703  5.33 
# 22 1     Varied     600     6 4          4.60  0.711  4.58 
# 23 1     Varied     600     6 5          6.48  0.623  6.45 
# 24 1     Varied     600     6 6          6.82  0.658  6.80 
# 25 1     Varied     600     6 7          5.88  0.718  5.87 
# 26 1     Varied     600     6 8          4.65  0.860  4.65 
# 27 1     Varied     600     6 9          5.65  0.744  5.65 
# 28 1     Varied     800     2 1          0.526 1.67   0.654
# 29 1     Varied     800     2 2          3.85  1.14   3.85 
# 30 1     Varied     800     2 3          6.59  0.891  6.56 
# 31 1     Varied     800     2 4          6.15  0.898  6.12 
# 32 1     Varied     800     2 5          8.21  0.827  8.17 
# 33 1     Varied     800     2 6          8.81  0.855  8.78 
# 34 1     Varied     800     2 7          8.09  0.903  8.06 
# 35 1     Varied     800     2 8          7.02  1.02   7.00 
# 36 1     Varied     800     2 9          8.22  0.922  8.21 
# 37 1     Varied    1000     1 1          0.944 1.79   1.05 
# 38 1     Varied    1000     1 2          4.88  1.31   4.85 
# 39 1     Varied    1000     1 3          7.85  1.09   7.79 
# 40 1     Varied    1000     1 4          7.71  1.10   7.66 
# 41 1     Varied    1000     1 5          9.95  1.04   9.90 
# 42 1     Varied    1000     1 6         10.8   1.07  10.8  
# 43 1     Varied    1000     1 7         10.3   1.11  10.3  
# 44 1     Varied    1000     1 8          9.39  1.20   9.35 
# 45 1     Varied    1000     1 9         10.8   1.12  10.8  
# 46 1     Varied    1200     3 1          1.36  1.94   1.46 
# 47 1     Varied    1200     3 2          5.90  1.50   5.86 
# 48 1     Varied    1200     3 3          9.10  1.31   9.03 
# 49 1     Varied    1200     3 4          9.26  1.32   9.20 
# 50 1     Varied    1200     3 5         11.7   1.26  11.6  
# 51 1     Varied    1200     3 6         12.8   1.28  12.7  
# 52 1     Varied    1200     3 7         12.5   1.32  12.5  
# 53 1     Varied    1200     3 8         11.8   1.39  11.7  
# 54 1     Varied    1200     3 9         13.4   1.33  13.3  



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


