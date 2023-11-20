

pacman::p_load(dplyr,tidyr,purrr,here,future,furrr)
purrr::walk(here::here(c("Functions/fit_funs.R", "Functions/alm_core.R","Functions/misc_model_funs.R")),source)

ds <- readRDS(here::here("data/e1_md_11-06-23.rds")) |> mutate(sbj=id) |> relocate(sbj,.after=id)




analyze_model_performance <- function(combined_data) {
  
  combined_data$scenario <- factor(combined_data$scenario, 
                                   levels = c("Default", "Plus 3 Nodes", "Plus 12 Nodes", "Plus 30 Nodes"),
                                   ordered = TRUE)
  
  summary_stats <- combined_data %>%
    group_by(scenario, condit) %>%
    summarise(
      mean_actual = mean(y),
      sd_actual = sd(y),
      mean_predicted = mean(pred),
      sd_predicted = sd(pred),
      mean_deviation = mean(dev),
      sd_deviation = sd(dev),
      .groups = 'drop'
    )
  print(summary_stats)
  
  # Detailed Statistics by Scenario, Condition, and Node
  detailed_stats <- combined_data %>%
    group_by(scenario, condit, x) %>%
    summarise(
      mean_actual = mean(y),
      sd_actual = sd(y),
      mean_predicted = mean(pred),
      sd_predicted = sd(pred),
      mean_deviation = mean(dev),
      sd_deviation = sd(dev),
      .groups = 'drop'
    )
  print(detailed_stats)
  

  overall_performance <- combined_data %>%
    group_by(condit, scenario) |>
    mutate(Average=mean(dev)) |>
    group_by(condit,scenario,x,Average) %>%
    summarise(
       dev = mean(dev),
      .groups = 'drop'
    ) |> pivot_wider(values_from = dev, names_from = x)
  print(overall_performance)
  
  
  overall_performance <- combined_data %>%
    group_by(scenario) %>%
    summarise(
      total_mean_deviation = mean(dev),
      total_sd_deviation = sd(dev),
      .groups = 'drop'
    )
  print(overall_performance)
  
  # ANOVA Test (Optional)
  anova_test <- aov(dev ~ scenario, data = combined_data)
  print(summary(anova_test))
  
  # Plotting
  p1 <- ggplot(combined_data, aes(x = y, y = pred, color = condit)) +
    geom_point() +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    facet_wrap(~ scenario + condit) +
    theme_minimal() +
    labs(title = "Model Fits Across Different Conditions",
         x = "Actual Values",
         y = "Predicted Values",
         color = "Condition")
  print(p1)
  
  p2 <- ggplot(combined_data, aes(x = scenario, y = dev, fill = condit)) +
    geom_boxplot() +
    facet_wrap(~ condit) +
    theme_minimal() +
    labs(title = "Distribution of Deviations Across Scenarios and Conditions",
         x = "Scenario",
         y = "Deviation")
  print(p2)
}






model_params = tibble::lst(input.layer = c(100,350,600,800,1000,1200), 
                           output.layer = c(100,350,600,800,1000,1200), 
                           trainVec=c(0,800,1000,1200))

initial_params <- c(c = 0.01, lr = .5, sigma = 200) 
#any(!is.finite(initial_params))


split_data <- split(ds, ds$id)
split_data <- ds %>% split(.$id) |> head(20)


opt.m="Nelder-Mead"
fit_params=list(pred_dat="test_avg", pred_fun="exam.response", loss_fun="nll", loss_data="test_error")

nm_ex_te_avg<- future_map(split_data, ~ optim_fit_id(.x, func="fit_lr_c0", fit_params, model_params,opt.m, initial_params))
nm_ex_te_avg[[1]]$test
nm_ex_te_avg[[1]]$Fit




generate_model_params <- function(additional_input_nodes = NULL, additional_output_nodes = NULL) {
  basis_nodes <- c(100, 350, 600, 800, 1000, 1200)
  input_layer <- unique(c(basis_nodes, additional_input_nodes)) |> sort()
  output_layer <- unique(c(basis_nodes, additional_output_nodes)) |> sort()
  tibble::lst(
    input.layer = input_layer, 
    output.layer = output_layer, 
    trainVec = c(0, 800, 1000, 1200)
  )
}

fit_model <- function(model_params, split_data, opt.m, initial_params) {
  fit_params <- list(pred_dat="test_data", pred_fun="exam.response", loss_fun="nll2", loss_data="test_error+train_error")
  future_map(split_data, ~ optim_fit_id(.x, func="fit_lr_c0", fit_params, model_params, opt.m=opt.m, initial_params))
}

fit_and_summarize <- function(split_data, opt.m, initial_params, additional_input_nodes = NULL, additional_output_nodes = NULL) {
  model_params <- generate_model_params(additional_input_nodes, additional_output_nodes)
  fit_results <- fit_model(model_params, split_data, opt.m, initial_params)
  summarize_results(fit_results)
}

summarize_results <- function(model_fit) {
  model_fit |> map_dfr("test") |> 
    group_by(id, condit, x) |> 
    summarise(y=mean(y), pred=mean(pred), dev=y-pred)
}


additional_output_nodes_1 <- c(300, 550, 1400)  # Scenario 2
additional_output_nodes_2 <- seq(100, 1200, length.out = 18)[-c(1, 4, 7, 10, 13, 14)]  # 8 additional nodes for Scenario 3
additional_output_nodes_3 <- seq(100, 1200, length.out = 36)[-c(1, 5, 9, 13, 17, 18)]  # 12 additional nodes for Scenario 4


plan(multisession)
default_ol_summary <- fit_and_summarize(split_data, opt.m, initial_params)
default_ol_summary$scenario <- "Default"

plan(multisession)
large_ol_summary_1 <- fit_and_summarize(split_data, opt.m, initial_params, 
                                        additional_input_nodes = additional_output_nodes_1,
                                        additional_output_nodes = additional_output_nodes_1)
large_ol_summary_1$scenario <- "Plus 3 Nodes"

plan(multisession)
large_ol_summary_2 <- fit_and_summarize(split_data, opt.m, initial_params, 
                                        additional_input_nodes = additional_output_nodes_2,
                                        additional_output_nodes = additional_output_nodes_2)
large_ol_summary_2$scenario <- "Plus 12 Nodes"

plan(multisession)
large_ol_summary_3 <- fit_and_summarize(split_data, opt.m, initial_params, 
                                        additional_input_nodes = additional_output_nodes_3,
                                        additional_output_nodes = additional_output_nodes_3)
large_ol_summary_3$scenario <- "Plus 30 Nodes"


combined_data <- rbind(default_ol_summary, large_ol_summary_1, large_ol_summary_2, large_ol_summary_3)


analyze_model_performance(combined_data)



ggplot(combined_data, aes(x = y, y = pred, color = condit)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +  # Adds a y=x reference line
  facet_wrap(~ scenario + condit) +
  theme_minimal() +
  labs(title = "Model Fits Across Different Conditions",
       x = "Actual Values",
       y = "Predicted Values",
       color = "Condition")






###############




summarize_results <- function(model_fit) {
  model_fit |> map_dfr("test") |> 
    group_by(id, condit, x) |> 
    summarise(y=mean(y), pred=mean(pred), dev=y-pred)
}



model_params_default <- tibble::lst(input.layer = c(100,350,600,800,1000,1200), 
                                    output.layer = c(100,350,600,800,1000,1200), 
                                    trainVec=c(0,800,1000,1200))

default_ol_fit <- fit_model(model_params_default, split_data, opt.m, initial_params)
default_ol_summary <- summarize_results(default_ol_fit)



model_params_large <- tibble::lst(input.layer = c(100,350,600,800,1000,1200), 
                                  output.layer = c(100,300,350,550,600,800,1000,1200,1400), 
                                  trainVec=c(0,800,1000,1200))


large_ol_fit <- fit_model(model_params_large, split_data, opt.m, initial_params)
large_ol_summary <- summarize_results(large_ol_fit)


default_ol_summary |> group_by(x, condit) |> summarise(y=mean(y), pred=mean(pred), dev=y-pred)
large_ol_summary |> group_by(x, condit) |> summarise(y=mean(y), pred=mean(pred), dev=y-pred)



default_ol_summary$layer_size <- "Default"
large_ol_summary$layer_size <- "Large"

combined_data <- rbind(default_ol_summary, large_ol_summary)


ggplot(combined_data, aes(x = y, y = pred, color = condit)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +  # Adds a y=x reference line
  facet_wrap(~ layer_size + condit) +
  theme_minimal() +
  labs(title = "Model Fits Across Different Conditions",
       x = "Actual Values",
       y = "Predicted Values",
       color = "Condition")





###########
model_params = tibble::lst(input.layer = c(100,350,600,800,1000,1200), 
                           output.layer = c(100,350,600,800,1000,1200), 
                           trainVec=c(0,800,1000,1200))

plan(multisession)
opt.m="Nelder-Mead"
fit_params=list(pred_dat="test_data", pred_fun="exam.response", loss_fun="nll", loss_data="test_error")
nm_ex_te_raw <- future_map(split_data, ~ optim_fit_id(.x, func="fit_lr_c0", fit_params, model_params,opt.m, initial_params))



nm_ex_te_raw[[1]]$test
nm_ex_te_raw[[1]]$Fit
nm_ex_te_raw[[1]]$test |> group_by(x) |> summarise(y=mean(y),pred=mean(pred))



model_params = tibble::lst(input.layer = c(100,350,600,800,1000,1200), 
                           output.layer = c(100,300, 350,550, 600,800,1000,1200, 1400), 
                           trainVec=c(0,800,1000,1200))


fit_params=list(pred_dat="test_data", pred_fun="exam.response", loss_fun="nll", loss_data="test_error")
nm_ex_te_raw2 <- future_map(split_data, ~ optim_fit_id(.x, func="fit_lr_c0", fit_params, model_params,opt.m, initial_params))


nm_ex_te_raw2[[1]]$test
nm_ex_te_raw2[[1]]$Fit
nm_ex_te_raw2[[1]]$test |> group_by(x) |> summarise(y=mean(y),pred=mean(pred))