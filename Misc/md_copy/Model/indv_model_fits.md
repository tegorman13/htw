

```r

runif(10)

```


```{r}

runif(10)

```


```r


pacman::p_load(tidyverse,patchwork,here, pander, latex2exp, flextable)
purrr::walk(here::here(c("Functions/Display_Functions.R", "Functions/alm_core.R","Functions/misc_model_funs.R")),source)
select <- dplyr::select; mutate <- dplyr::mutate 




ind_ex_te[[1]]$test
ind_ex_te[[1]]$Fit
ind_ex_te[[1]]$id



ind_ex_te=as.list(readRDS((here::here(paste0("data/model_cache/indv_nll_de2_ex_te15_26_16.rds")))))



# full combination of files, Model_Names, and Fit_Method - lined by by index. 
file_list <- list("indv_nll_de2_ex_te15_26_16", "indv_nll_de2_ex_tetr15_31_48", "indv_exam_de2_ex_tr15_37_27", "indv_nll_de2_alm_te15_42_48","indv_nll_de2_alm_tetr15_48_08", "indv_nll_de2_alm_tr15_53_32")
Model_Names <- list("EXAM", "EXAM", "EXAM","ALM","ALM","ALM")
Fit_Method <- list("Test Only", "Test & Train", "Train Only", "Test Only", "Test & Train", "Train Only")



compile_fit_info <- function(file_name, model_name, fit_method) {
  file_path <- here::here(paste0("data/model_cache/", file_name, ".rds"))
  data_list <- as.list(readRDS(file_path))

  map_dfr(data_list, ~ {
    test_df <- .x$test |> slice(1)
    fit_info <- .x$Fit |> slice(1)

    tibble(
      id = test_df$sbj,  
      condit = test_df$condit,  
      Model = model_name,
      Fit_Method = fit_method,
      fit_info,  
      train_error = .x$train_error,
      test_error = .x$test_error
    )
  })
}

indv_fit_info <- map2_dfr(file_list, Model_Names, ~ compile_fit_info(.x, .y, Fit_Method[[which(file_list == .x)]])) |> arrange(id,Fit_Method,Model)

indv_fit_info <- indv_fit_info %>% 
  group_by(id, Fit_Method) %>% 
  mutate(best_model = Model[which.min(value)]) %>%
  ungroup()







map_dfr(ind_ex_te,~tibble( pluck(.x$"test"))) |> select(id=sbj,condit)
map_dfr(ind_ex_te,~tibble( pluck(.x$"test"))) |> pull(sbj)


# example creating one
# ind_ex_te <- map_dfr(as.list(readRDS((here::here(paste0("data/model_cache/indv_nll_de2_ex_te15_26_16.rds"))))), ~ tibble(id = .x$id,condist=.x$condit, pluck(.x$"Fit"), train_error=.x$train_error, test_error=.x$test_error))


# params <- map_dfr(ind_ex_te,~pluck(.x$"Fit"))
# ids <- map_dfr(ind_ex_te, ~ tibble(id = .x$id))

i47 <- ind_ex_te[[47]]
i47$
```


```{r}

indv_fit_info |> ggplot(aes(x=Model,y=value,group=condit,fill=condit)) +geom_bar(position="dodge",stat="identity") + facet_wrap(~Fit_Method, scales="free_y")


summary_data <- indv_fit_info %>%
  group_by(Fit_Method, condit, best_model) %>%
  summarize(Count = n(), .groups = 'drop')

# Create the plot
ggplot(summary_data, aes(x = Fit_Method, y = Count, fill = best_model)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_wrap(~condit) +
  labs(title = "Best Model Fit by Condition and Fit Method",
       x = "Fit Method",
       y = "Number of Subjects",
       fill = "Best Model") +
  theme_minimal()

```



```{r}

compile_test_info <- function(file_name, model_name, fit_method) {
  file_path <- here::here(paste0("data/model_cache/", file_name, ".rds"))
  data_list <- as.list(readRDS(file_path))
   test <- map_dfr(data_list,~tibble( pluck(.x$"test"), Model = model_name,Fit_Method = fit_method))

}




indv_test <- map2_dfr(file_list, Model_Names, ~ compile_test_info(.x, .y, Fit_Method[[which(file_list == .x)]])) |> rename(id=sbj)

head(indv_test)

indv_test |> pivot_longer(c(y,pred), names_to="Resp", values_to = "vx") |>
ggplot(aes(x=x, y=vx,fill=condit)) +
    stat_summary(geom = "bar", position=position_dodge(), fun = mean) +
    stat_summary(geom = "errorbar", position=position_dodge(), fun.data = mean_se, alpha = .7) +  
    facet_wrap(Model~Resp~Fit_Method, scales="free_y")


 indv_test |> pivot_longer(c(y,pred), names_to="Resp", values_to = "vx")   |> filter(Fit_Method=="Test & Train") |>
 ggplot(aes(x=x, y=vx,fill=condit)) +
    stat_summary(geom = "bar", position=position_dodge(), fun = mean) +
    stat_summary(geom = "errorbar", position=position_dodge(), fun.data = mean_se, alpha = .7) +  
    facet_wrap(Model~Resp, scales="free_y")

 indv_test |> pivot_longer(c(y,pred), names_to="Resp", values_to = "vx")   |> filter(Fit_Method=="Test & Train") |>
 ggplot(aes(x=x, y=vx,fill=Resp)) +
    stat_summary(geom = "bar", position=position_dodge(), fun = mean) +
    stat_summary(geom = "errorbar", position=position_dodge(), fun.data = mean_se, alpha = .7) +  
    facet_wrap(Model~condit, scales="free_y")


 indv_test |> pivot_longer(c(y,pred), names_to="Resp", values_to = "vx")   |> filter(Fit_Method=="Train Only") |>
 ggplot(aes(x=x, y=vx,fill=Resp)) +
    stat_summary(geom = "bar", position=position_dodge(), fun = mean) +
    stat_summary(geom = "errorbar", position=position_dodge(), fun.data = mean_se, alpha = .7) +  
    facet_wrap(Model~condit, scales="free_y")

 indv_test |> pivot_longer(c(y,pred), names_to="Resp", values_to = "vx")   |> filter(Fit_Method=="Test Only") |>
 ggplot(aes(x=x, y=vx,fill=Resp)) +
    stat_summary(geom = "bar", position=position_dodge(), fun = mean) +
    stat_summary(geom = "errorbar", position=position_dodge(), fun.data = mean_se, alpha = .7) +  
    facet_wrap(Model~condit, scales="free_y")



#test <- map_dfr(ind_ex_te, ~pluck(.x$"test"))
#ind_ex_te[[1]]$test



```



```{r}


indv_fit_info |> filter(condit=="Constant",Fit_Method=="Test & Train") |> arrange(value)

idC <- indv_fit_info |> filter(condit=="Constant",Fit_Method=="Test & Train", Model=="EXAM") |> 
  arrange(test_error) |> slice(1:3)

 indv_test |> pivot_longer(c(y,pred), names_to="Resp", values_to = "vx") |> 
   filter(id %in% idC$id, Model=="EXAM",Fit_Method=="Test & Train") |>  
    ggplot(aes(x=x, y=vx,fill=Resp)) +
    stat_summary(geom = "bar", position=position_dodge(), fun = mean) +
    stat_summary(geom = "errorbar", position=position_dodge(), fun.data = mean_se, alpha = .7) +  
    facet_wrap(Fit_Method~id, scales="free_y")

idV <- indv_fit_info |> filter(condit=="Varied",Fit_Method=="Test Only", Model=="EXAM") |> 
  arrange(value) |> slice(1:3)

 indv_test |> pivot_longer(c(y,pred), names_to="Resp", values_to = "vx") |> 
   filter(id %in% idV$id, Model=="EXAM",Fit_Method=="Test Only") |>  
    ggplot(aes(x=x, y=vx,fill=Resp)) +
    stat_summary(geom = "bar", position=position_dodge(), fun = mean) +
    stat_summary(geom = "errorbar", position=position_dodge(), fun.data = mean_se, alpha = .7) +  
    facet_wrap(Fit_Method~id, scales="free_y")


```



```{r}

# find all files in model_cache dir that contain strings "_nm_" AND "_nll_", save paths in "model_files" list 

compile_fit_info <- function(file_name) {
  file_path <- here::here("data/model_cache", file_name)
  map_dfr(readRDS(file_path), ~ tibble(
    id = .x$test$sbj[1],  
    condit = .x$test$condit[1],  
    Fit_Method = .x$loss_dat,
    Model = .x$pred_fun,
     .x$Fit,  
    train_error = .x$train_error,
    test_error = .x$test_error
  )) %>% 
  mutate(
    Fit_Method = case_when(
      Fit_Method == "train_error"            ~ "Train Only",
      Fit_Method == "test_error"             ~ "Test Only",
      Fit_Method == "test_error+train_error" ~ "Test & Train"
    ),
    Model = case_when(
      Model == "exam.response"    ~ "EXAM",
      Model == "alm.responseOnly" ~ "ALM"
    )
  )
}

model_files <- list.files(here::here("data/model_cache"), pattern = "_nm_")
indv_fit_info_nm <- map_dfr(model_files, compile_fit_info) |> arrange(id)

model_files <- list.files(here::here("data/model_cache"), pattern = "_bfgs_")
indv_fit_info_bfgs <- map_dfr(model_files, compile_fit_info) |> arrange(id)

model_files <- list.files(here::here("data/model_cache"), pattern = "_sann_")
indv_fit_info_sann <- map_dfr(model_files, compile_fit_info) |> arrange(id)

model_files <- list.files(here::here("data/model_cache"), pattern = "_de2_", full.names = TRUE) %>%
    file.info() %>%
    arrange(desc(mtime)) %>%
    rownames() %>%
    head(6) |> basename()

indv_fit_info_de2 <- map_dfr(model_files, compile_fit_info) |> rename(Value=value) |> arrange(id)




combined_data <- bind_rows(
  indv_fit_info_nm %>% mutate(Optimization = "Nelder-Mead"),
  indv_fit_info_bfgs %>% mutate(Optimization = "BFGS"),
  indv_fit_info_sann %>% mutate(Optimization = "SANN"),
  indv_fit_info_de2 %>% mutate(Optimization = "DEoptim")
)



head(indv_fit_info_bfgs)
head(indv_fit_info_sann)
head(indv_fit_info_de2)





```


```{r}


# Grouped Bar Chart for Train and Test Errors
ggplot(combined_data, aes(x = Optimization, y = train_error, fill = Optimization)) +
  geom_bar(stat = "summary", fun = "mean", position = position_dodge()) +
  geom_bar(aes(y = test_error), stat = "summary", fun = "mean", position = position_dodge(), color = "blue") +
  labs(title = "Comparison of Train and Test Errors Across Optimization Methods", 
       y = "Error", 
       x = "Optimization Method") +
  theme_minimal()


long_data <- combined_data %>%
  pivot_longer(cols = c(Value, test_error, train_error), names_to = "Metric", values_to = "Value")

long_data |> pandoc.table()


# Create the plot
ggplot(long_data, aes(x = Optimization, y = Value, fill = Optimization)) +
  geom_bar(stat = "summary", fun = "median", position = position_dodge()) +
  facet_wrap(Fit_Method~Metric, scales = "free_y") +
  labs(title = "Comparison of Optimization Methods Across Different Metrics",
       y = "Value",
       x = "Optimization Method") +
  theme_minimal() +
  theme(legend.position = "bottom")

long_data <- combined_data %>%
  pivot_longer(cols = c(Value, test_error, train_error), names_to = "Metric", values_to = "Value") %>%
  group_by(Optimization, Metric) %>%
  summarize(
    Mean = mean(Value, na.rm = TRUE),
    Median = median(Value, na.rm = TRUE)
  ) %>%
  pivot_longer(cols = c(Mean, Median), names_to = "Statistic", values_to = "Value") %>%
  ungroup()

# Create the plot
ggplot(long_data, aes(x = Optimization, y = Value, fill = Optimization)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_grid(Statistic ~ Metric, scales = "free_y") +
  labs(title = "Comparison of Optimization Methods Across Different Metrics",
       y = "Value",
       x = "Optimization Method") +
  theme_minimal() +
  theme(legend.position = "bottom")




reshaped_df <- combined_data %>%
  select(Optimization, Fit_Method, Model, train_error, test_error) %>%
  pivot_longer(cols = c(train_error, test_error), names_to = "Parameter") %>%
  unite("Group", Optimization, Parameter, sep = "_", na.rm = TRUE) %>%  # Ensure no NAs in combined names
  pivot_wider(
    names_from = Group, 
    values_from = value, 
    values_fn = list(value = ~ mean(.x, na.rm = TRUE))  # Handle NAs in mean calculation
  ) %>%
  filter(complete.cases(.))  # Optional: Remove rows with any NA values

ft <- flextable(reshaped_df)

reshaped_df <- combined_data %>%
  select(Optimization, Fit_Method, Model, train_error, test_error) %>%
  pivot_longer(cols = c(train_error, test_error), names_to = "Parameter") %>%
  group_by(Optimization, Fit_Method, Model, Parameter) %>%
  summarize(Value = mean(value, na.rm = TRUE), .groups = 'drop') %>%
  arrange(Fit_Method, Model, Parameter, Optimization)

# Check the structure of reshaped_df
str(reshaped_df)

# Create a flextable from the reshaped data
ft <- flextable(reshaped_df)
ft





```

