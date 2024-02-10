

pacman::p_load(dplyr,purrr,tidyr,ggplot2, data.table, here, patchwork, conflicted, stringr,future,furrr, knitr, reactable)
conflict_prefer_all("dplyr", quiet = TRUE)
options(scipen = 999)
walk(c("Display_Functions","fun_alm","fun_indv_fit"), ~ source(here::here(paste0("Functions/", .x, ".R"))))
ds <- readRDS(here::here("data/e1_md_11-06-23.rds"))  |> as.data.table()


# should consider binning velocoties in 50s or 100s. 
# simpler version that just orders the extrapolation items
# compute bin widths as a function of SD or SE of each position


get_order_pattern <- function(df, pred_col = "pred",band_tol=50) {

  df <- df[order(df$x), ]
  # Calculate the standard deviation and determine the pattern
  pred_sd <- sd(df[[pred_col]]) / 2.8
  pattern <- ""
  # Iterate through the rows to build the pattern
  for (i in 1:(nrow(df) - 1)) {
    if (abs(df[[pred_col]][i] - df[[pred_col]][i + 1]) <= band_tol) {
      pattern <- paste0(pattern, "=")
    } else if (df[[pred_col]][i] < df[[pred_col]][i + 1]) {
      pattern <- paste0(pattern, "<")
    } else {
      pattern <- paste0(pattern, ">")
    }
  }

  return(pattern)
}


ind_ds <- ds |> 
  filter(expMode2 == "Test") |> 
  group_by(id, condit, x, expMode2) |> 
  summarise(y = mean(y), .groups = "keep") 

ind_pattern <- ind_ds |> 
  group_by(id, condit, expMode2) |>
  do(order_pattern = get_order_pattern(., pred_col = "y",band_tol=200)) 
  
ind_ds <- ind_ds |>
  left_join(ind_pattern, by=join_by(id,condit,expMode2)) |> 
  mutate(order_pattern = unlist(order_pattern))

lorder <- ind_ds |>
  group_by(order_pattern, x) |>
  summarize(y = mean(y), n_id = n_distinct(id)) %>% 
  pull(n_id) %>% 
  unique() %>% 
  sort(decreasing = TRUE)


ord_pattern <- ind_fits_df  |> filter(rank<=1) |>
  filter(expMode2=="Test") |> 
  group_by(condit,rank,Model,Fit_Method) %>% 
  do(order_pattern = unlist(get_order_pattern(.))) 





combo_pred <- ind_fits_df |>  filter(rank<=1) |>
  left_join(ord_pattern, by=join_by(condit,rank,Model,Fit_Method)) |>
  mutate(order_pattern = unlist(order_pattern))

head(combo_pred)
length(unique(combo_pred$order_pattern))


pattern_freq <- combo_pred %>%
  count(order_pattern) %>%
  arrange(desc(n))

# Print the frequency table
print(pattern_freq)

# Plot the frequency distribution
ggplot(pattern_freq, aes(x = order_pattern, y = n)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Frequency Distribution of Order Patterns",
       x = "Order Pattern",
       y = "Frequency")

combo_pred %>% group_by(condit,Model,Fit_Method) |> count(order_pattern) |> arrange(desc(n)) |> ggplot(aes(x = order_pattern, y = n)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Frequency Distribution of Order Patterns",
       x = "Order Pattern",
       y = "Frequency") +
  ggh4x::facet_nested_wrap(condit~Model)


combo_pred %>%
  group_by(Model, Fit_Method) %>%
  mutate(n_pattern = n_distinct(order_pattern)) %>%
  group_by(Model, Fit_Method, order_pattern) %>%
  mutate(pat_count = n()) %>%
  ungroup() %>%
  group_by(Model, Fit_Method) %>%
  arrange(Model, Fit_Method, desc(pat_count), order_pattern) %>%
  mutate(pat_rank = dense_rank(desc(pat_count))) |> 
  group_by(Model, Fit_Method, order_pattern) %>%
  # distinct(Model, Fit_Method, order_pattern, .keep_all = TRUE) %>%
  # arrange(Model, Fit_Method, pat_rank)
  filter(pat_rank<=5) |>
  filter(condit=="Varied")  |> 
 ggplot(aes(x=c,y=lr,color=as.factor(order_pattern))) + geom_point() + facet_wrap(Fit_Method~Model,scales="free") 

  
ordered_labels <- ind_ds |>
  group_by(order_pattern, x,condit) |>
  summarize(y = mean(y), n_id = n_distinct(id)) |>
  arrange(-n_id, order_pattern) 


ordered_labels$factor_label <- factor(ordered_labels$factor_label, levels =ordered_labels |>  pull(factor_label) %>% unique() )
levels(ordered_labels$factor_label)


ordered_labels |>
  ggplot(aes(x = x, y = y),fill=wes_palette("Darjeeling1")[1]) +
  stat_summary(fun = mean, geom = "bar") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.25) +
  facet_wrap(~factor_label, scales = "free_x")

k |> filter(condit=="Varied") |>
  ggplot(aes(x = x, y = y,fill=condit),fill=wes_palette("Darjeeling1")[1]) +
  stat_summary(fun = mean, geom = "bar", position = position_dodge(width = 0.8)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", position = position_dodge(width = 0.8), width = 0.25) +
  facet_wrap(~factor_label, scales = "free_x")

k |> filter(condit=="Constant") |>
  ggplot(aes(x = x, y = y),fill=wes_palette("Darjeeling1")[2]) +
  stat_summary(fun = mean, geom = "bar", position = position_dodge(width = 0.8)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", position = position_dodge(width = 0.8), width = 0.25) +
  facet_wrap(~factor_label, scales = "free_x")
