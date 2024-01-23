

abc_tables <- function(post_dat_l,ids=NULL){

  if(is.null(ids)) ids <- c(1,66,36)
  
et_sum <- post_dat_l |>
  group_by(id,condit, Fit_Method, Resp) |>
  summarise(val = mean(val), .groups = 'drop') |>
  pivot_wider(
    names_from = Resp,
    values_from = val,
    values_fill = list(val = NA)
  ) |>
  mutate(
    ALM_error = round(abs(ALM - Observed),1),
    EXAM_error = round(abs(EXAM - Observed),1),
    Best_Model = case_when(
      ALM_error < EXAM_error ~ "ALM",
      EXAM_error < ALM_error ~ "EXAM",
      TRUE ~ NA_character_  # In case of a tie or missing data
    )
  ) |>
  group_by(condit, Fit_Method) %>%
  summarise(
    Avg_ALM_error = round(mean(ALM_error, na.rm = TRUE),1),
    Avg_EXAM_error = round(mean(EXAM_error, na.rm = TRUE),1),
    N_Best_ALM = sum(Best_Model == "ALM", na.rm = TRUE),
    N_Best_EXAM = sum(Best_Model == "EXAM", na.rm = TRUE)
  ) %>%
  mutate(
    Best_Model = case_when(
      Avg_ALM_error < Avg_EXAM_error ~ "ALM",
      Avg_EXAM_error < Avg_ALM_error ~ "EXAM",
      TRUE ~ "Tie"  # In case of a tie or missing data
    )
  )

et_sum_x <- post_dat_l |>
  group_by(condit, Fit_Method, Resp, x) |>
  summarise(val = mean(val), .groups = 'drop') |>
  pivot_wider(
    names_from = Resp,
    values_from = val,
    values_fill = list(val = NA)
  ) |>
  group_by(condit, Fit_Method, x) |>
  transmute(
    ALM_error = round(abs(ALM - Observed),1),
    EXAM_error = round(abs(EXAM - Observed),1),
    Best_Model = case_when(
      ALM_error < EXAM_error ~ "ALM",
      EXAM_error < ALM_error ~ "EXAM",
      TRUE ~ NA_character_  # In case of a tie or missing data
    )
  ) 


et_sum_x_indv <- post_dat_l |> filter(id %in% ids) |>
  group_by(id,condit, Fit_Method, Resp, x) |>
  summarise(val = mean(val), .groups = 'drop') |>
  pivot_wider(
    names_from = Resp,
    values_from = val,
    values_fill = list(val = NA)
  ) |>
  group_by(id,condit, Fit_Method, x) |>
  transmute(
    ALM_error = round(abs(ALM - Observed),1),
    EXAM_error = round(abs(EXAM - Observed),1),
    Best_Model = case_when(
      ALM_error < EXAM_error ~ "ALM",
      EXAM_error < ALM_error ~ "EXAM",
      TRUE ~ NA_character_  # In case of a tie or missing data
    )
  ) 

  
  return(tibble::lst(et_sum,et_sum_x, et_sum_x_indv))


}





group_predictive_plots <- function(post_dat_l)
{

agg_posterior <- post_dat_l |> group_by(id,condit,Fit_Method, Resp,x) |>
  summarise(val=mean(val)) |>
  ggplot(aes(x = x, y = val, fill=Resp)) + 
  stat_bar + 
  facet_wrap(condit~Fit_Method, scales="free") + 
  labs(title = "Full Posterior Averages")

best_posterior <- post_dat_l |> group_by(id,condit,Fit_Method, Resp,x) |>
  filter(rank==1) |>
  summarise(val=mean(val)) |>
  ggplot(aes(x = x, y = val, fill=Resp)) + 
  stat_bar + 
  facet_wrap(condit~Fit_Method, scales="free") + 
  labs(title = "Predictions from Best Parameters")

agg_posterior / best_posterior

}


group_best_plots <- function(post_dat_l)
{

et <- post_dat_l |>
  group_by(id, condit, Fit_Method, Resp, x) |>
  summarise(val = mean(val), .groups = 'drop') |>
  pivot_wider(
    names_from = Resp,
    values_from = val,
    values_fill = list(val = NA)
  ) |>
  mutate(
    ALM_error = abs(ALM - Observed),
    EXAM_error = abs(EXAM - Observed),
    Best_Model = case_when(
      ALM_error < EXAM_error ~ "ALM",
      EXAM_error < ALM_error ~ "EXAM",
      TRUE ~ NA_character_  # In case of a tie or missing data
    )
  )

  fpb <- ggplot(et, aes(x = factor(x), fill = Best_Model)) +
  geom_bar(position = "dodge") +
  facet_wrap(condit~Fit_Method, scales = "free") +
  labs(
    title = "Frequency of Best Model - Full Posterior",
    x = "x Value",
    y = "Frequency",
    fill = "Best Model"
  ) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


et <- post_dat_l |> filter(rank==1) |>
  group_by(id, condit, Fit_Method, Resp, x) |>
  summarise(val = mean(val), .groups = 'drop') |>
  pivot_wider(
    names_from = Resp,
    values_from = val,
    values_fill = list(val = NA)
  ) |>
  mutate(
    ALM_error = abs(ALM - Observed),
    EXAM_error = abs(EXAM - Observed),
    Best_Model = case_when(
      ALM_error < EXAM_error ~ "ALM",
      EXAM_error < ALM_error ~ "EXAM",
      TRUE ~ NA_character_  # In case of a tie or missing data
    )
  )

bp<- ggplot(et, aes(x = factor(x), fill = Best_Model)) +
  geom_bar(position = "dodge") +
  facet_wrap(condit~Fit_Method, scales = "free") +
  labs(
    title = "Frequency of Best Model - Best Parameters",
    x = "x Value",
    y = "Frequency",
    fill = "Best Model"
  ) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

fpb / bp

}



indv_best_plots <- function(post_dat_l)
{

et <- post_dat_l |>
  group_by(id, condit, Fit_Method, Resp, x) |>
  summarise(val = mean(val), .groups = 'drop') |>
  pivot_wider(
    names_from = Resp,
    values_from = val,
    values_fill = list(val = NA)
  ) |>
  mutate(
    ALM_error = abs(ALM - Observed),
    EXAM_error = abs(EXAM - Observed),
    Best_Model = case_when(
      ALM_error < EXAM_error ~ "ALM",
      EXAM_error < ALM_error ~ "EXAM",
      TRUE ~ NA_character_  # In case of a tie or missing data
    )
  )

  fpb <- ggplot(et, aes( x = factor(x), y = id, fill = Best_Model)) +
  geom_tile(color = "white") +
  facet_wrap(condit~Fit_Method,scales="free") +
  labs(title = "Best Model for Each ID and X Value - Full Posterior",
       x = "X",
       y = "ID",
       fill = "Best Model") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


et <- post_dat_l |> filter(rank==1) |>
  group_by(id, condit, Fit_Method, Resp, x) |>
  summarise(val = mean(val), .groups = 'drop') |>
  pivot_wider(
    names_from = Resp,
    values_from = val,
    values_fill = list(val = NA)
  ) |>
  mutate(
    ALM_error = abs(ALM - Observed),
    EXAM_error = abs(EXAM - Observed),
    Best_Model = case_when(
      ALM_error < EXAM_error ~ "ALM",
      EXAM_error < ALM_error ~ "EXAM",
      TRUE ~ NA_character_  # In case of a tie or missing data
    )
  )

bp<- ggplot(et, aes( x = factor(x), y = id, fill = Best_Model)) +
  geom_tile(color = "white") +
  facet_wrap(condit~Fit_Method,scales="free") +
  labs(title = "Best Model for Each ID and X Value - Best Parameters",
       x = "X",
       y = "ID",
       fill = "Best Model") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

fpb / bp

}


indv_predictive_plots <- function(post_dat_l, ids=NULL)
{
if(is.null(ids)) ids <- c(1,66,36)

r_all <- post_dat_l |> filter(id %in% ids) |> mutate(dist="All")
r1 <- post_dat_l |> filter(id %in% ids, rank==1) |> mutate(dist="best")
r_10 <- post_dat_l |> filter(id %in% ids, rank<=10) |> mutate(dist="top10")

d <- bind_rows(r_all,r1)

d |> ggplot(aes(x = x, y = val, fill=Resp)) + 
  stat_bar + 
  ggh4x::facet_nested_wrap(id~Fit_Method~dist, scales="free",ncol=6) + 
  labs(title = "Individual Predictive Distributions")

}


plot_sampled_posterior <- function(ind_fits)
{
  sampled_posterior <- ind_fits$runInfo$kde_results |>
  imap_dfr(~map_dfr(.x, "kde_samples", .id = "Fit_Method"), .id = "mod_condit") |>
  filter(mod_condit != "abc_altv", mod_condit != "abc_altc") |>
  mutate(
    Model = case_when(
      str_detect(mod_condit, "e") ~ "Exam",
      str_detect(mod_condit, "alm") ~ "ALM",
      TRUE ~ as.character(NA) # in case there are other unexpected values
    ),
    condit = case_when(
      str_detect(mod_condit, "v") ~ "Varied",
      str_detect(mod_condit, "c") ~ "Constant",
      TRUE ~ as.character(NA) 
    )
  ) |> 
  select(-mod_condit) 
#sampled_posterior |> group_by(Model, condit) |> summarise(n=n())
pc <- sampled_posterior |> 
  ggplot(aes(x=c)) + geom_density(aes(fill=condit), alpha=.5) + 
  facet_wrap(Fit_Method~Model,scales="free",ncol=2) + 
  theme_bw() + theme(legend.position = "bottom") + labs(title = "Sampled Posterior of c") 

plr <- sampled_posterior |> 
  ggplot(aes(x=lr)) + geom_density(aes(fill=condit), alpha=.5) + 
  facet_wrap(Fit_Method~Model,scales="free",ncol=2) + 
  theme_bw() + theme(legend.position = "bottom") + labs(title = "Sampled Posterior  of lr") 

pc + plr
  
}


plot_indv_posterior <- function(ind_df)
{
plt_c_gi <- ind_df |> 
  ggplot(aes(x=c)) + geom_density(aes(fill=Group), alpha=.5) + 
  facet_wrap(Fit_Method~Model,scales="free",ncol=2) + 
  theme_bw() + theme(legend.position = "bottom") + labs(title = "Posterior distribution of c") 

plt_lr_gi <- ind_df |> 
  ggplot(aes(x=lr)) + geom_density(aes(fill=Group), alpha=.5) + 
  facet_wrap(Fit_Method~Model,scales="free",ncol=2) + 
  theme_bw() + theme(legend.position = "bottom") + labs(title = "Posterior distribution of lr") 
  plt_c_gi + plt_lr_gi
}



extract_info <- function(raw_names) {
  full_name <- raw_names
  n_samp <- str_extract(raw_names, "(?<=_)\\d+(?=_)")
  ng_value <- str_extract(raw_names, "(?<=ng)\\d+")
  buf_value <- str_replace(str_extract(raw_names, "(?<=buf)[0-9p]+"), "p", ".")
  run_type <- ifelse(str_detect(raw_names, "ss"), "ss", "trial")
  
  data.frame(full_name, n_samp, ng_value, buf_value, run_type)
}

repair_names <- function(names) {
  dupes <- names[duplicated(names)]
  if ("rank" %in% dupes) {
    first_rank <- which(names == "rank")[1]
    names <- make.unique(names, sep = "_")
    names[first_rank] <- "rank"
  }
  names
}