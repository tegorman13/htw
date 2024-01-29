

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


indv_predictive_dist <- function(post_dat_l, ind_fits_df, sbj=1)
{
  
  p1 <- indv_predictive_plots(post_dat_l, sbj) 
  
  p2 <-  ind_fits_df |> filter(id %in% sbj) |> 
    group_by(id, Fit_Method, Model) |>
    arrange(mean_error) |> mutate(rank=row_number()) |>
    arrange(id,rank, Fit_Method, Model) |>
    group_by(id, Fit_Method, Model) |>
    mutate(cBest=round(first(c),5), 
           flab = paste0(Fit_Method, "\n Best c: ", cBest)) |>
    ggplot(aes(x=c, Fill=Model))+geom_density(alpha=.5) + 
    ggh4x::facet_nested_wrap(~Model+flab, scales="free",ncol=3) +
    labs(title = "Posterior distribution of c")  
  
  p3 <-  ind_fits_df |> filter(id %in% sbj) |> 
    group_by(id, Fit_Method, Model) |>
    arrange(mean_error) |> mutate(rank=row_number()) |>
    arrange(id,rank, Fit_Method, Model) |>
    group_by(id, Fit_Method, Model) |>
    mutate(lrBest=round(first(lr),4), 
           flab = paste0(Fit_Method, "\n Best lr: ", lrBest)) |>
    ggplot(aes(x=lr, Fill=Model))+geom_density(alpha=.5) + 
    ggh4x::facet_nested_wrap(~Model+flab, scales="free",ncol=3) +
    labs(title = "Posterior distribution of lr") 
  
  p1 / p2 / p3  
  
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



####################



samp_priors <- function(n,cMean=-5,cSig=1.5,lrSig=1) {
  prior_samples <- tibble(
    c = rlnorm(n,cMean,sdlog=cSig),
    lr = extraDistr::rhnorm(n,sigma=lrSig),
  )
  return(prior_samples)
}

reject_abc <- function(simulation_function, prior_samples, data, num_iterations = 5000, n_try=500, return_dat="test_data") {
  input_layer =  c(100,350,600,800,1000,1200)
  output_layer = input_layer
  data <- data |> as.data.table()
  target_data <- case_when(
    return_dat == "test_data" ~ list(target_data_test <- data[data$expMode2 == "Test", ]), 
    return_dat == "train_data" ~ list( target_data_train <- data[data$expMode2 == "Train", ]),
    return_dat == "train_data, test_data" ~ list( target_data_train_test <- data[expMode2 %in% c("Test", "Train"), ])
    ) |> pluck(1)

  tol = target_data |> group_by(x) |> summarise(m=mean(y),sd=sd(y)) |> summarise(tol=mean(sd),.groups="drop") *tolM
  abc <- list()
  try_count=0;
  t1=system.time({
  for(j in 1:num_iterations) {
    #print(i)
    #print(paste0("sbj: ",data$id[1] ,"Particle:", j))

    found=0;
    
    inc_count=0;
    while(found==0) {
    try_count=try_count+1;
    current_theta <- prior_samples[sample(1:nrow(prior_samples), 1), ]
    sim_data <- simulation_function(data, current_theta$c, current_theta$lr, input_layer = input_layer, output_layer = output_layer, return_dat = return_dat) 
    dist_sd <- target_data |> mutate(pred=sim_data,error=abs(y-pred)) |> group_by(id,condit,x) |> 
      summarise(mean_error=mean(error),sd_error=sd(error),.groups="keep") |> 
      group_by(id,condit) |> 
      summarise(mean_error=mean(mean_error),sd_error=mean(sd_error),.groups="keep") 
    if(dist_sd$mean_error< tol) {
      #print(current_theta)
      abc$dist_sd[[j]] <- cbind(current_theta,dist_sd,tol,inc_count)
      found=1
     # try_count=0;
      }
      if (try_count>=n_try && j<2){
        message(print(paste0("increase tol for subject", data$id[1])))
        tol=tol*1.05
        inc_count=inc_count+1;
        try_count=0;
      }
    }
  }  
  })

  abc <- rbindlist(abc$dist_sd) |> arrange(mean_error) |> relocate(id,condit)
  best <- abc |> head(1) |> round_tibble(7)
  message((paste0( data$id[1], " completed in: ", t1[3])))
  message((paste0("Best c: ", best$c, " Best lr: ", best$lr, " Best error: ", best$mean_error)))
  return(abc)
}
