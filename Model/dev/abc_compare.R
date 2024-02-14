libraries <- c("dplyr", "purrr", "tidyr", "ggplot2", "data.table", 
               "conflicted", "here", "patchwork")
invisible(lapply(libraries, function(lib) {
  suppressPackageStartupMessages(require(lib, character.only = TRUE, quietly = TRUE))
}))


conflict_prefer_all("dplyr", quiet = TRUE)

ds <- readRDS(here::here("data/e1_md_11-06-23.rds"))  |> as.data.table()
dsv <- ds |> filter(condit=="Varied")  
dsc <- ds |> filter(condit=="Constant") 
tMax=84
avg_dsv <- ds |> filter(condit=="Varied",expMode2=="Train") |> group_by(tr) %>%
  mutate(bandInt2 = sample(rep(c(800, 1000, 1200), each = tMax / 3), tMax, replace = FALSE)[tr]) %>%
  filter(bandInt2 == x) |> select(-bandInt2) |> group_by(tr,condit,x,expMode2) |> summarise(y=mean(y),.groups="keep") |>
  rbind(dsv |> filter(expMode2=="Test") |> group_by(condit,x,expMode2) |> summarise(y=mean(y),tr=1,.groups="keep") ) |> setDT()


avg_dsc <- ds |> filter(condit=="Constant",expMode2=="Train",tr<=tMax) |> group_by(tr, condit,x,expMode2) |> 
  summarise(y=mean(y),.groups="keep") |> rbind(dsc |> filter(expMode2=="Test") |> group_by(condit,x,expMode2) |> summarise(y=mean(y),tr=1,.groups="keep") ) |> setDT()
 

# names(abc_500k)
# names(abc_500k[[1]])
# str(abc_500k[1:2])
#k = teter_combined [1,]$sim_dat

#abc_500k <- readRDS(here::here("data/abc_results_500k.rds"))
#abc_1M <- readRDS(here::here("data/abc_results_1M.rds"))




datasets <- list(
  abc_1M_p001 = readRDS(here::here("data/abc_1M_rmse_p001.rds")),
  abc_1M_p0001 = readRDS(here::here("data/abc_1M_rmse_p0001.rds")),
  abc_2M_p001 = readRDS(here::here("data/abc_2M_rmse_p001.rds")),
  abc_2M_p0001 = readRDS(here::here("data/abc_2M_rmse_p0001.rds")),
   abc_1p5M_p001 = readRDS(here::here("data/abc_1p5M_rmse_p001.rds")),
  abc_1p5M_p0001 = readRDS(here::here("data/abc_1p5M_rmse_p0001.rds"))
)


#######################


abc_1M <- readRDS(here::here("data/abc_1M_rmse_p001.rds"))
abc_1M <- readRDS(here::here("data/abc_1M_rmse_p0001.rds"))
abc_1M <- readRDS(here::here("data/abc_2M_rmse_p001.rds"))
abc_1M <- readRDS(here::here("data/abc_2M_rmse_p0001.rds"))
abc_1M <- readRDS(here::here("data/abc_1p5M_rmse_p001.rds"))
abc_1M <- readRDS(here::here("data/abc_1p5M_rmse_p0001.rds"))


teter_combined <- map_dfr(abc_1M, "teter_results", .id = "run_name") 
te_combined <- map_dfr(abc_1M, "te_results", .id = "run_name")
tr_combined <- map_dfr(abc_1M, "tr_results", .id = "run_name") 





teter_combined |> ggplot(aes(x=c,y=distance)) + geom_line(aes(color=Model)) +
  facet_wrap(~Group,scales = "free_y") 




# OLDER ######
############

abc_list <- readRDS('data/model_cache/abc_group_04_21_00.rds')
names(abc_list)


abc_v_exam <- abc_list[[1]]
names(abc_v_exam)

# extract unadjusted values from each item, and put with name of item
abc_v_exam_unadj <- map2_dfr(names(abc_v_exam), abc_v_exam, ~ .y$unadj.values |> as.data.frame() |> mutate(item = .x)) 
head(abc_v_exam_unadj)

# Process and combine all models into a single dataframe with pruned names
all_models_df <- purrr::imap_dfr(abc_list, ~ purrr::map2_dfr(
  names(.x), 
  .x, 
  ~ .y$unadj.values |> 
    as.data.frame() |> 
    mutate(Fit_Method = gsub("abc_", "", .x))
) |> 
  mutate(Model = gsub("abc_results_", "", .y))
)

all_models_df <- purrr::imap_dfr(abc_list, ~ {
  model_name <- gsub("abc_results_", "", .y)  # Prune the model name here
  purrr::map2_dfr(
    names(.x), 
    .x, 
    ~ {
      data_frame <- cbind(.y$unadj.values,dist=.y$dist[.y$region]) |> as.data.frame() 
      mutate(data_frame,
             Group = sub("_.*", "", model_name),  # Extracts everything before the underscore
             Model = sub("^[^_]*_", "", model_name),  # Extracts everything after the underscore
             Fit_Method = gsub("abc_", "", .x),
             .before="c") |> 
        mutate(
          n_param = .y$numparam,
          numstat=.y$numstat, 
          nsamples=length(.y$dist),
          method = .y$method
        ) |> group_by(Group,Model, Fit_Method) |>
        # find mind dist, and the c and lr values at that min dist
        mutate(min_dist = min(dist), 
               c_at_min_dist = c[which.min(dist)],
               lr_at_min_dist = lr[which.min(dist)], .before="method")
    }
  )
})
head(all_models_df)


all_models_df |> filter(Group=="varied",Fit_Method=="test") |> 
  pivot_longer(c:lr, names_to="param", values_to="value") |> 
  ggplot(aes(x=value)) + geom_density() + 
  facet_wrap(Model~param, scales="free") + theme_minimal() + 
  labs(x="Value", y="Density", title="Posterior Density Plots")

all_models_df |> filter(Group=="varied",Fit_Method=="train_test") |> 
  pivot_longer(c:lr, names_to="param", values_to="value") |> 
  ggplot(aes(x=value)) + geom_density() + 
  facet_wrap(Model~param, scales="free") + theme_minimal() + 
  labs(x="Value", y="Density", title="Posterior Density Plots")

all_models_df |> filter(Group=="constant",Fit_Method=="test") |> 
  pivot_longer(c:lr, names_to="param", values_to="value") |> 
  ggplot(aes(x=value)) + geom_density() + 
  facet_wrap(Model~param, scales="free") + theme_minimal() + 
  labs(x="Value", y="Density", title="Posterior Density Plots")

all_models_df |> filter(Group=="constant",Fit_Method=="train_test") |> 
  pivot_longer(c:lr, names_to="param", values_to="value") |> 
  ggplot(aes(x=value)) + geom_density() + 
  facet_wrap(Model~param, scales="free") + theme_minimal() + 
  labs(x="Value", y="Density", title="Posterior Density Plots")

find_mode <- function(x) {
  dens <- density(x)
  max_density_index <- which.max(dens$y)
  mode_value <- dens$x[max_density_index]
  return(mode_value)
}

calculate_HDI <- function(data, prob = 0.95) {
  hdi <- HDInterval::hdi(data, prob = prob)
  return(hdi)
}

# Apply the function to each combination
c_sum <- all_models_df |> 
  group_by(Group, Model, Fit_Method) |> 
  reframe(
    #HDI_c = HDInterval::hdi(c, .95),
    Median = median(c),
    mean=mean(c),
    mode=find_mode(c),
    IQR = IQR(c),
    sd=sd(c),
    Best_fit = min(dist),
    c_min= c[which.min(dist)], 
  )
c_sum
results <- all_models_df |> filter(min_dist==dist)


ggplot(all_models_df, aes(x=c, fill=Fit_Method)) +
  geom_density(alpha=0.5) +
  facet_grid(Model~Group) +
  theme_minimal() +
  labs(title="Density Plot of c by Group")





#######


abc_100k <- readRDS('data/model_cache/abc_group100k_03_32_46.rds')

all_models_df_100k <- purrr::imap_dfr(abc_100k, ~ {
  model_name <- gsub("abc_", "", .y)  # Prune the model name here
  relevant_lists <- .x[names(.x) %in% c("abc_train_test", "abc_test", "abc_train")]  # Filter relevant lists
  
  purrr::map2_dfr(
    names(relevant_lists), 
    relevant_lists, 
    ~ {
      data_frame <- cbind(.y$unadj.values, dist = .y$dist[.y$region]) |> as.data.frame() 
      mutate(data_frame,
             Group = sub("_.*", "", model_name),  # Extracts everything before the underscore
             Model = sub("^[^_]*_", "", model_name),  # Extracts everything after the underscore
             Fit_Method = gsub("abc_", "", .x)) |> 
        mutate(
          n_param = .y$numparam,
          numstat = .y$numstat, 
          nsamples = length(.y$dist),
          method = .y$method
        ) |> group_by(Group, Model, Fit_Method) |>
        # find min dist, and the c and lr values at that min dist
        mutate(min_dist = min(dist, na.rm = TRUE), 
               c_at_min_dist = c[which.min(dist)],
               lr_at_min_dist = lr[which.min(dist)], .before = "method")
    }
  )
})

head(all_models_df_100k)


all_models_df_100k |> filter(Group=="v",Fit_Method=="test") |> 
  pivot_longer(c:lr, names_to="param", values_to="value") |> 
  ggplot(aes(x=value)) + geom_density() + 
  facet_wrap(Model~param, scales="free") + theme_minimal() + 
  labs(x="Value", y="Density", title="Posterior Density Plots")

all_models_df_100k |> filter(Group=="varied",Fit_Method=="train_test") |> 
  pivot_longer(c:lr, names_to="param", values_to="value") |> 
  ggplot(aes(x=value)) + geom_density() + 
  facet_wrap(Model~param, scales="free") + theme_minimal() + 
  labs(x="Value", y="Density", title="Posterior Density Plots")

all_models_df_100k |> filter(Group=="constant",Fit_Method=="test") |> 
  pivot_longer(c:lr, names_to="param", values_to="value") |> 
  ggplot(aes(x=value)) + geom_density() + 
  facet_wrap(Model~param, scales="free") + theme_minimal() + 
  labs(x="Value", y="Density", title="Posterior Density Plots")

all_models_df_100k |> filter(Group=="c",Fit_Method=="train_test") |> 
  pivot_longer(c:lr, names_to="param", values_to="value") |> 
  ggplot(aes(x=value)) + geom_density() + 
  facet_wrap(Model~param, scales="free") + theme_minimal() + 
  labs(x="Value", y="Density", title="Posterior Density Plots")

c_sum100k <- all_models_df_100k |> 
  group_by(Group, Model, Fit_Method) |> 
  reframe(
    #HDI_c = HDInterval::hdi(c, .95),
    Median = median(c),
    mean=mean(c),
    mode=find_mode(c),
    IQR = IQR(c),
    sd=sd(c),
    Best_fit = min(dist),
    c_min= c[which.min(dist)], 
  )
c_sum100k



##############################

abc_500k <- readRDS('data/model_cache/abc_group500k_06_17_54.rds')

all_models_df_500k <- purrr::imap_dfr(abc_500k, ~ {
  model_name <- gsub("abc_", "", .y)  # Prune the model name here
  relevant_lists <- .x[names(.x) %in% c("abc_train_test", "abc_test", "abc_train")]  # Filter relevant lists
  
  purrr::map2_dfr(
    names(relevant_lists), 
    relevant_lists, 
    ~ {
      data_frame <- cbind(.y$unadj.values, dist = .y$dist[.y$region]) |> as.data.frame() 
      mutate(data_frame,
             Group = sub("_.*", "", model_name),  # Extracts everything before the underscore
             Model = sub("^[^_]*_", "", model_name),  # Extracts everything after the underscore
             Fit_Method = gsub("abc_", "", .x)) |> 
        mutate(
          n_param = .y$numparam,
          numstat = .y$numstat, 
          nsamples = length(.y$dist),
          method = .y$method
        ) |> group_by(Group, Model, Fit_Method) |>
        # find min dist, and the c and lr values at that min dist
        mutate(min_dist = min(dist, na.rm = TRUE), 
               c_at_min_dist = c[which.min(dist)],
               lr_at_min_dist = lr[which.min(dist)], .before = "method")
    }
  )
})

head(all_models_df_500k)


all_models_df_500k |> filter(Group=="v",Fit_Method=="test") |> 
  pivot_longer(c:lr, names_to="param", values_to="value") |> 
  ggplot(aes(x=value)) + geom_density() + 
  facet_wrap(Model~param, scales="free") + theme_minimal() + 
  labs(x="Value", y="Density", title="Posterior Density Plots")

all_models_df_500k |> filter(Group=="v",Fit_Method=="train_test") |> 
  pivot_longer(c:lr, names_to="param", values_to="value") |> 
  ggplot(aes(x=value)) + geom_density() + 
  facet_wrap(Model~param, scales="free") + theme_minimal() + 
  labs(x="Value", y="Density", title="Posterior Density Plots")

all_models_df_500k |> filter(Group=="c",Fit_Method=="test") |> 
  pivot_longer(c:lr, names_to="param", values_to="value") |> 
  ggplot(aes(x=value)) + geom_density() + 
  facet_wrap(Model~param, scales="free") + theme_minimal() + 
  labs(x="Value", y="Density", title="Posterior Density Plots")

all_models_df_500k |> filter(Group=="c",Fit_Method=="train_test") |> 
  pivot_longer(c:lr, names_to="param", values_to="value") |> 
  ggplot(aes(x=value)) + geom_density() + 
  facet_wrap(Model~param, scales="free") + theme_minimal() + 
  labs(x="Value", y="Density", title="Posterior Density Plots")

c_sum500k <- all_models_df_500k |> 
  group_by(Group, Model, Fit_Method) |> 
  reframe(
    #HDI_c = HDInterval::hdi(c, .95),
    Median = median(c),
    mean=mean(c),
    mode=find_mode(c),
    IQR = IQR(c),
    sd=sd(c),
    Best_fit = min(dist),
    c_min= c[which.min(dist)], 
  )
c_sum500k








update_columns <- function(df) {
  df <- df %>%
    mutate(
      Group = case_when(
        Model %in% c("abc_ev", "abc_almv", "abc_altv") ~ "Varied",
        Model %in% c("abc_ec", "abc_almc", "abc_altc") ~ "Constant",
        TRUE                                           ~ "Unknown"  
      ),
      Model = case_when(
        Model == "abc_ev"   ~ "EXAM",
        Model == "abc_almv" ~ "ALM",
        Model == "abc_altv" ~ "Alt_Exam",
        Model == "abc_ec"   ~ "EXAM",
        Model == "abc_almc" ~ "ALM",
        Model == "abc_altc" ~ "Alt_Exam",
        TRUE                ~ Model  
      )
    )
  return(df)
}