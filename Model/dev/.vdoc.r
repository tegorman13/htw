#
#
#
#

# Load required libraries
pacman::p_load(tidyverse, data.table, abc, future, furrr, here, patchwork, conflicted)



# Define function to process dataset
process_dataset <- function(name, dataset) {
  sample_size <- str_extract(name, "\\d+p?\\d*?M")
  posterior_cutoff <- str_extract(name, "p\\d+$")

  teter_combined <- map_dfr(dataset, "teter_results", .id = "run_name")
  te_combined <- map_dfr(dataset, "te_results", .id = "run_name")
  tr_combined <- map_dfr(dataset, "tr_results", .id = "run_name")

  plot_title_suffix <- sprintf("%s - %s", sample_size, posterior_cutoff)

  return(list(teter_combined = teter_combined,
              te_combined = te_combined,
              tr_combined = tr_combined,
              plot_title_suffix = plot_title_suffix))
}



create_density_plots <- function(processed_data) {
  teter_combined <- processed_data$teter_combined
  te_combined <- processed_data$te_combined
  tr_combined <- processed_data$tr_combined
  plot_title_suffix <- processed_data$plot_title_suffix

  # Base plot for density
  base_plot <- function(data, x_var, title, show_legend = FALSE) {
    ggplot(data, aes_string(x=x_var, color="Model")) +
      geom_density() +
      facet_wrap(~Group, scales = "free") +
      ggtitle(title) +
      theme(legend.position=ifelse(show_legend, "top", "none"))  # Toggle legend display
  }

  # Create plots, but only show legend on the first plot
  plot1 <- base_plot(teter_combined, "c", sprintf("c posterior - Test & Train - %s", plot_title_suffix),TRUE )
  plot2 <- base_plot(teter_combined, "lr", sprintf("lr posterior - Test & Train - %s", plot_title_suffix))
  plot3 <- base_plot(te_combined, "c", sprintf("c posterior - Test Only - %s", plot_title_suffix))
  plot4 <- base_plot(te_combined, "lr", sprintf("lr posterior - Test Only - %s", plot_title_suffix))
  plot5 <- base_plot(tr_combined, "c", sprintf("c posterior - Train Only - %s", plot_title_suffix))
  plot6 <- base_plot(tr_combined, "lr", sprintf("lr posterior - Train Only - %s", plot_title_suffix))

  # Combine plots and place the one with the legend (plot1) at the end so that the legend is not cut off
  density_combined <- plot1 + plot2 + plot3 + plot4 + plot5 + plot6 + 
    plot_layout(ncol = 2) 
  return(density_combined)
}


create_distance_plots <- function(processed_data) {
  teter_combined <- processed_data$teter_combined
  te_combined <- processed_data$te_combined
  tr_combined <- processed_data$tr_combined
  plot_title_suffix <- processed_data$plot_title_suffix

  # Base plot for distance
  base_dist_plot <- function(data, title, show_legend = FALSE) {
    ggplot(data, aes(x=Group, y=distance, fill=Model)) +
      stat_summary(fun=mean, geom="bar", position=position_dodge(width = 0.8)) +
      stat_summary(fun.data=mean_se, geom="errorbar", position=position_dodge(width = 0.8), width=0.25) +
      ggtitle(title) +
      theme(legend.position=ifelse(show_legend, "top", "none")) # Control the display of the legend
  }

  # Create distance plots
  dist_plot1 <- base_dist_plot(teter_combined, sprintf("Test & Train - %s", plot_title_suffix),TRUE)
  dist_plot2 <- base_dist_plot(te_combined, sprintf("Test Only - %s", plot_title_suffix),)
  dist_plot3 <- base_dist_plot(tr_combined, sprintf("Train Only - %s", plot_title_suffix))

  # Combine plots and place the one with the legend (dist_plot1) at the end
  distance_combined <- dist_plot1 + dist_plot2 + dist_plot3 + 
    plot_layout(ncol = 2) 
  return(distance_combined)
}

# Define function to save plots
save_plots <- function(density_plot, distance_plot, name) {
  sample_size <- str_extract(name, "\\d+p?\\d*?M")
  posterior_cutoff <- str_extract(name, "p\\d+$")

  density_filename <- paste0("assets/tmp_plots/density_plots_combo_abc_", sample_size, "_rmse_", posterior_cutoff, ".png")
  distance_filename <- paste0("assets/tmp_plots/distance_plots_combo_abc_", sample_size, "_rmse_", posterior_cutoff, ".png")

  ggsave(filename = here::here(density_filename), plot = density_plot, width = 10, height = 6)
  ggsave(filename = here::here(distance_filename), plot = distance_plot, width = 10, height = 6)
}

# Process datasets
processed_datasets <- lapply(names(datasets), function(name) process_dataset(name, datasets[[name]]))

# Create density plots
density_plots <- invisible(lapply(processed_datasets, create_density_plots))

# Create distance plots
distance_plots <- lapply(processed_datasets, create_distance_plots)

density_plots
distance_plots


# Save plots
invisible(lapply(seq_along(names(datasets)), function(i) save_plots(density_plots[[i]], distance_plots[[i]], names(datasets)[i])))

#
#
#
#
#

pacman::p_load(tidyverse, data.table, abc, future, furrr, here, patchwork, conflicted)

# Function to process dataset
process_dataset <- function(name, dataset) {
  sample_size <- str_extract(name, "\\d+p?\\d*?M")
  posterior_cutoff <- str_extract(name, "p\\d+$")

  teter_combined <- map_dfr(dataset, "teter_results", .id = "run_name")
  te_combined <- map_dfr(dataset, "te_results", .id = "run_name")
  tr_combined <- map_dfr(dataset, "tr_results", .id = "run_name")

  list(teter_combined, te_combined, tr_combined)
}

# Function to create and save density plots
create_density_plots <- function(data_combined, name) {
  sample_size <- str_extract(name, "\\d+p?\\d*?M")
  posterior_cutoff <- str_extract(name, "p\\d+$")
  plot_title_suffix <- sprintf("%s - %s", sample_size, posterior_cutoff)

  plot_list <- lapply(data_combined, function(data) {
    ggplot(data, aes(x = c)) + 
      geom_density(aes(color = Model)) +
      facet_wrap(~Group, scales = "free") + 
      ggtitle(sprintf("c posterior - %s - %s", data$context[1], plot_title_suffix))
  })

  density_combined <- wrap_plots(plot_list, ncol = 2)
  density_filename <- sprintf("assets/tmp_plots/density_plots_combo_abc_%srmse%s.png", sample_size, posterior_cutoff)
  ggsave(filename = here::here(density_filename), plot = density_combined, width = 10, height = 6)
}

# Function to create and save distance plots
create_distance_plots <- function(data_combined, name) {
  sample_size <- str_extract(name, "\\d+p?\\d*?M")
  posterior_cutoff <- str_extract(name, "p\\d+$")
  plot_title_suffix <- sprintf("%s - %s", sample_size, posterior_cutoff)

  dist_plot_list <- lapply(data_combined, function(data) {
    ggplot(data, aes(x = Group, y = distance, fill = Model)) +
      stat_summary(fun = mean, geom = "bar", position = position_dodge()) +
      stat_summary(fun.data = mean_se, geom = "errorbar", position = position_dodge()) +
      ggtitle(sprintf("%s - %s", data$context[1], plot_title_suffix))
  })

  distance_combined <- wrap_plots(dist_plot_list, ncol = 2)
  distance_filename <- sprintf("assets/tmp_plots/distance_plots_combo_abc_%srmse%s.png", sample_size, posterior_cutoff)
  ggsave(filename = here::here(distance_filename), plot = distance_combined, width = 10, height = 6)
}

# Load datasets
datasets <- list(
  abc_1M_p001 = readRDS(here::here("data/abc_1M_rmse_p001.rds")),
  abc_1M_p0001 = readRDS(here::here("data/abc_1M_rmse_p0001.rds")),
  abc_2M_p001 = readRDS(here::here("data/abc_2M_rmse_p001.rds")),
  abc_2M_p0001 = readRDS(here::here("data/abc_2M_rmse_p0001.rds")),
  abc_1p5M_p001 = readRDS(here::here("data/abc_1p5M_rmse_p001.rds")),
  abc_1p5M_p0001 = readRDS(here::here("data/abc_1p5M_rmse_p0001.rds"))
)

# Process and plot for each dataset
invisible(lapply(names(datasets), function(name) {
  processed_data <- process_dataset(name, datasets[[name]])
  create_density_plots(processed_data, name)
  create_distance_plots(processed_data, name)
}))

#
#
#
#
#
#
#
#
#
pacman::p_load(dplyr,purrr,tidyr,ggplot2,data.table,future,furrr,here,patchwork, conflicted)
conflict_prefer_all("dplyr", quiet = TRUE)
walk(c("fun_alm","fun_model"), ~ source(here::here(paste0("Functions/", .x, ".R"))))

ds <- readRDS(here::here("data/e1_md_11-06-23.rds"))  |> as.data.table()
dsv <- ds |> filter(condit=="Varied")  
dsc <- ds |> filter(condit=="Constant") 

## prepare group aggregated data
tMax=84
avg_dsv <- ds |> filter(condit=="Varied",expMode2=="Train") |> group_by(tr) %>%
  mutate(bandInt2 = sample(rep(c(800, 1000, 1200), each = tMax / 3), tMax, replace = FALSE)[tr]) %>%
  filter(bandInt2 == x) |> select(-bandInt2) |> group_by(tr,condit,x,expMode2) |> summarise(y=mean(y),.groups="keep") |>
  rbind(dsv |> filter(expMode2=="Test") |> group_by(condit,x,expMode2) |> summarise(y=mean(y),tr=1,.groups="keep") ) |> setDT()

avg_dsc <- ds |> filter(condit=="Constant",expMode2=="Train",tr<=tMax) |> group_by(tr, condit,x,expMode2) |> 
  summarise(y=mean(y),.groups="keep") |> rbind(dsc |> filter(expMode2=="Test") |> group_by(condit,x,expMode2) |> 
  summarise(y=mean(y),tr=1,.groups="keep") ) |> setDT()



# Cognitive Model helper functions 
## Organize train and test data, simulate ALM learning, and then test with given model. 

full_sim_exam <- function(data, c, lr,pred_fun=exam.response, input_layer, output_layer,return_dat="test_data",mode="sim") {
  train_data <- data[expMode2=="Train", c("condit","tr","expMode2", "x","y")] 
  test_data <- data[expMode2=="Test", c("condit","tr","expMode2", "x","y")] 
  trainVec=sort(unique(train_data$x))
  if (train_data$condit[1] != "Varied") {
    trainVec <- c(0, trainVec)
  }
  
  train_results <- alm.sim(train_data, c, lr, input_layer, output_layer)
  
  test_prediction <- map_dbl(test_data$x, ~ pred_fun(.x, c, input_layer, 
                                                     output_layer, train_results$wm,  trainVec=trainVec))
  
  train_data$pred <- train_results$d$almResp
  test_data$pred <- test_prediction
  
  fd = eval(parse(text=paste0("rbind(",return_dat,")")))
  if(mode=="sim"){return(fd$pred)
  }else {return(fd)}
  
}

full_sim_alm <- function(data, c, lr,pred_fun=alm.responseOnly, input_layer, output_layer,return_dat="test_data",mode="sim") {
  train_data <- data[expMode2=="Train", c("condit","tr","expMode2", "x","y")] 
  test_data <- data[expMode2=="Test", c("condit","tr","expMode2", "x","y")] 
  trainVec=sort(unique(train_data$x))
  
  train_results <- alm.sim(train_data, c, lr, input_layer, output_layer)
  test_prediction <- map_dbl(test_data$x, ~ alm.responseOnly(.x, c, input_layer, output_layer, train_results$wm,  trainVec=trainVec))
  
  train_data$pred <- train_results$d$almResp
  test_data$pred <- test_prediction
  
  #fd <- rbind(train_data,test_data)
  fd = eval(parse(text=paste0("rbind(",return_dat,")")))
  if(mode=="sim"){return(fd$pred)
  }else {return(fd)}
  
}

full_sim_alt_exam <- function(data, c, lr, pred_fun=alt_exam, input_layer, output_layer,return_dat="test_data",mode="sim") {
  train_data <- data[expMode2=="Train", c("condit","tr","expMode2", "x","y")] 
  test_data <- data[expMode2=="Test", c("condit","tr","expMode2", "x","y")] 
  trainVec=sort(unique(train_data$x))
  trainVecY=train_data$y
  
  train_results <- alm.sim(train_data, c, lr, input_layer, output_layer)
  test_prediction <- map_dbl(test_data$x, ~ alt_exam(.x, c, input_layer, output_layer, train_results$wm,  trainVecX=trainVec, trainVecY=trainVecY))
  
  train_data$pred <- train_results$d$almResp
  test_data$pred <- test_prediction
  
  fd = eval(parse(text=paste0("rbind(",return_dat,")")))
  if(mode=="sim"){return(fd$pred)
  }else {return(fd)}
  
}

# generate simulated dataset. 
sim_data_gen <- function(data, input_layer, output_layer, simulation_function, prior_samples, return_dat) {
  simulation_function <- match.fun(simulation_function)
  plan(multisession)
  suppressMessages(simulation_results <- future_map_dfc(seq_len(nrow(prior_samples)), function(idx) {
    params <- prior_samples[idx, ]
    simulation_function(data=as.data.table(data), c=params$c, lr=params$lr, input_layer=input_layer, output_layer=output_layer, return_dat = return_dat)
  }, .options = furrr_options(seed = TRUE)))
}


generate_prior_c_lr <- function(n) {
  prior_samples <- tibble(
    c = runif(n, 0.000000001, 3.5),
    lr = runif(n, 0.000001, 7.5),
  )
  return(prior_samples)
}

input_layer =  c(100,350,600,800,1000,1200)
output_layer = input_layer
n_prior_samples=1500000
prior_samples <- generate_prior_c_lr(n_prior_samples)
return_dat="train_data, test_data"

# Create a list of arguments for each sim_data_gen call
args_list <- list(
  Exam_Varied=list(avg_dsv, input_layer, output_layer, full_sim_exam, prior_samples, return_dat),
  Exam_Constant=list(avg_dsc, input_layer, output_layer, full_sim_exam, prior_samples, return_dat),
  ALM_Varied=list(avg_dsv, input_layer, output_layer, full_sim_alm, prior_samples, return_dat),
  ALM_Constant=list(avg_dsc, input_layer, output_layer, full_sim_alm, prior_samples, return_dat),
  Alt_Varied=list(avg_dsv, input_layer, output_layer, full_sim_alt_exam, prior_samples, return_dat),
  Alt_Constant=list(avg_dsc, input_layer, output_layer, full_sim_alt_exam, prior_samples, return_dat)
)

sim_data_wrapper <- function(args) {
  sim_data_gen(args[[1]], args[[2]], args[[3]], simulation_function = args[[4]], args[[5]], return_dat = args[[6]])
}


t1=system.time({
plan(multisession)
sim_dataAll <- future_map(args_list, sim_data_wrapper, .progress = TRUE,.options = furrr_options(seed = TRUE))
})
t1

# save simulated datasets
saveRDS(tibble::lst(sim_dataAll,prior_samples,args_list,time=t1[3]),
        file = here::here(paste0("data/sim_data/","sim_data_", n_prior_samples,"_",format(Sys.time(),"%H_%M_%OS"), ".rds")))



# load simulated datasets
sd <- readRDS(here::here("/Users/thomasgorman/Documents/Backup/htw_local/data/sim_data/sim_data_2M_12_21.rds"))


sd$sim_dataAll <- map(sd$sim_dataAll, setDT)
names(sd)
# "sim_dataAll"   "prior_samples" "args_list"     "time"   
map(sd$sim_dataAll, class)
# $Exam_Varied
# [1] "tbl_df"     "tbl"        "data.frame"
# $Exam_Constant
# [1] "tbl_df"     "tbl"        "data.frame"
# $ALM_Varied
# [1] "tbl_df"     "tbl"        "data.frame"
# $ALM_Constant
# [1] "tbl_df"     "tbl"        "data.frame"
# $Alt_Varied
# [1] "tbl_df"     "tbl"        "data.frame"
# $Alt_Constant
# [1] "tbl_df"     "tbl"        "data.frame"


# discrepacy function - computes distance between simulated and observed data
dist_rmse <- function(simulated, observed) {
  return(sqrt(mean((simulated - observed)^2))) #RMSE
}

# Core ABC function
## computes distance between simulated and observed data, using discrepancy function set by dist_fun
## filters to obtain the posterior distribution of the model parameters (c and lr)
### if pct_keep is provided, posterior is the top pct_keep of the parameter sets with smallest discrepancy
### if pct_keep is NULL, posterior is all parameter sets with discrepancy less than the tolerance (tol)
## distance computation, and filtering performed separately for each of the 3 Fit methods (te, tr, teter)
### te is test only; tr is train only, teter is test and train (weighted evenly)
run_abc_fits <- function(data,sim_data, prior_samples, dist_fun="dist_rmse",Model,Group,pct_keep) {
  tol = 500
  target_data_train_test <- data[expMode2 %in% c("Test", "Train"), ]$y
  target_data_test <- data[expMode2 == "Test", ]$y
  target_data_train <- data[expMode2 == "Train", ]$y
  
  fn_name=dist_fun
  dist_fun=get(dist_fun)
  
  te_distances <- purrr::map_dbl(sim_data[85:90, ], dist_fun, observed = target_data_test)
  tr_distances <- purrr::map_dbl(sim_data[1:84, ], dist_fun, observed = target_data_train)
  teter_distances <- 0.5 * te_distances + 0.5 * tr_distances
  
  teter_results <- tibble(distance = teter_distances, c = prior_samples$c, 
    lr = prior_samples$lr, sim_index= seq_along(teter_distances)) |> 
    arrange(distance) |> 
    filter(if (!is.null(pct_keep)) row_number() <= n() * pct_keep else distance <= tol) |>
    mutate(rank=row_number(),Model,Group,Fit_Method="Test & Train") |>
    rowwise() |>
    mutate(sim_dat = list(tibble(Model,Fit_Method,c,lr,distance,rank,data,setNames(sim_data[, .SD, .SDcols = sim_index], "pred")) |> 
    mutate(resid=y-pred)))
    
 # teter_pp <- sim_data[,.SD,.SDcols=teter_results$sim_index] |> bind_rows()
#  k = teter_results[1,]$sim_dat

  te_results <- tibble(distance = te_distances, c = prior_samples$c, 
  lr = prior_samples$lr,  sim_index= seq_along(te_distances)) |> 
    arrange(distance) |> 
    filter(if (!is.null(pct_keep)) row_number() <= n() * pct_keep else distance <= tol) |>
    mutate(rank=row_number(),Model,Group,Fit_Method="Test Only") |>
    rowwise() |>
    mutate(sim_dat = list(tibble(Model,Fit_Method,c,lr,distance,rank,data,setNames(sim_data[, .SD, .SDcols = sim_index], "pred")) |> 
    mutate(resid=y-pred)))
  
  tr_results <- tibble(distance = tr_distances, c = prior_samples$c, 
    lr = prior_samples$lr,  sim_index= seq_along(tr_distances)) |>
    arrange(distance) |> 
    filter(if (!is.null(pct_keep)) row_number() <= n() * pct_keep else distance <= tol) |>
    mutate(rank=row_number(),Model,Group,Fit_Method="Train Only") |>
    rowwise() |>
    mutate(sim_dat = list(tibble(Model,Fit_Method,c,lr,distance,rank,data,setNames(sim_data[, .SD, .SDcols = sim_index], "pred")) |> 
    mutate(resid=y-pred)))
  
  tibble::lst(teter_results = teter_results, te_results = te_results, tr_results = tr_results, Model, Group, pct_keep, fn_name)
}


pct_keep=.001
abc_ev <- run_abc_fits(avg_dsv, sim_data= sd$sim_dataAll$Exam_Varied ,sd$prior_samples, dist_fun="dist_rmse",Model="EXAM", Group="Varied",pct_keep)
abc_almv <- run_abc_fits(avg_dsv, sim_data=sd$sim_dataAll$ALM_Varied,sd$prior_samples, dist_fun="dist_rmse", Model="ALM", Group="Varied", pct_keep)
abc_altv <- run_abc_fits(avg_dsv, sim_data=sd$sim_dataAll$Alt_Varied,sd$prior_samples, dist_fun="dist_rmse",Model="Alt_EXAM", Group="Varied",pct_keep)

abc_ec <- run_abc_fits(avg_dsc, sim_data=sd$sim_dataAll$Exam_Constant,sd$prior_samples, dist_fun="dist_rmse",Model="EXAM", Group="Constant",pct_keep)
abc_almc <- run_abc_fits(avg_dsc, sim_data=sd$sim_dataAll$ALM_Constant,sd$prior_samples, dist_fun="dist_rmse",Model="ALM", Group="Constant",pct_keep)
abc_altc <- run_abc_fits(avg_dsc, sim_data=sd$sim_dataAll$Alt_Constant,sd$prior_samples, dist_fun="dist_rmse",Model="Alt_EXAM", Group="Constant",pct_keep)

saveRDS(tibble::lst(abc_ev,abc_almv,abc_altv,abc_ec,abc_almc,abc_altc),
  here::here(paste0("data/abc_1p5M_rmse_",sub("^0", "", gsub("\\.", "p", pct_keep)),".rds")))



# load group posteriors
group_prior=abc_2M_p001 = readRDS(here::here("data/abc_2M_rmse_p001.rds"))
group_prior=abc_2M_p0001 = readRDS(here::here("data/abc_2M_rmse_p0001.rds"))
names(group_prior)
# "abc_ev"   "abc_almv" "abc_altv" "abc_ec"   "abc_almc" "abc_altc"

# preview structure
str(group_prior$abc_ev[[1]] |> select(-sim_dat))
# rowws_df [2,000 × 8] (S3: rowwise_df/tbl_df/tbl/data.frame)
#  $ distance  : Named num [1:2000] 222 222 223 225 229 ...
#   ..- attr(*, "names")= chr [1:2000] "...1441560" "...184544" "...407068" "...428330" ...
#  $ c         : num [1:2000] 2.26e-05 2.53e-05 1.35e-04 2.27e-05 2.11e-04 ...
#  $ lr        : num [1:2000] 0.108 0.143 0.129 0.579 0.6 ...
#  $ sim_index : int [1:2000] 1441560 184544 407068 428330 597281 433816 112705 774776 1556947 544842 ...
#  $ rank      : int [1:2000] 1 2 3 4 5 6 7 8 9 10 ...
#  $ Model     : chr [1:2000] "EXAM" "EXAM" "EXAM" "EXAM" ...
#  $ Group     : chr [1:2000] "Varied" "Varied" "Varied" "Varied" ...
#  $ Fit_Method: chr [1:2000] "Test & Train" "Test & Train" "Test & Train" "Test & Train" ...
#  - attr(*, "groups")= tibble [2,000 × 1] (S3: tbl_df/tbl/data.frame)
#   ..$ .rows: list<int> [1:2000] 
#   .. ..$ : int 1

# create separate data frames for each fit method
teter <- group_prior |> map_dfr(~tibble(pluck(.x$teter_results))) 
te <- group_prior |> map_dfr(~tibble(pluck(.x$te_results)))
tr <- group_prior |> map_dfr(~tibble(pluck(.x$tr_results)))

head(teter,2)
# A tibble: 2 × 9
#   distance         c    lr sim_index  rank Model Group  Fit_Method   sim_dat           
#      <dbl>     <dbl> <dbl>     <int> <int> <chr> <chr>  <chr>        <list>            
# 1     222. 0.0000226 0.108   1441560     1 EXAM  Varied Test & Train <tibble [90 × 13]>
# 2     222. 0.0000253 0.143    184544     2 EXAM  Varied Test & Train <tibble [90 × 13]>



teter |> select(sim_dat) |> unnest(sim_dat)  %>%  filter(expMode2=="Test") |>
  mutate(facet_label = paste0("rank: ", rank, "\n", "c: ", round(c, 4), "\n", "lr: ", round(lr, 4), "\n", "distance: ", round(distance, 1))) |>
  ggplot(aes(x = x, y = pred, fill=condit)) + 
  stat_summary(fun=mean, geom="bar", position=position_dodge()) +
  stat_summary(fun.data=mean_se, geom="errorbar", position=position_dodge()) +
   facet_wrap(~Model)


tab_teter <- teter |> select(sim_dat) |> unnest(sim_dat)  %>%  filter(expMode2=="Test", rank<=3) |>
  group_by(Model, condit, x) |> 
  summarize(mean_pred = mean(pred),mean_y = mean(y), mean_distance = mean(distance), mean_c = mean(c), mean_lr = mean(lr)) |> arrange(condit, Model,x)
tab_teter 



te |> select(sim_dat) |> unnest(sim_dat)  %>%  filter(expMode2=="Test") |>
  group_by(Model, condit, x) |> 
  summarize(mean_pred = mean(pred),mean_y = mean(y), mean_distance = mean(distance), mean_c = mean(c), mean_lr = mean(lr)) 


#
#
#
