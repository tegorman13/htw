#
#
#
#
#
#
#
#
#
#
#
#
#
#
pacman::p_load(tidyverse,data.table,abc,future,furrr,here,patchwork, conflicted)
conflict_prefer_all("dplyr", quiet = TRUE)
walk(c("fun_alm","fun_model", "Display_Functions"), ~ source(here::here(paste0("Functions/", .x, ".R"))))
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



input_layer =  c(100,350,600,800,1000,1200)
output_layer = input_layer


#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#

generate_prior_c_lr <- function(n) {
  prior_samples <- tibble(
    c = runif(n, 0.000001, 7),
    lr = runif(n, 0.000001, 7),
  )
  return(prior_samples)
}




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
  
  #fd <- rbind(train_data,test_data)
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
  
  #fd <- rbind(train_data,test_data)
  fd = eval(parse(text=paste0("rbind(",return_dat,")")))
  if(mode=="sim"){return(fd$pred)
    }else {return(fd)}
  
}



abc_loss_plot <- function(abc_result,prior_samples)
{
  pgrid <- cbind(prior_samples,dist=abc_result$dist) |> mutate(dist=sqrt(dist))
  upper_limit <- quantile(pgrid$dist, 0.95)  # 95th percentile
  pgrid <- pgrid |> filter(dist <= upper_limit)
  gg_prior <- ggplot(pgrid, aes(x = c, y = lr, color = dist)) +
  geom_point() + 
  #scale_color_viridis(limits = c(NA, upper_limit)) +
  scale_color_gradient(low = "blue", high = "red",limits = c(NA, upper_limit)) +  
  geom_contour(aes(x = c, y = lr, z = dist), color = 'white') +
  labs(title = "Point Plot of pgrid", x = "c", y = "lr") +
  theme_minimal()+ ggtitle("Prior Loss")
  
  
  post_grid <- cbind(prior_samples[abc_result$region,],dist=abc_result$dist[abc_result$region],k=abc_result$unadj.values)
  gg_post <- ggplot(post_grid, aes(x = c, y = lr, color = dist)) +
  geom_point() + 
  scale_color_gradient(low = "blue", high = "red",limits = c(NA, upper_limit)) +  
  geom_contour(aes(x = c, y = lr, z = dist), color = 'white') +
  labs(title = "Point Plot of pgrid", x = "c", y = "lr") +
  theme_minimal() + ggtitle("Posterior Loss")
  gg_prior+gg_post
  
}



# Define the function
abc_plot <- function(abc_result, data) {
  
  pgrid <- cbind(prior_samples,dist=abc_result$dist) |> mutate(dist=sqrt(dist))

  # Extracting posterior samples
  posterior_samples <- abc_result$unadj.values
  colnames(posterior_samples) <- c("c", "lr")
  posterior_samples_long <- tidyr::pivot_longer(as.data.frame(posterior_samples), everything())

  # Generate posterior density plots
  postV <- ggplot(posterior_samples_long, aes(x=value)) +
    geom_density() +
    facet_wrap(~name, scales="free") +
    theme_minimal() +
    labs(x="Value", y="Density", title="Posterior Density Plots")
  

  # Compute summary statistics
  sum_stats <- data.frame(
    mean = apply(posterior_samples, 2, mean),
    median = apply(posterior_samples, 2, median),
    mode = apply(posterior_samples, 2, Mode) # Ensure the 'Mode' function is defined
  )
  sum_stats <- tibble::rownames_to_column(sum_stats, var="parameter")
  print(sum_stats)

  # Calculate mode for parameters c and lr
  density_c <- density(posterior_samples[, "c"])
  mode_c <- density_c$x[which.max(density_c$y)]
  density_lr <- density(posterior_samples[, "lr"])
  mode_lr <- density_lr$x[which.max(density_lr$y)]
  
  print(paste0("Mode of c: ", mode_c, ", Mode of lr: ", mode_lr))

  # Run simulation with full_sim_exam function
  fsv <- full_sim_exam(as.data.table(data), c=mode_c, lr=mode_lr, exam.response, input_layer, output_layer, mode="ret")

  # Generate plot for simulation results
  plot_fsv <- fsv |> 
    tidyr::pivot_longer(y:pred, names_to="Resp", values_to="vx") |> 
    ggplot(aes(x, vx, fill=Resp, group=Resp)) +
    stat_summary(geom="bar", fun="mean", position=position_dodge()) +
    scale_fill_manual(values=col_themes$wes2) +
    scale_x_continuous(breaks=sort(unique(data$x)), labels=sort(unique(fsv$x))) +
    theme(legend.title = element_blank(), legend.position="top") +
    ggtitle("Fit to Test Only")
 
  postV + plot_fsv
}




n_prior_samples <- 5000
prior_samples <- generate_prior_c_lr(n_prior_samples)


# 
# data=avg_dsc
# return_dat="train_data,test_data"
# return_dat="train_data"
# 
# eval(parse(text=paste0("rbind(",return_dat,")")))
# 
# paste0(rbind(parse(text=return_dat)))
# 

#
#
#
#
#
#
#
#
#

# General function for ABC fits with model flexibility
run_abc_fits <- function(data, input_layer, output_layer, simulation_function, n_prior_samples = 1000, tol = 0.01, return_dat = "train_data,test_data") {
  
  input_layer =  c(100,350,600,800,1000,1200)
  output_layer = input_layer

  # Generate prior samples
  prior_samples <- generate_prior_c_lr(n_prior_samples)
  
    # Convert the string to a function
  if (is.character(simulation_function)) {
    simulation_function <- get(simulation_function, mode = "function")
  } else {
    simulation_function <- match.fun(simulation_function)
  }

  plan(multisession)
  # Run simulations for each prior sample
  simulation_results <- future_map_dfc(seq_len(nrow(prior_samples)), function(idx) {
    params <- prior_samples[idx, ]
    simulation_function(data=as.data.table(data), c=params$c, lr=params$lr, input_layer=input_layer, output_layer=output_layer, return_dat = return_dat)
  }, .options = furrr_options(seed = TRUE))

  # Extract target data
  target_data_train_test <- data[expMode2 %in% c("Test", "Train"), ]$y
  target_data_test <- data[expMode2 == "Test", ]$y
  target_data_train <- data[expMode2 == "Train", ]$y

  # ABC for train and test data
  abc_train_test <- abc(
    target = target_data_train_test,
    param = prior_samples,
    sumstat = do.call(rbind, simulation_results),
    tol = tol,
    method = "rejection"
  )

  # ABC for test data only
  abc_test <- abc(
    target = target_data_test,
    param = prior_samples,
    sumstat = do.call(rbind, simulation_results[85:90, ]),
    tol = tol,
    method = "rejection"
  )
  
    # ABC for test data only
  abc_train <- abc(
    target = target_data_train,
    param = prior_samples,
    sumstat = do.call(rbind, simulation_results[1:84, ]),
    tol = tol,
    method = "rejection"
  )

  # Return results
  tibble::lst(abc_train_test = abc_train_test, abc_test = abc_test, abc_train,data=data,n_prior_samples,prior_samples,tol, simulation_function, targets=tibble::lst(target_data_train_test, target_data_test))
}

# Example usage for different models
n_priors=200
abc_results_varied_exam <- run_abc_fits(avg_dsv, input_layer, output_layer, full_sim_exam, n_prior_samples = n_priors)
abc_results_varied_alt_exam <- run_abc_fits(avg_dsv, input_layer, output_layer, full_sim_alt_exam, n_prior_samples = n_priors)
abc_results_varied_alm <- run_abc_fits(avg_dsv, input_layer, output_layer, full_sim_alm, n_prior_samples = n_priors)

abc_results_constant_alm <- run_abc_fits(avg_dsc, input_layer, output_layer, full_sim_alm, n_prior_samples = n_priors)
abc_results_constant_exam <- run_abc_fits(avg_dsc, input_layer, output_layer, full_sim_exam, n_prior_samples = n_priors)
abc_results_constant_alt_exam <- run_abc_fits(avg_dsc, input_layer, output_layer, full_sim_alt_exam, n_prior_samples = n_priors)



# store all models in list and save to data/model_cache
model_cache <- list(
  abc_results_varied_exam = abc_results_varied_exam,
  abc_results_varied_alt_exam = abc_results_varied_alt_exam,
  abc_results_varied_alm = abc_results_varied_alm,
  abc_results_constant_alm = abc_results_constant_alm,
  abc_results_constant_exam = abc_results_constant_exam,
  abc_results_constant_alt_exam = abc_results_constant_alt_exam
)

saveRDS(model_cache,here::here(paste0("data/model_cache/abc_group_200k_",format(Sys.time(), "%H_%M_%OS"),".rds")))



#
#
#
#
#
#
#

n_priors=20000
args_list <- tibble::lst(
  abc_v_exam=list(data = avg_dsv, input_layer = input_layer, output_layer = output_layer, simulation_function = full_sim_exam, n_prior_samples = n_priors),
  abc_v_alt_exam=list(data = avg_dsv, input_layer = input_layer, output_layer = output_layer, simulation_function = full_sim_alt_exam, n_prior_samples = n_priors),
  abc_v_alm=list(data = avg_dsv, input_layer = input_layer, output_layer = output_layer, simulation_function = full_sim_alm, n_prior_samples = n_priors),
  abc_c_alm=list(data = avg_dsc, input_layer = input_layer, output_layer = output_layer, simulation_function = full_sim_alm, n_prior_samples = n_priors),
  abc_c_exam=list(data = avg_dsc, input_layer = input_layer, output_layer = output_layer, simulation_function = full_sim_exam, n_prior_samples = n_priors),
  abc_c_alt_exam=list(data = avg_dsc, input_layer = input_layer, output_layer = output_layer, simulation_function = full_sim_alt_exam, n_prior_samples = n_priors)
)


results <- future_map(args_list, ~do.call(run_abc_fits, .x))
saveRDS(results,here::here(paste0("data/model_cache/abc_group12_",format(Sys.time(), "%H_%M_%OS"),".rds")))



#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#

k=full_sim_exam(avg_dsv, .5, .2,exam.response,input_layer, output_layer,return_dat = "train_data,test_data")


n_prior_samples <- 5000
prior_samples <- generate_prior_c_lr(n_prior_samples)


plan(multisession)
sd_v_tetr <- future_map_dfc(seq_len(nrow(prior_samples)), function(idx) {
  params <- prior_samples[idx, ]
  full_sim_exam(as.data.table(avg_dsv), params$c, params$lr,exam.response, input_layer, output_layer,return_dat = "train_data,test_data")
},.options = furrr_options(seed = T))


dst_v_tetr <- avg_dsv[expMode2=="Test" | expMode2=="Train",]$y
abc_v_tetr <- abc(
  target = dst_v_tetr,
  param = prior_samples,
  sumstat =do.call(rbind, sd_v_tetr),
  tol = .01,
  method = "rejection",
  names=colnames(dst_v_tetr)
)

dst_v_te <- avg_dsv[expMode2=="Test",]$y
sd_v_te <- sd_v_tetr[85:90,]
abc_v_te <- abc(
  target = dst_v_te,
  param = prior_samples,
  sumstat =do.call(rbind, sd_v_te),
  tol = .01,
  method = "rejection",
  names=colnames(dst_v_te)
)




abc_plot(abc_v_tetr, avg_dsv) / { abc_loss_plot(abc_v_tetr,prior_samples) }
abc_plot(abc_v_te, avg_dsv) / { abc_loss_plot(abc_v_te,prior_samples) }



#
#
#
#
#
#

n_prior_samples <- 5000
prior_samples <- generate_prior_c_lr(n_prior_samples)


plan(multisession)

sd_c_tetr <- future_map_dfc(seq_len(nrow(prior_samples)), function(idx) {
  params <- prior_samples[idx, ]
  full_sim_exam(as.data.table(avg_dsc), params$c, params$lr,exam.response, input_layer, output_layer,return_dat = "train_data,test_data")
},.options = furrr_options(seed = T))


dst_c_tetr <- avg_dsc[expMode2=="Test" | expMode2=="Train",]$y

abc_c_tetr<- abc(
  target = dst_c_tetr,
  param = prior_samples,
  sumstat =do.call(rbind, sd_c_tetr),
  tol = .01,
  method = "rejection",
  names=colnames(dst_c_tetr)
)

abc_plot(abc_c_tetr, avg_dsc)



dst_c_te <- avg_dsc[expMode2=="Test",]$y
sd_c_te <- sd_c_tetr[85:90,]

abc_c_te <- abc(
  target = dst_c_te,
  param = prior_samples,
  sumstat =do.call(rbind, sd_c_te),
  tol = .01,
  method = "rejection",
  names=colnames(dst_c_te)
)


abc_plot(abc_c_tetr, avg_dsc) / { abc_loss_plot(abc_c_tetr,prior_samples) }
abc_plot(abc_c_te, avg_dsc) / { abc_loss_plot(abc_c_te,prior_samples) }


#
#
#
#
#
#


#
#
#
#
#
#
#
#
#
#

lin <- abc(target=dst_obs, param=prior_samples, sumstat=do.call(rbind, simulated_data), tol=.08, hcorr =FALSE, method = "loclinear")


posterior_samples <- lin$unadj.values
colnames(posterior_samples) <- c("c", "lr")
head(posterior_samples)
posterior_samples_long <- tidyr::pivot_longer(as.data.frame(posterior_samples), everything())

postV <- ggplot(posterior_samples_long, aes(x=value)) +
  geom_density() + facet_wrap(~name, scales="free") +
  labs(x="Value", y="Density", title="Posterior Density Plots")
postV


sum_stats <- data.frame(mean = apply(posterior_samples, 2, mean),
  median = apply(posterior_samples, 2, median), mode=apply(posterior_samples, 2, Mode))
( sum_stats=rownames_to_column(sum_stats,var="parameter") )



summary(lin)

plot(lin,prior_samples,subsample=100)
hist(lin)




cv.lin <- cv4abc(param=prior_samples, sumstat=do.call(rbind, simulated_data),nval=15,tols=c(.01,.05,.1,.2,.3,.5), method="rejection", statistic="mode")
summary(cv.lin)

cv.lin2 <- cv4abc(param=prior_samples, sumstat=do.call(rbind, simulated_data),nval=15,tols=c(.01,.05,.1,.2,.3,.5),abc.out=abc_result, method="rejection", statistic="mode")

summary(cv.lin2)





nn <- abc(target=dst_obs, param=prior_samples, sumstat=do.call(rbind, simulated_data), tol=.01, hcorr =FALSE, method = "neuralnet")


#
#
#
#
#
#


postV <- ggplot(posterior_samples_long, aes(x=value)) +
  geom_density() +
  facet_wrap(~name, scales="free") +
  theme_minimal() +
  labs(x="Value", y="Density", title="Posterior Density Plots")
postV


sum_stats <- data.frame(mean = apply(posterior_samples, 2, mean),
  median = apply(posterior_samples, 2, median), mode=apply(posterior_samples,2,Mode))
( sum_stats=rownames_to_column(sum_stats,var="parameter") )

#
#
#
#
#
#

id1 <- ds |> filter(id==3)
n_prior_samples <- 500
prior_samples <- generate_prior_c_lr(n_prior_samples)

plan(multisession)
simulated_data <- future_map_dfc(seq_len(nrow(prior_samples)), function(idx) {
  params <- prior_samples[idx, ]
  full_sim_alt_exam(as.data.table(id1), params$c, params$lr,alt_exam, input_layer, output_layer,return_dat = "train_data,test_data")
},.options = furrr_options(seed = T))




dst_obs <- id1[expMode2=="Test",]$y
dst_obs <- id1[expMode2=="Test" | expMode2=="Train",]$y

tolerance <- 0.1 * sd(dst_obs)

abc_result <- abc(
  target = dst_obs,
  param = prior_samples,
  sumstat =do.call(rbind, simulated_data),
  tol = .1,
  method = "rejection",
  names=colnames(dst_obs)
)



posterior_samples <- abc_result$unadj.values
colnames(posterior_samples) <- c("c", "lr")

posterior_samples_long <- tidyr::pivot_longer(as.data.frame(posterior_samples), everything())


postV <- ggplot(posterior_samples_long, aes(x=value)) +
  geom_density() +
  facet_wrap(~name, scales="free") +
  theme_minimal() +
  labs(x="Value", y="Density", title="Posterior Density Plots")
postV


sum_stats <- data.frame(mean = apply(posterior_samples, 2, mean),
  median = apply(posterior_samples, 2, median), mode=apply(posterior_samples,2,Mode))
( sum_stats=rownames_to_column(sum_stats,var="parameter") )


density_c <- density(posterior_samples[, "c"])
mode_c <- density_c$x[which.max(density_c$y)]
density_lr <- density(posterior_samples[, "lr"])
mode_lr <- density_lr$x[which.max(density_lr$y)]




fsv <- full_sim_alt_exam(as.data.table(id1), c=mode_c, lr=mode_lr,exam.response, input_layer, output_layer,mode="ret")

fsv |> pivot_longer(y:pred,names_to="Resp",values_to="vx") |> 
  ggplot(aes(x,vx,fill=Resp, group=Resp)) +stat_summary(geom="bar",fun="mean",position=position_dodge())+
  scale_fill_manual(values=col_themes$wes2)+
  scale_x_continuous(breaks=sort(unique(ds$x)), labels=sort(unique(fsv$x))) +
  theme(legend.title = element_blank(), legend.position="top") +ggtitle("Fit to Test Only")





#
#
#
#
#
