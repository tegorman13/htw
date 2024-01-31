


save_abc_test_results <- function(fit_results, model, fit_method, ri_func, subjects_data, ids,save_folder,t1) {

  ri <- ri_func() %>% append(., t1[3])
  run_save <- tibble::lst(fit_results, Model = model, Fit_Method = fit_method, prior_samples, cMean, cSig, lrSig, tolM, min_accept_rate, ri)
  
  file_name <- paste0("data/abc_reject/", save_folder, "/", "reject_", run_save$Model, "_", run_save$Fit_Method, "_",
                      num_iterations, "_", n_try, "_", format(Sys.time(), "%M%OS"), ".rds")
  saveRDS(run_save, file = here::here(file_name))
}




ri_reject_indv <- function() {
  # Get all object names in the environment
  # all_objects <- ls(envir = globalenv())
  # function_names <- all_objects[sapply(all_objects, function(x) is.function(get(x)))]
  # functions_list <- setNames(lapply(function_names, function(x) get(x)), function_names)

  runInfo <- tibble::lst(
    #functions = functions_list,
    runScripts = list(
      indv_fit_script = readLines(here::here("Model/fit_reject_indv.R")),
     # model_script = readLines(here::here("Functions/fun_model.R")),
      alm_script = readLines(here::here("Functions/fun_alm.R"))
    ),
    path = getwd(),
    Computer = Sys.info()["nodename"],
    systemInfo = Sys.info(),
    sessionInfo = sessionInfo(),
    timeInit = Sys.time()
  )
  return(runInfo)
}



ri_pda_indv <- function() {
  runInfo <- tibble::lst(
    runScripts = list(
      indv_fit_script = readLines(here::here("Model/fit_pda_indv.R")),
     # model_script = readLines(here::here("Functions/fun_model.R")),
      alm_script = readLines(here::here("Functions/fun_alm.R"))
    ),
    path = getwd(),
    Computer = Sys.info()["nodename"],
    systemInfo = Sys.info(),
    sessionInfo = sessionInfo(),
    timeInit = Sys.time()
  )

  return(runInfo)
}


sample_from_kde <- function(kde_result, n_samples) {
  # Flatten the z matrix and sample indices based on these probabilities
  sampled_indices <- sample(length(kde_result$z), size = n_samples, replace = TRUE, prob = c(kde_result$z))
  
  # Convert indices to matrix row and column numbers
  rows <- ((sampled_indices - 1) %% nrow(kde_result$z)) + 1
  cols <- ((sampled_indices - 1) %/% nrow(kde_result$z)) + 1
  
  # Map back to x and y values and create a tibble
  samples <- tibble(c = kde_result$x[rows], lr = kde_result$y[cols])
  return(samples)
}

compute_kde <- function(c, lr, ngrid = 100, nsamples = 100, lim_buffer = 0.15) {
  iqr_c = IQR(c)
  iqr_lr = IQR(lr)
  
  # Calculate limits ensuring they are not below the specified minimum
  min_value = .00000001
  c_lims = c(max(min(c) - iqr_c * lim_buffer, min_value), max(c) + iqr_c * lim_buffer)
  lr_lims = c(max(min(lr) - iqr_lr * lim_buffer, min_value), max(lr) + iqr_lr * lim_buffer)
  
  kde = MASS::kde2d(c, lr, n = ngrid, lims = c(c_lims, lr_lims))
  kde_samples = sample_from_kde(kde, nsamples)
  
  return(list(kde = kde, kde_samples = kde_samples))
}


dist_rmse <- function(simulated, observed) {
  return(sqrt(mean((simulated - observed)^2))) #RMSE
}

dist_mse <- function(simulated, observed) {
  return(mean((simulated - observed)^2)) #MSE
}

rho=function(x,y) abs(sum(x)-sum(y))/length(x)			# rho function



# nll2 <- function(obsv,pred,sigma)
# {
#   nll= -sum(dnorm(obsv,mean=pred,sd=sigma,log=TRUE)) 
#   if (is.nan(nll)) {
#     nll <- 1e4 # Large penalty
#   }
#   return(nll)
# }


Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

weight_plot <- function(weight.mat){
  tibbleFromMat <-
    weight.mat %>%
    tibble::as_tibble() %>%
    tibble::rownames_to_column("Var1") %>%
    tidyr::pivot_longer(-Var1, names_to = "Var2", values_to = "value") %>%
    dplyr::mutate(
      Output_Layer = factor(Var1, levels = 1:dim(weight.mat)[1]),
      Input_Layer = factor(gsub("V", "", Var2), levels = 1:dim(weight.mat)[2])
    )
  ggplot(tibbleFromMat, aes(Output_Layer, Input_Layer)) +
    geom_tile(aes(fill = value)) +
    geom_text(aes(label = round(value, 3))) +
    scale_fill_gradient(low = "white", high = "red")
}


round_tibble <- function(tbl, rn) {
  tbl %>% 
    mutate(across(where(is.numeric), ~round(., rn)))
}

strip_list_notation <- function(str) {
  # Remove 'list(' at the beginning
  str <- gsub("^list\\(", "", str)
  # Remove ')' at the end
  str <- gsub("\\)$", "", str)
  return(str)
}


adjust_layer <- function(input.layer, k){
  if(k == 1) return(input.layer)
  new.layer <- c()
  for(i in 1:(length(input.layer) - 1)){
    new.layer <- c(new.layer, seq(input.layer[i], input.layer[i+1], length.out = k + 1))
  }
  new.layer <- unique(new.layer)
  return(new.layer)
}


### Loss Functions

#  mse <- mean((predictions - test_data$y)^2)

nll <- function(obsv,pred,sigma)
{
 -sum(dnorm(obsv,mean=pred,sd=sigma,log=TRUE))
}


RMSE <- function(x,y){
  sqrt(mean((x-y)^2, na.rm=TRUE))
}

# RMSE.blocked <- function(x,y,blocks=6){
#   #print("rmseBlocked")
#   data.table(x=x,y=y,t=seq(1,length(x))) %>% 
#     .[, `:=`(fitBins = cut(t, breaks = ..blocks, labels = c(1:..blocks)))] %>%
#     .[, .(predMean = mean(x), obsMean = mean(y)), keyby = .(fitBins)] %>%
#     .[, RMSE(predMean,obsMean)] %>% as.numeric()
# }

MAE <- function(x, y) {
  mean(abs(x - y))
}

# MAPE <- function(x, y) {
#   mean(abs((x - y) / y)) * 100
# }

# MedAE <- function(x, y) {
#   median(abs(x - y))
# }

# HuberLoss <- function(x, y, delta = 1) {
#   error <- x - y
#   abs_error <- abs(error)
#   loss <- ifelse(abs_error <= delta, 0.5 * error^2, delta * (abs_error - 0.5 * delta))
#   mean(loss)
# }


