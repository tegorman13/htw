pacman::p_load(dplyr,purrr,tidyr, data.table, here, conflicted, future, furrr)
conflict_prefer_all("dplyr", quiet = TRUE)
walk(c("fun_alm","fun_model","fun_indv_fit", "prep_mpar"), ~ source(here::here(paste0("Functions/", .x, ".R"))))

ds <- readRDS(here::here("data/e1_md_11-06-23.rds"))  |> as.data.table()
watch_ids <<- c(1)

seed <- round(runif(1,min=1,max=10),0)
set.seed(seed)
print(paste("seed:",seed))

####################################

samp_priors <- function(n=1) {
  prior_samples <- tibble(
    c = rlnorm(n,cMean,sdlog=cSig),
    lr = extraDistr::rhnorm(n,sigma=lrSig),
  )
  return(prior_samples)
}

reject_abc <- function(simulation_function, data, num_iterations = 5000, n_try=500, return_dat="test_data") {
  input_layer =  c(100,350,600,800,1000,1200)
  output_layer = input_layer
  data <- data |> as.data.table()
  target_data <- case_when(
    return_dat == "test_data" ~ list(target_data_test <- data[data$expMode2 == "Test", ]), 
    return_dat == "train_data" ~ list( target_data_train <- data[data$expMode2 == "Train", ]),
    return_dat == "train_data, test_data" ~ list( target_data_train_test <- data[expMode2 %in% c("Test", "Train"), ])
    ) |> pluck(1)

  tol <- target_data[, .(m = mean(y), sd = sd(y)), by = .(x)][, .(tol = mean(sd) * tolM)]

  start_tol = round(tol,3); 
  abc <- list(); abc$dist_sd <- vector("list", num_iterations)
  try_count=0;
  inc_count=0;
  cur_tol_success=0;
  success_rate=0;
  closest_mean_errors <- numeric(0)

  t1=system.time({
  for(j in 1:num_iterations) {

    found=0;

    while(found==0) {
    try_count=try_count+1;
    #current_theta <- prior_samples[sample(1:nrow(prior_samples), 1), ]
    current_theta <- samp_priors()
    sim_data <- simulation_function(data, current_theta$c, current_theta$lr, input_layer = input_layer, output_layer = output_layer, return_dat = return_dat) 

    # dist_sd <- target_data |> mutate(pred=sim_data,error=abs(y-pred)) |> group_by(id,condit,x) |> 
    #   summarise(mean_error=mean(error),.groups="keep") |> 
    #   group_by(id,condit) |> 
    #   summarise(mean_error=mean(mean_error),.groups="keep") 

    target_data[, pred := sim_data]
    target_data[, error := abs(y - pred)]
    intermediate_result <- target_data[, .(mean_error = mean(error)), by = .(id, condit, x)]
    dist_sd <- intermediate_result[, .(mean_error = mean(mean_error)), by = .(id, condit)]

    closest_mean_errors <- sort(c(closest_mean_errors, dist_sd$mean_error))
    if (length(closest_mean_errors) > 5) {
      closest_mean_errors <- head(closest_mean_errors, 5)  
    }

    if(dist_sd$mean_error< tol) {
      abc$dist_sd[[j]] <- cbind(current_theta,dist_sd,tol,inc_count)
      found=1
      try_count=try_count+1;
      cur_tol_success = cur_tol_success+1;
      success_rate = cur_tol_success/try_count;
      } else if (try_count > n_try && success_rate < min_accept_rate){
       
        average_closest_error <- mean(closest_mean_errors)
        tol <- (tol*tolInc) + (abs(average_closest_error-tol)/2) #* tolInc  # Adjust tolerance
        #tol=tol*tolInc
        inc_count=inc_count+1;
        try_count=0;
        cur_tol_success=0;

       if (data$id[1] %in% watch_ids){
        message(paste0("increase tol(",round(tol,3),") for subject", data$id[1]," current iteration: ",j,"\n",
        "cur accept rate: ",round(success_rate,4),"\n",
         "avg closest error: ",round(average_closest_error,2), " new tol: ",round(tol,2),"\n"))
        }

      }
    }
  }  
  })

  abc <- rbindlist(abc$dist_sd) |> arrange(mean_error) |> relocate(id,condit)
  best <- abc |> head(1) |> round_tibble(7)
  message((paste0("\n", data$id[1], " completed in: ", round(t1[3],2))))
  message((paste0("inc count: ",inc_count," Start tol: ",start_tol ," end tol: ", round(tol,2)," success rate: ",round(success_rate,4))))
  message((paste0("Best c: ", best$c, " Best lr: ", round(best$lr,3), " Best error: ", round(best$mean_error,2),"\n\n\n")))
  return(abc)
}

####################################
####################################

run_abc_tests <- function(simulation_function, data_list, return_dat, ids) {
  p_abc()
  cat("\nRunning ABC Test: ",as.character(substitute(simulation_function))," ",return_dat, "\n", "Parallel Execution\n")

  t1 <- system.time({
    results <- future_map(data_list, 
                          ~reject_abc(simulation_function = simulation_function, 
                                      #prior_samples = prior_samples, 
                                      data = .x, 
                                      num_iterations = num_iterations, 
                                      n_try = n_try,
                                      return_dat = return_dat), 
                          .options = furrr_options(seed = TRUE),
                          future.stdout = NA) %>% 
              setNames(ids)
  })
  
  print(t1[3])
  return(results)
}

####################################
####################################

args <- commandArgs(trailingOnly = TRUE)
num_iterations = ifelse(length(args) > 0, as.numeric(args[1]), 500)
n_try = ifelse(length(args) > 1, as.numeric(args[2]), 3500)
tolM <<- ifelse(length(args) > 2, as.numeric(args[3]), .85)
tolInc <<- ifelse(length(args) > 3, as.numeric(args[4]), 1.01)
min_accept_rate <<- ifelse(length(args) > 4, as.numeric(args[5]), .0004)

cMean <<- -6; 
cSig <<- 3.5; 
lrSig <<- 4.0

# uniform dist between -4.5 and -6.0
# cMean <<- runif(1, min = -4.5, max = -6.0)
# cSig <<- runif(1, min = 1.5, max = 3.5)
# lrSig <<- runif(1, min = 1.0, max = 3.0)

prior_samples <- samp_priors(n=13000) 

p_abc <- function(){
  message(paste0("\n","cMean: ",cMean," cSig: ",cSig," lrSig: ",lrSig,"\n", "tolM: ",
  tolM, " tolInc: ",tolInc," accept_rate: ",min_accept_rate, "\n",
  "num_iterations: ",num_iterations," n_try: ",n_try,"\n",
  "median-c=",round(median(prior_samples$c),4)," median-lr=",round(median(prior_samples$lr),4),"\n"))
  }
 
####################################
# ids1 <- 1
# ids1 <- c(1,33,66)
# ids1 <- as.numeric(levels(ds$id))[1:8]
ids1 <- as.numeric(levels(ds$id))
subjects_data <-  ds |> filter(id %in% ids1)  %>% with(split(.,f =c(id), drop=TRUE))


save_folder <- paste0("n_iter_",num_iterations,"_ntry_",n_try,"_",format(Sys.time(),"%M%OS"))
dir.create(paste0("data/abc_reject/",save_folder))

parallel <<- 1
#parallel <<- runif(1) <.5

if (parallel==1) {
  print(nc <- future::availableCores())
  future::plan(cluster, workers = nc)
  message("Running in parallel\n")
} else if(parallel==2){
spec <- make_spec()
pc <- parallel::makeCluster(type='PSOCK', master=macpro, spec)
plan(cluster, workers = pc)
} else {
  message("Running in serial\n")
}


# parallel::stopCluster(pc)

run_function <- ifelse(parallel>0,run_abc_tests, run_abc_tests_serial)



t1<- system.time( exam_test <- run_function(full_sim_exam, subjects_data, "test_data", ids1) )
save_abc_test_results(exam_test, "EXAM", "Test", ri_reject_indv, subjects_data, ids1,save_folder, t1)
rm(exam_test); gc()

t1<- system.time( alm_test <- run_function(full_sim_alm, subjects_data, "test_data", ids1) )
save_abc_test_results(alm_test, "ALM", "Test", ri_reject_indv, subjects_data, ids1,save_folder, t1)
rm(alm_test); gc()

t1<- system.time( exam_test_train <- run_function(full_sim_exam, subjects_data, "train_data, test_data", ids1) )
save_abc_test_results(exam_test_train, "EXAM", "Test_Train", ri_reject_indv, subjects_data, ids1,save_folder, t1)
rm(exam_test_train); gc()

t1<- system.time( alm_test_train <- run_function(full_sim_alm, subjects_data, "train_data, test_data", ids1) )
save_abc_test_results(alm_test_train, "ALM", "Test_Train", ri_reject_indv, subjects_data, ids1,save_folder, t1)
rm(alm_test_train); gc()

t1<- system.time( exam_train <- run_function(full_sim_exam, subjects_data, "train_data", ids1) )
save_abc_test_results(exam_train, "EXAM", "Train", ri_reject_indv, subjects_data, ids1,save_folder, t1)
rm(exam_train); gc()

t1<- system.time( alm_train <- run_function(full_sim_alm, subjects_data, "train_data", ids1) )
save_abc_test_results(alm_train, "ALM", "Train", ri_reject_indv, subjects_data, ids1,save_folder, t1)



# cat("\nend")
# p_abc()

# cat("\n EXAM Test: \n")
# print(knitr::kable(exam_test[[1]] |> head(3),format="markdown"))

# cat("\n ALM Test: \n")
# print(knitr::kable(alm_test[[1]] |> head(3),format="markdown"))



# t1<- ( exam_test <- run_abc_tests_serial(full_sim_exam, subjects_data, "test_data", ids1) )
