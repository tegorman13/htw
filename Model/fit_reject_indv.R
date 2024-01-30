pacman::p_load(dplyr,purrr,tidyr, data.table, here, conflicted, future, furrr)
conflict_prefer_all("dplyr", quiet = TRUE)
walk(c("fun_alm","fun_model"), ~ source(here::here(paste0("Functions/", .x, ".R"))))
set.seed(123)
ds <- readRDS(here::here("data/e1_md_11-06-23.rds"))  |> as.data.table()
watch_ids <<- c(1,33,66,20,76)

samp_priors <- function(n,cMean=-5,cSig=1.5,lrSig=1) {
  prior_samples <- tibble(
    c = rlnorm(n,cMean,sdlog=cSig),
    lr = extraDistr::rhnorm(n,sigma=lrSig),
    #lr = fdrtool::rhalfnorm(n,lrSig)
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
  start_tol = round(tol,3); 
  abc <- list()
  try_count=0;
   inc_count=0;
   cur_tol_success=0;
   success_rate=0;
  t1=system.time({
  for(j in 1:num_iterations) {
    #print(i)
    #print(paste0("sbj: ",data$id[1] ,"Particle:", j))

    found=0;
   
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
      try_count=try_count+1;
      cur_tol_success = cur_tol_success+1;
      success_rate = cur_tol_success/try_count;
      }
      if (try_count > n_try && success_rate < min_accept_rate){
        if (data$id[1] %in% watch_ids){
        message((paste0("increase tol for subject", data$id[1]," current iteration: ",j," cur accept rate: ",round(success_rate,4))))
        }
        tol=tol*tolInc
        inc_count=inc_count+1;
        try_count=0;
        cur_tol_success=0;
      }
    }
  }  
  })

  abc <- rbindlist(abc$dist_sd) |> arrange(mean_error) |> relocate(id,condit)
  best <- abc |> head(1) |> round_tibble(7)
  message((paste0("\n", data$id[1], " completed in: ", t1[3])))
  message((paste0("inc count: ",inc_count," Start tol: ",start_tol ," end tol: ", round(tol,2)," success rate: ",round(success_rate,4))))
  message((paste0("Best c: ", best$c, " Best lr: ", best$lr, " Best error: ", best$mean_error,"\n")))
  return(abc)
}

####################################
####################################

# ids1 <- 1
#ids1 <- c(1,33,66)
#ids1 <- as.numeric(levels(ds$id))[1:8]
ids1 <- as.numeric(levels(ds$id))
#ids1 <- c(49)

cMean <<- -5.5; cSig <<- 3.5; lrSig <<- 2.5
prior_samples <- samp_priors(n=150000, cMean=cMean, cSig=cSig, lrSig=lrSig) 


subjects_data <-  ds |> filter(id %in% ids1)  %>% with(split(.,f =c(id), drop=TRUE))




args <- commandArgs(trailingOnly = TRUE)
# print(args)
# print(args[1])

# if args[1] is defined, set to num_iterations, otherwise 1000
num_iterations = ifelse(length(args) > 0, as.numeric(args[1]), 100)
n_try = ifelse(length(args) > 1, as.numeric(args[2]), 50)
tolM <<- ifelse(length(args) > 2, as.numeric(args[3]), .8)
tolInc <<- ifelse(length(args) > 3, as.numeric(args[4]), 1.1)
min_accept_rate <<- ifelse(length(args) > 4, as.numeric(args[5]), .01)


print(num_iterations)
print(n_try)
p_abc <- function(){message(paste0("cMean: ",cMean," cSig: ",cSig," lrSig: ",lrSig,"\n", "tolM: ",tolM, " tolInc: ",tolInc," accept_rate: ",min_accept_rate))}
p_abc()

save_folder <- paste0("n_iter_",num_iterations,"_ntry_",n_try,"_",format(Sys.time(),"%M%OS"))
dir.create(paste0("data/abc_reject/",save_folder))


(nc <- future::availableCores())
future::plan(multisession, workers = nc)

### EXAM Test 

cat("\nEXAM Test\n")
t1=system.time({
exam_test <- future_map(subjects_data, ~reject_abc(simulation_function = full_sim_exam, 
                                                    prior_samples = prior_samples, 
                                                    data = .x, 
                                                    num_iterations = num_iterations, 
                                                    n_try = n_try,
                                                    return_dat="test_data"), .options = furrr_options(seed = TRUE),
                                                    future.stdout = NA) %>% 
                                                    setNames(ids1)

                                                    })

#map(exam_test, ~{.x |> arrange(mean_error) |> head(5)})


ri=ri_reject_indv() %>% append(.,t1[3])
run_save <- tibble::lst(exam_test,Model="EXAM",Fit_Method="Test",prior_samples,cMean,cSig,lrSig,tolM,ri)
file_name <- paste0("data/abc_reject/",save_folder,"/","reject_",run_save$Model,"_",run_save$Fit_Method,"_",
                    num_iterations,"_",n_try,"_",format(Sys.time(),"%M%OS"), ".rds")
saveRDS(run_save,file=here::here(file_name))




#### ALM Test
cat("\nALM Test\n")
p_abc()
rm(run_save)

t1=system.time({
  alm_test <- future_map(subjects_data, ~reject_abc(simulation_function = full_sim_alm, 
                                                  prior_samples = prior_samples, 
                                                  data = .x, 
                                                  num_iterations = num_iterations, 
                                                  n_try = n_try,
                                                  return_dat="test_data"), 
                                                  .options = furrr_options(seed = TRUE),
                                                  future.stdout = NA) %>% 
    setNames(ids1)
  
})
print(t1[3])

ri=ri_pda_indv() %>% append(.,t1[3])
run_save <- tibble::lst(alm_test,Model="ALM",Fit_Method="Test",prior_samples,cMean,cSig,lrSig,tolM,ri)
file_name <- paste0("data/abc_reject/",save_folder,"/","reject_",run_save$Model,"_",run_save$Fit_Method,"_",
                    num_iterations,"_",n_try,"_",format(Sys.time(),"%M%OS"), ".rds")
saveRDS(run_save,file=here::here(file_name))





#### EXAM Test & Train
cat("\nEXAM Test Train\n")
p_abc()
rm(run_save)

t1=system.time({
  exam_test_train <- future_map(subjects_data, ~reject_abc(simulation_function = full_sim_exam, 
                                                  prior_samples = prior_samples, 
                                                  data = .x, 
                                                  num_iterations = num_iterations, 
                                                  n_try = n_try,
                                                  return_dat="train_data, test_data"), 
                                                  .options = furrr_options(seed = TRUE),
                                                  future.stdout = NA) %>% 
    setNames(ids1)
  
})
print(t1[3])

ri=ri_pda_indv() %>% append(.,t1[3])
run_save <- tibble::lst(exam_test_train,Model="EXAM",Fit_Method="Test_Train",prior_samples,cMean,cSig,lrSig,tolM,ri)
file_name <- paste0("data/abc_reject/",save_folder,"/","reject_",run_save$Model,"_",run_save$Fit_Method,"_",
                    num_iterations,"_",n_try,"_",format(Sys.time(),"%M%OS"), ".rds")
saveRDS(run_save,file=here::here(file_name))




#### ALM Test & Train
cat("\nALM Test Train\n")
p_abc()
rm(run_save)

t1=system.time({
  alm_test_train <- future_map(subjects_data, ~reject_abc(simulation_function = full_sim_alm, 
                                                  prior_samples = prior_samples, 
                                                  data = .x, 
                                                  num_iterations = num_iterations, 
                                                  n_try = n_try,
                                                  return_dat="train_data, test_data"), 
                                                  .options = furrr_options(seed = TRUE),
                                                  future.stdout = NA) %>% 
    setNames(ids1)
  
})
print(t1[3])

ri=ri_pda_indv() %>% append(.,t1[3])
run_save <- tibble::lst(alm_test_train,Model="ALM",Fit_Method="Test_Train",prior_samples,cMean,cSig,lrSig,tolM,ri)
file_name <- paste0("data/abc_reject/",save_folder,"/","reject_",run_save$Model,"_",run_save$Fit_Method,"_",
                    num_iterations,"_",n_try,"_",format(Sys.time(),"%M%OS"), ".rds")
saveRDS(run_save,file=here::here(file_name))




#### EXAM Train
cat("\nExam Train\n")
p_abc()
rm(run_save)

t1=system.time({
  exam_train <- future_map(subjects_data, ~reject_abc(simulation_function = full_sim_exam, 
                                                  prior_samples = prior_samples, 
                                                  data = .x, 
                                                  num_iterations = num_iterations, 
                                                  n_try = n_try,
                                                  return_dat="train_data"), 
                                                  .options = furrr_options(seed = TRUE),
                                                  future.stdout = NA) %>% 
    setNames(ids1)
  
})
print(t1[3])
print(paste0("cMean: ",cMean," cSig: ",cSig," lrSig: ",lrSig," tolM: ",tolM))
print(paste0("cMean: ",cMean," cSig: ",cSig," lrSig: ",lrSig," tolM: ",tolM, "tolInc: ",tolInc))


ri=ri_pda_indv() %>% append(.,t1[3])
run_save <- tibble::lst(exam_train,Model="EXAM",Fit_Method="Train",prior_samples,cMean,cSig,lrSig,tolM,ri)
file_name <- paste0("data/abc_reject/",save_folder,"/","reject_",run_save$Model,"_",run_save$Fit_Method,"_",
                    num_iterations,"_",n_try,"_",format(Sys.time(),"%M%OS"), ".rds")
saveRDS(run_save,file=here::here(file_name))





#### ALM Train
cat("\nALM Train\n")
p_abc()
rm(run_save)

t1=system.time({
  alm_train <- future_map(subjects_data, ~reject_abc(simulation_function = full_sim_alm, 
                                                       prior_samples = prior_samples, 
                                                       data = .x, 
                                                       num_iterations = num_iterations, 
                                                       n_try = n_try,
                                                       return_dat="train_data"), 
                                                       .options = furrr_options(seed = TRUE),
                                                       future.stdout = NA) %>% 
    setNames(ids1)
  
})
print(t1[3])

ri=ri_pda_indv() %>% append(.,t1[3])
run_save <- tibble::lst(alm_train,Model="ALM",Fit_Method="Train",prior_samples,cMean,cSig,lrSig,tolM,ri)
file_name <- paste0("data/abc_reject/",save_folder,"/","reject_",run_save$Model,"_",run_save$Fit_Method,"_",
                    num_iterations,"_",n_try,"_",format(Sys.time(),"%M%OS"), ".rds")
saveRDS(run_save,file=here::here(file_name))

cat("\nend")
p_abc()

cat("\n EXAM Test: \n")
print(knitr::kable(exam_test[[1]] |> head(3),format="markdown"))

cat("\n ALM Test: \n")
print(knitr::kable(alm_test[[1]] |> head(3),format="markdown"))
