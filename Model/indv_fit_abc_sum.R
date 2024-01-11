pacman::p_load(dplyr,purrr,tidyr,ggplot2, data.table, here, patchwork, conflicted)
conflict_prefer_all("dplyr", quiet = TRUE)
walk(c("fun_alm","fun_model"), ~ source(here::here(paste0("Functions/", .x, ".R"))))
set.seed(123)

ds <- readRDS(here::here("data/e1_md_11-06-23.rds"))  |> as.data.table()
dsv <- ds |> filter(condit=="Varied")  
dsc <- ds |> filter(condit=="Constant") 


group_posterior_all <- readRDS(here::here("data/abc_2M_rmse_p001.rds"))  


teter <- group_posterior_all |> map_dfr(~tibble(pluck(.x$teter_results))) 
te <- group_posterior_all |> map_dfr(~tibble(pluck(.x$te_results))) 
tr <- group_posterior_all |> map_dfr(~tibble(pluck(.x$tr_results))) 


# kernel density estimate of group posteriors
n_prior_samples=100; ngrid=100; buf = .25
kde_results <- purrr::map(group_posterior_all, ~
                            purrr::map_if(.x, is.list, ~compute_kde(.x$c, .x$lr, ngrid = ngrid, nsamples=n_prior_samples, lim_buffer=buf)))



fit_indv <- function(sbj_id, simulation_function, prior_samples,Model, Group) {

  data <- filter(ds,id==sbj_id) |> as.data.table()
  input_layer =  c(100,350,600,800,1000,1200)
  output_layer = input_layer
  return_dat="train_data, test_data"
  train_idx <- which(data$expMode2=="Train")
  test_idx <- which(data$expMode2=="Test")
  
  target_data_train_test <- data[expMode2 %in% c("Test", "Train"), ]$y
  target_data_test <- data[expMode2 == "Test", ]$y
  target_data_train <- data[expMode2 == "Train", ]$y
  
  sim_data_gen_s <- function(data, input_layer, output_layer, simulation_function, prior_samples, return_dat) {
    simulation_function <- match.fun(simulation_function)
    suppressMessages(simulation_results <- map_dfc(seq_len(nrow(prior_samples)), function(idx) {
      params <- prior_samples[idx, ]
      simulation_function(data=as.data.table(data), c=params$c, lr=params$lr, input_layer=input_layer, output_layer=output_layer, return_dat = return_dat)    
      }))
    }
  
    
    pct_keep=.8
    prior_samples_teter <- prior_samples$teter_results$kde_samples
    sim_data <- sim_data_gen_s(data, input_layer, output_layer,simulation_function, prior_samples_teter,return_dat) |> as.data.table()

    te_distances <- purrr::map_dbl(sim_data[test_idx, ], dist_rmse, observed = target_data_test)
    tr_distances <- purrr::map_dbl(sim_data[train_idx, ], dist_rmse, observed = target_data_train)
    teter_distances <- 0.5 * te_distances + 0.5 * tr_distances
    
    
    teter_results <- tibble(distance = teter_distances, c = prior_samples_teter$c, 
                            lr = prior_samples_teter$lr, sim_index= seq_along(teter_distances)) |> 
      arrange(distance) |> 
      filter(if (!is.null(pct_keep)) row_number() <= n() * pct_keep else distance <= tol) |>
      mutate(rank=row_number(),Model,Group,Fit_Method="Test & Train") |>
      rowwise() |>
      mutate(sim_dat = list(tibble(Model,Fit_Method,c,lr,distance,rank,data,pred=(sim_data[, .SD, .SDcols = sim_index])) |> 
                              mutate(resid=y-pred)))
    
    
    prior_samples_te <- prior_samples$te_results$kde_samples
    sim_data <- sim_data_gen_s(data, input_layer, output_layer,simulation_function, prior_samples_te,return_dat) |> as.data.table()
    
    
    te_distances <- purrr::map_dbl(sim_data[test_idx, ], dist_rmse, observed = target_data_test)

    te_results <- tibble(distance = te_distances, c = prior_samples_te$c, 
                         lr = prior_samples_te$lr,  sim_index= seq_along(te_distances)) |> 
      arrange(distance) |> 
      filter(if (!is.null(pct_keep)) row_number() <= n() * pct_keep else distance <= tol) |>
      mutate(rank=row_number(),Model,Group,Fit_Method="Test Only") |>
      rowwise() |>
      mutate(sim_dat = list(tibble(Model,Fit_Method,c,lr,distance,rank,data,pred=(sim_data[, .SD, .SDcols = sim_index])) |> 
                              mutate(resid=y-pred
                                     )))
    
    
    prior_samples_tr <- prior_samples$tr_results$kde_samples
    sim_data <- sim_data_gen_s(data, input_layer, output_layer,simulation_function,prior_samples_tr,return_dat) |> as.data.table()
    tr_distances <- purrr::map_dbl(sim_data[train_idx, ], dist_rmse, observed = target_data_train)
    
    tr_results <- tibble(distance = tr_distances, c = prior_samples_tr$c, 
                         lr = prior_samples_tr$lr,  sim_index= seq_along(tr_distances)) |>
      arrange(distance) |> 
      filter(if (!is.null(pct_keep)) row_number() <= n() * pct_keep else distance <= tol) |>
      mutate(rank=row_number(),Model,Group,Fit_Method="Train Only") |>
      rowwise() |>
      mutate(sim_dat = list(tibble(Model,Fit_Method,c,lr,distance,rank,data,pred=(sim_data[, .SD, .SDcols = sim_index])) |> 
                              mutate(resid=y-pred)))
    
    tibble::lst(teter_results = teter_results, te_results = te_results, tr_results = tr_results,id = as.numeric(sbj_id),
                         Model, Group, pct_keep, fn_name=dist_rmse) 
}


#id1 <- fit_indv(sbj_id=1,simulation_function=full_sim_exam, prior_samples=kde_results$abc_ev$teter_results$kde_samples)
#unique(dsv$id)

id_fits_exam_varied <- map(unique(dsv$id), ~ fit_indv(sbj_id=.x,
                                          simulation_function=full_sim_exam, 
                                          prior_samples=kde_results$abc_ev, 
                                          Model="EXAM", 
                                          Group="Varied"))
#names(id_fits_exam_varied) <- unique(dsv$id)
  
  
id_fits_alm_varied <- map(unique(dsv$id), ~ fit_indv(sbj_id=.x,
                                                simulation_function=full_sim_alm, 
                                                prior_samples=kde_results$abc_almv, 
                                                Model="ALM", 
                                                Group="Varied"))  



id_fits_exam_constant <- map(unique(dsc$id), ~ fit_indv(sbj_id=.x,
                                                      simulation_function=full_sim_exam, 
                                                      prior_samples=kde_results$abc_ec, 
                                                      Model="EXAM", 
                                                      Group="Constant"))



id_fits_alm_constant <- map(unique(dsc$id), ~ fit_indv(sbj_id=.x,
                                                     simulation_function=full_sim_alm, 
                                                     prior_samples=kde_results$abc_almc, 
                                                     Model="ALM", 
                                                     Group="Constant"))  
  



# save prior samples in sim instead of full kde results?

runInfo=tibble::lst(kde_results,
                        # data=ds,
                         functions=tibble::lst(fit_indv,dist_rmse),
                 runScripts=list(indv_fit_script=readLines(here::here("Model/indv_fit.R")),
                                 model_script=readLines(here::here("Functions/fun_model.R")),
                                 alm_script=readLines(here::here("Functions/fun_alm.R"))),
                 path=getwd(),
                 Computer=Sys.info()["nodename"],
                 systemInfo=Sys.info(),
                 sessionInfo=sessionInfo(),timeInit=Sys.time())


indv_fit_results <- tibble::lst(id_fits_alm_varied,id_fits_exam_varied,id_fits_alm_constant,id_fits_exam_constant, runInfo) 

saveRDS(indv_fit_results,
        file=here::here(paste0("data/indv_sim/","ind_abc_", n_prior_samples,"_",format(Sys.time(),"%H_%M_%OS"), ".rds")))





# id1 <- ds |> filter(id==1) 
# Model="EXAM"; Group="Varied"
# data=id1
# 
# sim <- full_sim_exam(data=as.data.table(id1), c=.02, lr=2.9, input_layer=input_layer, output_layer=output_layer, return_dat = return_dat)
# id1p<- id1 |> mutate(pred=sim1, resid=y-pred)
test1 <- id1 |> filter(expMode2=="Test") 
sd(test1$y)
id1 |> filter(expMode2=="Test") |> group_by(x) |> summarise(m=mean(y),sd=sd(y),n=n())
id1 |> filter(expMode2=="Train") |> group_by(x) |> summarise(m=mean(y),sd=sd(y),n=n()) |> ungroup() |> summarise(m=mean(m), sd=mean(sd))
 


