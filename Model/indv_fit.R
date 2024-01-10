pacman::p_load(dplyr,purrr,tidyr,ggplot2, data.table, here, patchwork, conflicted)
conflict_prefer_all("dplyr", quiet = TRUE)

walk(c("fun_alm","fun_model"), ~ source(here::here(paste0("Functions/", .x, ".R"))))

ds <- readRDS(here::here("data/e1_md_11-06-23.rds"))  |> as.data.table()
dsv <- ds |> filter(condit=="Varied")  
dsc <- ds |> filter(condit=="Constant") 


group_posterior_all <- readRDS(here::here("data/abc_2M_rmse_p001.rds"))  


teter <- group_posterior_all |> map_dfr(~tibble(pluck(.x$teter_results))) 
te <- group_posterior_all |> map_dfr(~tibble(pluck(.x$te_results))) 
tr <- group_posterior_all |> map_dfr(~tibble(pluck(.x$tr_results))) 




kde_results <- purrr::map(group_posterior_all, ~purrr::map_if(.x, is.list, ~compute_kde(.x$c, .x$lr, ngrid = 100, nsamples=100)))


n_prior_samples=500

# function that fits model to data



fit_indv <- function(sbj_id, simulation_function, prior_samples,Model, Group) {

  data <- filter(ds,id==sbj_id) |> as.data.table()
  input_layer =  c(100,350,600,800,1000,1200)
  output_layer = input_layer
  return_dat="train_data, test_data"
  train_idx <- which(data$expMode2=="Train")
  test_idx <- which(data$expMode2=="Test")
  
  sim_data_gen_s <- function(data, input_layer, output_layer, simulation_function, prior_samples, return_dat) {
    simulation_function <- match.fun(simulation_function)
    suppressMessages(simulation_results <- map_dfc(seq_len(nrow(prior_samples)), function(idx) {
      params <- prior_samples[idx, ]
      simulation_function(data=as.data.table(data), c=params$c, lr=params$lr, input_layer=input_layer, output_layer=output_layer, return_dat = return_dat)    
      }))
  }
  
  t1=system.time({
    sim_data <- sim_data_gen_s(data, input_layer, output_layer,simulation_function, prior_samples,return_dat) |> as.data.table()
  })

    
    target_data_train_test <- data[expMode2 %in% c("Test", "Train"), ]$y
    target_data_test <- data[expMode2 == "Test", ]$y
    target_data_train <- data[expMode2 == "Train", ]$y
    
    te_distances <- purrr::map_dbl(sim_data[test_idx, ], dist_rmse, observed = target_data_test)
    tr_distances <- purrr::map_dbl(sim_data[train_idx, ], dist_rmse, observed = target_data_train)
    teter_distances <- 0.5 * te_distances + 0.5 * tr_distances
    
    
    pct_keep=.1
    
    teter_results <- tibble(distance = teter_distances, c = prior_samples$c, 
                            lr = prior_samples$lr, sim_index= seq_along(teter_distances)) |> 
      arrange(distance) |> 
      filter(if (!is.null(pct_keep)) row_number() <= n() * pct_keep else distance <= tol) |>
      mutate(rank=row_number(),Model,Group,Fit_Method="Test & Train") |>
      rowwise() |>
      mutate(sim_dat = list(tibble(Model,Fit_Method,c,lr,distance,rank,data,pred=(sim_data[, .SD, .SDcols = sim_index])) |> 
                              mutate(resid=y-pred)))
    
    te_results <- tibble(distance = te_distances, c = prior_samples$c, 
                         lr = prior_samples$lr,  sim_index= seq_along(te_distances)) |> 
      arrange(distance) |> 
      filter(if (!is.null(pct_keep)) row_number() <= n() * pct_keep else distance <= tol) |>
      mutate(rank=row_number(),Model,Group,Fit_Method="Test Only") |>
      rowwise() |>
      mutate(sim_dat = list(tibble(Model,Fit_Method,c,lr,distance,rank,data,pred=(sim_data[, .SD, .SDcols = sim_index])) |> 
                              mutate(resid=y-pred
                                     )))
    
    tr_results <- tibble(distance = tr_distances, c = prior_samples$c, 
                         lr = prior_samples$lr,  sim_index= seq_along(tr_distances)) |>
      arrange(distance) |> 
      filter(if (!is.null(pct_keep)) row_number() <= n() * pct_keep else distance <= tol) |>
      mutate(rank=row_number(),Model,Group,Fit_Method="Train Only") |>
      rowwise() |>
      mutate(sim_dat = list(tibble(Model,Fit_Method,c,lr,distance,rank,data,pred=(sim_data[, .SD, .SDcols = sim_index])) |> 
                              mutate(resid=y-pred)))
    
    tibble::lst(teter_results = teter_results, te_results = te_results, tr_results = tr_results,sbj_id, Model, Group, pct_keep, fn_name=dist_rmse)
}


#id1 <- fit_indv(sbj_id=1,simulation_function=full_sim_exam, prior_samples=kde_results$abc_ev$teter_results$kde_samples)
#unique(dsv$id)

id_fits_exam_varied <- map(unique(dsv$id), ~ fit_indv(sbj_id=.x,
                                          simulation_function=full_sim_exam, 
                                          prior_samples=kde_results$abc_ev$teter_results$kde_samples, 
                                          Model="EXAM", 
                                          Group="Varied"))
  
  
  
id_fits_alm_varied <- map(unique(dsv$id), ~ fit_indv(sbj_id=.x,
                                                simulation_function=full_sim_alm, 
                                                prior_samples=kde_results$abc_ev$teter_results$kde_samples, 
                                                Model="ALM", 
                                                Group="Varied"))  



id_fits_exam_constant <- map(unique(dsc$id), ~ fit_indv(sbj_id=.x,
                                                      simulation_function=full_sim_exam, 
                                                      prior_samples=kde_results$abc_ev$teter_results$kde_samples, 
                                                      Model="EXAM", 
                                                      Group="Constant"))



id_fits_alm_constant <- map(unique(dsc$id), ~ fit_indv(sbj_id=.x,
                                                     simulation_function=full_sim_alm, 
                                                     prior_samples=kde_results$abc_ev$teter_results$kde_samples, 
                                                     Model="ALM", 
                                                     Group="Constant"))  
  






# id1 <- ds |> filter(id==1) 
# Model="EXAM"; Group="Varied"
# data=id1
# 
# sim <- full_sim_exam(data=as.data.table(id1), c=.02, lr=2.9, input_layer=input_layer, output_layer=output_layer, return_dat = return_dat)
# id1p<- id1 |> mutate(pred=sim1, resid=y-pred)

