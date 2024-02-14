
pacman::p_load(dplyr,purrr,tidyr,ggplot2, data.table, here, patchwork, conflicted,coda)
conflict_prefer_all("dplyr", quiet = TRUE)
walk(c("fun_alm","fun_model"), ~ source(here::here(paste0("Functions/", .x, ".R"))))
set.seed(123)

ds <- readRDS(here::here("data/e1_md_11-06-23.rds"))  |> as.data.table()
dsv <- ds |> filter(condit=="Varied")  
dsc <- ds |> filter(condit=="Constant") 

dsId <- ds |> select(id,condit) |> unique()
ids <- c(1,2,4,5,6,7,8, 10,11,12,13)


lg_generate_prior_c_lr <- function(n,cMean, cSig=2,lrSig=1) {
  prior_samples <- tibble(
    c = rlnorm(n,-5,sdlog=cSig),
    lr = extraDistr::rhnorm(n,sigma=lrSig),
  )
  return(prior_samples)
}

n=5000

cMean <<- -7.0; cSig <<- 1.5; lrSig <<- 2.0
prior_samples <- lg_generate_prior_c_lr(n=6000, cMean=cMean, cSig=cSig, lrSig=lrSig) 

mean(prior_samples$c)
median(prior_samples$c)
min(prior_samples$c)
max(prior_samples$c)
quantile(prior_samples$c)
plot(density(prior_samples$c))

mean(prior_samples$lr)
median(prior_samples$lr)
min(prior_samples$lr)
max(prior_samples$lr)
plot(density(prior_samples$lr))



list.files("data/abc_pda/")
list.files("data/abc_pda/n_iter_3000_nc_6_130909/")
list.files("data/abc_pda/n_iter_3000_nc_6_130909/",pattern="EXAM_Test")
list.files("data/abc_pda/n_iter_3000_nc_6_130909/",pattern="EXAM_Test")


exam_test <- readRDS(here(grep("Train", 
                                list.files("data/abc_pda/n_iter_3000_nc_6_130909/", 
                                           pattern="EXAM_Test",full.names = TRUE), 
                                invert=TRUE, value=TRUE))) |> pluck("exam_test")

exam_test <- readRDS(here("data/abc_pda/pda_EXAM_Test_8000_6_212144.rds")) |> pluck("exam_test")
exam_test <- readRDS(here("data/abc_pda/pda_EXAM_Test_12000_3_002125.rds")) |> pluck("exam_test")


k = rbindlist(exam_test)|> left_join(dsId, join_by(id))
k |> filter(condit=="Varied") |> ggplot(aes(x=lr, y=c, color=condit)) + geom_point() 
#tidyselect:::select("1","5","10","33","66") 
#tidyselect:::select("1","33","66") 

map(exam_test |> tidyselect:::select("1","33","66")  , possibly(~{
  trim = .x |> group_by(chain) |> filter(row_number() < nrow(.x)/1.1) |> ungroup()  
  post <- split(trim |> select(c,lr), f = trim$chain)
  mcmc_list <- lapply(post, as.mcmc, start=nr/2,thin=4)
  (gelman_result <- gelman.diag(mcmc_list))
} ))

map(exam_test |> tidyselect:::select("1","33","66") , ~{
  post <- split(.x |> select(c,lr), f = .x$chain)
  mcmc_list <- lapply(post, as.mcmc)
  plot(mcmc.list(mcmc_list), sub=.x$id[1]) 
} )

map(exam_test |> tidyselect:::select("1","5","10"), ~{
  post <- split(.x |> select(c,lr), f = .x$chain)
  mcmc_list <- lapply(post, as.mcmc)
  plot(mcmc.list(mcmc_list)) 
} )

library(tidybayes)

library(ggmcmc)

k1 = exam_test |> pluck("1")
kmc <- as.mcmc(split(k1 |> select(c,lr), f = k1$chain))
kmc <- lapply(split(k1 |> select(c,lr), f = k1$chain),as.mcmc)
as.mcmc.list(kmc)
ggs(kmc)

gelman.plot(mcmc.list(kmc))
plot(mcmc.list(kmc))

k1 |> group_by(chain) |> mutate(iteration=1:n(), chain=as.factor(chain)) |> 
  filter(iteration>10) |> 
  ggplot(aes(x=iteration, y=c, color=chain)) + geom_line() 


# data <- dsv |> filter(id == 1)
# chain_list <- pda_abc(full_sim_exam, prior_samples, data, num_iterations = 4000, num_chains = 2)


# chains_df <- imap_dfr(chain_list, ~mutate(.x, chain = .y))
# chain_dfs2 <- split(chains_df |> select(-chain), f = chains_df$chain)
# mcmc_list <- lapply(chain_dfs2, as.mcmc)
# (gelman_result <- gelman.diag(mcmc_list))
# plot(mcmc.list(mcmc_list))

# simulation_function <- full_sim_exam
# prior_samples <- lg_generate_prior_c_lr(n=1000) 
# num_iterations = 100
# num_chains = 4

# posterior_predictive_data <- generate_posterior_predictive(full_sim_exam, chains_by_sbj[[2]], data, num_samples = 100)

# ids1 <- c(1,36,66, 76,101,192)


ids1 <- as.numeric(levels(ds$id))
prior_samples <- lg_generate_prior_c_lr(n=5000) 
subjects_data <-  ds |> filter(id %in% ids1)  %>% split(f =c(.$id), drop=FALSE)

# chains_by_sbj <- map(subjects_data, ~pda_abc(simulation_function = full_sim_exam, 
#                                                  prior_samples = prior_samples, 
#                                                  data = .x, 
#                                                  num_iterations = 10, 
#                                                  num_chains = 1)) %>% # 
#                                                  setNames(ids1)
                                                


# 3070  
# multicore doesn't work 
# 3857s (1 hour 4 min) for 5000 iter; 4 chains on windows_3070. 
# 2238 for 3000 iter, 4 chain on 3070
# 2132; 3082s teter; 2270 alm teter for 4000 iter, 4 chains on 3070 (cluster plan - 14 cores)
# 1497s 1060s; 1521; 1087; 1500-1068  -for 2000 iter, 4 chains (multicore 15 cores)

# tg_m1
# 887 secs for 2K iterations, 4 chains on tg_m1
# 1997s; 1420s; 2066s teter; 1512s; 2189; 1538 for 3000 iter, 6 chain on tg_m1 (cluster plan 8 cores)

# M1 iMacs
# 3559-2529-3675-2676-2539 (m1l4) for 4000 iter, 4 chains (cluster plan - 6 cores)
# 3506-2496-3630-2683-2577 (m1l1 ) 
# 3277 (m1l4) 4000 iter, 4 chains (multicore plan - 8 cores)
# 3318 (m1l1) 4000 iter, 4 chains (multicore plan - 8 cores)

# 3476 on m1r1 4000 iter, 4 chains (multicore plan - 6 cores)
# 4509; 3208; 4671; 3430; 4320-3102 for m1l2 for 3000 iter, 6 chain (multicore plan - 6 cores)


# (m1r2) 1812-1294-1884-1374-1740-1245 - 2000 iter, 4 chains (multicore 7 cores?)
# (m1l3) 1805-1295-1871-1361-1750 1258 - 2000 iter, 4 chains (multicore 7 cores?)
# (m1l5) 1638-1174-1696-1240 2000 iter, 4 chains (multicore 8 cores)






map(chains_by_sbj, ~{
  chain_dfs2 <- split(.x |> select(c,lr), f = .x$chain)
  mcmc_list <- lapply(chain_dfs2, as.mcmc)
  plot(mcmc.list(mcmc_list))
} )


map(chains_by_sbj, possibly(~{
  nr = nrow(.x)
  trim = .x |> group_by(chain) |> filter(row_number() < nr/2) |> ungroup()  
  chain_dfs2 <- split(trim |> select(c,lr), f = trim$chain)
  mcmc_list <- lapply(chain_dfs2, as.mcmc, start=nr/2,thin=4)
  (gelman_result <- gelman.diag(mcmc_list))
} ))


map(chains_by_sbj, ~{
  .x |> select(c,lr) |> head(2)
} )


map(chains_by_sbj, ~{
  nr = nrow(.x)
  trim=.x[nr/2:nr, ]
  paste0("c: ", mean(trim$c), " lr: ", mean(trim$lr))
} )


map(chains_by_sbj, ~{
  unique(.x$chain) %>% print()
  nr = nrow(.x)
  trim = .x |> group_by(chain) |> filter(row_number() < nr/2) |> ungroup()
  
  unique(trim$chain) %>% print()
  { trim |> ggplot(aes(x=c,fill=as.factor(chain))) +geom_density() } +
  trim |> ggplot(aes(x=lr,fill=as.factor(chain))) +geom_density() 
  
} )



map(chains_by_sbj, ~{
  unique(.x$chain) %>% print()
  nr = nrow(.x)
  trim = .x |> group_by(chain) |> filter(row_number() < nr/2) |> ungroup()  
  posterior_predictive_data <- generate_posterior_predictive(full_sim_exam, .x, data, num_samples = 100)
  head(posterior_predictive_data$sim_rank)
  head(posterior_predictive_data$sim_vx_agg)
} )

chains_by_sbj = exam_test
chains_by_sbj = alm_test
chains_by_sbj |> tidyselect:::select("1") %>% map(., ~{
  trim = .x |> group_by(chain) |> filter(row_number() < (nrow(.x))/1.5) |> ungroup() 
  data <- ds |> filter(id == trim$id[1])
  posterior_predictive_data <- generate_posterior_predictive(full_sim_exam, trim, data, num_samples = 100)
   print({ trim |> ggplot(aes(x=c,fill=as.factor(chain))) +geom_density() } +
  trim |> ggplot(aes(x=lr,fill=as.factor(chain))) +geom_density(alpha=.5) )
  print(head(posterior_predictive_data$sim_rank))
  head(posterior_predictive_data$sim_vx_agg)
} )




k = rbindlist(chains_by_sbj)|> left_join(dsId, join_by(id))


k |> filter(c<.1) |> ggplot(aes(x=c,fill=condit)) + geom_density()
k  |> ggplot(aes(x=lr,fill=condit)) + geom_density()

k  |> filter(rank<10) |>  ggplot(aes(x=condit,y=lr)) + geom_col()

posterior_predictive_data <- generate_posterior_predictive(full_sim_exam, chains_by_sbj[[2]], data, num_samples = 100)




mean(chains_by_sbj[[1]]$c)
mean(chains_by_sbj[[2]]$c)


map(chains_by_sbj, ~{
  head(.x,5)

} )


# Processing the chains
chain_dfs <- lapply(chain_list, function(chain) {
  chain_df <- do.call(rbind, chain) |> as.data.frame()
  colnames(chain_df) <- c("c", "lr")
  chain_df[100:nrow(chain_df), ] # trimming the first 1000 iterations
})
combined_chain_df <- do.call(rbind, chain_dfs) |> as.data.frame()

#do.call(rbind, chain_list[[1]]) |> as.data.frame()


library(coda)
mcmc_chains <- lapply(chain_list, function(chain) {
  chain_df <- do.call(rbind, chain) |> as.data.frame()
  colnames(chain_df) <- c("c", "lr")
  chain_df[1000:nrow(chain_df), ] # trimming the first 1000 iterations
  mcmc(as.matrix(chain_df))
})

# Create an mcmc.list object for diagnostics
mcmc_list <- mcmc.list(mcmc_chains)

# Gelman-Rubin Diagnostic
gelman_diag <- gelman.diag(mcmc_list)
print(gelman_diag)

plot(mcmc_list)

for (param in colnames(combined_chain_df)) {
    print(param)
 g1<- ggplot(data = bind_rows(chain_dfs, .id = 'Chain'), aes_string(x = param, group = 'Chain', color = 'Chain')) +
    geom_density() +
    labs(title = paste("Density Plot for", param))
print(g1)
}



summary_stats <- combined_chain_df %>% 
  summarise(across(everything(), list(mean = mean, sd = sd, quantile25 = ~quantile(., probs = 0.25), median = median, quantile75 = ~quantile(., probs = 0.75))))
print(summary_stats)


summary_stats_c <- summary(combined_chain_df$c)
summary_stats_lr <- summary(combined_chain_df$lr)

# Printing summary statistics
print(summary_stats_c)
print(summary_stats_lr)



generate_posterior_predictive <- function(simulation_function, combined_chain_df, data, num_samples = 5) {a

  input_layer =  c(100,350,600,800,1000,1200)
  output_layer = input_layer
  return_dat = "test_data"
  test_idx <- which(data$expMode2 == "Test")
  target_data_test <- data[data$expMode2 == "Test", ]$y
  # Randomly select parameter sets from the posterior
  sampled_params <- combined_chain_df[sample(1:nrow(combined_chain_df), num_samples), ]
  
  # Generate simulated data for each sampled parameter set
  simulated_data_list <- lapply(1:nrow(sampled_params), function(i) {
    params <- sampled_params[i, ]
    sim_data <- simulation_function(data, params$c, params$lr, input_layer = input_layer, output_layer = output_layer, return_dat = "test_data") |> as_tibble() |>
        mutate(y=target_data_test,x=data[expMode2 == "Test", ]$x,c=params$c, lr=params$lr, sim=i ) |> rename(pred="value")
  })

  # Combine simulated data
  sim_trial_data <- bind_rows(simulated_data_list)

  sim_vx_agg <- sim_trial_data |> group_by(sim,x,c,lr) |> summarise(pred=mean(pred), y=mean(y), error=abs(y-pred)) |> ungroup() |> group_by(sim,c,lr) |>
    mutate(meanError=mean(error)) |> ungroup() |>
    mutate(rank = dense_rank(meanError)) |>
    arrange(rank)
   

  sim_rank <- sim_trial_data |> group_by(sim,c,lr) |> mutate(error=abs(y-pred)) |> summarise(mae=mean(error)) |> arrange(mae)

  tibble::lst(sim_vx_agg,sim_rank)

}

# Generate posterior predictive data
posterior_predictive_data <- generate_posterior_predictive(full_sim_exam, combined_chain_df, data, num_samples = 100)


