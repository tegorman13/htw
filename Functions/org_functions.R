

gen_trials <- function(blocks=3,numinputs=3,shuffle=FALSE) {
  if (shuffle) {
    examples <- as.vector(apply(replicate(blocks,seq(1, numinputs)), 2,
                                sample, numinputs))
  } else{
    examples <- as.vector(replicate(blocks, seq(1, numinputs)))
  }
}



normal_predictive_distribution <- function(mu_samples,
                                           sigma_samples,
                                           N_obs) {
  map2_dfr(mu_samples, sigma_samples, function(mu, sigma) {
    tibble(
      trialn = seq_len(N_obs),
      t_pred = rnorm(N_obs, mu, sigma)
    )
  }, .id = "iter") %>%
    # .id is always a string and
    # needs to be converted to a number
    mutate(iter = as.numeric(iter))
}



# pathFN <- paste0(path,'/', modelName, ".rds")
# saveRDS(fit, file = pathFN)

#Wrapper for brm models such that it saves the full model the first time it is run, otherwise it loads it from disk
# run_model <- function(expr, modelName, path=here::here("data/model_cache"), reuse = TRUE, chains=4,parallel=TRUE) {
#   
#   if (parallel) {
#     withr::local_options(list(mc.cores =  parallel::detectCores()))
#     if (chains >  parallel::detectCores()) {
#       chains <-  parallel::detectCores()
#     }
#   }
#   
#   
#   pathFN <- paste0(path,'/', modelName, ".rda")
#   if (reuse) {
#     fit <- suppressWarnings(try(readRDS(pathFN), silent = TRUE))
#   }
#   if (is(fit, "try-error")) {
#     fit <- eval(expr)
#     saveRDS(fit, file = pathFN, compress = "xz")
#   }
#   fit
# }



