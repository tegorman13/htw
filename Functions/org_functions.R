

gen_trials <- function(blocks=3,numinputs=3,shuffle=FALSE) {
  if (shuffle) {
    examples <- as.vector(apply(replicate(blocks,seq(1, numinputs)), 2,
                                sample, numinputs))
  } else{
    examples <- as.vector(replicate(blocks, seq(1, numinputs)))
  }
}

custom_scale <- function(x) {
  if (sd(x) == 0) {
    return(as.vector(rep(0, length(x))))
  } else {
    return(as.vector(scale(x)))
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






brmsfit_to_df <- function(brmsfit, predictor_to_label, odds_ratios=FALSE){
  
  # Get fixed effects ------------------------------------------------------####
  
  fixed_effects <- as.data.frame(fixef(brmsfit, summary = TRUE))
  
  # Get posterior probabilities of Est > 0, for null model -----------------####
  
  print("Posterior probabilities:")
  
  post_samples <- posterior_samples(brmsfit)
  
  fixed_effects[,"P(Est. $>$ 0)"] <- NA
  
  for (fixed_effect in rownames(fixed_effects)){
    print(fixed_effect)
    post_prob <- round(mean(post_samples[,paste0("b_",fixed_effect)] > 0), 2)
    print(post_prob)
    fixed_effects[fixed_effect,"P(Est. $>$ 0)"] <- post_prob
  }
  
  # Convert results to latex table -----------------------------------------####
  
  # add predictors in a column
  fixed_effects$Predictor <- rownames(fixed_effects)
  
  if (odds_ratios) {
    # Convert log-odds estimates to odds ratio's by exponentiating
    fixed_effects[,c("Estimate","Q2.5","Q97.5")] <- 
      exp(fixed_effects[,c("Estimate","Q2.5","Q97.5")])
  }
  
  # round estimates
  fixed_effects[,c("Estimate","Q2.5","Q97.5")] <- 
    round(fixed_effects[,c("Estimate","Q2.5","Q97.5")], digits = 2)
  
  # combine Q2.5 and Q97.5 into Credible Interval
  
  fixed_effects[,"95\\% CrI"] <- 
    paste0('[',sprintf("%.2f",fixed_effects$Q2.5), ", 
           ", sprintf("%.2f",fixed_effects$Q97.5),']')
  
  # Drop Estimated Error column
  fixed_effects <- fixed_effects[,c("Predictor", "Estimate","95\\% CrI",
                                    "P(Est. $>$ 0)")]
  
  if (odds_ratios) {
    colnames(fixed_effects)[colnames(fixed_effects) == "Estimate"] <- 
      "Odds Ratio"
  }
  
  # Rename predictors
  for (i in 1:nrow(predictor_to_label)){
    fixed_effects[
      fixed_effects$Predictor == predictor_to_label[i,"predictor"],c("Predictor")] <- 
      predictor_to_label[i,"label"]
    # print(predictor_to_label[i,"predictor"])
    # print(predictor_to_label[i,"label"])
  }
  
  return(fixed_effects)
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



