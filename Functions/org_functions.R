

#walk(paste0("Functions/", c("Display_Functions.R", "org_functions.R")), ~source(here::here(.)))


# options(digits = 2, width = 120,
#         dplyr.summarise.inform = FALSE,
#         knitr.kable.NA = "")


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



# ks2 = model_parameters(e1_testDistRF2_0,effects="random",keep="^r_id*") 
# ks2 <- ks2 %>%
#   mutate(
#     id = str_extract(Parameter, "(?<=\\[)\\d+"), # Extract digits after the opening square bracket
#     vb = str_extract(Parameter, "\\d+M\\d+") %>%
#       str_remove("vb") %>%
#       str_replace("M", "-")
#   )

GetIndvFits <- function(model) {
ks2 = parameters::model_parameters(model,effects="random",keep="^r_id*") 
ks2 <- ks2 |>
  mutate(
    id = factor(gsub(".*\\[([0-9]+),.*", "\\1", Parameter)),
    vb = factor(gsub(".*vb([0-9]+)M([0-9]+).*", "\\1-\\2", Parameter), levels = levels(test$vb))
  ) |>
  left_join(select(model$data,id,condit) |> distinct(),by=join_by("id"),keep=FALSE) |>
  select(-Group) |> relocate(id,condit,vb)
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



# p1 <- GetModelStats(e1_testDistRF2_0)
# 
# kable(p1,escape=F, booktabs=T) |> column_spec(1:8,width="1.5em") |> 
#   kable_styling(full_width = F)

GetModelStats <- function(model, type="brms") {
  
    # Get current model stats for brms model
    if (type == "brms") {
      m1 <- as.data.frame(describe_posterior(model, centrality = "Median"))
      m2 <- fixef(model)
      df <- cbind(m1[, c(1,2)], m2[, 2], m1[, c(4,5, 6, 11, 12)])
      # colnames(df) <- c("Term", "\\(\\beta_{Median}\\)", "SD",
      #                   "95% CrI \nLower", "95% CrI \nUpper", "pd", "\\(\\widehat R\\)", "ESS")
      
      colnames(df) <- c("Term", "Estimate", "Est.Error",
                        "95% CrI Lower", "95% CrI Upper", "pd", "Rhat", "ESS")
      
      # Add model name and re-order columns
      #df$Model <- mnames[n]
      #df <- df[, c(9, 1:8)]
      
      # Get current model stats for lmer model
    } else {
      df <- data.frame(summary(model$coefficients))
      df$Parameter <- rownames(df)
      df <- df[, c(7, 6, 1:5)]
    }

  
  df <- df |> 
    mutate(across(where(is.numeric), \(x) round(x, 2))) |>
    mutate(Rhat = round(Rhat, 3)) |>
    tibble::remove_rownames() |> 
    # Replace "bandInt" with "Band" wherever it occurs in the Term column
    mutate(Term = stringr::str_replace_all(Term, "bandInt", "Band")) |>
    # Remove "b_" from the start of the Term column if it's present
    mutate(Term = stringr::str_replace_all(Term, "^b_", ""))
  
  return(df)
}


GetBrmsModelStats <- function(model1, model2 = NA, BF = NA) {
  
  cat("Calculating stats for model1....\n")
  LL1 <- logLik(model1) %>% apply(., 2, map_estimate) %>% sum()
  cat(paste("Model1 logLikelihood:", LL1, "\n"))
  Deviance1 <- -2 * LL1
  cat(paste("Model1 Deviance:", Deviance1, "\n"))
  
  if (class(model2) == "logical") {
    loo1 <- loo(model1)$loo
    cat(paste("Model1 loo:", loo1, "\n"))
    waic1 <- waic(model1)$waic
    cat(paste("Model1 WAIC:", waic1, "\n\n"))
    return(c("logLik" = LL1, "Deviance" = Deviance1, "loo" = loo1, "WAIC" = waic1))
  }
  
  
  cat("\nCalculating stats for model2....\n")
  LL2 <- logLik(model2) %>% apply(., 2, map_estimate) %>% sum()
  cat(paste("Model2 logLikelihood:", LL2, "\n"))
  Deviance2 <- -2 * LL2
  cat(paste("Model2 Deviance:", Deviance2, "\n"))
  
  modelsLoo <- loo_compare(add_criterion(model1, "loo"), add_criterion(model2, "loo")) %>% data.frame()
  modelsWaic <- loo_compare(add_criterion(model1, "waic"), add_criterion(model2, "waic"), criterion = "waic") %>% data.frame()
  
  loo1 <- modelsLoo["add_criterion(model1, \"loo\")", c("looic", "se_looic")]
  loo2 <- modelsLoo["add_criterion(model2, \"loo\")", c("looic", "se_looic")]
  cat(paste("\nModel1 loo:", loo1[1], "SE:", loo1[2], "\n"))
  cat(paste("\nModel2 loo:", loo2[1], "SE:", loo2[2], "\n"))
  
  waic1 <- modelsWaic["add_criterion(model1, \"waic\")",  c("waic", "se_waic")]
  waic2 <- modelsWaic["add_criterion(model2, \"waic\")", c("waic", "se_waic")]
  cat(paste("\nModel1 waic:", waic1[1], "SE:", waic1[2], "\n"))
  cat(paste("\nModel2 waic:", waic2[1], "SE:", waic2[2], "\n"))
  
  model1 = c("logLik" = LL1, "Deviance" = Deviance1, "loo" = loo1, "WAIC" = waic1)    
  model2 = c("logLik" = LL2, "Deviance" = Deviance2, "loo" = loo2, "WAIC" = waic2)
  
  cat(paste("\nWaic difference:", -2 * modelsWaic[2, 1], "SE:", 2 * modelsWaic[2, 2],  "\n"))
  cat(paste("\nLoo difference:", -2 * modelsLoo[2, 1], "SE:", 2 * modelsLoo[2, 2],  "\n"))
  
  if (!class(BF) == "logical") {
    cat("\nExtracting BF stats...\n")
    BF = BF$bf
    cat(paste("Bayes Factor:", BF, "\n"))
    pProb = BF/(1+BF)
    cat(paste("Posterior probability in favour of model1:", pProb, "\n"))
    Comparisons = c("loo_diff" = c(-2 * modelsLoo[2, 1], 2 * modelsLoo[2, 2]), "waic_diff" = c(-2 * modelsWaic[2, 1], 2 * modelsWaic[2, 2]), "BF" = BF, "pProb" = pProb)
    
  } else {Comparisons = c("loo_diff" = c(-2 * modelsLoo[2, 1], 2 * modelsLoo[2, 2]), "waic_diff" = c(-2 * modelsWaic[2, 1], 2 * modelsWaic[2, 2]))}
  
  return(list("Model1" = model1, "Model2" = model2, "Comparisons" = Comparisons))
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



