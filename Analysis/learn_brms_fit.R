pacman::p_load(dplyr,purrr,tidyr,brms,bayestestR,stringr, here,conflicted)
walk(c("brms","dplyr"), conflict_prefer_all, quiet = TRUE)
defaults <- list(nchain=2, niter=1000,adapt_delta=.94, max_treedepth=11)
args <- commandArgs(trailingOnly = TRUE)
args_numeric <- map(args, as.numeric)
names_to_use <- names(defaults)[1:length(args_numeric)]
provided_values <- set_names(args_numeric, names_to_use)
final_values <- modifyList(defaults, provided_values)
invisible(list2env(final_values, envir = .GlobalEnv))
print(paste(names(final_values), ":", unlist(final_values), collapse = "; "))

options(brms.backend="cmdstanr",mc.cores=nchain)

e1 <- readRDS(here("data/e1_08-21-23.rds")) 
e2 <- readRDS(here("data/e2_08-04-23.rds")) 
e3 <- readRDS(here("data/e3_08-04-23.rds")) 
d <- rbind(e1,e2,e3)




e1_e2_form <-  "condit * bandType + (1 + bandInt|id)"
e3_form <- "condit * bandType * bandOrder + (1 + bandInt|id)"
combo_form <- "condit * bandType * bandOrder * fb + (1 + bandInt|id)"

nlform <- bf(
  dist ~ a + (b - a) * exp(-c * gt.train),
  a ~ 0 + condit + bandInt + (1|id), 
  b ~ 0 + condit + bandInt + (1|id), 
  c ~ 0 + condit + (1|id), 
  nl = TRUE 
)





nlform <- bf(
  dist ~ a + (b - a) * exp(-c * gt.train),
  a ~ 1 + condit + bandInt,
  b ~ 1 + condit + bandInt,
  c ~ 1 + condit + bandInt,
  nl = TRUE 
)

prior <- c(
  prior(normal(200, 250), lb=0, nlpar = "a"),
  prior(normal(800, 400), lb=0,  nlpar = "b"),
  prior(normal(.3, .3), lb=0,  nlpar = "c")
)



drun <- "e1"; mstage <- "Train_"; iv2="conBand"; dv="dist"
modelName_suffix <- paste0(drun, "_", mstage, str_to_title(dv),iv2, "_", paste(unlist(final_values), collapse = "_"))
modelName <- paste0(modelName_suffix, "_", format(Sys.time(), "%M%OS"), ".rds")
file_path <- paste0(here::here("data/model_cache", modelName))


fit <- brm(
  formula = nlform,
  data = e1 |> filter(expMode2=="Train"),
  family = gaussian(), # Assuming 'dist' is continuous; adjust as needed
  prior = prior,
  iter = final_values$niter,
  chains = final_values$nchain,
  silent = 0,
  control = list(adapt_delta = final_values$adapt_delta, max_treedepth = final_values$max_treedepth),
  file = file_path
)

print(summary(fit))
print(bayestestR::describe_posterior(fit) )


  modelNameTxt <- sub("\\.rds$", ".txt", modelName) # Replace .rds with .txt
  file_path_txt <- paste0(here::here("data/model_cache/brms", modelNameTxt))
  sink(file_path_txt)
  print(summary(fit))
  print(bayestestR::describe_posterior(fit))
  sink()



print("Start E1 RF learn fits")

nlform <- bf(
  dist ~ a + (b - a) * exp(-c * gt.train),
  a ~ 0 + condit + (1|bandInt) + (1|id), 
  b ~ 0 + condit + (1|bandInt)  + (1|id), 
  c ~ 0 + condit + (1|bandInt) + (1|id), 
  nl = TRUE 
)



drun <- "e1"; mstage <- "Train_"; iv2="conBandRF"; dv="dist"
modelName_suffix <- paste0(drun, "_", mstage, str_to_title(dv),iv2, "_", paste(unlist(final_values), collapse = "_"))
modelName <- paste0(modelName_suffix, "_", format(Sys.time(), "%M%OS"), ".rds")
file_path <- paste0(here::here("data/model_cache", modelName))


fit <- brm(
  formula = nlform,
  data = e1 |> filter(expMode2=="Train"),
  family = gaussian(), # Assuming 'dist' is continuous; adjust as needed
  prior = prior,
  iter = final_values$niter,
  chains = final_values$nchain,
  silent = 1,
  control = list(adapt_delta = final_values$adapt_delta, max_treedepth = final_values$max_treedepth),
  file = file_path
)

print(summary(fit))
print(bayestestR::describe_posterior(fit) )


  modelNameTxt <- sub("\\.rds$", ".txt", modelName) # Replace .rds with .txt
  file_path_txt <- paste0(here::here("data/model_cache/brms", modelNameTxt))
  sink(file_path_txt)
  print(summary(fit))
  print(bayestestR::describe_posterior(fit))
  sink()



##### E2
print("start E2")



nlform <- bf(
  dist ~ a + (b - a) * exp(-c * gt.train),
  a ~ 1 + condit + bandInt,
  b ~ 1 + condit + bandInt,
  c ~ 1 + condit + bandInt,
  nl = TRUE 
)

drun <- "e2"; mstage <- "Train_"; iv2="conBand"; dv="dist"
modelName_suffix <- paste0(drun, "_", mstage, str_to_title(dv),iv2, "_", paste(unlist(final_values), collapse = "_"))
modelName <- paste0(modelName_suffix, "_", format(Sys.time(), "%M%OS"), ".rds")
file_path <- paste0(here::here("data/model_cache", modelName))


fit <- brm(
  formula = nlform,
  data = e2 |> filter(expMode2=="Train"),
  family = gaussian(), # Assuming 'dist' is continuous; adjust as needed
  prior = prior,
  iter = final_values$niter,
  chains = final_values$nchain,
  silent = 0,
  control = list(adapt_delta = final_values$adapt_delta, max_treedepth = final_values$max_treedepth),
  file = file_path
)

print(summary(fit))
print(bayestestR::describe_posterior(fit) )

  modelNameTxt <- sub("\\.rds$", ".txt", modelName) # Replace .rds with .txt
  file_path_txt <- paste0(here::here("data/model_cache/brms", modelNameTxt))
  sink(file_path_txt)
  print(summary(fit))
  print(bayestestR::describe_posterior(fit))
  sink()

print("Start E2 RF learn fits")

nlform <- bf(
  dist ~ a + (b - a) * exp(-c * gt.train),
  a ~ 0 + condit + (1|bandInt) + (1|id), 
  b ~ 0 + condit + (1|bandInt)  + (1|id), 
  c ~ 0 + condit + (1|bandInt) + (1|id), 
  nl = TRUE 
)



drun <- "e2"; mstage <- "Train_"; iv2="conBandRF"; dv="dist"
modelName_suffix <- paste0(drun, "_", mstage, str_to_title(dv),iv2, "_", paste(unlist(final_values), collapse = "_"))
modelName <- paste0(modelName_suffix, "_", format(Sys.time(), "%M%OS"), ".rds")
file_path <- paste0(here::here("data/model_cache", modelName))


fit <- brm(
  formula = nlform,
  data = e2 |> filter(expMode2=="Train"),
  family = gaussian(), # Assuming 'dist' is continuous; adjust as needed
  prior = prior,
  iter = final_values$niter,
  chains = final_values$nchain,
  silent = 1,
  control = list(adapt_delta = final_values$adapt_delta, max_treedepth = final_values$max_treedepth),
  file = file_path
)

print(summary(fit))
print(bayestestR::describe_posterior(fit) )



  modelNameTxt <- sub("\\.rds$", ".txt", modelName) # Replace .rds with .txt
  file_path_txt <- paste0(here::here("data/model_cache/brms", modelNameTxt))
  sink(file_path_txt)
  print(summary(fit))
  print(bayestestR::describe_posterior(fit))
  sink()



##### E3
print("start E3")


nlform <- bf(
  dist ~ a + (b - a) * exp(-c * gt.train),
  a ~ 1 + condit *bandOrder + bandInt,
  b ~ 1 + condit *bandOrder + bandInt,
  c ~ 1 + condit *bandOrder + bandInt,
  nl = TRUE 
)


drun <- "e3"; mstage <- "Train_"; iv2="conBand"; dv="dist"
modelName_suffix <- paste0(drun, "_", mstage, str_to_title(dv),iv2, "_", paste(unlist(final_values), collapse = "_"))
modelName <- paste0(modelName_suffix, "_", format(Sys.time(), "%M%OS"), ".rds")
file_path <- paste0(here::here("data/model_cache", modelName))


fit <- brm(
  formula = nlform,
  data = e3 |> filter(expMode2=="Train"),
  family = gaussian(), # Assuming 'dist' is continuous; adjust as needed
  prior = prior,
  iter = final_values$niter,
  chains = final_values$nchain,
  silent = 0,
  control = list(adapt_delta = final_values$adapt_delta, max_treedepth = final_values$max_treedepth),
  file = file_path
)

print(summary(fit))
print(bayestestR::describe_posterior(fit) )


  modelNameTxt <- sub("\\.rds$", ".txt", modelName) # Replace .rds with .txt
  file_path_txt <- paste0(here::here("data/model_cache/brms", modelNameTxt))
  sink(file_path_txt)
  print(summary(fit))
  print(bayestestR::describe_posterior(fit))
  sink()

print("Start E3 RF learn fits")

nlform <- bf(
  dist ~ a + (b - a) * exp(-c * gt.train),
  a ~ 0 + condit *bandOrder + (1|bandInt) + (1|id), 
  b ~ 0 + condit *bandOrder + (1|bandInt)  + (1|id), 
  c ~ 0 + condit *bandOrder + (1|bandInt) + (1|id), 
  nl = TRUE 
)



drun <- "e3"; mstage <- "Train_"; iv2="conBandRF"; dv="dist"
modelName_suffix <- paste0(drun, "_", mstage, str_to_title(dv),iv2, "_", paste(unlist(final_values), collapse = "_"))
modelName <- paste0(modelName_suffix, "_", format(Sys.time(), "%M%OS"), ".rds")
file_path <- paste0(here::here("data/model_cache", modelName))


fit <- brm(
  formula = nlform,
  data = e3 |> filter(expMode2=="Train"),
  family = gaussian(), # Assuming 'dist' is continuous; adjust as needed
  prior = prior,
  iter = final_values$niter,
  chains = final_values$nchain,
  silent = 1,
  control = list(adapt_delta = final_values$adapt_delta, max_treedepth = final_values$max_treedepth),
  file = file_path
)

print(summary(fit))
print(bayestestR::describe_posterior(fit) )


  modelNameTxt <- sub("\\.rds$", ".txt", modelName) # Replace .rds with .txt
  file_path_txt <- paste0(here::here("data/model_cache/brms", modelNameTxt))
  sink(file_path_txt)
  print(summary(fit))
  print(bayestestR::describe_posterior(fit))
  sink()