
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
testAll <- d |> filter(expMode2 == "Test")
testE1 <- e1 |> filter(expMode2 == "Test")
testE2 <- e2 |> filter(expMode2 == "Test")
testE3 <- e3 |> filter(expMode2 == "Test")

# testAll <- d |> filter(expMode == "test-Nf")
# testE1 <- e1 |> filter(expMode == "test-Nf") 
# testE2 <- e2 |> filter(expMode == "test-Nf") 
# testE3 <- e3 |> filter(expMode == "test-Nf") 


#prior = c(prior(normal(30,200),lb=0,class= Intercept))

modelName <- "combo_testVxBand_RF_5K.rds"
param_str <- paste(unlist(final_values), sep="_",collapse = "_")

iv2 <- "Bt_testNf"
iv2 <<- "Bt_test6_"

e1_e2_form <-  "condit * bandInt + (1 + bandInt|id)"
e3_form <- "condit * bandInt * bandOrder + (1 + bandInt|id)"
combo_form <- "condit * bandInt * bandOrder * fb + (1 + bandInt|id)"

# e1_e2_form <-  "condit * bandType + (1 + bandInt|id)"
# e3_form <- "condit * bandType * bandOrder + (1 + bandInt|id)"
# combo_form <- "condit * bandType * bandOrder * fb + (1 + bandInt|id)"


e1_e2_form <-  "condit * bandType + (1 + bandType|id)"
e3_form <- "condit * bandType * bandOrder + (1 + bandType|id)"
combo_form <- "condit * bandType * bandOrder * fb + (1 + bandType|id)"




run_model <- function(data, dv, formula_str, modelName_suffix, final_values) {
  formula <- as.formula(paste(dv, "~", formula_str))
  modelName <- paste0(modelName_suffix, "_",iv2, format(Sys.time(), "%M%OS"), ".rds")
  file_path <- paste0(here::here("data/model_cache", modelName))
  brm(formula,
      data = data,
      file = file_path,
      iter = final_values$niter,
      chains = final_values$nchain,
      silent = 1,
      #prior = prior,
      control = list(adapt_delta = final_values$adapt_delta, max_treedepth = final_values$max_treedepth))
}


run_e1_model <- function(dv, formula_str, final_values) {
  drun <- "e1"; mstage <- "Te"
  modelName_suffix <- paste0(drun, "_", mstage, str_to_title(dv),iv2, "_", paste(unlist(final_values), collapse = "_"))
  run_model(data = testE1, dv = dv, formula_str = formula_str, modelName_suffix = modelName_suffix, final_values = final_values)
}

run_e2_model <- function(dv, formula_str, final_values) {
  drun <- "e2"; mstage <- "Te"
  modelName_suffix <- paste0(drun, "_", mstage, str_to_title(dv),iv2, "_", paste(unlist(final_values), collapse = "_"))
  run_model(data = testE2, dv = dv, formula_str = formula_str, modelName_suffix = modelName_suffix, final_values = final_values)
}

run_e3_model <- function(dv, formula_str, final_values) {
  drun <- "e3"; mstage <- "Te"
  modelName_suffix <- paste0(drun,"_", mstage, str_to_title(dv),iv2, "_", paste(unlist(final_values), collapse = "_"))
  run_model(data = testE3, dv = dv, formula_str = formula_str, modelName_suffix = modelName_suffix, final_values = final_values)
}

run_combo_model <- function(dv, formula_str, final_values) {
  drun <- "combo"; mstage <- "Te"
  modelName_suffix <- paste0(drun, "_", mstage, str_to_title(dv),iv2, "_", paste(unlist(final_values), collapse = "_"))
  run_model(data = testAll, dv = dv, formula_str = formula_str, modelName_suffix = modelName_suffix, final_values = final_values)
}


# run_e1_model("vx", e1_e2_form, final_values)
# run_combo_model("vx", combo_form, final_values)

model_flags <- if(length(args) >= 5) strsplit(args[5], "")[[1]] else c("0", "0", "0", "0")
model_flags <- as.numeric(model_flags)  # Convert to numeric to enable logical operations

# Check if there are enough digits, pad with zeros if necessary
model_flags <- c(model_flags, rep(0, max(0, 4 - length(model_flags))))

print(model_flags)

if(model_flags[1] == 1) {
  print(paste("Running e1 model for vx and dist"))
  run_e1_model("vx", e1_e2_form, final_values)
  run_e1_model("dist", e1_e2_form, final_values)
}

if(model_flags[2] == 1) {
  print(paste("Running e2 model for vx and dist"))
  run_e2_model("vx", e1_e2_form, final_values)
  run_e2_model("dist", e1_e2_form, final_values)
}

if(model_flags[3] == 1) {
  print(paste("Running e3 model for vx and dist"))
  run_e3_model("vx", e3_form, final_values)
  run_e3_model("dist", e3_form, final_values)
}

if(model_flags[4] == 1) {
  print(paste("Running combo model for vx and dist"))
  run_combo_model("vx", combo_form, final_values)
  run_combo_model("dist", combo_form, final_values)
}



e1_e2_form <-  "condit * bandInt + (1 + bandInt|id)"
e3_form <- "condit * bandInt * bandOrder + (1 + bandInt|id)"
combo_form <- "condit * bandInt * bandOrder * fb + (1 + bandInt|id)"

iv2 <<- "Band_test6_"



if(model_flags[1] == 1) {
  print(paste("Running e1 model for vx and dist"))
  run_e1_model("vx", e1_e2_form, final_values)
  run_e1_model("dist", e1_e2_form, final_values)
}

if(model_flags[2] == 1) {
  print(paste("Running e2 model for vx and dist"))
  run_e2_model("vx", e1_e2_form, final_values)
  run_e2_model("dist", e1_e2_form, final_values)
}

if(model_flags[3] == 1) {
  print(paste("Running e3 model for vx and dist"))
  run_e3_model("vx", e3_form, final_values)
  run_e3_model("dist", e3_form, final_values)
}

if(model_flags[4] == 1) {
  print(paste("Running combo model for vx and dist"))
  run_combo_model("vx", combo_form, final_values)
  run_combo_model("dist", combo_form, final_values)
}