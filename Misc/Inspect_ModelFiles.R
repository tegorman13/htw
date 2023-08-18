source(here::here("Functions", "packages.R"))

library(tidyverse)
library(brms)


# file = file_list[1]
# result <- readRDS(file_list[1])

extract_info <- function(file) {
  result <- try(readRDS(file), silent = TRUE)
  if (inherits(result, "try-error")) {
    return(NULL)
  } else {
    file_name <- basename(file)
    md <- strftime(file.info(file)$mtime, format = "%Y-%m-%d %H:%M")
    model_formula <- deparse(result$formula$formula)
    model_dv <- result$formula$resp
    m_fixed <- tryCatch(paste0(result$formula$formula[[3]][[2]], collapse=" "), error = function(e) NA)
    m_rand <- tryCatch(paste0(result$formula$formula[[3]][[3]],collapse=" "), error = function(e) NA)
    m_family <- result$family$family
    
    bayesR2 <- tryCatch(bayes_R2(result)[1], error = function(e) NA)
    waic <- tryCatch(waic(result)$estimates["waic","Estimate"], error = function(e) NA)
    looic <- tryCatch(loo(result)$estimates["looic", "Estimate"], error = function(e) NA)

    n_param <- length(get_variables(result))
    
    n_id <- tryCatch(nrow(unique(result$data["id"])), error = function(e) NA)
   
    ms <- summary(result)
    n_observe <- tryCatch(ms$nobs, error = function(e) NA)
    n_chains <- ms$chains
    n_iter <- ms$iter
    n_diverge <- rstan::get_num_divergent(result$fit)
    adapt_delta <- tryCatch(control_params(result)[["adapt_delta"]], error = function(e) NA)
    max_treedepth <- tryCatch(control_params(result)[["max_treedepth"]], error = function(e) NA)

    elapsed_time <- rstan::get_elapsed_time(result$fit)
    total_time <- sum(elapsed_time[, "warmup"]) + sum(elapsed_time[, "sample"])
    avgChainMin <- (total_time / n_chains) / 60
    
    df <- data.frame(file_name, md, model_formula, model_dv, m_fixed, m_rand, m_family, bayesR2, waic, looic, n_param, n_observe, n_id,
                     n_diverge, n_chains, n_iter, avgChainMin, adapt_delta, max_treedepth) %>%
      mutate(across(where(is.numeric), \(x) round(x, 3)))

    return(df)
  }
}


#saveRDS(model_inventory, here::here("data", "model_inventory.rds"))

# Path to model_inventory.rds
model_inventory_path <- here::here("data", "model_inventory.rds")

# Load existing model_inventory if it exists
if (file.exists(model_inventory_path)) {
  model_inventory <- readRDS(model_inventory_path)
} else {
  model_inventory <- data.frame()
}

# List model files
file_list <- list.files(here::here("data/model_cache"), full.names = TRUE, recursive = FALSE)
file_list <- file_list[file.info(file_list)$isdir == FALSE]
file_list <- file_list[order(file.info(file_list)$mtime, decreasing = TRUE)] %>% head(55)


existing_files <- basename(file_list)
new_files_base <- base::setdiff(existing_files, model_inventory$file_name)
new_files <- file_list[existing_files %in% new_files_base]
new_model_inventory <- map_dfr(new_files, ~ extract_info(.x))

model_inventory_updated <- bind_rows(model_inventory, new_model_inventory)
# Check if changes occurred, and if so, write updated model_inventory
if (nrow(model_inventory) != nrow(model_inventory_updated)) {
  saveRDS(model_inventory_updated, model_inventory_path)
}




######### Inspecting/comparing models

model_inventory <- readRDS(here::here("data", "model_inventory.rds"))
vxModels <- model_inventory |> filter(model_dv=="vx") |> arrange(desc(bayesR2))

distModels <- model_inventory |> filter(model_dv=="dist") |> 
  arrange(desc(bayesR2)) |> # filter for models that have bandInt in the forumula
  filter(grepl("bandInt", model_formula)) 






#################


###### OLD




# Handle any error files if needed
error_files <- setdiff(file_list, model_inventory$file_name)
if (length(error_files) > 0) {
  message("Error reading the following files:")
  message(error_files)
} else {
  message("All files read successfully.")
}

# Print the model_inventory dataframe
print(model_inventory)
#This code will create the model_inventory dataframe with the columns "file_name," "file_modified_date," and "model_formula" as required.














list.files(here::here("data/model_cache"), full.names = TRUE) %>% 
  file.info() %>% 
  arrange(desc(mtime)) %>% head(3)



# modList <- list()
# for (i in seq_along(file_list)) {
#   modList[[i]] <- readRDS(file_list[i])
# }



# print the model formulas for each model
for (i in seq_along(modList)) {
  cat("\n\nModel", i, "formula:\n")
  print(modList[[i]]$formula$formula) # Adjust this line based on the structure of the object
}


tmod <- c()
for (i in seq_along(modList)) {
  formula_as_string <- deparse(modList[[i]]$formula$formula)
  if (grepl("trunc", formula_as_string)) {
    tmod <- c(tmod, paste0("Model ", i))
  }
}

# identify models which have "dist" as the dependent variable
dmod <- c()
for (i in seq_along(modList)) {
  formula_as_string <- deparse(modList[[i]]$formula$formula)
  if (grepl("dist", formula_as_string)) {
    dmod <- c(dmod, paste0("Model ", i))
  }
}

# identify models which have both vb and condit as predictors
# vmod <- c()
# for (i in seq_along(modList)) {
#   formula_as_string <- deparse(modList[[i]]$formula$formula)
#   if (grepl("vb", formula_as_string) & grepl("condit", formula_as_string)) {
#     vmod <- c(vmod, paste0("Model ", i))
#   }
# }

vmod <- map(seq_along(modList), ~ {
  formula_as_string <- deparse(modList[[.x]]$formula$formula)
  if (grepl("vb", formula_as_string) & grepl("condit", formula_as_string)) {
    return(paste0("Model ", .x))
  } else {
    return(NULL)
  }
}) %>% compact() %>% unlist()



# identify models which have bandInt as a predictor
bmod <- c()
for (i in seq_along(modList)) {
  formula_as_string <- deparse(modList[[i]]$formula$formula)
  if (grepl("bandInt", formula_as_string)) {
    bmod <- c(bmod, paste0("Model ", i))
  }
}

bmod <- map(seq_along(modList), ~ {
  formula_as_string <- deparse(modList[[.x]]$formula$formula)
  if (grepl("bandInt", formula_as_string)) {
    return(paste0("Model ", .x))
  } else {
    return(NULL)
  }
}) %>% compact() %>% unlist()


# identify models with gaussian family
gmod <- map(seq_along(modList), ~ {
  if (!is.null(modList[[.x]]$family) && 
      !is.null(modList[[.x]]$family$family) && 
      modList[[.x]]$family$family == "gaussian") {
    return(paste0("Model ", .x))
  } else {
    return(NULL)
  }
}) %>% compact() %>% unlist |> unique()



# identify models with weibull family
wmod <- c()
for (i in seq_along(modList)) {
  if (!is.null(modList[[i]]$family) && 
      !is.null(modList[[i]]$family$family) && 
      modList[[i]]$family$family == "weibull") {
    wmod <- c(wmod, paste0("Model ", i))
  }
}




# print summaries for all models with bandInt as predictor
for (i in seq_along(modList)) {
  formula_as_string <- deparse(modList[[i]]$formula$formula)
  if (grepl("bandInt", formula_as_string)) {
    cat("\n\nModel", i, "summary:\n")
    print(summary(modList[[i]]))
  }
}

# print summaries for all models with vb and condit as predictors, but which don't have truncation
for (i in seq_along(modList)) {
  formula_as_string <- deparse(modList[[i]]$formula$formula)
  if (grepl("vb", formula_as_string) & grepl("condit", formula_as_string) & !grepl("trunc", formula_as_string)) {
    cat("\n\nModel", i, "summary:\n")
    print(summary(modList[[i]]))
  }
}

