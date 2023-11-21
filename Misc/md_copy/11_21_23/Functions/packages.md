#### INSTALL AND LOAD PACKAGES ==========================================================
# install pacman package if not installed -----------------------------------------------
suppressWarnings(if (!require("pacman")) install.packages("pacman"))
#remotes::install_github("mattcowgill/ggannotate")
# load packages and install if not installed --------------------------------------------
pacman::p_load(tidyverse,knitr,conflicted, tidybayes,brms, lme4, bayesplot,bayestestR,parameters,marginaleffects,
                emmeans, equatiomatic, here, pacman,  broom,
               broom.mixed,lme4,emmeans,gt,
                wesanderson,glue, ggdist,ggforce,patchwork,gghalves,ggh4x,
                install = TRUE,
                update = FALSE
               )

               
# show loaded packages ------------------------------------------------------------------
# cat("loaded packages\n")
# print(pacman::p_loaded())
conflict_prefer_all("dplyr", quiet = TRUE)
# select <- dplyr::select
# mutate <- dplyr::mutate
# filter <- dplyr::filter
# map <- purrr::map
# walk(c(here("Functions/Display_Functions.R"), here("Functions/org_functions.R"), 
#        here("Functions/Table_Functions.R")), source)


# purrr::walk(c("con_group_exam_fits", "var_group_exam_fits", "hybrid_group_exam_fits"), 
#             ~ list2env(readRDS(here::here(paste0("data/model_cache/", .x, ".rds"))), 
#             envir = .GlobalEnv))

# load function scripts with succinct purrr
walk(c("Display_Functions", "org_functions", "Table_Functions"), ~ source(here::here(paste0("Functions/", .x, ".R"))))