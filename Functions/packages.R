#### INSTALL AND LOAD PACKAGES ==========================================================
# install pacman package if not installed -----------------------------------------------
suppressWarnings(if (!require("pacman")) install.packages("pacman"))
#remotes::install_github("mattcowgill/ggannotate")
# load packages and install if not installed --------------------------------------------
pacman::p_load(tidyverse,tidybayes,brms, lme4, bayesplot,bayestestR,parameters,marginaleffects,
                emmeans, equatiomatic, here, pacman,  broom.mixed,
               broom.mixed,lme4,emmeans,here,knitr,kableExtra,gt,
                wesanderson,glue, ggdist,ggforce,patchwork,gghalves,
                install = TRUE,
                update = FALSE
               )

# show loaded packages ------------------------------------------------------------------
cat("loaded packages\n")
print(pacman::p_loaded())

select <- dplyr::select
mutate <- dplyr::mutate
filter <- dplyr::filter
map <- purrr::map
walk(c(here("Functions/Display_Functions.R"), here("Functions/org_functions.R")), source)