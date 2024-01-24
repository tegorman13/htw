# https://www.digitalocean.com/community/tutorials/how-to-configure-ssh-key-based-authentication-on-a-linux-server#step-2-copying-an-ssh-public-key-to-your-server
# https://kb.iu.edu/d/aews
#https://www.mm218.dev/posts/2023-03-03-parallel/
#https://win-vector.com/2016/01/22/running-r-jobs-quickly-on-many-machines/
#m1r1 10.147.17.45
# m1l2 10.147.17.69

#m1r2 10.147.17.211

# tg_m1 10.147.17.202



library(parallel)

cl <- parallel::makePSOCKcluster(
  c("10.147.17.45", "10.147.17.202"),
  master = "10.147.17.202",
  user = "thomasgorman"
)
future::plan(future::cluster, workers = cl)


primary <- '10.147.17.202'

machineAddresses <- list(
  list(host=primary,user='thomasgorman',
       ncore=4),
  list(host='10.147.17.45',user='thomasgorman',
       ncore=4)
)

spec <- lapply(machineAddresses,
               function(machine) {
                 rep(list(list(host=machine$host,
                               user=machine$user)),
                     machine$ncore)
               })
spec <- unlist(spec,recursive=FALSE)

parallelCluster <- parallel::makeCluster(type='PSOCK',
                                         master=primary,
                                         spec=spec)
print(parallelCluster)
#print(parallelCluster)
# socket cluster with 8 nodes on hosts ‘10.147.17.202’, ‘10.147.17.45’



machineAddresses <- list(
  list(host=primary,user='thomasgorman',
       ncore=4),
  list(host='10.147.17.45',user='thomasgorman',
       ncore=4), 
  list(host='10.147.17.69',user='thomasgorman',
       ncore=4)
)
spec <- lapply(machineAddresses,
               function(machine) {
                 rep(list(list(host=machine$host,
                               user=machine$user)),
                     machine$ncore)
               })
spec <- unlist(spec,recursive=FALSE)
parallelCluster <- parallel::makeCluster(type='PSOCK',
                                         master=primary,
                                         spec=spec)
print(parallelCluster)
# socket cluster with 12 nodes on hosts ‘10.147.17.202’, ‘10.147.17.45’, ‘10.147.17.69’




# https://mvuorre.github.io/posts/parallel-multiverse/
generate_data <- function(seed = NA, n = 1e5) {
  if (!is.na(seed)) set.seed(seed)
  dat <- tibble(
    x1 = rnorm(n),
    x2 = rnorm(n),
    y1 = rnorm(n) + x1*.1,
    y2 = rnorm(n) + x1*.2,
    c1 = rnorm(n) + x1*.3,
    c2 = rnorm(n),
    group = sample(c("a", "b", "c", "d"), n, replace = TRUE)
  )
}
dat <- generate_data(9, n=1e6)

library(specr)
library(kableExtra)
library(scales)
library(ggthemes)
library(tictoc)
library(tidyverse)
library(furrr)
library(broom)


expand_covariate <- function(covariate) {
  list(
    "1",
    do.call(
      "c",
      map(
        seq_along(covariate), 
        ~combn(covariate, .x, FUN = list))
    ) %>%
      map(~paste(.x, collapse = " + "))
  ) %>%
    unlist
}

lm2 <- function(formula, data) {
  fit <- lm(formula = formula, data = data)
  out <- tidy(fit, conf.int = TRUE) # Tidy table of parameters
  out <- slice(out, 2) # Second row (slope parameter)
  bind_cols(out, n = nobs(fit))
}
glm2 <- function(formula, data) {
  fit <- glm(formula = formula, data = data)
  out <- tidy(fit, conf.int = TRUE)
  out <- slice(out, 2)
  bind_cols(out, n = nobs(fit))
}


specs <- expand_grid(
  x = c("x1", "x2"),
  y = c("y1", "y2"),
  covariate = expand_covariate(c("c1", "c2")),
  model = c("lm", "glm"),
  distinct(dat, group)
) %>% 
  mutate(formula = str_glue("{y} ~ {x} + {covariate}"))

specs <- mutate(specs, model = paste0(model, "2"))


tic()
results_dplyr <- specs %>% 
  mutate(
    out = pmap(
      list(model, formula, group), 
      ~do.call(
        ..1, 
        list(
          formula = ..2, 
          data = filter(dat, group == ..3)
        )
      )
    )
  )
toc()
# 13.293 sec elapsed - 1e5
# 171.405 sec elapsed 1e6


future::plan(future::cluster, workers = parallelCluster)
plan(multisession, workers = 2)
opts <- furrr_options(
  globals = list(dat = dat, lm2 = lm2, glm2 = glm2),
  packages = c("dplyr", "broom")
)

tic()
results_furrr <- specs %>% 
  mutate(
    out = future_pmap(
      list(model, formula, group), 
      ~do.call(
        what = ..1, 
        args = list(
          formula = ..2, 
          data = filter(dat, group == ..3)
        )
      ), 
      .options = opts
    )
  ) %>% 
  unnest(out)
toc()
#p2 18.137 sec elapsed 1e5
#p2 125.214 sec elapsed 1e6


# Shutdown cluster neatly
if(!is.null(parallelCluster)) {
  parallel::stopCluster(parallelCluster)
  parallelCluster <- c()
}
