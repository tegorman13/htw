


pacman::p_load(tidyverse, data.table, here, patchwork, conflicted)
conflict_prefer_all("dplyr", quiet = TRUE)


walk(c("fun_alm","fun_model", "Display_Functions"), ~ source(here::here(paste0("Functions/", .x, ".R"))))
ds <- readRDS(here::here("data/e1_md_11-06-23.rds"))  |> as.data.table()
dsv <- ds |> filter(condit=="Varied")  
dsc <- ds |> filter(condit=="Constant") 
ind_ds <- ds |> 
  filter(expMode2 == "Test") |> 
  group_by(id, condit, x, expMode2) |> 
  summarise(y = mean(y), .groups = "keep") 



tMax=84
avg_dsv <- ds |> filter(condit=="Varied",expMode2=="Train") |> group_by(tr) %>%
  mutate(bandInt2 = sample(rep(c(800, 1000, 1200), each = tMax / 3), tMax, replace = FALSE)[tr]) %>%
  filter(bandInt2 == x) |> select(-bandInt2) |> group_by(tr,condit,x,expMode2) |> summarise(y=mean(y),.groups="keep") |>
  rbind(dsv |> filter(expMode2=="Test") |> group_by(condit,x,expMode2) |> summarise(y=mean(y),tr=1,.groups="keep") ) |> setDT()

avg_dsc <- ds |> filter(condit=="Constant",expMode2=="Train",tr<=tMax) |> group_by(tr, condit,x,expMode2) |> 
  summarise(y=mean(y),.groups="keep") |> 
  rbind(dsc |> filter(expMode2=="Test") |> 
          group_by(condit,x,expMode2) |> 
          summarise(y=mean(y),tr=1,.groups="keep") ) |> 
  setDT()

input_layer =  c(100,350,600,800,1000,1200)
output_layer = input_layer

group_posterior_all <- readRDS(here::here("data/abc_2M_rmse_p001.rds"))  

Model_Names <- list("EXAM", "Alt_EXAM", "ALM")
Fit_Methods <- list("Test & Train", "Test Only","Train Only")
condits <- list("Constant", "Varied")


# posteriors for each Fit_Method 
teter <- group_posterior_all |> map_dfr(~tibble(pluck(.x$teter_results))) 
te <- group_posterior_all |> map_dfr(~tibble(pluck(.x$te_results))) 
tr <- group_posterior_all |> map_dfr(~tibble(pluck(.x$tr_results))) 

# posterior distribution from fit to entire varied group
gp_teter_v <- teter |>  filter(Group=="Varied",Model=="EXAM") |> select(c,lr)
head(gp_teter_v)

# individual data
dsv |> filter(condit=="Varied") |> head()
dsv %>% pull(id) %>% unique() %>% length()


kde_result <- MASS::kde2d(gp_teter_v$c, gp_teter_v$lr, n = 100)
contour(kde_result)
image(kde_result)
persp(kde_result)





density_c <- apply(kde_result$z, 1, sum)
density_c <- density_c / sum(density_c)

density_lr <- apply(kde_result$z, 2, sum)
density_lr <- density_lr / sum(density_lr)

plot(kde_result$x, density_c, type = 'l', xlab = 'c', ylab = 'Density')
plot(kde_result$y, density_lr, type = 'l', xlab = 'lr', ylab = 'Density')
ggplot(data.frame(x=gp_teter_v$c), aes(x=x)) +
      geom_density() + geom_rug()+
ggplot(data.frame(x=gp_teter_v$lr), aes(x=x)) +
      geom_density() +geom_rug()



group_posterior_all <- readRDS(here::here("data/abc_2M_rmse_p001.rds"))  

Model_Names <- list("EXAM", "Alt_EXAM", "ALM")
Fit_Methods <- list("Test & Train", "Test Only","Train Only")
condits <- list("Constant", "Varied")

sample_from_kde <- function(kde_result, n_samples) {
  # Flatten the z matrix and sample indices based on these probabilities
  sampled_indices <- sample(length(kde_result$z), size = n_samples, replace = TRUE, prob = c(kde_result$z))
  
  # Convert indices to matrix row and column numbers
  rows <- ((sampled_indices - 1) %% nrow(kde_result$z)) + 1
  cols <- ((sampled_indices - 1) %/% nrow(kde_result$z)) + 1
  
  # Map back to x and y values and create a tibble
  samples <- tibble(c = kde_result$x[rows], lr = kde_result$y[cols])
  return(samples)
}

compute_kde <- function(c, lr, ngrid = 100, nsamples = 500) {
  kde=MASS::kde2d(c, lr, n = ngrid)
  kde_samples=sample_from_kde(kde, nsamples)
  return(list(kde=kde, kde_samples=kde_samples))
}

kde_results <- purrr::map(group_posterior_all, ~purrr::map_if(.x, is.list, ~compute_kde(.x$c, .x$lr, ngrid = 100, nsamples=500)))








compute_kde <- function(c, lr, ngrid = 100) {
  MASS::kde2d(c, lr, n = ngrid)
}

kde_results <- purrr::map(group_posterior_all, ~purrr::map_if(.x, is.list, ~compute_kde(.x$c, .x$lr, ngrid = 100)))







# Write a function to  Obtain kde for each combo of Model, Fit_Method, and condit. Use MASS::kde2d
# 
# group_posterior_all has 6 lists in its first level ("abc_ev"   "abc_almv" "abc_altv" "abc_ec"   "abc_almc" "abc_altc"), 1 for each combo of condit (c=Constant, v=Varied) and Model (e=EXAM, alm=ALM, alt=Alt_Exam). Each of those 6 lists contains 3 further lists for each Fit_Method (teter_results, te_results, tr_results), the fit method lists contain the posterior parameters for c and lr. 
# 
# Use purrr functions to map through each of the first and 2nd level list items. 
# 
# pacman::p_load(dplyr,purrr,tidyr,ggplot2, data.table, here, patchwork, conflicted)
# conflict_prefer_all("dplyr", quiet = TRUE)
# 
# 
# # load group level abc posterior - posterior parameters and predicted data for each combo of condit, Fit_Method, and Model
# group_posterior_all <- readRDS(here::here("data/abc_2M_rmse_p001.rds"))  
# 
# # In Model names, e is short hand for EXAM,alm for ALM, and alt for Alt_EXAM, c is for constant, v is for varied, teter is for test & train, te is for test only, tr is for train only
# Model_Names <- list("EXAM", "Alt_EXAM", "ALM")
# Fit_Methods <- list("Test & Train", "Test Only","Train Only")
# condits <- list("Constant", "Varied")
# 
# ## Example of getting to group level posterior data for one of the combos: 
# names(group_posterior_all)
# [1] "abc_ev"   "abc_almv" "abc_altv" "abc_ec"   "abc_almc" "abc_altc"
# > names(group_posterior_all$abc_ev)
# [1] "teter_results" "te_results"    "tr_results"    "Model"         "Group"         "pct_keep"      "fn_name"      
# > group_posterior_all$abc_ev$te_results |> select(c,lr)
# # A tibble: 2,000 × 2
# # Rowwise: 
# c    lr
# <dbl> <dbl>
#   1 0.268   2.43
# 2 0.268   2.42


###
hist(kde_result$x, prob=TRUE)
rug(kde_result$x)


hist(gp_teter_v$c, prob=TRUE)
rug(gp_teter_v$c)

kern_dens <- function(x, h, m = 512) {
  rg <- range(x)
  # xx is equivalent to grid points in 'density()'
  xx <- seq(rg[1] - 3 * h, rg[2] + 3 * h, length.out = m)
  y <- numeric(m) # The evaluations, initialized as a vector of zeros
  # The actual computation is done using nested for-loops. The outer loop
  # is over the grid points, and the inner loop is over the data points.
  for (i in seq_along(xx))
    for (j in seq_along(x))
      y[i] <- y[i] + exp(- (xx[i] - x[j])^2 / (2 * h^2))  
  y <- y / (sqrt(2 * pi) * h * length(x))
  list(x = xx, y = y)
}
f_hat <- kern_dens(kde_result$x, 0.00005, m=1024)
plot(f_hat, type = "l", lwd = 4, xlab = "x", ylab = "Density")