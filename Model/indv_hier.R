

pacman::p_load(dplyr,purrr,tidyr,ggplot2, data.table, here, patchwork, conflicted)
conflict_prefer_all("dplyr", quiet = TRUE)
walk(c("fun_alm","fun_model"), ~ source(here::here(paste0("Functions/", .x, ".R"))))
set.seed(123)

ds <- readRDS(here::here("data/e1_md_11-06-23.rds"))  |> as.data.table()
dsv <- ds |> filter(condit=="Varied")  
dsc <- ds |> filter(condit=="Constant") 

ids <- c(1,2,4,5,6,7,8, 10,11,12,13)
ids2 <- c(1,66,36)


group_posterior_all <- readRDS(here::here("data/abc_2M_rmse_p001.rds"))  

teter <- group_posterior_all |> map_dfr(~tibble(pluck(.x$teter_results))) 
te <- group_posterior_all |> map_dfr(~tibble(pluck(.x$te_results))) 
tr <- group_posterior_all |> map_dfr(~tibble(pluck(.x$tr_results))) 


# Function to perform Gibbs sampling for the subject-level parameters
sample_subject_parameters <- function(y, hyperparameters, tolerance) {
  # Assuming hyperparameters is a list with mu and sigma squared
  mu <- hyperparameters$mu
  sigma_squared <- hyperparameters$sigma_squared
  
  # Initialize subject parameters
  subject_params <- exp(rnorm(length(y), mu, sqrt(sigma_squared)))
  
  # Perform rejection sampling for each subject
  for (i in seq_along(y)) {
    repeat {
      proposed_param <- exp(rnorm(1, mu, sqrt(sigma_squared)))
      simulated_data <- simulate_poisson_data(proposed_param, length(y[[i]]))
      if (rho(simulated_data, y[[i]]) < tolerance) {
        subject_params[i] <- proposed_param
        break
      }
    }
  }
  
  subject_params
}

# Function to sample new hyperparameters from their conditional distributions
sample_hyperparameters <- function(subject_params, hyper_prior_params) {
  # Assuming hyper_prior_params is a list with mu and sigma squared for the Gaussian prior
  mu_prior <- hyper_prior_params$mu
  sigma_squared_prior <- hyper_prior_params$sigma_squared
  
  # Log-transform subject parameters
  log_params <- log(subject_params)
  
  # Calculate the posterior parameters for mu
  n <- length(log_params)
  sigma_squared_post <- 1 / (1 / sigma_squared_prior + n / hyper_prior_params$sigma_squared)
  mu_post <- sigma_squared_post * (mu_prior / sigma_squared_prior + sum(log_params) / hyper_prior_params$sigma_squared)
  
  # Sample new mu from the normal distribution
  mu_sampled <- rnorm(1, mean = mu_post, sd = sqrt(sigma_squared_post))
  
  # Calculate the posterior parameters for sigma squared
  alpha_post <- hyper_prior_params$alpha + n / 2
  beta_post <- hyper_prior_params$beta + 0.5 * sum((log_params - mu_sampled)^2)
  
  # Sample new sigma squared from the inverse gamma distribution
  sigma_squared_sampled <- 1 / rgamma(1, shape = alpha_post, rate = beta_post)
  
  list(mu = mu_sampled, sigma_squared = sigma_squared_sampled)
}


# Main Gibbs ABC function for individual fits
fit_indv_gibbs_abc <- function(sbj_id, simulation_function, hyper_prior_params, tolerance, iterations, Model, Group) {
  data <- filter(ds, id == sbj_id) |> as.data.table()
  target_data <- data$y
  
  # Initialize hyperparameters
  hyper_samples <- list(mu = hyper_prior_params$mu, sigma_squared = hyper_prior_params$sigma_squared)
  
  # Initialize subject-level parameters
  subject_samples <- lapply(target_data, mean) # Just a simple initialization
  
  for (i in 2:iterations) {
    # Sample hyperparameters given subject-level parameters
    hyper_samples <- sample_hyperparameters(unlist(subject_samples), hyper_prior_params)
    
    # Sample subject-level parameters given hyperparameters
    subject_samples <- sample_subject_parameters(target_data, hyper_samples, tolerance, simulation_function)
  }
  
  list(hyper_samples = hyper_samples, subject_samples = subject_samples, Model = Model, Group = Group)
}

# ... [rest of the code] ...

# Example usage for a single subject
hyper_prior_params <- list(mu = log(0.5), sigma_squared = 1) # Example hyperprior parameters
tolerance <- 0.01
iterations <- 100

# Fit the EXAM model using hierarchical Gibbs ABC for a single subject
id_fit_exam_varied <- fit_indv_gibbs_abc(sbj_id = 1,
                                          simulation_function = full_sim_exam, 
                                          hyper_prior_params = hyper_prior_params, 
                                          tolerance = tolerance, 
                                          iterations = iterations, 
                                          Model = "EXAM", 
                                          Group = "Varied")

# ... [rest of the code] ...




J = 2
T = 100
theta = rexp(J,lambda) # generate parameters for each subject
lambda = .3
Y = matrix(NA,J,T)
for(j in 1:J) { Y[j,] = rpois(T,theta[j])} #generate subject data

eps=1e-10# tolerance threshold
N=50 # number of particles
p.lambda=c(.01,.01) #prior settings(Gamma distribution)
rho=function(x,y)abs(sum(x)-sum(y))/T #rho function

abc=NULL
abc$theta=matrix(NA,N,J) # declare a matrix for storage
abc$lambda=numeric(N) # declare a matrix for storage

for(i in 1:N){ #loop over particles

for(j in 1:J){ #loop over subjects

d=eps+1 #initialized to be greater than eps

#continue proposal generation until condition is satisfied
while(d>eps){
  print(paste0("Subject ",j, "Particle:", i))

theta.1=rnorm(1,mean(Y[j,1]),1)#sample from proposal dist.

x=rpois(T,theta.1) # simulate data
d=ifelse(theta.1>0,rho(Y[j,1],x),eps+1) #compute distance
print(d)
}
abc$theta[i,j]=theta.1 #store the accepted value
}

#sample from conditional distribution of lambda
abc$lambda[i]=rgamma(1,J+p.lambda[1],sum(abc$theta[i,])+p.lambda[2])
} 

abc$lambda
abc$theta