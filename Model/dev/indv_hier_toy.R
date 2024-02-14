

pacman::p_load(dplyr,purrr,tidyr,ggplot2, data.table, here, patchwork, conflicted)
conflict_prefer_all("dplyr", quiet = TRUE)
walk(c("fun_alm","fun_model"), ~ source(here::here(paste0("Functions/", .x, ".R"))))
set.seed(123)

ds <- readRDS(here::here("data/e1_md_11-06-23.rds"))  |> as.data.table()
dsv <- ds |> filter(condit=="Varied")  
dsc <- ds |> filter(condit=="Constant") 






samp_priors <- function(n,cMean=-5,cSig=2,lrSig=1) {
  prior_samples <- tibble(
    c = rlnorm(n,cMean,sdlog=cSig),
    lr = extraDistr::rhnorm(n,sigma=lrSig),
  )
  return(prior_samples)
}

cMean <<- -5.5; cSig <<- 2.0; lrSig <<- 2.0
prior_samples <- samp_priors(n=300000, cMean=cMean, cSig=cSig, lrSig=lrSig) 
return_dat="test_data"
input_layer =  c(100,350,600,800,1000,1200)
output_layer = input_layer
simulation_function = full_sim_exam

n_particle = 1000
ntry=300; 
ids1 <- 1
ids1 <- c(1,33,66)
sbjList <- list()

for(i in 1:length(ids1)) {
    abc = list()
    abc$theta <- vector("list", n_particle)
    abc$dist_sd <-  vector("list", n_particle)
    data <- ds[id==ids1[i],]
    target_data <- case_when(
    return_dat == "test_data" ~ list(target_data_test <- data[data$expMode2 == "Test", ]), 
    return_dat == "train_data" ~ list( target_data_train <- data[data$expMode2 == "Train", ]),
    return_dat == "train_data, test_data" ~ list( target_data_train_test <- data[expMode2 %in% c("Test", "Train"), ])
    ) |> pluck(1)

    tol = target_data |> group_by(x) |> summarise(m=mean(y),sd=sd(y)) |> summarise(tol=mean(sd)) *.9
  for(j in 1:n_particle) {
    #print(i)
    print(paste0("sbj: ",data$id[1] ,"Particle:", j))


    found=0;
    try_count=0;
    while(found==0) {
    try_count=try_count+1;
    current_theta <- prior_samples[sample(1:nrow(prior_samples), 1), ]
    sim_data <- simulation_function(data, current_theta$c, current_theta$lr, input_layer = input_layer, output_layer = output_layer, return_dat = return_dat) 
    dist_sd <- target_data |> mutate(pred=sim_data,error=abs(y-pred)) |> group_by(id,x) |> 
      summarise(mean_error=mean(error),sd_error=sd(error),.groups="keep") |> 
      group_by(id) |> 
      summarise(mean_error=mean(mean_error),sd_error=mean(sd_error)) 
    if(dist_sd$mean_error< tol) {
      #print(current_theta)
      abc$theta[[j]] <- current_theta 
      abc$dist_sd[[j]] <- cbind(current_theta,dist_sd)
      found=1
      try_count=0;
      
      }
      if (try_count>=ntry){
        print(paste0("increase tol for subject", data$id[1]))
        tol=tol*1.1
        try_count=0;
      }
    }

  }
  sbjList[[i]]<- abc
}


map(sbjList, ~{.x$dist_sd |> rbindlist() |> arrange(mean_error) |> head(5)})



sbjList


target_data |> mutate(pred=sim_data,error=abs(y-pred)) |> group_by(x) |> summarise(mean_error=mean(error),sd_error=sd(error), sdY=sd(y)) |> 
  summarise(mean_error=mean(mean_error),sd_error=mean(sd_error), sdY=mean(sdY))



J = 4
T = 100
lambda = 5 # This might represent the mean of the normal distribution for simplicity
sigma = 1 # Assume a known standard deviation for the normal distribution subject-level parameters
mu = rnorm(J, lambda, sigma) # generate means for each subject

Y = matrix(NA, J, T)
for(j in 1:J) {
  Y[j,] = rnorm(T, mu[j], sigma) # generate subject data from a normal distribution
}

eps = 1 # tolerance threshold
N = 1000 # number of particles
p.mu = c(0, 1) # prior settings for mu (normally distributed)
rho = function(x, y) mean(abs(x-y)) # rho function

abc = list()
abc$mu = matrix(NA, N, J) # declare a matrix for storage
abc$sigma = numeric(N) # declare a vector for storage (assuming sigma is known and fixed)

for(i in 1:N) { # loop over particles
  print(paste0("Particle:", i))
  for(j in 1:J) { # loop over subjects
    d = eps + 1 # initialized to be greater than eps

    # continue proposal generation until condition is satisfied
    while(d > eps) {
      mu.1 = rnorm(1, mean(Y[j,]), sigma) # sample from proposal distribution
      x = rnorm(T, mu.1, sigma) # simulate data
      print(d)
      d = ifelse(mu.1 > 0, rho(Y[j,], x), eps + 1) # compute distance
    }
    abc$mu[i,j] = mu.1 # store the accepted value
  }

  # sample from conditional distribution of lambda (mean), assuming a normal prior for mu
  abc$lambda[i] = rnorm(1, mean(abc$mu[i,]), sqrt(sigma^2 / J))
}

# Note that sigma is assumed to be known and fixed, so it's not updated in the ABC algorithm.


colMeans(abc$mu) 
# median of posterior theta's
apply(abc$mu, 2, median)
#true mu
mu

# mean of posterior lambda
mean(abc$lambda)
# median of posterior lambda
median(abc$lambda)

#true lambda
lambda
