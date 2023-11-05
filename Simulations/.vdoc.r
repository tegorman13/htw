#
#
#
#
#
#
#
#
#
#
#
#
#
#
# load and view data
pacman::p_load(tidyverse,lme4,future,furrr,patchwork,here)
purrr::walk(here(c("Functions/Display_Functions.R", "Functions/alm_core.R")),source)
select <- dplyr::select; mutate <- dplyr::mutate 

d <- readRDS(here("data/e1_08-21-23.rds"))
levels(d$condit)
dtest <- d %>% filter(expMode %in% c("test-Nf","test-train-nf")) %>% group_by(id,lowBound) %>% 
  mutate(nBand=n(),band=bandInt,id=factor(id)) %>% group_by(id) %>% mutate(nd=n_distinct(lowBound))
dtest <- dtest %>% group_by(id,lowBound) %>% filter(nBand>=5 & nd==6)
dtest <- dtest %>% group_by(id) %>% filter(!id %in% unique(dtest$id[dtest$nBand<5]))

ds <- d %>% filter(expMode2 %in% c("Train","Test")) |> 
  filter(!id %in% unique(dtest$id[dtest$nBand<5])) |>   
  group_by(id,condit,expMode2) |> 
  # create tr as trial number within id, condit, expMode2, 
  mutate(input=bandInt,x=bandInt, y=vx, tr= row_number()) |>
  select(id,condit,expMode2,tr,x,y) 
  #select(id,condit,expMode2,trial,gt.train,vb,band,bandInt,input,vx,dist,vxb) 
  
# dst <- ds %>% filter(expMode2=="Train")
# dst <- dst %>%
#   group_by(id, vb) %>%
#   mutate(trial_band = row_number())

# create dsAvg, averaging over id, within tr
dsAvg <- ds |> group_by(condit,expMode2,tr, x) |> 
  summarise(y=mean(y),.groups="keep") 

head(dsAvg) |> pander::pandoc.table(style = "simple")

#
#
#
#

inNodes <- 7; outNodes = 14; 
input.layer <- seq(50,1600,length.out=inNodes*1)  
output.layer <- seq(50,1600,length.out=outNodes*1) 

input.layer <- c(100,350,600,800,1000,1200)
output.layer <- c(100,350,600,800,1000,1200)

weights = matrix(.0000,nrow=length(output.layer),ncol=length(input.layer))


adjust_layer <- function(input.layer, k){
  if(k == 1) return(input.layer)
  new.layer <- c()
  for(i in 1:(length(input.layer) - 1)){
    new.layer <- c(new.layer, seq(input.layer[i], input.layer[i+1], length.out = k + 1))
  }
  new.layer <- unique(new.layer)
  return(new.layer)
}

default.layer <- c(100,350,600,800,1000,1200)
k = 2
input.layer <- adjust_layer(default.layer, k)
output.layer <- adjust_layer(default.layer, k)


#
#
#
#

optimize_params <- function(dsAvg, c_values, lr_values, input.layer, output.layer){
    best_error <- Inf
    best_params <- list(c = NA, lr = NA)
    
    # Split the data into training and testing sets
    train_data <- dsAvg[dsAvg$expMode2 == "Train", ]

    test_data <- dsAvg[dsAvg$expMode2 == "Test", ] |> group_by(x) |> summarise(y=mean(y))
    
    # Grid search over c and lr values
    for (c in c_values) {
        for (lr in lr_values) {
            # Train ALM with current parameters
            sim_results <- alm.sim(train_data, c, lr, input.layer, output.layer)
            weight.mat <- sim_results$wm
            trainVec <- c(0,unique(train_data$x))
            # Test on the testing data
            predictions <- sapply(test_data$x, function(input) {
                exam.response(input, c,trainVec, input.layer, output.layer, weight.mat)
            })
            (test_data <- cbind(test_data, pred=predictions))
            # Compute mean squared error
            mse <- mean((predictions - test_data$y)^2)
        
        # Update best parameters if current error is lower
            if (!is.na(mse) && mse < best_error) {
                best_error <- mse
                best_params$c <- c
                best_params$lr <- lr }
        }
    }
    print(best_params)
    print(best_error)
    
    return(best_params)
  }

# Example usage
c_values <- seq(0.000001, 0.01, length.out=40)
lr_values <- seq(0.01, 1.5, length.out=20)
best_paramsV <- optimize_params(dsAvg |> filter(condit=="Varied"), c_values, lr_values, input.layer, output.layer)
c=best_paramsV$c; lr=best_paramsV$lr



best_paramsC <- optimize_params(dsAvg |> filter(condit=="Constant"), c_values, lr_values, input.layer, output.layer)
c=best_paramsC$c; lr=best_paramsC$lr

#
#
#
#
#
#

optimize_params <- function(dsAvg, c_values, lr_values, input.layer, output.layer){
    best_error <- Inf
    best_params <- list(c = NA, lr = NA)
    
    # Split the data into training and testing sets
    train_data <- dsAvg[dsAvg$expMode2 == "Train", ]
    test_data <- dsAvg[dsAvg$expMode2 == "Test", ] |> group_by(x) |> summarise(y = mean(y))

    # Determine the second half of the training data
    train_data_half <- train_data[(round(nrow(train_data)/2) + 1):nrow(train_data), ]
    
    # Grid search over c and lr values
    for (c in c_values) {
        for (lr in lr_values) {
            # Train ALM with current parameters
            sim_results <- alm.sim(train_data, c, lr, input.layer, output.layer)
            weight.mat <- sim_results$wm
            trainVec <- c(unique(train_data$x))

            # Test on the testing data
            predictions_test <- sapply(test_data$x, function(input) {
                #exam.response(input, c, trainVec, input.layer, output.layer, weight.mat)
                exam.response2(input, c, train_data, input.layer, output.layer, weight.mat)

            })
           
            # Test on the latter half of the training data
            predictions_train <- sapply(train_data_half$x, function(input) {
                alm.response(input, c, input.layer, output.layer, weight.mat)$mean.response
            })

            # Compute mean squared error for test and latter half of training data
            mse_test <- mean((predictions_test - test_data$y)^2)
            mse_train <- mean((predictions_train - train_data_half$y)^2)

            # Use the average of the two MSEs as the final error
            final_error <- mean(c(mse_test, mse_train), na.rm = TRUE)

            # Update best parameters if current error is lower
            if (!is.na(final_error) && final_error < best_error) {
                best_error <- final_error
                best_params$c <- c
                best_params$lr <- lr
            }
        }
    }
     (test_data <- cbind(test_data, pred = predictions_test))
    print(best_params)
    print(best_error)
    
    return(best_params)
}


# Example usage
c_values <- seq(0.000001, 0.01, length.out=40)
lr_values <- seq(0.01, 1.5, length.out=20)
best_paramsV <- optimize_params(dsAvg |> filter(condit=="Varied"), c_values, lr_values, input.layer, output.layer)
c = best_paramsV$c; lr = best_paramsV$lr
viewPreds(dsAvg |> filter(condit=="Varied"), c_values, lr_values, input.layer, output.layer)


best_paramsC <- optimize_params(dsAvg |> filter(condit=="Constant"), c_values, lr_values, input.layer, output.layer)
c = best_paramsC$c; lr = best_paramsC$lr
viewPreds(dsAvg |> filter(condit=="Constant"), c_values, lr_values, input.layer, output.layer)


#
#
#
#
#
#
#
#

dsAvg <- ds |> group_by(condit,expMode2,tr, x) |> summarise(y=mean(y),.groups="keep") 
#dsAvg <- dsAvg |> filter(condit=="Varied")

optimize_params <- function(dsAvg, c_values, lr_values, input.layer, output.layer){
    best_error <- Inf
    best_params <- list(c = NA, lr = NA)
    
    # Split the data into training and testing sets
    train_data <- dsAvg[dsAvg$expMode2 == "Train", ]
    test_data <- dsAvg[dsAvg$expMode2 == "Test", ] |> group_by(x) |> summarise(y = mean(y))

    # Determine the second half of the training data
    train_data_half <- train_data[(round(nrow(train_data)/2) + 1):nrow(train_data), ]
    
    # Grid search over c and lr values
    for (c in c_values) {
        for (lr in lr_values) {
            # Train ALM with current parameters
            sim_results <- alm.sim(train_data, c, lr, input.layer, output.layer)
            weight.mat <- sim_results$wm
            trainVec <- c(unique(train_data$x))


            # predictions_train <- sim_results$d[(round(nrow(sim_results$d)/2) + 1):nrow(sim_results$d), ] |> mutate(error=(y-almResp))

             predictions_train <- sim_results$d[2:nrow(sim_results$d),] |> mutate(error=(y-almResp))


            # Test on the testing data
            predictions_test <- sapply(test_data$x, function(input) {
                #exam.response(input, c, trainVec, input.layer, output.layer, weight.mat)
                alm.response(input, c, input.layer, output.layer, weight.mat)$mean.response

            })

            # Compute mean squared error for test and latter half of training data
            mse_test <- mean((predictions_test - test_data$y)^2, na.rm=TRUE)
            mse_train <- mean(predictions_train$error^2)

            # Use the average of the two MSEs as the final error
            final_error <- mean(c( mse_train), na.rm = TRUE)

            # Update best parameters if current error is lower
            if (!is.na(final_error) && final_error < best_error) {
                best_error <- final_error
                best_params$c <- c
                best_params$lr <- lr
            }
        }
    }
     (test_data <- cbind(test_data, pred = predictions_test))
    print(best_params)
    print(best_error)
    
    return(best_params)
}


# Example usage
c_values <- seq(0.000001, 0.01, length.out=50)
lr_values <- seq(0.01, 1.5, length.out=30)
best_paramsV <- optimize_params(dsAvg |> filter(condit=="Varied"), c_values, lr_values, input.layer, output.layer)
c = best_paramsV$c; lr = best_paramsV$lr
viewPreds(dsAvg |> filter(condit=="Varied"), c, lr, input.layer, output.layer)


best_paramsC <- optimize_params(dsAvg |> filter(condit=="Constant"), c_values, lr_values, input.layer, output.layer)
c = best_paramsC$c; lr = best_paramsC$lr
viewPreds(dsAvg |> filter(condit=="Constant"), c, lr, input.layer, output.layer)




#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#

# ALM-only model prediction
predict_alm_only <- function(input, c, trainData,input.layer, output.layer, weight.mat) {
  alm.response(input, c, input.layer, output.layer, weight.mat)$mean.response
}

# ALM-EXAM model prediction
predict_alm_exam <- function(input, c, trainData, input.layer, output.layer, weight.mat) {
  exam.response2(input, c, trainData, input.layer, output.layer, weight.mat)
}

# ALM-EXAM hybrid model prediction
predict_alm_exam_hybrid <- function(input, c, trainData, input.layer, output.layer, weight.mat) {
  alm_pred <- predict_alm_only(input, c,trainData, input.layer, output.layer, weight.mat)
  exam_pred <- predict_alm_exam(input, c, trainData, input.layer, output.layer, weight.mat)
  
  # Average the predictions from ALM and EXAM
  (alm_pred + exam_pred) / 2
}


optimize_params <- function(dsAvg, c_values, lr_values, input.layer, output.layer, predict_function) {
    best_error <- Inf
    best_params <- list(c = NA, lr = NA)
    
    # Split the data into training and testing sets
    train_data <- dsAvg[dsAvg$expMode2 == "Train", ]
    test_data <- dsAvg[dsAvg$expMode2 == "Test", ] |> group_by(x) |> summarise(y=mean(y))
    
    # Grid search over c and lr values
    for (c in c_values) {
        for (lr in lr_values) {
            # Train ALM with current parameters
            sim_results <- alm.sim(train_data, c, lr, input.layer, output.layer)
            weight.mat <- sim_results$wm
            
            # Test on the testing data using the provided prediction function
            predictions <- sapply(test_data$x, function(input) {
                predict_function(input, c, train_data, input.layer, output.layer, weight.mat)
            })
            
            # Compute mean squared error
            mse <- mean((predictions - test_data$y)^2)
            
            # Update best parameters if current error is lower
            if (!is.na(mse) && mse < best_error) {
                best_error <- mse
                best_params$c <- c
                best_params$lr <- lr
            }
        }
    }
    
    list(best_params=best_params, best_error=best_error)
}


c_values <- seq(0.000001, 0.01, length.out=50)
lr_values <- seq(0.01, 1.5, length.out=50)

result_alm_only <- optimize_params(dsAvg |> filter(condit=="Varied"), c_values, lr_values, input.layer, output.layer, predict_alm_only)


result_alm_exam <- optimize_params(dsAvg |> filter(condit=="Varied"), c_values, lr_values, input.layer, output.layer, predict_alm_exam)


result_alm_exam_hybrid <- optimize_params(dsAvg |> filter(condit=="Varied"), c_values, lr_values, input.layer, output.layer, predict_alm_exam_hybrid)


map(list(result_alm_only, result_alm_exam, result_alm_exam_hybrid), 
    ~ viewPreds(dsAvg |> filter(condit=="Varied"), .x$best_params$c, .x$best_params$lr, input.layer,output.layer))


#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
predict_alm_exam_weighted_hybrid <- function(input, c, weight_exam, trainData, input.layer, output.layer, weight.mat) {

  alm_pred <- predict_alm_only(input, c, trainData, input.layer, output.layer, weight.mat)
    exam_pred <- predict_alm_exam(input, c, trainData, input.layer, output.layer, weight.mat)

  
  # Weighted combination of the predictions from ALM and EXAM
  (1 - weight_exam) * alm_pred + weight_exam * exam_pred
}

optimize_params_weighted <- function(dsAvg, c_values, lr_values, weight_exam_values, input.layer, output.layer) {
    best_error <- Inf
    best_params <- list(c = NA, lr = NA, weight_exam = NA)
    error_data <- data.frame() # Empty data frame to store errors

    # Split the data into training and testing sets
    train_data <- dsAvg[dsAvg$expMode2 == "Train", ]
    test_data <- dsAvg[dsAvg$expMode2 == "Test", ] |> group_by(x) |> summarise(y=mean(y))
    
    # Grid search over c, lr, and weight_exam values
    for (c in c_values) {
        for (lr in lr_values) {
            for (weight_exam in weight_exam_values) {
                # Train ALM with current parameters
                sim_results <- alm.sim(train_data, c, lr, input.layer, output.layer)
                weight.mat <- sim_results$wm
                
                # Test on the testing data using the weighted hybrid prediction function
                predictions <- sapply(test_data$x, function(input) {
                    predict_alm_exam_weighted_hybrid(input, c, weight_exam, train_data, input.layer, output.layer, weight.mat)
                })
                
                # Compute mean squared error
                mse <- mean((predictions - test_data$y)^2)
                error_data <- rbind(error_data, data.frame(c=c, lr=lr, weight_exam=weight_exam, mse=mse))
                
                # Update best parameters if current error is lower
                if (!is.na(mse) && mse < best_error) {
                    best_error <- mse
                    best_params$c <- c
                    best_params$lr <- lr
                    best_params$weight_exam <- weight_exam
                }
            }
        }
    }
    
    list(best_params=best_params, best_error=best_error, error_data=error_data)
}

# Example usage
weight_exam_values <- seq(.2, .8, length.out=6) # or any other desired sequence
c_values <- seq(0.0001, 0.01, length.out=40)
lr_values <- seq(0.01, 1.5, length.out=30)
result_weighted_hybrid <- optimize_params_weighted(dsAvg |> filter(condit=="Varied"), c_values, lr_values, weight_exam_values, input.layer, output.layer)
c = result_weighted_hybrid$best_params$c; lr = result_weighted_hybrid$best_params$lr; weight_exam = result_weighted_hybrid$best_params$weight_exam;
viewPreds(dsAvg |> filter(condit=="Varied"),c, lr, input.layer, output.layer, weight_exam)


error_data <- result_weighted_hybrid$error_data

# Visualize the error surface
p <- ggplot(error_data, aes(x=c, y=lr, z=mse)) +
    geom_line(aes(col=mse)) +
    scale_fill_gradient(low="blue", high="red") +
    theme_void() +
    # use viridis color palette
    scale_color_viridis_c(option="plasma") +
    labs(x="c", y="lr", z="mse") +
    guides(color=FALSE) +
    theme(legend.position="none")


plot_ly(error_data, x = ~c, y = ~lr, z = ~mse, type = "scatter3d", mode = "markers", color = ~mse, colors = c("blue", "red")) %>%
  layout(scene = list(xaxis = list(title = "c"),
                      yaxis = list(title = "lr"),
                      zaxis = list(title = "mse")))

#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
weight_exam_values <- seq(0, 1, length.out=10) # or any other desired sequence
c_values <- seq(0.00001, 0.01, length.out=50)
lr_values <- seq(0.00001, 1.5, length.out=40)
result_weighted_hybridC <- optimize_params_weighted(dsAvg |> filter(condit=="Constant"), c_values, lr_values, weight_exam_values, input.layer, output.layer)

best_params <- result_weighted_hybridC$best_params

c = best_params$c; lr = best_params$lr; weight_exam = best_params$weight_exam;
result_weighted_hybridC$best_error
viewPreds(dsAvg |> filter(condit=="Constant"),c, lr, input.layer, output.layer, weight_exam)



#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#


optimize_params_weighted_individual <- function(ds, c_values, lr_values, weight_exam_values, input.layer, output.layer) {
    all_results <- list()
    
    # Loop through each unique id
    for (individual in unique(ds$id)) {
        indiv_data <- ds[ds$id == individual, ]
        
        # Run the optimization function for the individual's data
        result <- optimize_params_weighted(indiv_data, c_values, lr_values, weight_exam_values, input.layer, output.layer)
        
        all_results[[as.character(individual)]] <- result
    }
    
    all_results
}

dss <- ds |> filter(id %in% c(1,2))


weight_exam_values <- seq(0, 1, length.out=15) # or any other desired sequence
c_values <- seq(0.000001, 0.5, length.out=50)
lr_values <- seq(0.000001, 1.5, length.out=50)

all_results_weighted_hybrid <- optimize_params_weighted_individual(ds, c_values, lr_values, weight_exam_values, input.layer, output.layer)

saveRDS(all_results_weighted_hybrid,file=paste0(here::here("data/model_cache/indv_hybrid_fits2.rds")))


all_results_weighted_hybrid <- readRDS(here::here('data/model_cache/indv_hybrid_fits.rds'))
ma <- map(all_results_weighted_hybrid, "best_params") |> map("c")


data.frame(id=names(ma),c=as.numeric(ma))
ma = map(all_results_weighted_hybrid, "best_params") |> map(magrittr::extract,c("c","lr")


map(all_results_weighted_hybrid,"best_params") |> pluck("c")
all_results_weighted_hybrid[["1"]]$b

 map(~ map(.x$best_params, pluck, "c"))

 map_df(~ map_df(.x$train, pluck, "d"), .id = "density")
#
#
#
#
#



viewPreds <- function(dsAvg, c, lr, input.layer, output.layer,weightExam=.5)
{
    train_data <- dsAvg[dsAvg$expMode2 == "Train", ]
    test_data <- dsAvg[dsAvg$expMode2 == "Test", ] |> group_by(x) |> summarise(y = mean(y))
     sim_results <- alm.sim(train_data, c, lr, input.layer, output.layer)
      weight.mat <- sim_results$wm
      trainVec <- c(unique(train_data$x))

    # predictions_test <- sapply(test_data$x, function(input) {
    #       exam.response(input, c, trainVec, input.layer, output.layer, weight.mat)
    #         })
    # (test_data <- cbind(test_data, exam.pred = predictions_test))  
    # 
    #  trainVec <- c(100,unique(train_data$x))
    # predictions_test <- sapply(test_data$x, function(input) {
    #       exam.response(input, c, trainVec, input.layer, output.layer, weight.mat)
    #         })
    # (test_data <- cbind(test_data, exam.pred0 = predictions_test))  
      
      predictions_test <- sapply(test_data$x, function(input) {
          exam.response2(input, c, train_data, input.layer, output.layer, weight.mat)
            })
    (test_data <- cbind(test_data, exam2.pred = predictions_test))

      predictions_test <- sapply(test_data$x, function(input) {
            #exam.response(input, c, trainVec, input.layer, output.layer, weight.mat)
            alm.response(input, c, input.layer, output.layer, weight.mat)$mean.response
            })

    (test_data <- cbind(test_data, alm.pred = predictions_test))
  predictions_test <- sapply(test_data$x, function(input) {
            #exam.response(input, c, trainVec, input.layer, output.layer, weight.mat)
             predict_alm_exam_hybrid(input, c, train_data,input.layer, output.layer, weight.mat)
            })
    (test_data <- cbind(test_data, hybrid.pred = predictions_test))


    (test_data <- cbind(test_data, alm.pred = predictions_test))
  predictions_test <- sapply(test_data$x, function(input) {
            #exam.response(input, c, trainVec, input.layer, output.layer, weight.mat)
             predict_alm_exam_weighted_hybrid(input, c, weightExam, train_data,input.layer, output.layer, weight.mat)
            })
    (test_data <- cbind(test_data, hybridFit.pred = predictions_test))
    

}

#
#
#
#
#
#
i=1
dat=train_data
input=800; corRep=975
c=.001; lr=.5
 weight.mat <- matrix(0.01, nrow = length(output.layer), ncol = length(input.layer))
  xt <- dat$x
  n <- nrow(dat)
  st <- numeric(n) # Initialize the vector to store mean responses
  for(i in 1:n) {
    trial <- alm.trial(dat$x[i], dat$y[i], c, lr, input.layer, output.layer, weight.mat)
    weight.mat <- trial$weight.mat
    st[i] <- trial$mean.response
  }
  dat <- cbind(dat,almResp=st)

alm_resp <- alm.response(input, c, input.layer,output.layer, weight.mat)
  

input.activation <- exp(-c * (input.layer - input)^2) / sum(exp(-c * (input.layer - input)^2))
output.activation <- weight.mat %*% input.activation
output.probability <- output.activation / sum(output.activation)
mean.response <- sum(output.layer * output.probability)
#
#
#
#
#
