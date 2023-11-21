


generate_simulated_data <- function(n_agents = 100, n_train_trials = 151,trainVec=c(4,5,6), input_layer, output_layer) {
  simulated_data <- list()

  for (i in 1:n_agents) {
    c_sim <- sample(c_range, 1)  # Randomly sample c parameter
    lr_sim <- sample(lr_range, 1)  # Randomly sample learning rate

    # Initialize a weight matrix for the ALM model
    weight.mat <- matrix(0.000001, 
                         nrow = length(output_layer), 
                         ncol = length(input_layer))
    #diag(weight.mat) <- 2

    unique_inputs <- trainVec
    repetitions <- n_train_trials / length(unique_inputs)
    training_inputs <- rep(unique_inputs, each = repetitions)
    correct_outputs <- training_inputs * 2

    # Pre-create agent_data tibble for both training and testing trials
    train_data <- tibble(
      x = training_inputs,
      y = correct_outputs,
      expMode = rep('train', n_train_trials),
      modelType = rep('ALM', n_train_trials),
      model_resp = rep(NA_real_, n_train_trials)
    )
    # randomize order of train_data
    train_data <- train_data[sample(nrow(train_data)), ]

    test_inputs <- rep(input_layer, each = 4)
    test_data <- tibble(
      x = test_inputs,
      y = test_inputs*2,  # Placeholder, as correct outputs are not relevant in testing
      expMode = rep('test', length(test_inputs)),
      modelType = rep(c('ALM', 'EXAM',"exam2","exam0"), length(input_layer)), 
      model_resp = rep(NA_real_, length(test_inputs))
    )

    agent_data <- bind_rows(train_data, test_data) |> mutate(agent=i,c=c_sim, lr=lr_sim) |>
        relocate(agent, c, lr,.before=x)

    # Training loop
    for (t in 1:n_train_trials) {
      trial_result <- alm.trial(agent_data$x[t], agent_data$y[t], c_sim, lr_sim, input_layer, output_layer, weight.mat)
      agent_data$model_resp[t] <- trial_result$mean.response
      weight.mat <- trial_result$weight.mat
    }
    trainVecY <- agent_data$model_resp[agent_data$expMode == 'train']
    # Testing phase: Generate responses using map_dbl

  
    alm_preds <- map_dbl(input_layer, ~alm.responseOnly(.x, c_sim, input_layer, output_layer, weight.mat))
    exam_preds <- map_dbl(input_layer, ~exam.response(.x, c_sim, input_layer, output_layer, weight.mat, trainVec = trainVec ))
    exam_preds0 <- map_dbl(input_layer, ~exam.response(.x, c_sim, input_layer, output_layer, weight.mat,trainVec=c(0,trainVec)))
    
    exam2_preds <- map_dbl(input_layer, ~alt_exam(.x, c_sim, input_layer, output_layer, weight.mat,trainVecX = trainVec,trainVecY ))

    
    agent_data$model_resp[agent_data$modelType == 'ALM' & agent_data$expMode == 'test'] <- alm_preds
    agent_data$model_resp[agent_data$modelType == 'EXAM' & agent_data$expMode == 'test'] <- exam_preds
    agent_data$model_resp[agent_data$modelType == 'exam2' & agent_data$expMode == 'test'] <- exam2_preds
    agent_data$model_resp[agent_data$modelType == 'exam0' & agent_data$expMode == 'test'] <- exam_preds0


    # agent_data |> filter(expMode=="test")
    simulated_data[[i]] <- agent_data
  }

  return(simulated_data)
}

input_layer <- c(1,2,3,4,5,6)
output_layer <- input_layer*2

# Defining parameter ranges for simulation
c_range <- seq(from = 0.1, to = 1, by = 0.1)  # Example range for generalization parameter
lr_range <- seq(from = 0.5, to = 1, by = 0.5)  # Example range for learning rate


sd <- generate_simulated_data(n_agents = 100, n_train_trials = 90,trainVec=c(4), input_layer, output_layer)
d <- bind_rows(sd) |> mutate(dev=y-model_resp, abs_dev=abs(dev), lab=paste0(agent,"_c:",c,"_lr:",lr)) |> 
    relocate(agent, c, lr,lab,.before=x) |> group_by(agent,expMode,lab,modelType) |> 
    mutate(mean_dev=mean(dev)) |> # sort by mean_dev
    arrange(mean_dev,agent,expMode)




d |> filter(expMode=="train") |> 
    ggplot(aes(x=x,y=y))+
    stat_summary(geom="point",fun="mean")+ylim(c(0,25))+
    stat_summary(geom="bar",fun="mean", position=position_dodge(),aes(y=model_resp,fill=modelType)) 


d |> filter(expMode=="test") |> 
    ggplot(aes(x=x,y=y))+
    stat_summary(geom="point",fun="mean")+ylim(c(0,25))+
    stat_summary(geom="bar",fun="mean", position=position_dodge(),aes(y=model_resp,fill=modelType)) +
    stat_summary(geom="errorbar",fun.data="mean_se", position=position_dodge(),aes(x=x,y=model_resp,group=modelType)) 







d |> filter(expMode=="test", agent %in% unique(d$agent)[1:20]) |> 
    ggplot(aes(x=x,y=y))+
    geom_point()+ylim(c(0,15))+
    geom_line(aes(y=model_resp,color=modelType)) +
    facet_wrap(~lab) 


d |> filter(expMode=="train", agent %in% unique(d$agent)[1:50]) |> 
    mutate(trial = row_number()) |>
    ggplot(aes(x=trial,y=dev, col=factor(x)))+
    stat_summary(geom="line",fun="mean")+ylim(c(0,25))+
    stat_summary(geom="errorbar",fun.data="mean_se")+ylim(c(0,25))+
    facet_wrap(modelType~lab) +ylim(c(0,15))


d |> filter(expMode=="test", agent %in% unique(d$agent)[1:10]) |> 
    select(agent,x,y,modelType, model_resp,mean_dev) |> DT::datatable()


head(d)
head(sd)

head(sd[[1]])
sd[[1]] |> filter(expMode=="train")
sd[[1]] |> filter(expMode=="test")

