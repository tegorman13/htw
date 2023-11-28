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
pacman::p_load(tidyverse,data.table,abc,future,furrr,here,patchwork, conflicted)
conflict_prefer_all("dplyr", quiet = TRUE)
walk(c("fun_alm","fun_model", "Display_Functions"), ~ source(here::here(paste0("Functions/", .x, ".R"))))
ds <- readRDS(here::here("data/e1_md_11-06-23.rds"))  |> as.data.table()
dsv <- ds |> filter(condit=="Varied")  
dsc <- ds |> filter(condit=="Constant") 





tMax=84
avg_dsv <- ds |> filter(condit=="Varied",expMode2=="Train") |> group_by(tr) %>%
  mutate(bandInt2 = sample(rep(c(800, 1000, 1200), each = tMax / 3), tMax, replace = FALSE)[tr]) %>%
  filter(bandInt2 == x) |> select(-bandInt2) |> group_by(tr,condit,x,expMode2) |> summarise(y=mean(y),.groups="keep") |>
  rbind(dsv |> filter(expMode2=="Test") |> group_by(tr,condit,x,expMode2) |> summarise(y=mean(y),.groups="keep") ) |> as.data.table()

avg_dsc <- dsc |> group_by(condit,expMode2,tr,x) |> summarise(y=mean(y),.groups="keep") |> as.data.table()



# ds |> filter(expMode2=="Test") |> group_by(condit,expMode2,tr,x) |> summarise(y=mean(y),n=n(),.groups="keep") |> 
#   ggplot(aes(x=tr,y=y,fill=condit)) + stat_summary(geom="bar",fun=mean,position=position_dodge()) + facet_wrap(~x) + ggtitle("Mean per Trial")

# ds |> filter(expMode2=="Test") |> group_by(condit,expMode2,tr,x) |> summarise(y=median(y),n=n(),.groups="keep") |> 
#   ggplot(aes(x=tr,y=y,fill=condit)) + stat_summary(geom="bar",fun=mean,position=position_dodge()) + facet_wrap(~x) + ggtitle("Median per Trial")


# avg_dsv |> filter(expMode2=="Train") |> group_by(x,tr) |> summarise(y=mean(y),n=n(),.groups="keep") |> 
#   ggplot(aes(x=tr,y=y)) + stat_summary(geom="bar",fun=mean) + facet_wrap(~x) + ggtitle("N per trial")



input_layer =  c(100,350,600,800,1000,1200)
output_layer = input_layer



#
#
#
#
#
#

generate_prior_c_lr <- function(n) {
  prior_samples <- tibble(
    c = runif(n, 0.000001, 5),
    lr = runif(n, 0.000001, 5),
  )
  return(prior_samples)
}




full_sim_exam <- function(data, c, lr,pred_fun=exam.response, input_layer, output_layer,return_dat="test_data",mode="sim") {
  train_data <- data[expMode2=="Train", c("condit","tr","expMode2", "x","y")] 
  test_data <- data[expMode2=="Test", c("condit","tr","expMode2", "x","y")] 
  trainVec=sort(unique(train_data$x))
 if (train_data$condit[1] != "Varied") {
  trainVec <- c(0, trainVec)
}
  
  train_results <- alm.sim(train_data, c, lr, input_layer, output_layer)

  test_prediction <- map_dbl(test_data$x, ~ pred_fun(.x, c, input_layer, output_layer, train_results$wm,  trainVec=trainVec))

  train_data$pred <- train_results$d$almResp
  test_data$pred <- test_prediction
  
  #fd <- rbind(train_data,test_data)
  fd = eval(parse(text=paste0("rbind(",return_dat,")")))
  if(mode=="sim"){return(fd$pred)
    }else {return(fd)}
  
}


full_sim_alm <- function(data, c, lr,pred_fun=alm.responseOnly, input_layer, output_layer,return_dat="test_data",mode="sim") {
  train_data <- data[expMode2=="Train", c("condit","tr","expMode2", "x","y")] 
  test_data <- data[expMode2=="Test", c("condit","tr","expMode2", "x","y")] 
  trainVec=sort(unique(train_data$x))
  
  train_results <- alm.sim(train_data, c, lr, input_layer, output_layer)
  test_prediction <- map_dbl(test_data$x, ~ alm.responseOnly(.x, c, input_layer, output_layer, train_results$wm,  trainVec=trainVec))

  train_data$pred <- train_results$d$almResp
  test_data$pred <- test_prediction
  
  #fd <- rbind(train_data,test_data)
  fd = eval(parse(text=paste0("rbind(",return_dat,")")))
  if(mode=="sim"){return(fd$pred)
    }else {return(fd)}
  
}

full_sim_alt_exam <- function(data, c, lr, pred_fun=alt_exam, input_layer, output_layer,return_dat="test_data",mode="sim") {
  train_data <- data[expMode2=="Train", c("condit","tr","expMode2", "x","y")] 
  test_data <- data[expMode2=="Test", c("condit","tr","expMode2", "x","y")] 
  trainVec=sort(unique(train_data$x))
  trainVecY=train_data$y
  
  train_results <- alm.sim(train_data, c, lr, input_layer, output_layer)
  test_prediction <- map_dbl(test_data$x, ~ alt_exam(.x, c, input_layer, output_layer, train_results$wm,  trainVecX=trainVec, trainVecY=trainVecY))

  train_data$pred <- train_results$d$almResp
  test_data$pred <- test_prediction
  
  #fd <- rbind(train_data,test_data)
  fd = eval(parse(text=paste0("rbind(",return_dat,")")))
  if(mode=="sim"){return(fd$pred)
    }else {return(fd)}
  
}



n_prior_samples <- 5000
prior_samples <- generate_prior_c_lr(n_prior_samples)





# 
# data=avg_dsc
# return_dat="train_data,test_data"
# return_dat="train_data"
# 
# eval(parse(text=paste0("rbind(",return_dat,")")))
# 
# paste0(rbind(parse(text=return_dat)))
# 

#
#
#
#
#

k=full_sim_exam(avg_dsv, .5, .2,exam.response,input_layer, output_layer,return_dat = "train_data,test_data")


plan(multisession)
simulated_data <- future_map_dfc(seq_len(nrow(prior_samples)), function(idx) {
  params <- prior_samples[idx, ]
  full_sim_exam(as.data.table(avg_dsv), params$c, params$lr,exam.response, input_layer, output_layer,return_dat = "train_data,test_data")
},.options = furrr_options(seed = T))




dst_obs <- avg_dsv[expMode2=="Test",]$y
dst_obs <- avg_dsv[expMode2=="Test" | expMode2=="Train",]$y

tolerance <- 0.1 * sd(dst_obs)

abc_result <- abc(
  target = dst_obs,
  param = prior_samples,
  sumstat =do.call(rbind, simulated_data),
  tol = .03,
  method = "rejection",
  names=colnames(dst_obs)
)


posterior_samples <- abc_result$unadj.values
colnames(posterior_samples) <- c("c", "lr")
head(posterior_samples)


posterior_samples_long <- tidyr::pivot_longer(as.data.frame(posterior_samples), everything())




postV <- ggplot(posterior_samples_long, aes(x=value)) +
  geom_density() +
  facet_wrap(~name, scales="free") +
  theme_minimal() +
  labs(x="Value", y="Density", title="Posterior Density Plots")
postV


sum_stats <- data.frame(mean = apply(posterior_samples, 2, mean),
  median = apply(posterior_samples, 2, median), mode=apply(posterior_samples, 2, Mode))
( sum_stats=rownames_to_column(sum_stats,var="parameter") )


density_c <- density(posterior_samples[, "c"])
mode_c <- density_c$x[which.max(density_c$y)]
density_lr <- density(posterior_samples[, "lr"])
mode_lr <- density_lr$x[which.max(density_lr$y)]


fsv <- full_sim_exam(as.data.table(avg_dsv), c=mode_c, lr=mode_lr,exam.response, input_layer, output_layer,mode="ret")

#fsv <- full_sim(as.data.table(avg_dsv), c=.007, lr=.047,exam.response, input_layer, output_layer)

fsv |> pivot_longer(y:pred,names_to="Resp",values_to="vx") |> 
  ggplot(aes(x,vx,fill=Resp, group=Resp)) +stat_summary(geom="bar",fun="mean",position=position_dodge())+
  scale_fill_manual(values=col_themes$wes2)+
  scale_x_continuous(breaks=sort(unique(ds$x)), labels=sort(unique(fsv$x))) +
  theme(legend.title = element_blank(), legend.position="top") +ggtitle("Fit to Test Only")





pgrid <- cbind(prior_samples,dist=abc_result$dist) |> mutate(dist=sqrt(dist))

pgrid |> filter(dist<50) |> ggplot(aes(x=c,y=lr,col=dist)) + geom_point() + viridis::scale_color_viridis()


pgrid |> filter(dist<1700) |> ggplot(aes(x=c,y=dist)) + geom_line()
pgrid |> filter(dist<1700) |> ggplot(aes(x=c,y=dist)) + geom_smooth()



pgrid |> filter(dist<50) |> ggplot(aes(x = c, y = lr, col = factor(dist))) +
  geom_point() +
  viridis::scale_color_viridis(discrete = TRUE)


pgrid |> ggplot(aes(x=c,y=lr,fill=as.factor(dist))) + geom_point() + viridis::scale_fill_viridis()


pgrid |> ggplot(aes(x = c, y = lr, col = factor(dist))) +
  geom_point() +
  viridis::scale_color_viridis(discrete = TRUE)

# 
# dim(simulated_data$y)[2]
# dim(sim)
# dim(t(sim))
#
#
#
#
#
#

lin <- abc(target=dst_obs, param=prior_samples, sumstat=do.call(rbind, simulated_data), tol=.08, hcorr =FALSE, method = "loclinear")


posterior_samples <- lin$unadj.values
colnames(posterior_samples) <- c("c", "lr")
head(posterior_samples)
posterior_samples_long <- tidyr::pivot_longer(as.data.frame(posterior_samples), everything())

postV <- ggplot(posterior_samples_long, aes(x=value)) +
  geom_density() + facet_wrap(~name, scales="free") +
  labs(x="Value", y="Density", title="Posterior Density Plots")
postV


sum_stats <- data.frame(mean = apply(posterior_samples, 2, mean),
  median = apply(posterior_samples, 2, median), mode=apply(posterior_samples, 2, Mode))
( sum_stats=rownames_to_column(sum_stats,var="parameter") )



summary(lin)

plot(lin,prior_samples,subsample=100)
hist(lin)




cv.lin <- cv4abc(param=prior_samples, sumstat=do.call(rbind, simulated_data),nval=15,tols=c(.01,.05,.1,.2,.3,.5), method="rejection", statistic="mode")
summary(cv.lin)

cv.lin2 <- cv4abc(param=prior_samples, sumstat=do.call(rbind, simulated_data),nval=15,tols=c(.01,.05,.1,.2,.3,.5),abc.out=abc_result, method="rejection", statistic="mode")

summary(cv.lin2)


#
#
#
#
#
#


postV <- ggplot(posterior_samples_long, aes(x=value)) +
  geom_density() +
  facet_wrap(~name, scales="free") +
  theme_minimal() +
  labs(x="Value", y="Density", title="Posterior Density Plots")
postV


sum_stats <- data.frame(mean = apply(posterior_samples, 2, mean),
  median = apply(posterior_samples, 2, median), mode=apply(posterior_samples,2,Mode))
( sum_stats=rownames_to_column(sum_stats,var="parameter") )

#
#
#
#
#
#

id1 <- ds |> filter(id==3)
n_prior_samples <- 500
prior_samples <- generate_prior_c_lr(n_prior_samples)

plan(multisession)
simulated_data <- future_map_dfc(seq_len(nrow(prior_samples)), function(idx) {
  params <- prior_samples[idx, ]
  full_sim_alt_exam(as.data.table(id1), params$c, params$lr,alt_exam, input_layer, output_layer,return_dat = "train_data,test_data")
},.options = furrr_options(seed = T))




dst_obs <- id1[expMode2=="Test",]$y
dst_obs <- id1[expMode2=="Test" | expMode2=="Train",]$y

tolerance <- 0.1 * sd(dst_obs)

abc_result <- abc(
  target = dst_obs,
  param = prior_samples,
  sumstat =do.call(rbind, simulated_data),
  tol = .1,
  method = "rejection",
  names=colnames(dst_obs)
)



posterior_samples <- abc_result$unadj.values
colnames(posterior_samples) <- c("c", "lr")

posterior_samples_long <- tidyr::pivot_longer(as.data.frame(posterior_samples), everything())


postV <- ggplot(posterior_samples_long, aes(x=value)) +
  geom_density() +
  facet_wrap(~name, scales="free") +
  theme_minimal() +
  labs(x="Value", y="Density", title="Posterior Density Plots")
postV


sum_stats <- data.frame(mean = apply(posterior_samples, 2, mean),
  median = apply(posterior_samples, 2, median), mode=apply(posterior_samples,2,Mode))
( sum_stats=rownames_to_column(sum_stats,var="parameter") )


density_c <- density(posterior_samples[, "c"])
mode_c <- density_c$x[which.max(density_c$y)]
density_lr <- density(posterior_samples[, "lr"])
mode_lr <- density_lr$x[which.max(density_lr$y)]




fsv <- full_sim_alt_exam(as.data.table(id1), c=mode_c, lr=mode_lr,exam.response, input_layer, output_layer,mode="ret")

fsv |> pivot_longer(y:pred,names_to="Resp",values_to="vx") |> 
  ggplot(aes(x,vx,fill=Resp, group=Resp)) +stat_summary(geom="bar",fun="mean",position=position_dodge())+
  scale_fill_manual(values=col_themes$wes2)+
  scale_x_continuous(breaks=sort(unique(ds$x)), labels=sort(unique(fsv$x))) +
  theme(legend.title = element_blank(), legend.position="top") +ggtitle("Fit to Test Only")





#
#
#
