


load_init <- function(){
  # Load libraries
pacman::p_load(dplyr,purrr,tidyr,ggplot2, data.table, here, patchwork, conflicted, 
               stringr,future,furrr, knitr, reactable, flextable,ggstance, htmltools,ggdist)
#conflict_prefer_all("dplyr", quiet = TRUE)
walk(c("flextable","dplyr"), conflict_prefer_all, quiet = TRUE)

options(digits=2, scipen=999, dplyr.summarise.inform=FALSE)
walk(c("Display_Functions","fun_alm","fun_indv_fit","fun_model"), ~ source(here::here(paste0("Functions/", .x, ".R"))))

}

# Only call load_init if not already loaded
if (!all(c("dplyr", "purrr", "conflicted", "here") %in% .packages())) {
    load_init()
    print("Loaded libraries and functions.")
}
#  && exists("dplyr.summarise.inform")



generate_data <- function(Model, post_samples, data, num_samples = 1, return_dat = "train_data, test_data") {
  # Filter data for the specific id without invalidating selfref
  sbj_data <- copy(data[id == post_samples$id[1]])
  simulation_function <- ifelse(Model == "EXAM", full_sim_exam, full_sim_alm)

  target_data <- switch(return_dat,
                        "test_data" = copy(sbj_data[expMode2 == "Test"]),
                        "train_data" = copy(sbj_data[expMode2 == "Train"]),
                        "train_data, test_data" = copy(sbj_data[expMode2 %in% c("Test", "Train")]))
  
  post_samples <- post_samples[order(mean_error)][1:num_samples, .(c, lr, mean_error, rank = .I)]

  simulated_data_list <- lapply(1:nrow(post_samples), function(i) {
    params <- post_samples[i]
    sim_data <- simulation_function(sbj_data, params$c, params$lr, input_layer = input_layer, 
                                    output_layer = output_layer, return_dat = return_dat)
    sim_data_dt <- data.table(id = sbj_data$id[1], condit = sbj_data$condit[1], 
                              expMode2 = target_data$expMode2, Model = Model,tr=target_data$tr,
                              y = target_data$y, x = target_data$x, c = params$c, 
                              lr = params$lr, mean_error = params$mean_error, rank = i,
                              pred = sim_data)
    return(sim_data_dt)
  })
  
  result_dt <- rbindlist(simulated_data_list)
  setcolorder(result_dt, c("id", "condit", "expMode2","tr", "c", "lr", "x", "y", "pred"))
  return(result_dt)
}




load_sbj_data <- function(){
e1 <- readRDS(here("data/e1_08-21-23.rds")) 
e2 <- readRDS(here("data/e2_08-04-23.rds")) 
e3 <- readRDS(here("data/e3_08-04-23.rds")) 

# combine all 3 experiments
d <- rbind(e1,e2,e3)

d <- d |> 
    mutate(trainCon=case_when(
    bandOrder=="Original" ~ "800",
    bandOrder=="Reverse" ~ "600",
    TRUE ~ NA_character_
    ), trainCon=as.numeric(trainCon)) 


nbins=5
train <-  d |> filter(expMode2=="Train") |> group_by(id,condit,fb,bandOrder, vb) |> 
    mutate(Trial_Bin = cut( gt.train, breaks = seq(1, max(gt.train),length.out=nbins+1),include.lowest = TRUE, labels=FALSE)) 
train_max <- train |> filter(Trial_Bin == nbins, bandInt==trainCon)
train_avg <- train_max |> group_by(id,condit,fb,bandOrder) |> summarise(train_end = mean(dist))
train_avg2 <- train_max |> select(id,condit,fb,bandOrder,expMode2,vb,bandInt,dist,vx)

test2 <- d |> filter(expMode2=="Test") |> 
  select(id,condit,fb,bandOrder,expMode2,vb,bandInt,dist,vx) |> 
  rbind(train_avg2) |>
  left_join(train_avg, by=c("id","condit", "fb", "bandOrder")) 

test <- d |> filter(expMode2=="Test") |> left_join(train_avg, by=c("id","condit", "fb", "bandOrder")) |>
  select(id,condit,bandType,bandInt,vb,vx,dist,train_end,fb,bandOrder)

testE1 <- e1 <-  e1 |> filter(expMode2 == "Test") 
testAvgE1 <- testE1 %>% group_by(id, condit, vb, bandInt,bandType,tOrder) %>%
  summarise(nHits=sum(dist==0),vxAvg=mean(vx),distAvg=mean(dist),sdist=mean(sdist),n=n(),Percent_Hit=nHits/n)


return(tibble::lst(e1=e1, e2=e2, e3=e3, d=d, train=train, testAvgE1))

}




load_e1 <- function(n_posterior=50, n_post_train=50){


e1 <- readRDS(here::here("data/e1_md_11-06-23.rds"))  |> as.data.table()
nbins <- 3
e1Sbjs <- e1 |> group_by(id,condit) |> summarise(n=n())
fd <- readRDS(here("data/e1_08-21-23.rds"))
testE1 <-  fd |> filter(expMode2 == "Test") 
testE1Avg <- testE1 %>% group_by(id, condit, vb, bandInt,bandType,tOrder) %>%
  summarise(nHits=sum(dist==0),vxAvg=mean(vx),distAvg=mean(dist),sdist=mean(sdist),n=n(),Percent_Hit=nHits/n)

trainE1Avg <- fd |> filter(expMode2 == "Train") |> group_by(id) |> 
  mutate(tr=trial,x=vb,Block=case_when(expMode2=="Train" ~ cut(tr,breaks=seq(1,max(tr), length.out=nbins+1),include.lowest=TRUE,labels=FALSE),
                                         expMode2=="Test" ~ 4)) |> 
  group_by(id,condit,vb,x,Block) |> 
  summarise(dist=mean(dist),y=mean(vx))

input_layer <<- output_layer <<-  c(100,350,600,800,1000,1200)
ids2 <- c(1,66,36)
file_name <- "n_iter_200_ntry_300_5354"
#file_name <- "n_iter_400_ntry_100_2944"


ind_fits <- map(list.files(here(paste0('data/abc_reject/'),file_name),full.names=TRUE), readRDS)
ind_fits_df <- ind_fits |> map(~list(dat=.x[[1]], Model = .x[["Model"]], Fit_Method=.x[["Fit_Method"]]))
ind_fits_df <- ind_fits_df |> map(~rbindlist(.x$dat) |> mutate(Model = .x$Model, Fit_Method = .x$Fit_Method)) |> rbindlist() 




#future::plan(multisession)

nestSbjModelFit <- ind_fits_df %>% nest(.by=c(id,Model,Fit_Method))

# organize test data predictions
# post_dat <- nestSbjModelFit |> mutate(pp=furrr::future_pmap(list(id,Model,Fit_Method,data), ~{
#    generate_data(..2, ..4 |> mutate(id=..1), ds, num_samples = n_posterior, return_dat="test_data")
#    })) |> 
#   select(Fit_Method,pp,-data) |>  
#   unnest(pp) |>  filter(expMode2=="Test") |> as.data.table()
# 
# saveRDS(post_dat, here("data/model_cache/post_dat.rds"))

post_dat <- readRDS(here("data/model_cache/post_dat.rds"))

post_dat_avg <- post_dat |> group_by(id, condit, Model, Fit_Method, x, c, lr, rank) |> 
  mutate(error2 = y - pred) |>
  summarise(y = mean(y), pred = mean(pred), error = y - pred, error2=mean(error2)) |> as.data.table()

setorder(post_dat_avg, id, x, rank)
post_dat_l <- melt(post_dat_avg, id.vars = c("id", "condit", "Model", "Fit_Method", "x", "c", "lr", "rank","error"),
                   measure.vars = c("pred", "y"), variable.name = "Resp", value.name = "val")
post_dat_l[, Resp := fifelse(Resp == "y", "Observed",
                             fifelse(Model == "ALM", "ALM", "EXAM"))]
setorder(post_dat_l, id, Resp)
#rm(post_dat_avg)

post_dat_l <- post_dat_l |> mutate(dist = case_when(
    val >= x & val <= x + 200 ~ 0,                 
    val < x ~ abs(x - val),                       
    val > x + 200 ~ abs(val - (x + 200)),           
    TRUE ~ NA_real_                                 
  ), signed_dist = case_when(
    val >= x & val <= x + 200 ~ 0,                 
    val < x ~ x - val,                             
    val > x + 200 ~ val - (x + 200),               
    TRUE ~ NA_real_                                 
  ))

post_dat <- post_dat |>  mutate(dist = case_when(
    y >= x & y <= x + 200 ~ 0,                 
    y < x ~ abs(x - y),                       
    y > x + 200 ~ abs(y - (x + 200)),           
    TRUE ~ NA_real_                                 
  ), pred_dist = case_when(
    pred >= x & pred <= x + 200 ~ 0,                 
    pred < x ~ abs(x - pred),                       
    pred > x + 200 ~ abs(pred - (x + 200)),           
    TRUE ~ NA_real_                                 
  ))

post_dat <- post_dat |> 
  left_join(testAvgE1 |> 
              select(id,condit,bandInt,bandType,vb,bandInt), 
            by=join_by(id,condit,x==bandInt))

post_dat_l <- post_dat_l |> 
  left_join(testAvgE1 |> 
              select(id,condit,bandInt,bandType,vb,bandInt), 
            by=join_by(id,condit,x==bandInt))

# organize training data predictions
# pd_train <- nestSbjModelFit |> mutate(pp=furrr::future_pmap(list(id,Model,Fit_Method,data), ~{
#    generate_data(..2, ..4 |> mutate(id=..1), ds, num_samples = n_post_train, return_dat="train_data")
#    })) |>
#   select(Fit_Method,pp,-data) |>
#   unnest(pp) |> as.data.table() |> filter(expMode2=="Train")

#saveRDS(pd_train, here("data/model_cache/pd_train.rds"))

pd_train <- readRDS(here("data/model_cache/pd_train.rds"))

nbins <- 3
pd_train <- pd_train |> group_by(id,condit,Model,Fit_Method) |>
  mutate(Block=cut(tr,breaks=seq(1,max(tr), length.out=nbins+1),include.lowest=TRUE,labels=FALSE))
setorder(pd_train, id, x,Block, rank)

pd_train_l <- melt(pd_train, id.vars = c("id", "condit", "Model","Block", "Fit_Method", "x", "c", "lr", "rank"),
                   measure.vars = c("pred", "y"), variable.name = "Resp", value.name = "val") |> as.data.table()
pd_train_l[, Resp := fifelse(Resp == "y", "Observed",
                             fifelse(Model == "ALM", "ALM", "EXAM"))] 
setorder(pd_train_l, id,Block, Resp) 

pd_train_l <- pd_train_l  |>
  mutate(dist = case_when(
    val >= x & val <= x + 200 ~ 0,                 
    val < x ~ abs(x - val),                       
    val > x + 200 ~ abs(val - (x + 200)),           
    TRUE ~ NA_real_                                 
  ))

#plan(sequential)

return(e1_model <- tibble::lst(post_dat=post_dat, post_dat_l=post_dat_l, 
            pd_train=pd_train, pd_train_l=pd_train_l, post_dat_avg=post_dat_avg))

}







load_e2 <- function(n_posterior=50, n_post_train=50){



ds <- readRDS(here::here("data/e2_md_02-23-24.rds"))  |> as.data.table()
nbins <- 3

fd <- readRDS(here("data/e2_08-21-23.rds"))
testE2 <- fd |> filter(expMode2 == "Test") 
testAvgE2 <- testE2 %>% group_by(id, condit, vb, bandInt,bandType,tOrder) %>%
  summarise(nHits=sum(dist==0),vx=mean(vx),dist=mean(dist),sdist=mean(sdist),n=n(),Percent_Hit=nHits/n)

trainAvg <- fd |> filter(expMode2 == "Train") |> group_by(id) |> 
  mutate(tr=trial,x=vb,Block=case_when(expMode2=="Train" ~ cut(tr,breaks=seq(1,max(tr), length.out=nbins+1),include.lowest=TRUE,labels=FALSE),
                                         expMode2=="Test" ~ 4)) |> 
  group_by(id,condit,vb,x,Block) |> 
  summarise(dist=mean(dist),y=mean(vx))

input_layer <<- output_layer <<-  c(100,350,600,800,1000,1200)
ids2 <- c(1,66,36)

#file_name <- "e2_n_iter_50_ntry_200_2506"
#file_name <- "n_iter_400_ntry_100_2944"
#file_name <- "e2_n_iter_100_ntry_200_3436"
file_name <- "e2_n_iter_200_ntry_300_2344"

ind_fits <- map(list.files(here(paste0('data/abc_reject/'),file_name),full.names=TRUE), readRDS)
ind_fits_df <- ind_fits |> map(~list(dat=.x[[1]], Model = .x[["Model"]], Fit_Method=.x[["Fit_Method"]]))
ind_fits_df <- ind_fits_df |> map(~rbindlist(.x$dat) |> mutate(Model = .x$Model, Fit_Method = .x$Fit_Method)) |> rbindlist() 




generate_data <- function(Model, post_samples, data, num_samples = 1, return_dat = "train_data, test_data") {
  # Filter data for the specific id without invalidating selfref
  sbj_data <- copy(data[id == post_samples$id[1]])
  simulation_function <- ifelse(Model == "EXAM", full_sim_exam, full_sim_alm)

  target_data <- switch(return_dat,
                        "test_data" = copy(sbj_data[expMode2 == "Test"]),
                        "train_data" = copy(sbj_data[expMode2 == "Train"]),
                        "train_data, test_data" = copy(sbj_data[expMode2 %in% c("Test", "Train")]))
  
  post_samples <- post_samples[order(mean_error)][1:num_samples, .(c, lr, mean_error, rank = .I)]

  simulated_data_list <- lapply(1:nrow(post_samples), function(i) {
    params <- post_samples[i]
    sim_data <- simulation_function(sbj_data, params$c, params$lr, input_layer = input_layer, 
                                    output_layer = output_layer, return_dat = return_dat)
    sim_data_dt <- data.table(id = sbj_data$id[1], condit = sbj_data$condit[1], 
                              expMode2 = target_data$expMode2, Model = Model,tr=target_data$tr,
                              y = target_data$y, x = target_data$x, c = params$c, 
                              lr = params$lr, mean_error = params$mean_error, rank = i,
                              pred = sim_data)
    return(sim_data_dt)
  })
  
  result_dt <- rbindlist(simulated_data_list)
  setcolorder(result_dt, c("id", "condit", "expMode2","tr", "c", "lr", "x", "y", "pred"))
  return(result_dt)
}

#future::plan(multisession)

nestSbjModelFit <- ind_fits_df %>% nest(.by=c(id,Model,Fit_Method))


#  post_dat <- nestSbjModelFit |> mutate(pp=furrr::future_pmap(list(id,Model,Fit_Method,data), ~{
#     generate_data(..2, ..4 |> mutate(id=..1), ds, num_samples = 50, return_dat="test_data")
#     })) |> 
#    select(Fit_Method,pp,-data) |>  
#    unnest(pp) |>  filter(expMode2=="Test") |> as.data.table()

# saveRDS(post_dat, here("data/model_cache/post_dat_e2.rds"))

post_dat <- readRDS(here("data/model_cache/post_dat_e2.rds"))



post_dat_avg <- post_dat |> group_by(id, condit, Model, Fit_Method, x, c, lr, rank) |> 
  mutate(error2 = y - pred) |>
  summarise(y = mean(y), pred = mean(pred), error = y - pred, error2=mean(error2)) |> as.data.table()

setorder(post_dat_avg, id, x, rank)
post_dat_l <- melt(post_dat_avg, id.vars = c("id", "condit", "Model", "Fit_Method", "x", "c", "lr", "rank","error"),
                   measure.vars = c("pred", "y"), variable.name = "Resp", value.name = "val")
post_dat_l[, Resp := fifelse(Resp == "y", "Observed",
                             fifelse(Model == "ALM", "ALM", "EXAM"))]
setorder(post_dat_l, id, Resp)
#rm(post_dat_avg)

post_dat_l <- post_dat_l |> mutate(dist = case_when(
    val >= x & val <= x + 200 ~ 0,                 
    val < x ~ abs(x - val),                       
    val > x + 200 ~ abs(val - (x + 200)),           
    TRUE ~ NA_real_                                 
  ))

post_dat <- post_dat |> 
  left_join(testAvgE2 |> 
              select(id,condit,bandInt,bandType,vb,bandInt), 
            by=join_by(id,condit,x==bandInt))

post_dat_l <- post_dat_l |> 
  left_join(testAvgE2 |> 
              select(id,condit,bandInt,bandType,vb,bandInt), 
            by=join_by(id,condit,x==bandInt))


# organize training data predictions
 pd_train <- nestSbjModelFit |> mutate(pp=furrr::future_pmap(list(id,Model,Fit_Method,data), ~{
   generate_data(..2, ..4 |> mutate(id=..1), ds, num_samples = 20, return_dat="train_data")
    })) |>
   select(Fit_Method,pp,-data) |>
  unnest(pp) |> as.data.table() |> filter(expMode2=="Train")

#saveRDS(pd_train, here("data/model_cache/pd_train2.rds"))

#pd_train <- readRDS(here("data/model_cache/pd_train2.rds"))

nbins <- 3
pd_train <- pd_train |> group_by(id,condit,Model,Fit_Method) |>
  mutate(Block=cut(tr,breaks=seq(1,max(tr), length.out=nbins+1),include.lowest=TRUE,labels=FALSE))
setorder(pd_train, id, x,Block, rank)

pd_train_l <- melt(pd_train, id.vars = c("id", "condit", "Model","Block", "Fit_Method", "x", "c", "lr", "rank"),
                   measure.vars = c("pred", "y"), variable.name = "Resp", value.name = "val") |> as.data.table()
pd_train_l[, Resp := fifelse(Resp == "y", "Observed",
                             fifelse(Model == "ALM", "ALM", "EXAM"))] 
setorder(pd_train_l, id,Block, Resp) 

pd_train_l <- pd_train_l  |>
  mutate(dist = case_when(
    val >= x & val <= x + 200 ~ 0,                 
    val < x ~ abs(x - val),                       
    val > x + 200 ~ abs(val - (x + 200)),           
    TRUE ~ NA_real_                                 
  ))

#plan(sequential)
### Accounting for individual patterns


return(e2_model <- tibble::lst(post_dat=post_dat, post_dat_l=post_dat_l, 
                               pd_train=pd_train, pd_train_l=pd_train_l, post_dat_avg=post_dat_avg))


}












load_e3 <- function(n_posterior=50, n_post_train=50){

#ds <- readRDS(here::here("data/e1_md_11-06-23.rds"))  |> as.data.table()
ds <- readRDS(here::here("data/e3_md_02-23-24.rds"))  |> as.data.table()
nbins <- 3

fd <- readRDS(here("data/e3_08-21-23.rds"))
group_ids <- fd |> distinct(id,condit,bandOrder,fb)
testE3 <- fd |> filter(expMode2 == "Test") 
testAvgE3 <- testE3 %>% group_by(id, condit, vb, bandInt,bandType,tOrder) %>%
  summarise(nHits=sum(dist==0),vx=mean(vx),dist=mean(dist),sdist=mean(sdist),n=n(),Percent_Hit=nHits/n)

trainAvg <- fd |> filter(expMode2 == "Train") |> group_by(id) |> 
  mutate(tr=trial,x=vb,Block=case_when(expMode2=="Train" ~ cut(tr,breaks=seq(1,max(tr), length.out=nbins+1),include.lowest=TRUE,labels=FALSE),
                                         expMode2=="Test" ~ 4)) |> 
  group_by(id,condit,vb,x,Block) |> 
  summarise(dist=mean(dist),y=mean(vx))

input_layer <<- output_layer <<-  c(100,350,600,800,1000,1200)
ids2 <- c(1,66,36)

#file_name <- "e2_n_iter_50_ntry_200_2506"
#file_name <- "n_iter_400_ntry_100_2944"
#file_name <- "e2_n_iter_100_ntry_200_3436"
#file_name <- "e3_n_iter_50_ntry_100_3244"
file_name <- "e3_n_iter_200_ntry_300_2033"

ind_fits <- map(list.files(here(paste0('data/abc_reject/'),file_name),full.names=TRUE), readRDS)
ind_fits_df <- ind_fits |> map(~list(dat=.x[[1]], Model = .x[["Model"]], Fit_Method=.x[["Fit_Method"]]))
ind_fits_df <- ind_fits_df |> map(~rbindlist(.x$dat) |> mutate(Model = .x$Model, Fit_Method = .x$Fit_Method)) |> rbindlist() 


generate_data <- function(Model, post_samples, data, num_samples = 1, return_dat = "train_data, test_data") {
  # Filter data for the specific id without invalidating selfref
  sbj_data <- copy(data[id == post_samples$id[1]])
  simulation_function <- ifelse(Model == "EXAM", full_sim_exam, full_sim_alm)

  target_data <- switch(return_dat,
                        "test_data" = copy(sbj_data[expMode2 == "Test"]),
                        "train_data" = copy(sbj_data[expMode2 == "Train"]),
                        "train_data, test_data" = copy(sbj_data[expMode2 %in% c("Test", "Train")]))
  
  post_samples <- post_samples[order(mean_error)][1:num_samples, .(c, lr, mean_error, rank = .I)]

  simulated_data_list <- lapply(1:nrow(post_samples), function(i) {
    params <- post_samples[i]
    sim_data <- simulation_function(sbj_data, params$c, params$lr, input_layer = input_layer, 
                                    output_layer = output_layer, return_dat = return_dat)
    sim_data_dt <- data.table(id = sbj_data$id[1], condit = sbj_data$condit[1], 
                              expMode2 = target_data$expMode2, Model = Model,tr=target_data$tr,
                              y = target_data$y, x = target_data$x, c = params$c, 
                              lr = params$lr, mean_error = params$mean_error, rank = i,
                              pred = sim_data)
    return(sim_data_dt)
  })
  
  result_dt <- rbindlist(simulated_data_list)
  setcolorder(result_dt, c("id", "condit", "expMode2","tr", "c", "lr", "x", "y", "pred"))
  return(result_dt)
}

#future::plan(multisession)

nestSbjModelFit <- ind_fits_df %>% nest(.by=c(id,Model,Fit_Method))

# organize test data predictions
#  post_dat <- nestSbjModelFit |> mutate(pp=furrr::future_pmap(list(id,Model,Fit_Method,data), ~{
#     generate_data(..2, ..4 |> mutate(id=..1), ds, num_samples = 50, return_dat="test_data")
#     })) |> 
#    select(Fit_Method,pp,-data) |>  
#    unnest(pp) |>  filter(expMode2=="Test") |> as.data.table()
# 
#saveRDS(post_dat, here("data/model_cache/post_dat_e3.rds"))

post_dat <- readRDS(here("data/model_cache/post_dat_e3.rds"))

post_dat_avg <- post_dat |> group_by(id, condit, Model, Fit_Method, x, c, lr, rank) |> 
  mutate(error2 = y - pred) |>
  summarise(y = mean(y), pred = mean(pred), error = y - pred, error2=mean(error2)) |> as.data.table() |>
  left_join(group_ids, by = join_by(id,condit))

setorder(post_dat_avg, id, x, rank)
post_dat_l <- melt(post_dat_avg, id.vars = c("id", "condit","bandOrder","fb", "Model", "Fit_Method", "x", "c", "lr", "rank","error"),
                   measure.vars = c("pred", "y"), variable.name = "Resp", value.name = "val")
post_dat_l[, Resp := fifelse(Resp == "y", "Observed",
                             fifelse(Model == "ALM", "ALM", "EXAM"))]
setorder(post_dat_l, id, Resp)
#rm(post_dat_avg)


post_dat_l <- post_dat_l |> mutate(dist = case_when(
    val >= x & val <= x + 200 ~ 0,                 
    val < x ~ abs(x - val),                       
    val > x + 200 ~ abs(val - (x + 200)),           
    TRUE ~ NA_real_                                 
  ))

post_dat <- post_dat |> 
  left_join(testAvgE3 |> 
              select(id,condit,bandInt,bandType,vb,bandInt), 
            by=join_by(id,condit,x==bandInt))

post_dat_l <- post_dat_l |> 
  left_join(testAvgE3 |> 
              select(id,condit,bandInt,bandType,vb,bandInt), 
            by=join_by(id,condit,x==bandInt))


# organize training data predictions
#  pd_train <- nestSbjModelFit |> mutate(pp=furrr::future_pmap(list(id,Model,Fit_Method,data), ~{
#    generate_data(..2, ..4 |> mutate(id=..1), ds, num_samples = 20, return_dat="train_data")
#     })) |>
#    select(Fit_Method,pp,-data) |>
#   unnest(pp) |> as.data.table() |> filter(expMode2=="Train")

#saveRDS(pd_train, here("data/model_cache/pd_train_e3.rds"))

pd_train <- readRDS(here("data/model_cache/pd_train_e3.rds"))

nbins <- 3
pd_train <- pd_train |> group_by(id,condit,Model,Fit_Method) |>
  mutate(Block=cut(tr,breaks=seq(1,max(tr), length.out=nbins+1),include.lowest=TRUE,labels=FALSE))
setorder(pd_train, id, x,Block, rank)

pd_train_l <- melt(pd_train, id.vars = c("id", "condit", "Model","Block", "Fit_Method", "x", "c", "lr", "rank"),
                   measure.vars = c("pred", "y"), variable.name = "Resp", value.name = "val") |> as.data.table()
pd_train_l[, Resp := fifelse(Resp == "y", "Observed",
                             fifelse(Model == "ALM", "ALM", "EXAM"))] 
setorder(pd_train_l, id,Block, Resp) 

pd_train_l <- pd_train_l  |>
  mutate(dist = case_when(
    val >= x & val <= x + 200 ~ 0,                 
    val < x ~ abs(x - val),                       
    val > x + 200 ~ abs(val - (x + 200)),           
    TRUE ~ NA_real_                                 
  ))

#plan(sequential)

return(e3_model <- tibble::lst(post_dat=post_dat, post_dat_l=post_dat_l, 
                               pd_train=pd_train, pd_train_l=pd_train_l, post_dat_avg=post_dat_avg))



}