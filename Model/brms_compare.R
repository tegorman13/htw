
pacman::p_load(dplyr,purrr,tidyr,ggplot2, data.table, here, patchwork, conflicted, 
               stringr,future,furrr, knitr, reactable, flextable,ggstance, htmltools,ggdist)
#conflict_prefer_all("dplyr", quiet = TRUE)
walk(c("flextable","dplyr"), conflict_prefer_all, quiet = TRUE)
#options(brms.backend="cmdstanr",mc.cores=4)
options(digits=2, scipen=999, dplyr.summarise.inform=FALSE)
walk(c("Display_Functions","fun_alm","fun_indv_fit","fun_model", "prep_model_data"), ~ source(here::here(paste0("Functions/", .x, ".R"))))

library(brms)
options(brms.backend="cmdstanr",mc.cores=1)


invisible(list2env(load_sbj_data(), envir = .GlobalEnv))
invisible(list2env(load_e1(), envir = .GlobalEnv))

e1Sbjs <- e1 |> group_by(id,condit) |> summarise(n=n())

pdl <- post_dat_l |> rename("bandInt"=x) |> left_join(testAvgE1,by=c("id","condit","bandInt")) |> 
    filter(rank<=1,Fit_Method=="Test_Train", !(Resp=="Observed")) |> mutate(aerror = abs(error))

pdl_all <- post_dat_l |> rename("bandInt"=x) |> left_join(testAvgE1,by=c("id","condit","bandInt")) |> 
    filter(rank<=1, !(Resp=="Observed")) |> mutate(aerror = abs(error))


pd <- post_dat |> rename("bandInt"=x) |> left_join(testAvgE1,by=c("id","condit","bandInt")) |> 
    filter(rank<=1,Fit_Method=="Test_Train") 






e1_ee_brm <- brm(data=pdl,
  error ~ 1 + (Model*condit)+bandType + (1|id) + (1|bandInt), 
  #file = paste0(here("data/model_cache/e1_ee1_1_rf2.rds")),
  chains=1,silent=0, iter=500, control=list(adapt_delta=0.92, max_treedepth=11))

summary(e1_ee_brm)
as.data.frame(bayestestR::describe_posterior(e1_ee_brm, centrality = "Mean"))[seq(1:6)]



e1_ee_abs <- brm(data=pdl,
  aerror ~ 1 + (Model*condit)+bandType + (1|id) + (1|bandInt), 
  #file = paste0(here("data/model_cache/e1_ee1_1_rf2.rds")),
  chains=1,silent=0, iter=500, control=list(adapt_delta=0.92, max_treedepth=11))

summary(e1_ee_abs)
as.data.frame(bayestestR::describe_posterior(e1_ee_abs, centrality = "Mean"))[seq(1:6)]


e1_ee_abs2 <- brm(data=pdl,
  aerror ~ 1 + Model*condit*bandType + (1|id) + (1|bandInt), 
  #file = paste0(here("data/model_cache/e1_ee1_1_rf2.rds")),
  chains=1,silent=0, iter=500, control=list(adapt_delta=0.92, max_treedepth=11))

summary(e1_ee_abs2)
as.data.frame(bayestestR::describe_posterior(e1_ee_abs2, centrality = "Mean"))[seq(1:6)]


e1_ee_abs3 <- brm(data=pdl,
  aerror ~ 1 + Model*condit*bandType + (1+Model|id), 
  #file = paste0(here("data/model_cache/e1_ee1_1_rf2.rds")),
  chains=1,silent=0, iter=500, control=list(adapt_delta=0.92, max_treedepth=11))

summary(e1_ee_abs3)
as.data.frame(bayestestR::describe_posterior(e1_ee_abs3, centrality = "Mean"))[seq(1:6)]





re <- data.frame(ranef(e1_ee_abs3, pars = "ModelEXAM")$id[, ,'ModelEXAM']) |> 
  tibble::rownames_to_column("id") |> 
  left_join(e1Sbjs,by="id")

re |> ggplot(aes(x=Estimate, fill=condit)) + geom_density()
  
re |>  mutate(id=reorder(id,Estimate)) |>
  ggplot(aes(y=id, x=Estimate,fill=condit,color=condit)) + 
    geom_pointrange(aes(xmin=Q2.5, xmax=Q97.5)) +
     ggh4x::facet_wrap2(~condit,axes="all",scales="free_y")



e1_ee_abs_id <- brm(data=pdl,
  aerror ~ 0 + id:Model, 
  chains=1,silent=0, iter=500, control=list(adapt_delta=0.92, max_treedepth=11))

summary(e1_ee_abs_id)
as.data.frame(bayestestR::describe_posterior(e1_ee_abs_id, centrality = "Mean"))[seq(1:6)]

b3_coef <- fixef(e1_ee_brm_ae_fm) |> as_tibble() %>% mutate(term=rownames(fixef(e1_ee_brm_ae_fm))) 
exam_coef <- b3_coef |> filter(term=="ModelEXAM") |> pull(Estimate)
head(b3_coef)
# extract id's from term column - Intercept is id 1, term=ModelEXAM doesn't have an id, so don't include
k <- b3_coef |> filter(term!="Intercept" & term!="ModelEXAM") |> mutate(id=as.factor(str_extract(term,"\\d+"))) |> left_join(testAvgE1 |> filter(bandInt==100) |> ungroup() |> select(-bandType,-vb,-bandInt),by="id") |> 
  mutate(exam_coef=exam_coef) 

# add new variable parType, which is set to "Model" if term contains "ModelEXAM" and "id" otherwise
k <- k |> mutate(parType=case_when(
  str_detect(term,"ModelEXAM") ~ "Model",
  TRUE ~ "id"
))
k |> ggplot(aes(x=parType,y=Estimate,fill=condit)) + geom_boxplot()

k |> filter(parType=="Model") |> mutate(id=reorder(id,Estimate), Estimate=Estimate+exam_coef) |> 
  ggplot(aes(x=id,y=Estimate,fill=condit)) + 
  geom_col() +  ggh4x::facet_grid2(~condit,axes="all",scales="free_y", independent = "y") + coord_flip()

# pull out ids in the form of b_id67:bandInt:ModelEXAM etc.
k |> filter(parType=="Model") |> pull(term) |> str_extract("b_id\\d+:bandInt:ModelEXAM") |> unique()


e1_ee_abs_id1 <- brm(data=pdl,
  aerror ~ 1 + id:Model, 
  chains=1,silent=0, iter=500, control=list(adapt_delta=0.92, max_treedepth=11))

as.data.frame(bayestestR::describe_posterior(e1_ee_abs_id1, centrality = "Mean"))[seq(1:6)]


e1_ee_abs_rf<- brm(data=pdl,
  aerror ~ 1 + condit*Model + (0+Model|id), 
  chains=1,silent=0, iter=500, control=list(adapt_delta=0.92, max_treedepth=11))

as.data.frame(bayestestR::describe_posterior(e1_ee_abs_rf, centrality = "Mean"))[seq(1:6)]





e1_ee_abs4 <- brm(data=pdl,
  aerror ~ 1 + Model*condit + (1+Model|id), 
  chains=1,silent=0, iter=500, control=list(adapt_delta=0.92, max_treedepth=11))
as.data.frame(bayestestR::describe_posterior(e1_ee_abs4, centrality = "Mean"))[seq(1:6)]

re <- data.frame(ranef(e1_ee_abs4, pars = "ModelEXAM")$id[, ,'ModelEXAM']) |> 
  tibble::rownames_to_column("id") |> 
  left_join(e1Sbjs,by="id")

re |> ggplot(aes(x=Estimate, fill=condit)) + geom_density()
  
re |>  mutate(id=reorder(id,Estimate)) |>
  ggplot(aes(y=id, x=Estimate,fill=condit,color=condit)) + 
    geom_pointrange(aes(xmin=Q2.5, xmax=Q97.5)) +
     ggh4x::facet_wrap2(~condit,axes="all",scales="free_y")


re |>  mutate(adjust=fixef(e1_ee_brm_ae_fm)[,1]["ModelEXAM"] + 
   fixef(e1_ee_abs4)[,1]["ModelEXAM:conditVaried"] * (condit=="Varied"), 
   Estimate=Estimate+adjust ) |> 
mutate(id=reorder(id,Estimate)) |>
  ggplot(aes(y=id, x=Estimate,fill=condit,color=condit)) + 
    geom_pointrange(aes(xmin=Q2.5, xmax=Q97.5)) +
     ggh4x::facet_wrap2(~condit,axes="all",scales="free_y")






fm_model_condit<- brm(data=pdl_all,
  aerror ~ 1 + condit*Model*Fit_Method*bandInt, 
  chains=1,silent=0, iter=500, control=list(adapt_delta=0.92, max_treedepth=11))

as.data.frame(bayestestR::describe_posterior(fm_model_condit, centrality = "Mean"))[seq(1:6)]





model_id = brm(data=pdl,
  aerror ~ Model + (1+Model|id), 
  chains=1,silent=0, iter=500, control=list(adapt_delta=0.92, max_treedepth=11))

as.data.frame(bayestestR::describe_posterior(model_id, centrality = "Mean"))[seq(1:6)]


re <- data.frame(ranef(model_id, pars = "ModelEXAM")$id[, ,'ModelEXAM']) |> 
  tibble::rownames_to_column("id") |> 
  left_join(e1Sbjs,by="id")




model_id = brm(data=pdl,
  aerror ~ Model + (1+Model|id), 
  chains=1,silent=0, iter=500, control=list(adapt_delta=0.92, max_treedepth=11))

as.data.frame(bayestestR::describe_posterior(model_id, centrality = "Mean"))[seq(1:6)]


re <- data.frame(ranef(model_id, pars = "ModelEXAM")$id[, ,'ModelEXAM']) |> 
  tibble::rownames_to_column("id") |> 
  left_join(e1Sbjs,by="id")


re |> ggplot(aes(x=Estimate, fill=condit)) + geom_density()
  
re |>  mutate(id=reorder(id,Estimate)) |>
  ggplot(aes(y=id, x=Estimate,fill=condit,color=condit)) + 
    geom_pointrange(aes(xmin=Q2.5, xmax=Q97.5)) +
     ggh4x::facet_wrap2(~condit,axes="all",scales="free_y")


re |>  mutate(adjust=fixef(model_id)[,1]["ModelEXAM"],
   #fixef(model_id)[,1]["ModelEXAM:conditVaried"] * (condit=="Varied"), 
   Estimate=Estimate+adjust ) |> 
mutate(id=reorder(id,Estimate)) |>
  ggplot(aes(y=id, x=Estimate,fill=condit,color=condit)) + 
    geom_pointrange(aes(xmin=Q2.5, xmax=Q97.5)) +
     ggh4x::facet_wrap2(~condit,axes="all",scales="free_y")







model_id = brm(data=pdl,
  aerror ~ Model*bandInt + (1+Model|id), 
  chains=1,silent=0, iter=500, control=list(adapt_delta=0.92, max_treedepth=11))

as.data.frame(bayestestR::describe_posterior(model_id, centrality = "Mean"))[seq(1:6)]


re <- data.frame(ranef(model_id, pars = "ModelEXAM")$id[, ,'ModelEXAM']) |> 
  tibble::rownames_to_column("id") |> 
  left_join(e1Sbjs,by="id")


re |> ggplot(aes(x=Estimate, fill=condit)) + geom_density()
  
re |>  mutate(id=reorder(id,Estimate)) |>
  ggplot(aes(y=id, x=Estimate,fill=condit,color=condit)) + 
    geom_pointrange(aes(xmin=Q2.5, xmax=Q97.5)) +
     ggh4x::facet_wrap2(~condit,axes="all",scales="free_y")


re |>  mutate(adjust=fixef(model_id)[,1]["ModelEXAM"] +
   fixef(model_id)[,1]["ModelEXAM:bandInt"], 
   Estimate=Estimate+adjust ) |> 
mutate(id=reorder(id,Estimate)) |>
  ggplot(aes(y=id, x=Estimate,fill=condit,color=condit)) + 
    geom_pointrange(aes(xmin=Q2.5, xmax=Q97.5)) +
     ggh4x::facet_wrap2(~condit,axes="all",scales="free_y")


