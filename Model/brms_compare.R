
pacman::p_load(dplyr,purrr,tidyr,ggplot2, data.table, here, patchwork, conflicted, 
               stringr,future,furrr, knitr, reactable, flextable,ggstance, htmltools,ggdist)
#conflict_prefer_all("dplyr", quiet = TRUE)
walk(c("flextable","dplyr"), conflict_prefer_all, quiet = TRUE)
#options(brms.backend="cmdstanr",mc.cores=4)
options(digits=2, scipen=999, dplyr.summarise.inform=FALSE)
walk(c("Display_Functions","fun_alm","fun_indv_fit","fun_model", "prep_model_data"), ~ source(here::here(paste0("Functions/", .x, ".R"))))

library(emmeans)
library(brms)
library(tidybayes)
options(brms.backend="cmdstanr",mc.cores=1)


invisible(list2env(load_sbj_data(), envir = .GlobalEnv))
invisible(list2env(load_e1(), envir = .GlobalEnv))

e1Sbjs <- e1 |> group_by(id,condit) |> summarise(n=n())

pdl <- post_dat_l |> rename("bandInt"=x) |> #left_join(testAvgE1,by=c("id","condit","bandInt")) |> 
    filter(rank<=1,Fit_Method=="Test_Train", !(Resp=="Observed")) |> mutate(aerror = abs(error))

pdl_all <- post_dat_l |> rename("bandInt"=x) |> #left_join(testAvgE1,by=c("id","condit","bandInt")) |> 
    filter(rank<=1, !(Resp=="Observed")) |> mutate(aerror = abs(error))


pd <- post_dat_avg |> rename("bandInt"=x) |> left_join(testAvgE1,by=c("id","condit","bandInt")) |> 
    filter(rank<=1,Fit_Method=="Test_Train")  |> mutate(aerror = abs(error))


full_e1_id_emmeans <- readRDS("~/Library/CloudStorage/GoogleDrive-tegorman13@gmail.com/My Drive/HTW/gl/full_e1_id_emmeans.rds") |> 
  pluck(contrast) |> as_tibble() 

k=full_e1_id_emmeans |> mutate(Model = str_extract(contrast, "^\\w+"), bandInt = as.numeric(str_extract(contrast, "(?<=bandInt)\\d+"))) |> 
  ggplot(aes(x=condit,y=estimate,fill=Model)) + stat_halfeye() + 
  #ggh4x::facet_grid2(~bandInt,axes="all",scales="free_y", independent = "y")
  facet_wrap(~bandInt)


full_e1_id_emmeans |> filter(id %in% levels(e1$id)[1:10]) |> 
  mutate(Model = str_extract(contrast, "^\\w+"), bandInt = as.numeric(str_extract(contrast, "(?<=bandInt)\\d+"))) |> 
  ggplot(aes(x=bandInt,y=estimate,col=Model)) + stat_halfeye() + 
  facet_wrap(~id)

testAvgE1 |> filter(id %in% levels(e1$id)[1:10]) |> ggplot(aes(x=bandInt,y=vxAvg,col=condit)) + 
  geom_point() + facet_wrap(~id)

post_dat_avg |> filter(id %in% levels(e1$id)[1:10], rank==1,Fit_Method=="Test_Train") |>
  ggplot(aes(x=x,y=pred,col=Model)) + 
  geom_point() + facet_wrap(~id)  + geom_point(aes(x=x,y=y),col="black",alpha=.2)


e1_ee_brm_ae <- brm(data=pdl,
                    aerror ~  Model * condit *dist+ (1+bandInt|id), 
                    #file = paste0(here("data/model_cache/e1_ae_modelCond_RFint.rds")),
                    chains=2,silent=1, iter=1000, control=list(adapt_delta=0.92, max_treedepth=11))

bayestestR::describe_posterior(e1_ee_brm_ae)
wrap_plots(plot(conditional_effects(e1_ee_brm_ae),points=FALSE,plot=FALSE))


e1_ee_brm_ae_vb <- brm(data=pdl,
                    aerror ~ Model * condit *vb + (1+vb+model|id), 
                    file = paste0(here("data/model_cache/e1_ee_brm_ae_vb2.rds")),
                    chains=2,silent=1, iter=1000, control=list(adapt_delta=0.92, max_treedepth=11))

bayestestR::describe_posterior(e1_ee_brm_ae_vb)
wrap_plots(plot(conditional_effects(e1_ee_brm_ae_vb),points=FALSE,plot=FALSE))


con8vb <- e1_ee_brm_ae_vb |>   emmeans(pairwise ~ Model | vb, by = "id",at=list(vb=levels(pdl$vb)[1:3]), re_formula = NULL)


# condit   id    contrast   estimate lower.HPD upper.HPD     n
# <chr>    <chr> <chr>         <dbl>     <dbl>     <dbl> <int>
#   1 Constant 3     ALM - EXAM    134.       85.0     188.     58

con8vb$contrasts |>  as_tibble() |>  left_join(e1Sbjs, by=c("id","condit")) |> 
  mutate(id=reorder(id,estimate)) |> 
  ggplot(aes(x=estimate,y=id)) + stat_halfeye() + 
  geom_vline(xintercept=0,linetype="dashed") + facet_wrap(~vb)






plot(conditional_effects(e1_ee_brm_ae, 
                         effects = "condit:Model", 
                         conditions=make_conditions(e1_ee_brm_ae,vars=c("dist","aerror"))),
     points=FALSE,plot=TRUE)

marginaleffects::plot_predictions(e1_ee_brm_ae, condition = list(
  "Model","condit","dist"
)
)



pd |> ggplot(aes(x=distAvg,y=error,color=condit)) + geom_point() + facet_wrap(~Model+vb)

brm_dist_ae <- pd %>% brm(data=.,aerror~distAvg*condit*bandInt*Model,chains=1)
bayestestR::describe_posterior(brm_dist_ae)
wrap_plots(plot(conditional_effects(brm_dist_ae),points=FALSE,plot=FALSE))

brm_dist_ae2 <- pd %>% brm(data=.,distAvg~c*condit*lr*Model+vb + (1 + Model|id),chains=1)
bayestestR::describe_posterior(brm_dist_ae2)
wrap_plots(plot(conditional_effects(brm_dist_ae2),points=FALSE,plot=FALSE))


plot(conditional_effects(brm_dist_ae2, 
                         effects = "condit:Model", 
                         conditions=make_conditions(brm_dist_ae2,vars=c("vb","aerror"))),
     points=FALSE,plot=TRUE)

marginaleffects::plot_predictions(brm_dist_ae, condition = list(
  "Model","condit","aerror",
  "bandInt" = c(100)
)
)

marginaleffects::plot_predictions(brm_dist_ae, condition = list(
  "Model","distAvg"
)
)




e1_c_brm1 <- pdl %>% brm(data=.,c~condit,chains=1)
bayestestR::describe_posterior(e1_c_brm1)

e1_c_brm2 <- pdl %>% brm(data=.,c~condit*Model*aerror,chains=1)
bayestestR::describe_posterior(e1_c_brm2)

# try again taking in to consideration that c is typically small, ~ .001, it cannot be zero, and it was generated from a lognormal distribution
#e1_c_brm3 <- pdl %>% brm(data=.,c~condit*Model*aerror,chains=1, prior = c(prior(normal(0,.1), class = "Intercept"), prior(normal(0,.1), class = "b"))


e1_c_brm3 <- pdl %>% 
  brm(data = ., c ~ condit * Model * aerror, chains = 1,
      prior = c(prior(normal(0, 0.1), class = "Intercept"),
                prior(normal(0, 0.1), class = "b"),
                prior(lognormal(0, 0.1), class = "sigma")),
      control = list(adapt_delta = 0.92, max_treedepth = 11))
bayestestR::describe_posterior(e1_c_brm3)









e1_ee_brm <- brm(data=pdl,
                 aerror ~  Model * condit*bandType + (1+bandInt|id), 
                 #file = paste0(here("data/model_cache/e1_ee1_1_rf2.rds")),
                 chains=1,silent=0, iter=500, control=list(adapt_delta=0.92, max_treedepth=11))

summary(e1_ee_brm)
as.data.frame(bayestestR::describe_posterior(e1_ee_brm, centrality = "Mean"))[seq(1:6)]


wrap_plots(plot(conditional_effects(e1_ee_brm),points=FALSE,plot=FALSE))


e1_ee_brm <- brm(data=pdl_all,
                 aerror ~  Model * condit*Fit_Method + (1+bandInt|id), 
                 #file = paste0(here("data/model_cache/e1_ee1_1_rf2.rds")),
                 chains=1,silent=0, iter=500, control=list(adapt_delta=0.92, max_treedepth=11))

as.data.frame(bayestestR::describe_posterior(e1_ee_brm, centrality = "Mean"))[seq(1:6)]

plot(conditional_effects(e1_ee_brm, 
                         effects = "Model:condit", 
                         conditions=make_conditions(e1_ee_brm,"Fit_Method" )),
     points=FALSE,plot=TRUE)


marginaleffects::plot_predictions(e1_ee_brm, condition = list(
  "Model","condit",
  "bandInt" = c(100,350,600,800,1000,1200)
)
)



### with Fit_Method

e1_ee_brm_ae_fm <- brm(data=pdl_all,
                       aerror ~ 1 + Model*condit*bandInt*Fit_Method  + (1+Model|id), 
                       file = paste0(here("data/model_cache/e1_ae1_1_fm_rf.rds")),
                       chains=2,silent=1, iter=1000, control=list(adapt_delta=0.92, max_treedepth=11))

wrap_plots(plot(conditional_effects(e1_ee_brm_ae_fm),points=FALSE,plot=FALSE))
as.data.frame(bayestestR::describe_posterior(e1_ee_brm_ae_fm, centrality = "Mean"))[seq(1:6)]



plot(conditional_effects(e1_ee_brm_ae_fm, 
                         effects = "Model:condit", 
                         conditions=make_conditions(e1_ee_brm,"Fit_Method" )),
     points=FALSE,plot=TRUE)


plot(conditional_effects(e1_ee_brm_ae_fm, 
                         effects = "Model:condit"),
     points=FALSE,plot=TRUE)

plot(conditional_effects(e1_ee_brm_ae_fm, 
                         effects = "Model:condit",
                         conditions=make_conditions(e1_ee_brm_ae_fm,"Fit_Method" )),
     points=FALSE,plot=TRUE)


plot(conditional_effects(e1_ee_brm_ae_fm, 
                         effects = "Fit_Method:bandInt",
                         conditions=make_conditions(e1_ee_brm_ae_fm,"Model" )),
     points=FALSE,plot=TRUE)

e1_ee_brm_ae_fm |> emmeans( ~condit *bandInt*Fit_Method*Model, 
                       at=list(bandInt=c(100,350,600,800,1000,1200) )) |>
  gather_emmeans_draws() |>
  ggplot(aes(x=bandInt,y=.value,color=condit)) + 
  stat_dist_pointinterval(position="dodge") +
  facet_wrap(~Model+Fit_Method)


e1_ee_brm_ae_fm |> emmeans( ~condit *Fit_Method*Model) |>
  gather_emmeans_draws() |>
  ggplot(aes(x=Model,y=.value,color=condit)) + 
  stat_dist_pointinterval(position="dodge") +
  facet_wrap(~Fit_Method)

e1_ee_brm_ae_fm |> emmeans( ~condit *Model) |>
  gather_emmeans_draws() |>
  ggplot(aes(x=Model,y=.value,color=condit)) + 
  stat_dist_pointinterval(position="dodge") 

marginaleffects::plot_predictions(e1_ee_brm_ae_fm, condition = list(
  "Model","condit","Fit_Method",
  "bandInt" = c(100,350,600)
)
)

marginaleffects::plot_predictions(e1_ee_brm_ae_fm, condition = list(
  "Model","condit","Fit_Method",
  "bandInt" = c(800,1000,1200)
)
)


conditional_effects(e1_ee_brm_ae_fm,effects = "Model:condit",conditions=make_conditions(e1_ee_brm,"bandInt" ))

marginaleffects::plot_predictions(e1_ee_brm_ae_fm, condition = list(
  "Model","condit","Fit_Method",
  "bandInt" = c(100,350,600,800,1000,1200)
)
)



e1_ee_brm_ae <- brm(data=pdl,
                    aerror ~ 1 + Model*condit*bandType + (1+Model*bandInt|id) + (1+Model|bandInt), 
                    file = paste0(here("data/model_cache/e1_ae1_1_rf6.rds")),
                    chains=2,silent=1, iter=1000, control=list(adapt_delta=0.92, max_treedepth=11))
summary(e1_ee_brm)

plot(conditional_effects(e1_ee_brm_ae,effects = "condit:bandType"), points=FALSE,plot=TRUE)




marginaleffects::plot_predictions(e1_ee_brm_ae, condition = list(
  "Model","condit","bandType"
)
)





marginaleffects::plot_predictions(e1_ee_brm_ae_fm, condition = list(
  "Model","condit","Fit_Method",
  "bandInt" = c(100,350,600,800,1000,1200)
)
)




  conditional_effects(e1_ee_brm_ae_fm, 
                      effects = "condit:bandInt",
                      int_conditions = list(bandInt=c(100,350, 600, 800, 1000, 1200)))

  conditional_effects(e1_ee_brm_ae_fm, 
                      effects = "Model:bandInt",
                      int_conditions = list(bandInt=c(100,350, 600, 800, 1000, 1200)))


  
 mc <-  make_conditions(e1_ee_brm,"bandInt" )  

 new_data_grid=map_dfr(1, ~data.frame(unique(testAvgE1[,c("bandInt")])))
 out <- list(bandInt=c(100,350,600,800,1000,1200))
 out <- list(100,350,600,800,1000,1200)
 out$cond__ <- brms:::rows2labels(new_data_grid)
 
 conditional_effects(e1_ee_brm_ae_fm,effects = "Model:condit",conditions=out)
 
 
 
 marginaleffects::plot_predictions(e1_ee_brm_ae_fm, condition = list(
    "Model","condit",
   "bandInt" = c(100,350,600,800,1000,1200)
 )
 )
 
 
 
 
 #make_condtions
 vars <- rev(as.character(vars))
 if (!is.data.frame(x) && "data" %in% names(x)) {
   x <- x$data
 }
 x <- as.data.frame(x)
 out <- brms:::named_list(vars)
 for (v in vars) {
   tmp <- get(v, x)
   if (brms:::is_like_factor(tmp)) {
     tmp <- levels(as.factor(tmp))
   }
   else {
     tmp <- mean(tmp, na.rm = TRUE) + (-1:1) * sd(tmp, 
                                                  na.rm = TRUE)
   }
   out[[v]] <- tmp
 }
 out <- rev(expand.grid(out))
 out$cond__ <- brms:::rows2labels(out, ...)
 out
 
 
 
 
 
 
 
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


