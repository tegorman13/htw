#
#
#
#
#
#
source(here::here("Functions", "packages.R"))
test <- readRDS(here("data/e1_08-21-23.rds")) |> filter(expMode2 == "Test") 
options(brms.backend="cmdstanr",mc.cores=4)
e1Sbjs <- test |> group_by(id,condit) |> summarise(n=n())
testAvg <- test %>% group_by(id, condit, vb, bandInt,bandType,tOrder) %>%
  summarise(nHits=sum(dist==0),vx=mean(vx),dist=mean(dist),sdist=mean(sdist),n=n(),Percent_Hit=nHits/n)


get_coef_details <- function(model, term_name) {
  broom.mixed::tidy(model) |>
    filter(term == term_name) |> 
    select(estimate, conf.low, conf.high) |> 
    mutate(across(where(is.numeric), \(x) round(x, 2)))
}


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
#| label: tbl-e1-test-nf-deviation
#| tbl-cap: "Testing Deviation - Empirical Summary"
#| tbl-subcap: ["Constant Testing - Deviation", "Varied Testing - Deviation"]
#| layout-ncol: 2
#| echo: fenced

result <- test_summary_table(test, "dist", mfun = list(mean = mean, median = median, sd = sd))
result$constant
result$varied
#
#
#
#| label: fig-e1-test-dev
#| fig-cap: E1. Deviations from target band during testing without feedback stage. 
test |>  ggplot(aes(x = vb, y = dist,fill=condit)) +
    stat_summary(geom = "bar", position=position_dodge(), fun = mean) +
    stat_summary(geom = "errorbar", position=position_dodge(.9), fun.data = mean_se, width = .4, alpha = .7) + 
  labs(x="Band", y="Deviation From Target")
#
#
#
#
#| label: tbl-e1-bmm-dist
#| tbl-cap: "Experiment 1. Bayesian Mixed Model predicting absolute deviation as a function of condition (Constant vs. Varied) and Velocity Band"

#contrasts(test$condit) 
#contrasts(test$vb)
modelName <- "e1_testDistBand_RF_5K"
e1_distBMM <- brm(dist ~ condit * bandInt + (1 + bandInt|id),
                      data=test,file=paste0(here::here("data/model_cache",modelName)),
                      iter=5000,chains=4)
GetModelStats(e1_distBMM) |> kable(escape=F,booktabs=T,caption="Model Coefficients")

e1_distBMM |> 
  emmeans("condit",by="bandInt",at=list(bandInt=c(100,350,600,800,1000,1200)),
          epred = TRUE, re_formula = NA) |> 
  pairs() |> gather_emmeans_draws()  |> 
   summarize(median_qi(.value),pd=sum(.value>0)/n()) |>
   select(contrast,Band=bandInt,value=y,lower=ymin,upper=ymax,pd) |> 
   mutate(across(where(is.numeric), \(x) round(x, 2)),
          pd=ifelse(value<0,1-pd,pd)) |>
   kbl(caption="Contrasts")


coef_details <- broom.mixed::tidy(e1_distBMM) |>
  filter(term == "conditVaried") |> 
  select(estimate, conf.low, conf.high) |> mutate(across(where(is.numeric), \(x) round(x, 2)))


#
#
#
#
#
#
#
#| label: fig-e1-bmm-dist
#| fig-cap: E1. Conditioinal Effect of Training Condition and Band. Ribbon indicated 95% Credible Intervals. 

condEffects <- function(m){
  m |> ggplot(aes(x = bandInt, y = .value, color = condit, fill = condit)) + 
  stat_dist_pointinterval() + stat_halfeye(alpha=.2) +
  stat_lineribbon(alpha = .25, size = 1, .width = c(.95)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5)) +
  ylab("Predicted X Velocity") + xlab("Band")
}

e1_distBMM |> emmeans( ~condit + bandInt, 
                       at = list(bandInt = c(100, 350, 600, 800, 1000, 1200))) |>
  gather_emmeans_draws() |>
 condEffects()+
  scale_x_continuous(breaks = c(100, 350, 600, 800, 1000, 1200), 
                     labels = levels(test$vb), 
                     limits = c(0, 1400)) 
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
#| label: fig-e1-test-vx
#| fig-cap: E1 testing x velocities. Translucent bands with dash lines indicate the correct range for each velocity band. 
#| fig-width: 11
#| fig-height: 9
test %>% group_by(id,vb,condit) |> plot_distByCondit()

#
#
#
#| label: tbl-e1-test-nf-vx
#| tbl-cap: "Testing vx - Empirical Summary"
#| tbl-subcap: ["Constant Testing - vx", "Varied Testing - vx"]
#| layout-ncol: 2
#| echo: fenced

result <- test_summary_table(test, "vx", mfun = list(mean = mean, median = median, sd = sd))
result$constant 
result$varied 

#
#
#
#| label: tbl-e1-bmm-vx
#| tbl-cap: "Experiment 1. Bayesian Mixed Model Predicting Vx as a function of condition (Constant vs. Varied) and Velocity Band"
#| tbl-subcap: ["Model fit to all 6 bands", "Model fit to 3 extrapolation bands"]
#| layout-ncol: 2
#| echo: fenced
e1_vxBMM <- brm(vx ~ condit * bandInt + (1 + bandInt|id),
                        data=test,file=paste0(here::here("data/model_cache", "e1_testVxBand_RF_5k")),
                        iter=5000,chains=4,silent=0,
                        control=list(adapt_delta=0.94, max_treedepth=13))
GetModelStats(e1_vxBMM ) |> kable(escape=F,booktabs=T, caption="Fit to all 6 bands")

cd1 <- get_coef_details(e1_vxBMM, "conditVaried")
sc1 <- get_coef_details(e1_vxBMM, "bandInt")
intCoef1 <- get_coef_details(e1_vxBMM, "conditVaried:bandInt")


modelName <- "e1_extrap_testVxBand"
e1_extrap_VxBMM <- brm(vx ~ condit * bandInt + (1 + bandInt|id),
                  data=test |>
                    filter(expMode=="test-Nf"),file=paste0(here::here("data/model_cache",modelName)),
                  iter=5000,chains=4)
GetModelStats(e1_extrap_VxBMM ) |> kable(escape=F,booktabs=T, caption="Fit to 3 extrapolation bands")


sc2 <- get_coef_details(e1_extrap_VxBMM, "bandInt")
intCoef2 <- get_coef_details(e1_extrap_VxBMM, "conditVaried:bandInt")

#
#
#
#
#
#
#
#| label: fig-e1-bmm-vx
#| fig-cap: Conditional effect of training condition and Band. Ribbons indicate 95% HDI. The steepness of the lines serves as an indicator of how well participants discriminated between velocity bands. 
#| fig-subcap: ["Model fit to all 6 bands", "Model fit to only 3 extrapolation bands"]
#| layout-ncol: 2

e1_vxBMM |> emmeans( ~condit + bandInt, 
                       at = list(bandInt = c(100, 350, 600, 800, 1000, 1200))) |>
  gather_emmeans_draws() |> 
  condEffects() +
  scale_x_continuous(breaks = c(100, 350, 600, 800, 1000, 1200), 
                     labels = levels(test$vb), 
                     limits = c(0, 1400))

e1_extrap_VxBMM |> emmeans( ~condit + bandInt, 
                       at = list(bandInt = c(100, 350, 600))) |>
  gather_emmeans_draws() |>
  condEffects() +
  scale_x_continuous(breaks = c(100, 350, 600), 
                     labels = levels(test$vb)[1:3], 
                     limits = c(0, 1000)) 

#
#
#
#
#
#| label: tbl-e1-slope-quartile
#| tbl-cap: "Slope coefficients by quartile, per condition"

new_data_grid=map_dfr(1, ~data.frame(unique(test[,c("id","condit","bandInt")]))) |> 
  dplyr::arrange(id,bandInt) |> 
  mutate(condit_dummy = ifelse(condit == "Varied", 1, 0)) 

indv_coefs <- as_tibble(coef(e1_vxBMM)$id, rownames="id")|> 
  select(id, starts_with("Est")) |>
  left_join(e1Sbjs, by=join_by(id) ) 


fixed_effects <- e1_vxBMM |> 
  spread_draws(`^b_.*`,regex=TRUE) |> arrange(.chain,.draw,.iteration)


random_effects <- e1_vxBMM |> 
  gather_draws(`^r_id.*$`, regex = TRUE, ndraws = 1500) |> 
  separate(.variable, into = c("effect", "id", "term"), sep = "\\[|,|\\]") |> 
  mutate(id = factor(id,levels=levels(test$id))) |> 
  pivot_wider(names_from = term, values_from = .value) |> arrange(id,.chain,.draw,.iteration)


 indvDraws <- left_join(random_effects, fixed_effects, by = join_by(".chain", ".iteration", ".draw")) |> 
  rename(bandInt_RF = bandInt,RF_Intercept=Intercept) |>
  right_join(new_data_grid, by = join_by("id")) |> 
  mutate(
    Slope = bandInt_RF+b_bandInt,
    Intercept= RF_Intercept + b_Intercept,
    estimate = (b_Intercept + RF_Intercept) + (bandInt*(b_bandInt+bandInt_RF)) + (bandInt * condit_dummy) * `b_conditVaried:bandInt`,
    SlopeInt = Slope + (`b_conditVaried:bandInt`*condit_dummy)
  ) 

  indvSlopes <- indvDraws |> group_by(id) |> median_qi(Slope,SlopeInt, Intercept,b_Intercept,b_bandInt) |>
  left_join(e1Sbjs, by=join_by(id)) |> group_by(condit) |>
    select(id,condit,Intercept,b_Intercept,starts_with("Slope"),b_bandInt, n) |>
  mutate(rankSlope=rank(Slope)) |> arrange(rankSlope)   |> ungroup()
 
  
  indvSlopes |> mutate(Condition=condit) |>  group_by(Condition) |> 
    reframe(enframe(quantile(SlopeInt, c(0.0,0.25, 0.5, 0.75,1)), "quantile", "SlopeInt")) |> 
  pivot_wider(names_from=quantile,values_from=SlopeInt,names_prefix="Q_") |>
  group_by(Condition) |>
  summarise(across(starts_with("Q"), list(mean = mean))) |> kbl()

#
#
#
#
#
#| label: fig-e1-bmm-bx2
#| fig-cap: Slope distributions between condition
#| fig-subcap: ["Slope estimates by participant - ordered from lowest to highest within each condition. ", "Destiny of slope coefficients by training group"]
#| layout-ncol: 2
#| fig-height: 9
#| fig-width: 10

  indvSlopes |> ggplot(aes(y=rankSlope, x=SlopeInt,fill=condit,color=condit)) + 
  geom_pointrange(aes(xmin=SlopeInt.lower , xmax=SlopeInt.upper)) + 
  labs(x="Estimated Slope", y="Participant")  + facet_wrap(~condit)

   ggplot(indvSlopes, aes(x = SlopeInt, color = condit)) + 
  geom_density() + labs(x="Slope Coefficient",y="Density")


#
#
#
#| label: fig-e1-indv-slopes
#| fig-cap: Subset of Varied and Constant Participants with the smallest and largest estimated slope values. Red lines represent the best fitting line for each participant, gray lines are 200 random samples from the posterior distribution. Colored points and intervals at each band represent the empirical median and 95% HDI. 
#| fig-subcap: ["subset with largest slopes", "subset with smallest slopes"]
#| fig-height: 11
#| fig-width: 12
nSbj <- 3
indvDraws  |> indv_model_plot(indvSlopes, testAvg, SlopeInt,rank_variable=Slope,n_sbj=nSbj,"max")
indvDraws |> indv_model_plot(indvSlopes, testAvg,SlopeInt, rank_variable=Slope,n_sbj=nSbj,"min")

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
