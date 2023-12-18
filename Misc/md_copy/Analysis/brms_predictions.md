
source(here::here("Functions", "packages.R"))
test <- readRDS(here("data/e1_08-04-23.rds")) |> 
  filter(expMode2 == "Test") 
options(brms.backend="cmdstanr",mc.cores=4)

################
### Setup #####
################
set.seed(12345)

e1Sbjs <- test |> group_by(id,condit) |> summarise(n=n())
testAvg <- test |> group_by(id,condit,vb,bandInt) |> 
  summarise(meanVx=mean(vx),medianVx=median(vx),sdVx=sd(vx))

                 
e1_vxBMM <- brm(bf(vx ~ condit * bandInt + (1 + bandInt|id),
                   sigma ~ condit * bandInt + (1+bandInt|id)),
                data=test,
                file=paste0(here::here("data/model_cache", "e1_testVxBand_RF_5k_Ml1")),
                iter=5000,chains=4,silent=0,
                control=list(adapt_delta=0.94, max_treedepth=13))


indv_coefs <- coef(e1_vxBMM)$id |> 
  as_tibble(rownames="id") |> 
  select(id, starts_with("Est")) |>
  left_join(e1Sbjs, by=join_by(id) ) |> 
  group_by(condit) |> 
  mutate(rank = rank(desc(Estimate.bandInt))) 

nTop <- 5
top_indv_coefs <- indv_coefs %>%
  filter(rank <=nTop) |>
  select(id,condit,rank,Est.BandInt=Estimate.bandInt, Est.Intercept=Estimate.Intercept) 

new_data_grid=map_dfr(1, ~data.frame(unique(test[,c("id","condit","bandInt")]))) |> 
  dplyr::arrange(id,bandInt) |> 
  mutate(condit_dummy = ifelse(condit == "Varied", 1, 0)) # |>
  #filter(id %in% top_indv_coefs$id)

fixed_effects <- fixef(e1_vxBMM)
random_effects <- coef(e1_vxBMM)$id


################


  




e1_vxBMM |> spread_draws(`^r_id.*$`,regex=TRUE)
e1_vxBMM |> gather_draws(`^r_id.*$`,regex=TRUE)





## Method 0 -linear prediction
linpred_manual <- new_data_grid |> 
  left_join(top_indv_coefs,by=join_by(id,condit)) |>
  mutate(mu = Est.Intercept + (Est.BandInt * bandInt) + fixed_effects[2] * condit_dummy) |>
  left_join(testAvg,by=join_by(id,condit,bandInt))


linpred_manual2 <- new_data_grid %>%
  left_join(top_indv_coefs, by=join_by(id,condit)) %>%
  mutate(mu = Est.Intercept +
           (Est.BandInt * bandInt) +
           fixed_effects[2] * condit_dummy +   # main effect of conditVaried
           fixed_effects[4] * bandInt * condit_dummy) |> # interaction term
  left_join(testAvg,by=join_by(id,condit,bandInt)) 



## Method 1 - predict
nd1 <- new_data_grid |> 
  mutate(pred=predict(e1_vxBMM, newdata = new_data))



## Method 2 - posterior predict
nd2 <- new_data_grid |> 
  mutate(pred=t(posterior_predict(e1_vxBMM, newdata = new_data_grid,ndraws=105))) |>
    mutate(pred=rowMeans(across(starts_with("pred")))) |>
  left_join(top_indv_coefs, by=join_by(id,condit)) |>
  left_join(testAvg,by=join_by(id,condit,bandInt))

sigma_draws <- spread_draws(e1_vxBMM, sigma)$sigma

postpred_manual <- new_data_grid %>%
  add_linpred_draws(e1_vxBMM, newdata = .) |> 
  mutate(condit_dummy = as.numeric(condit == "Varied")) |>
  mutate(y_new = round(rnorm(n(), mean = .linpred, sd = sigma_draws),1))

## 2B

# issue due to the fact that different draws as sampled for RF and FF, leading to issue when left_joining 
# try gathering both effects all at once?
random_effects <- e1_vxBMM |> 
  gather_draws(`^r_id.*$`, regex = TRUE, ndraws = 1000) |> 
  separate(.variable, into = c("effect", "id", "term"), sep = "\\[|,|\\]") |> 
  mutate(id = factor(id,levels=levels(test$id))) |> 
  pivot_wider(names_from = term, values_from = .value) |> arrange(id,.chain,.draw,.iteration)

fixed_effects <- e1_vxBMM |> 
  spread_draws(`^b_.*`,regex=TRUE) |> arrange(.chain,.draw,.iteration)

combined_df <- left_join(random_effects, fixed_effects, by = join_by(".chain", ".iteration", ".draw")) |> 
  rename(bandInt_RF = bandInt) |>
  right_join(new_data_grid, by = join_by("id")) |> 
  mutate(
    fixed_effects = b_Intercept +
      (bandInt * b_bandInt) +
      condit_dummy * b_conditVaried + 
      (bandInt * condit_dummy) * `b_conditVaried:bandInt` ,  # Note: Replace : with correct operator
    random_effects = Intercept + bandInt_RF,  # Assuming random effects for intercept and bandInt
    estimate = fixed_effects + random_effects
  )


# shouldn't need this as long as the entirety of fixed_effects is drawn
all_effects <- e1_vxBMM |> 
  gather_draws(b_Intercept, b_conditVaried, b_bandInt, `b_conditVaried:bandInt`, `^r_id.*$`, regex = TRUE, ndraws = 1) |> 
  separate(.variable, into = c("effect_type", "effect_details"), sep = "\\[", remove = FALSE, fill = "right", extra = "merge") |> 
  mutate(id = str_extract(effect_details, "\\d+"),
         term = str_extract(effect_details, "[A-Za-z]+")) |> 
  pivot_wider(names_from = c(effect_type, term), values_from = .value, names_glue = "{effect_type}_{term}")




## Method 3 - epred
nd3 <- e1_vxBMM %>%
  epred_draws(newdata = new_data_grid,ndraws = 1000) %>%
  as_tibble() |> 
  select(-.chain,-.iteration) |>
  left_join(top_indv_coefs, by=join_by(id,condit)) |>
  left_join(testAvg,by=join_by(id,condit,bandInt)) 
 



###### Marginal Effects

e1_distBMM |> emmeans( ~condit +bandInt, at=list(bandInt=c(100,350,600,800,1000,1200) ))
e1_distBMM |> emmeans( ~condit +bandInt, 
                       at=list(bandInt=c(100,350,600,800,1000,1200) ),
                       epred=TRUE, re_formula = NULL)

e1_distBMM |> emtrends(~bandInt+condit,var="bandInt",at=list(bandInt=c(100,350,600,800,1000,1200)))
e1_distBMM |> emtrends(~bandInt+condit,var="bandInt")
e1_distBMM |> marginaleffects(newdata = datagrid(bandInt=c(100,350,600,800,1000,1200)))

new_data_grid=map_dfr(1, ~data.frame(unique(test[,c("id","condit","bandInt")])))

e1_distBMM |> marginaleffects(newdata = datagrid(new_data_grid))
e1_distBMM |> plot_cme(variables="bandInt",condition="condit")


e1_distBMM |> 
  emmeans("condit",by="bandInt",at=list(bandInt=c(100,350,600,800,1000,1200)),epred = TRUE) |> 
  pairs() |> gather_emmeans_draws() |> mean_hdi(.width=.95)
  

cSamp <- e1_distBMM |> 
  emmeans("condit",by="bandInt",at=list(bandInt=c(100,350,600,800,1000,1200)),epred = TRUE) |> 
  pairs() |> gather_emmeans_draws()  |>
  group_by(contrast, .draw,bandInt) |> summarise(value=mean(.value), n=n())


 cSamp |> ggplot(aes(value))+geom_histogram() + 
  geom_vline(xintercept=0,alpha=.4)+
  facet_wrap(~bandInt,ncol=1) +
  ggtitle("Posterior of Constant - Varied contrasts")



 ameBand <- cSamp |> ggplot(aes(x=value,y="")) + 
  stat_halfeye() + 
  geom_vline(xintercept=0,alpha=.4)+
  facet_wrap(~bandInt,ncol=1) + labs(x="Marginal Effect (Constant - Varied)", y= NULL)+
  ggtitle("Average Marginal Effect")



bothConditRf <- e1_distBMM %>%
  epred_draws(newdata = new_data_grid,ndraws = 2000) |>
  ggplot(aes(x=.epred,y="Mean",fill=condit)) + stat_halfeye() +facet_wrap(~bandInt)

bothConditGM <- e1_distBMM %>%
  epred_draws(newdata = new_data_grid,ndraws = 2000, re_formula = NA) |>
  ggplot(aes(x=.epred,y="Mean",fill=condit)) + 
  stat_halfeye() +facet_wrap(~bandInt, ncol = 1) +
  labs(x="Predicted Deviation", y=NULL)+
  ggtitle("Grand Means") +theme(legend.position = "bottom")


bothConditGM | ameBand
(bothConditGM | ameBand) + plot_layout(widths=c(2,1.5))

e1_distBMM |> emmeans( ~condit +bandInt, 
                       at=list(bandInt=c(100,350,600,800,1000,1200) )) |>
  gather_emmeans_draws() |>
  ggplot(aes(x=bandInt,y=.value,color=condit)) + stat_lineribbon(alpha=.25,size=1)
  geom_half_violin()


# Need to investiage whey including RF's changes the results so much
  
  new_data_grid=map_dfr(1, ~data.frame(unique(test[,c("id","condit","bandInt")])))
  
  cSamp <- e1_distBMM |> 
    emmeans("condit",by="bandInt",at=list(bandInt=c(100,350,600,800,1000,1200)),
            epred = TRUE, re_formula = NULL) |> 
    pairs() |> gather_emmeans_draws()  |>
    group_by(contrast, .draw,bandInt) |> summarise(value=mean(.value), n=n())
  
  ameBand <- cSamp |> ggplot(aes(x=value,y="")) + 
    stat_halfeye() + 
    geom_vline(xintercept=0,alpha=.4)+
    facet_wrap(~bandInt,ncol=1) + labs(x="Marginal Effect (Constant - Varied)", y= NULL)+
    ggtitle("Average Marginal Effect")
  
  bothConditGM <- e1_distBMM %>%
    epred_draws(newdata = new_data_grid,ndraws = 2000, re_formula = NULL) |>
    ggplot(aes(x=.epred,y="Mean",fill=condit)) + 
    stat_halfeye() +facet_wrap(~bandInt, ncol = 1) +
    labs(x="Predicted Deviation", y=NULL)+
    ggtitle("Grand Means") +theme(legend.position = "bottom")
  
  (bothConditGM | ameBand) + plot_layout(widths=c(2,1.0))
  
  
  
  
  
  
  
##################
# Plot
ggplot(nd3, aes(x = vb, y = .epred)) +
  #stat_lineribbon(alpha = 0.3) +
  stat_halfeye(aes(fill=vb)) +
  labs(title = "Top 10 High Estimate.bandInt",
       x = "BandInt", y = "Predicted vx") +
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5)) +
  ggh4x::facet_nested (rank~condit)


ggplot(nd2, aes(x = vb, y = pred)) +
  #stat_lineribbon(alpha = 0.3) +
  stat_halfeye(aes(fill=vb)) +
  labs(title = "Top 10 High Estimate.bandInt",
       x = "BandInt", y = "Predicted vx") +
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5)) +
  ggh4x::facet_nested (rank~condit)


top_indv_coefs |> left_join(test,by=join_by(id,condit)) |> 
  ggplot(aes(x=bandInt,y=vx))+ stat_lineribbon(alpha = 0.3) +
  stat_halfeye(aes(fill=as.factor(bandInt))) +
  ggh4x::facet_nested (rank~condit)


plotSlope <- function(df,title="",colour=NULL){
  rectWidth=50
  df |> ggplot(aes(x = bandInt, y = vx))+ 
    geom_abline(
      aes(intercept = Est.Intercept, slope = Est.BandInt,color="black"),
      size = 1,alpha=.9
    ) +
    geom_point(aes(color=vb)) +
    geom_half_violin(aes(fill=vb),position=position_dodge(.5)) +
    geom_half_point(aes(fill=vb))+
    stat_halfeye(aes(fill=vb))+
    facet_wrap("id") +
    geom_rect(aes(xmin=bandInt-rectWidth,xmax=bandInt+rectWidth,ymin=bandInt,ymax=highBound,fill=vb),alpha=.1)+
    geom_segment(aes(x=bandInt-rectWidth,xend=bandInt+rectWidth,y=highBound,yend=highBound),alpha=1,linetype="dashed")+
    geom_segment(aes(x=bandInt - rectWidth,xend=bandInt+rectWidth,y=bandInt,yend=bandInt),alpha=1,linetype="dashed")+
    labs(x = "Velocity Band", y = "vxMean")+
    scale_x_continuous(labels=sort(unique(df$bandInt)),breaks=sort(unique(df$bandInt)))+
    ggtitle(title) + theme(legend.position = "none")+theme_classic()+guides(fill="none",color="none")
}


top_indv_coefs |> left_join(test,by=join_by(id,condit)) |> plotSlope()

testSlopeIndv %>% left_join(dtestAgg, by = c("id","condit")) %>% filter(condit=="Varied",rank<=6) %>% 
  plotSlope(.,colour="black",title="Largest Individually fit Varied Sbj. Slopes")