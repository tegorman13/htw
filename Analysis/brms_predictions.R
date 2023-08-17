
source(here::here("Functions", "packages.R"))
test <- readRDS(here("data/e1_08-04-23.rds")) |> 
  filter(expMode2 == "Test") 
options(brms.backend="cmdstanr",mc.cores=4)

set.seed(12345)

e1Sbjs <- test |> group_by(id,condit) |> summarise(n=n())
testAvg <- test |> group_by(id,condit,vb,bandInt) |> 
  summarise(meanVx=mean(vx),medianVx=median(vx),
            sdVx=sd(vx))


nested_settings <- strip_nested(
  background_x = list(element_rect(fill = "grey92"), NULL),
  by_layer_x = TRUE)


# e1_vxBMM <- brm(bf(vx ~ condit * bandInt + (1 + bandInt|id),
#                    sigma ~ condit * bandInt + (1+bandInt|id)),
#                 data=test,
#                 file=paste0(here::here("data/model_cache", "e1_testVxBand_RF_5k_Ml1")),
#                 iter=5000,chains=4,silent=0,
#                 control=list(adapt_delta=0.94, max_treedepth=13))


e1_vxBMM <- brm(vx ~ condit * bandInt + (1 + bandInt|id),
                        data=test,file=paste0(here::here("data/model_cache", "e1_testVxBand_RF_5k_Ml1")),
                        iter=5000,chains=4,silent=0,
                        control=list(adapt_delta=0.94, max_treedepth=13))
mt2 <-GetModelStats(e1_vxBMM ) |> kable(escape=F,booktabs=T)
mt2


indv_coefs <- coef(e1_vxBMM)$id |> 
  as_tibble(rownames="id") |> 
  select(id, starts_with("Est")) |>
  left_join(e1Sbjs, by=join_by(id) ) |> 
  # rank by Estimate.bandInt - within condit, higher is better
  group_by(condit) |> 
  mutate(rank = rank(desc(Estimate.bandInt))) 

top_indv_coefs <- indv_coefs %>%
  filter(rank <= 4) |>
  select(id,condit,rank,Estimate.bandInt, Estimate.Intercept) 
top_ids <- top_indv_coefs$id

fixed_effects <- fixef(e1_vxBMM)
random_effects <- coef(e1_vxBMM)$id



new_data_grid=map_dfr(1, ~data.frame(unique(test[,c("id","condit","bandInt")]))) |> 
  dplyr::arrange(id,bandInt) |> 
  mutate(condit_dummy = ifelse(condit == "Varied", 1, 0)) |>
  filter(id %in% top_ids) |>
  




e1_vxBMM |> spread_draws(`^r_id.*$`,regex=TRUE)
e1_vxBMM |> gather_draws(`^r_id.*$`,regex=TRUE)





## Method 0 -linear prediction
linpred_manual <- new_data_grid |> 
  left_join(top_indv_coefs,by=join_by(id,condit)) |>
  mutate(mu = Estimate.Intercept + (Estimate.bandInt * bandInt) + fixed_effects[2] * condit_dummy) |>
  left_join(testAvg,by=join_by(id,condit,bandInt))


linpred_manual2 <- new_data_grid %>%
  left_join(top_indv_coefs, by=join_by(id,condit)) %>%
  mutate(mu = Estimate.Intercept +
           (Estimate.bandInt * bandInt) +
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
  pivot_wider(names_from = term, values_from = .value)

# Combine with the fixed effects
fixed_effects <- e1_vxBMM |> 
  spread_draws(`^b_.*`,regex=TRUE,ndraws=1000) 

combined_df <- left_join(random_effects, fixed_effects, by = c(".chain", ".iteration", ".draw")) |> 
  rename(bandInt_RF = bandInt) |>
  left_join(new_data_grid, by = join_by("id")) |> 
  mutate(
    fixed_effects = b_Intercept +
      (bandInt * b_bandInt) +
      condit_dummy * b_conditVaried + 
      (bandInt * condit_dummy) * `b_conditVaried:bandInt` ,  # Note: Replace : with correct operator
    random_effects = Intercept + bandInt_RF,  # Assuming random effects for intercept and bandInt
    estimate = fixed_effects + random_effects
  )





## Method 3 - epred
nd3 <- e1_vxBMM %>%
  epred_draws(newdata = new_data_grid,ndraws = 1000) %>%
  as_tibble() |> 
  select(-.chain,-.iteration) |>
  left_join(top_indv_coefs, by=join_by(id,condit)) |>
  left_join(testAvg,by=join_by(id,condit,bandInt)) 
 



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
