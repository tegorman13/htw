source(here::here("Functions", "packages.R"))

test <- readRDS(here("data/e1_08-21-23.rds")) |>  filter(expMode2 == "Test") |>
  select(id,condit,bandInt,vb,vx,dist,sdist,bandType)



get_preds <- function(model) {
  vx_pred <- 
    posterior_predict(model, ndraws = 500) |> 
    array_branch(margin=1) |> 
    map_dfr( 
      function(yrep_iter) {
        test  |>
          mutate(vx_pred = yrep_iter)
      },
      .id = 'iter'
    ) |>
    mutate(iter = as.numeric(iter))
  
  
  vx_pred_agg <- vx_pred |> group_by(id,condit,bandInt) |> 
    mutate(resid=vx-vx_pred) |> summarise(vx=mean(vx),vx_pred=mean(vx_pred),resid=mean(resid))
  
  (vx_pred_agg |> group_by(condit) |> summarise(resid=mean(resid)) |> pandoc.table())
  vx_pred_agg |> group_by(condit) |> summarise(resid=mean(abs(resid))) |> pandoc.table()
  
  testAvgResid <- testAvg |> select(id,condit,vx,vb,bandInt) |> 
    left_join(vx_pred |> group_by(id,bandInt) |> median_qi(vx_pred), by=join_by(id,bandInt)) |>
    mutate(resid=vx-vx_pred) |> relocate(vx,resid,.before=vx_pred)
  
  (testAvgResid |> group_by(condit) |> summarise(resid=mean(resid)))
  return(list(trialPred=vx_pred_agg,avgPred=testAvgResid))
}

c.t.b <- get_preds(e1_vxBMM)

c.t.b |> map(head)

c.b.gr <- get_preds(e1_testVx_grRF8)
c.t.b.gr <- get_preds(e1_testVxBand_grRF)



vx_pred |> ggplot(aes(x=vx_pred,fill=condit)) + geom_density() + facet_wrap(~vb)
vx_pred |> ggplot(aes(x=vx_pred,fill=condit)) + geom_histogram() + facet_wrap(~vb)

vx_pred |> ggplot(aes(x=vx,fill=condit)) + geom_density() + facet_wrap(~vb)
vx_pred |> ggplot(aes(x=vx,fill=condit)) + geom_histogram() + facet_wrap(~vb)


vx_pred |> ggplot(aes(x=(vx-vx_pred),fill=condit)) + geom_density(alpha=.4) + facet_wrap(~vb)
vx_pred |> ggplot(aes(x=(vx-vx_pred),fill=condit)) + geom_histogram(alpha=.4) + facet_wrap(~vb)




vx_pred <- 
  posterior_predict(e1_vxBMM, ndraws = 500) |> 
  array_branch(margin=1) |> 
  map_dfr( 
    function(yrep_iter) {
      test  |>
        mutate(vx_pred = yrep_iter)
    },
    .id = 'iter'
  ) |>
  mutate(iter = as.numeric(iter))


vx_pred_agg <- vx_pred |> group_by(id,condit,bandInt) |> 
  mutate(resid=vx-vx_pred) |> summarise(vx=mean(vx),vx_pred=mean(vx_pred),resid=mean(resid))

vx_pred_agg |> group_by(condit) |> summarise(resid=mean(resid))



#####################
 
modelName <- "e1_testVxCondit_grRF2"
b.eq2 <- bf(vx ~ condit  + (bandInt||condit) + (1 + bandInt||id))
e1_gr2_vxBMM <- brm(b.eq2, data=test,file=paste0(here::here("data/model_cache", modelName)),
                      iter=2000,chains=3,silent=0, prior=prior, 
                      control=list(adapt_delta=0.90, max_treedepth=11))

new_data_grid=map_dfr(1, ~data.frame(unique(test[,c("id","condit")]))) |> 
  dplyr::arrange(id) |> 
  mutate(condit_dummy = ifelse(condit == "Varied", 1, 0)) 

random_effectsCondit <- e1_gr2_vxBMM |> 
  gather_draws(`^r_condit.*$`, regex = TRUE, ndraws = 1000) |> 
  separate(.variable, into = c("effect", "condit", "term"), sep = "\\[|,|\\]") |> 
  mutate(condit = factor(condit,levels=levels(test$condit))) |> 
  pivot_wider(names_from = term, values_from = .value) |> arrange(condit,.chain,.draw,.iteration)


random_effectsCondit |> group_by(condit) |> median_hdi(bandInt)  |> pandoc.table()
ranef(e1_gr2_vxBMM)$condit[, ,"bandInt"]





#### Indvidual Slopes
indv_coefs <- coef(e1_gr2_vxBMM)$id |> 
  as_tibble(rownames="id") |> 
  select(id, starts_with("Est")) |>
  left_join(e1Sbjs, by=join_by(id) ) 

fixed_effects <- e1_gr2_vxBMM |> 
  spread_draws(`^b_.*`,regex=TRUE) |> arrange(.chain,.draw,.iteration)

random_effects <- e1_gr2_vxBMM |> 
  gather_draws(`^r_id.*$`, regex = TRUE, ndraws = 1000) |> 
  separate(.variable, into = c("effect", "id", "term"), sep = "\\[|,|\\]") |> 
  mutate(id = factor(id,levels=levels(test$id))) |> 
  pivot_wider(names_from = term, values_from = .value) |> arrange(id,.chain,.draw,.iteration)


indvDraws <- left_join(random_effects, fixed_effects, by = join_by(".chain", ".iteration", ".draw")) |> 
  right_join(new_data_grid, by = join_by("id")) 

indvSlopes <- indvDraws |> group_by(id) |> median_qi(bandInt)  |> 
  left_join(e1Sbjs, by=join_by(id)) |> group_by(condit) |>
  mutate(rankSlope=rank(bandInt)) |> arrange(rankSlope)  

indvSlopes |> reframe(enframe(quantile(bandInt, c(0.0,0.25, 0.5, 0.75,1)), "quantile", "bandInt")) |> 
  pivot_wider(names_from=quantile,values_from=bandInt,names_prefix="Q_") |>
  group_by(condit) |>
  summarise(across(starts_with("Q"), list(mean = mean))) |> pandoc.table() 
  

indvDraws |> group_by(condit) |> median_hdi() |> pandoc.table() 



ggplot(indvSlopes, aes(x = bandInt, color = condit)) + 
  geom_density()



indvSlopes |> group_by(condit) |> 
  summarise(med_slope = median(bandInt), sd_slope = sd(bandInt) ) |> pandoc.table()

indvDraws |> group_by(id) |> median_hdi(bandInt)  |> 
  left_join(e1Sbjs, by=join_by(id)) |> group_by(condit) |>
  median_qi(bandInt) |> pandoc.table()


indvDraws |> group_by(condit) |> median_hdi(bandInt)


vx_pred <- 
  posterior_predict(e1_gr2_vxBMM, ndraws = 500) |> 
  array_branch(margin=1) |> 
  map_dfr( 
    function(yrep_iter) {
      test  |>
        mutate(vx_pred = yrep_iter)
    },
    .id = 'iter'
  ) |>
  mutate(iter = as.numeric(iter))


vx_pred  |>
  filter(iter < 100) %>%
  ggplot(aes(vx_pred, group = iter)) +
  geom_line(alpha = .03, stat = 'density', color = 'blue') +
  geom_density(data = test,
               aes(vx,col=vb),
               inherit.aes = FALSE,
               size = 0.7) + # 1
  facet_grid(condit ~ vb) +
  xlab('Vx')


bayes_R2(e1_gr2_vxBMM)



bayes_R2(e1_testVxCondit_grRF4)

indv_coefs <- coef(e1_testVxCondit_grRF4)$id |> 
  as_tibble(rownames="id") |> 
  select(id, starts_with("Est")) |>
  left_join(e1Sbjs, by=join_by(id) ) 

fixed_effects <- e1_testVxCondit_grRF4 |> 
  spread_draws(`^b_.*`,regex=TRUE) |> arrange(.chain,.draw,.iteration)

random_effects <- e1_testVxCondit_grRF4 |> 
  gather_draws(`^r_id.*$`, regex = TRUE, ndraws = 1000) |> 
  separate(.variable, into = c("effect", "id", "term"), sep = "\\[|,|\\]") |> 
  mutate(id = factor(id,levels=levels(test$id))) |> 
  pivot_wider(names_from = term, values_from = .value) |> arrange(id,.chain,.draw,.iteration)


indvDraws <- left_join(random_effects, fixed_effects, by = join_by(".chain", ".iteration", ".draw")) |> 
  right_join(new_data_grid, by = join_by("id")) 

indvSlopes <- indvDraws |> group_by(id) |> median_qi(bandInt)  |> 
  left_join(e1Sbjs, by=join_by(id)) |> group_by(condit) |>
  mutate(rankSlope=rank(bandInt)) |> arrange(rankSlope)  

indvSlopes |> reframe(enframe(quantile(bandInt, c(0.0,0.25, 0.5, 0.75,1)), "quantile", "bandInt")) |> 
  pivot_wider(names_from=quantile,values_from=bandInt,names_prefix="Q_") |>
  group_by(condit) |>
  summarise(across(starts_with("Q"), list(mean = mean))) |> pandoc.table() 

indvDraws |> group_by(id) |> median_hdi(bandInt)  |> 
  left_join(e1Sbjs, by=join_by(id)) |> group_by(condit) |>
  median_qi(bandInt) |> pandoc.table()





###########

bayes_R2(e1_testVxCondit_grRF6)

indv_coefs <- coef(e1_testVxCondit_grRF6)$id |> 
  as_tibble(rownames="id") |> 
  select(id, starts_with("Est")) |>
  left_join(e1Sbjs, by=join_by(id) ) 

fixed_effects <- e1_testVxCondit_grRF6 |> 
  spread_draws(`^b_.*`,regex=TRUE) |> arrange(.chain,.draw,.iteration)

random_effects <- e1_testVxCondit_grRF6 |> 
  gather_draws(`^r_id.*$`, regex = TRUE, ndraws = 1000) |> 
  separate(.variable, into = c("effect", "id", "term"), sep = "\\[|,|\\]") |> 
  mutate(id = factor(id,levels=levels(test$id))) |> 
  pivot_wider(names_from = term, values_from = .value) |> arrange(id,.chain,.draw,.iteration)


indvDraws <- left_join(random_effects, fixed_effects, by = join_by(".chain", ".iteration", ".draw")) |> 
  right_join(new_data_grid, by = join_by("id")) 

indvSlopes <- indvDraws |> group_by(id) |> median_qi(bandInt)  |> 
  left_join(e1Sbjs, by=join_by(id)) |> group_by(condit) |>
  mutate(rankSlope=rank(bandInt)) |> arrange(rankSlope)  

indvSlopes |> reframe(enframe(quantile(bandInt, c(0.0,0.25, 0.5, 0.75,1)), "quantile", "bandInt")) |> 
  pivot_wider(names_from=quantile,values_from=bandInt,names_prefix="Q_") |>
  group_by(condit) |>
  summarise(across(starts_with("Q"), list(mean = mean))) |> pandoc.table() 

indvDraws |> group_by(id) |> median_hdi(bandInt)  |> 
  left_join(e1Sbjs, by=join_by(id)) |> group_by(condit) |>
  median_qi(bandInt) |> pandoc.table()




###########

bayes_R2(e1_testVxCondit_grRF5)

indv_coefs <- coef(e1_testVxCondit_grRF5)$id |> 
  as_tibble(rownames="id") |> 
  select(id, starts_with("Est")) |>
  left_join(e1Sbjs, by=join_by(id) ) 

fixed_effects <- e1_testVxCondit_grRF5 |> 
  spread_draws(`^b_.*`,regex=TRUE) |> arrange(.chain,.draw,.iteration)

random_effects <- e1_testVxCondit_grRF5 |> 
  gather_draws(`^r_id.*$`, regex = TRUE, ndraws = 1000) |> 
  separate(.variable, into = c("effect", "id", "term"), sep = "\\[|,|\\]") |> 
  mutate(id = factor(id,levels=levels(test$id))) |> 
  pivot_wider(names_from = term, values_from = .value) |> arrange(id,.chain,.draw,.iteration)


indvDraws <- left_join(random_effects, fixed_effects, by = join_by(".chain", ".iteration", ".draw")) |> 
  right_join(new_data_grid, by = join_by("id")) 

indvSlopes <- indvDraws |> group_by(id) |> median_qi(bandInt)  |> 
  left_join(e1Sbjs, by=join_by(id)) |> group_by(condit) |>
  mutate(rankSlope=rank(bandInt)) |> arrange(rankSlope)  

indvSlopes |> reframe(enframe(quantile(bandInt, c(0.0,0.25, 0.5, 0.75,1)), "quantile", "bandInt")) |> 
  pivot_wider(names_from=quantile,values_from=bandInt,names_prefix="Q_") |>
  group_by(condit) |>
  summarise(across(starts_with("Q"), list(mean = mean))) |> pandoc.table() 

indvDraws |> group_by(id) |> median_hdi(bandInt)  |> 
  left_join(e1Sbjs, by=join_by(id)) |> group_by(condit) |>
  median_qi(bandInt) |> pandoc.table()
