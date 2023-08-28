
source(here::here("Functions", "packages.R"))
options(brms.backend="cmdstanr",mc.cores=4)


test <- readRDS(here("data/e1_08-21-23.rds")) |> filter(expMode2 == "Test") 
testE2 <- readRDS(here("data/e2_08-04-23.rds")) |> filter(expMode2 == "Test") 
testE3 <- readRDS(here("data/e3_08-04-23.rds")) |> filter(expMode2 == "Test") 

prior = c(prior(normal(30,200),lb=0,class= Intercept))


# e1Sbjs <- test |> group_by(id,condit) |> summarise(n=n())
# testAvg <- test %>% group_by(id, condit, vb, bandInt,bandType,tOrder) %>%
#   summarise(nHits=sum(dist==0),vx=mean(vx),dist=mean(dist),sdist=mean(sdist),n=n(),Percent_Hit=nHits/n)




######## Testing Models Fit to All 6 Bands

modelName <- "e1_testDistBand_RF_5K"
e1_distBMM <- brm(dist ~ condit * bandInt + (1 + bandInt|id),
                      data=test,file=paste0(here::here("data/model_cache",modelName)),
                      iter=5000,chains=4, prior=prior)


modelName <- "e1_testVxBand_RF_5k"
e1_vxBMM <- brm(vx ~ condit * bandInt + (1 + bandInt|id),
                        data=test,file=paste0(here::here("data/model_cache", modelName)),
                        iter=5000,chains=4,silent=0,, prior=prior, 
                        control=list(adapt_delta=0.94, max_treedepth=13))



modelName <- "e2_testDistBand_RF_5K"
e2_distBMM <- brm(dist ~ condit * bandInt + (1 + bandInt|id),
                      data=testE2,file=paste0(here::here("data/model_cache",modelName)),
                      iter=5000,chains=4)

e2_vxBMM <- brm(vx ~ condit * bandInt + (1 + bandInt|id),
                        data=testE2,file=paste0(here::here("data/model_cache", "e2_testVxBand_RF_5k")),
                        iter=5000,chains=4,silent=0,
                        control=list(adapt_delta=0.90, max_treedepth=11))



modelName <- "e3_testDistBand_RF_5K"
e3_distBMM <- brm(dist ~ condit * bandInt + (1 + bandInt|id),
                      data=testE3,file=paste0(here::here("data/model_cache",modelName)),
                      iter=5000,chains=4)


e3_vxBMM <- brm(vx ~ condit * bandInt + (1 + bandInt|id),
                        data=testE3,file=paste0(here::here("data/model_cache", "e3_testVxBand_RF_5k")),
                        iter=5000,chains=4,silent=0,
                        control=list(adapt_delta=0.90, max_treedepth=11))






######## Testing Models Fit to 3 Extrapolation Bands


modelName <- "e1_extrap_testDistBand"
e1_extrap_distBMM <- brm(dist ~ condit * bandInt + (1 + bandInt|id),
                  data=test |> filter(expMode=="test-Nf"),file=paste0(here::here("data/model_cache",modelName)),
                  iter=5000,chains=4, prior=prior)

modelName <- "e1_extrap_testVxBand"
e1_extrap_VxBMM <- brm(vx ~ condit * bandInt + (1 + bandInt|id),
                  data=test |> filter(expMode=="test-Nf"),file=paste0(here::here("data/model_cache",modelName)),
                  iter=5000,chains=4, prior=prior)




modelName <- "e2_extrap_testDistBand"
e2_extrap_distBMM <- brm(dist ~ condit * bandInt + (1 + bandInt|id),
                  data=testE2 |> filter(expMode=="test-Nf"),file=paste0(here::here("data/model_cache",modelName)),
                  iter=5000,chains=4)

modelName <- "e2_extrap_testVxBand"
e2_extrap_VxBMM <- brm(vx ~ condit * bandInt + (1 + bandInt|id),
                  data=testE2 |> filter(expMode=="test-Nf"),file=paste0(here::here("data/model_cache",modelName)),
                  iter=5000,chains=4)


modelName <- "e3_extrap_testDistBand"
e3_extrap_distBMM <- brm(dist ~ condit * bandInt + (1 + bandInt|id),
                  data=testE3 |> filter(expMode=="test-Nf"),file=paste0(here::here("data/model_cache",modelName)),
                  iter=5000,chains=4)

modelName <- "e3_extrap_testVxBand"
e3_extrap_VxBMM <- brm(vx ~ condit * bandInt + (1 + bandInt|id),
                  data=testE3 |> filter(expMode=="test-Nf"),file=paste0(here::here("data/model_cache",modelName)),
                  iter=5000,chains=4)






######## TRUNC Testing Models Fit to All 6 Bands


prior = c(prior(normal(30,200),lb=0,class= Intercept))
modelName <- "e1_trncInt_testVxBand"
e1_trncInt_vxBMM <-  brm(vx ~ condit * bandInt + (1 + bandInt|id),
                        data=test,
                        file=paste0(here::here("data/model_cache", modelName)),
                        iter=1000,chains=4,silent=0, prior=prior,
                        control=list(adapt_delta=0.91, max_treedepth=11))



###### No interaction

modelName <- "e1_conditPlusBand_RF"
e1_vxB_CpB<- brm(vx ~ condit + bandInt + (1 + bandInt|id),
                        data=test,file=paste0(here::here("data/model_cache", modelName)),
                        iter=2000,chains=4,silent=0, prior=prior, 
                        control=list(adapt_delta=0.92, max_treedepth=11))

modelName <- "e1_conditPlusBand_RF_gr"
e1_vxB_CpB<- brm(vx ~ condit + bandInt + (1 + bandInt|gr(id,by=condit)),
                        data=test,file=paste0(here::here("data/model_cache", modelName)),
                        iter=2000,chains=4,silent=0, prior=prior, 
                        control=list(adapt_delta=0.92, max_treedepth=11))



#### Grouped RF models
# - barely any difference from non-grouped models
#   - bf(vx ~ condit * bandInt + (1 + bandInt|gr(id, by = condit)))


modelName <- "e1_testVxBand_grRF"
b.eq <- bf(vx ~ condit * bandInt + (1 + bandInt|gr(id, by = condit)))
e1_gr_vxBMM <- brm(b.eq, data=test,file=paste0(here::here("data/model_cache", modelName)),
                   iter=5000,chains=4,silent=0,, prior=prior, 
                   control=list(adapt_delta=0.94, max_treedepth=13))


modelName <- "e1_testVxCondit_grRF2"
b.eq2 <- bf(vx ~ condit  + (bandInt||condit) + (1 + bandInt||id))
e1_gr2_vxBMM <- brm(b.eq2, data=test,file=paste0(here::here("data/model_cache", modelName)),
                   iter=2000,chains=3,silent=0, prior=prior, 
                   control=list(adapt_delta=0.90, max_treedepth=11))



modelName <- "e1_testVxCondit_grRF3"
b.eq3 <- bf(vx ~ condit  + (0 + bandInt|condit) + (0 + bandInt|gr(id, by = condit)))
e1_gr3_vxBMM <- brm(b.eq3, data=test,file=paste0(here::here("data/model_cache", modelName)),
                   iter=2000,chains=3,silent=0, prior=prior, 
                   control=list(adapt_delta=0.91, max_treedepth=12))



modelName <- "e1_testVxCondit_grRF4"
b.eq4 <- bf(vx ~ condit  + (1 + bandInt|gr(id, by = condit)))
e1_gr4_vxBMM <- brm(b.eq4, data=test,file=paste0(here::here("data/model_cache", modelName)),
                   iter=2000,chains=3,silent=0, prior=prior, 
                   control=list(adapt_delta=0.91, max_treedepth=12))



modelName <- "e1_testVxCondit_grRF5"
b.eq5 <- bf(vx ~ 0+condit*bandInt  + (1 + bandInt||gr(id, by = condit)))
e1_gr5_vxBMM <- brm(b.eq5, data=test,file=paste0(here::here("data/model_cache", modelName)),
                   iter=2000,chains=3,silent=0,  
                   control=list(adapt_delta=0.91, max_treedepth=12))


modelName <- "e1_testVxCondit_grRF6"
b.eq6 <- bf(vx ~ condit*bandInt  + (0 + bandInt||gr(id, by = condit)))
e1_gr6_vxBMM <- brm(b.eq6, data=test,file=paste0(here::here("data/model_cache", modelName)),
                   iter=2000,chains=3,silent=0,  
                   control=list(adapt_delta=0.91, max_treedepth=12))


modelName <- "e1_testVxCondit_grRF7"
b.eq7 <- bf(vx ~ condit*bandInt  + (bandInt|condit) + ( 1 + bandInt||gr(id, by = condit)))
e1_gr7_vxBMM <- brm(b.eq7, data=test,file=paste0(here::here("data/model_cache", modelName)),
                   iter=2000,chains=3,silent=0,  
                   control=list(adapt_delta=0.91, max_treedepth=12))



modelName <- "e1_testVx_grRF8"
b.eq8 <- bf(vx ~ bandInt + ( 1 + bandInt||gr(id, by = condit)))
e1_gr8_vxBMM <- brm(b.eq8, data=test,file=paste0(here::here("data/model_cache", modelName)),
                   iter=2000,chains=3,silent=0,  
                   control=list(adapt_delta=0.91, max_treedepth=12))






#################


modelName <- "e1__condRF_vx"
e1_condRF_vx <- brm(vx ~ condit + (1|bandInt) + (1 + bandInt||condit) + (1 + bandInt||id),
                  data=test,file=paste0(here::here("data/model_cache",modelName)),
                  iter=5000,chains=4)


modelName <- "e1__condRF_vx2"
e1_condRF_vx <- brm(vx ~ condit + (1+condit|bandInt) + (1 + bandInt||id),
                    data=test,file=paste0(here::here("data/model_cache",modelName)),
                    iter=2000,chains=3)

modelName <- "e1__condRF_vx3"
e1_condRF_vx3 <- brm(vx ~ condit + (1|bandInt) + (1 + bandInt||id),
                    data=test,file=paste0(here::here("data/model_cache",modelName)),
                    iter=2000,chains=3,silent=0)

modelName <- "e1__condRF_vx3_gr"
e1_condRF_vx3_gr <- brm(vx ~ condit + (1|bandInt) + (1 + bandInt||gr(id,by=condit)),
                     data=test,file=paste0(here::here("data/model_cache",modelName)),
                     iter=2000,chains=3,silent=0)





######## No Correlation - Testing Models Fit to 3 Extrapolation Bands


modelName <- "e1_NC_extrap_testDistBand"
e1_NC_extrap_distBMM <- brm(dist ~ condit * bandInt + (1 + bandInt||id),
                  data=test |> filter(expMode=="test-Nf"),file=paste0(here::here("data/model_cache",modelName)),
                  iter=5000,chains=4)

modelName <- "e1_NC_extrap_testVxBand"
e1_NC_extrap_VxBMM <- brm(vx ~ condit * bandInt + (1 + bandInt||id),
                  data=test |> filter(expMode=="test-Nf"),file=paste0(here::here("data/model_cache",modelName)),
                  iter=5000,chains=4)


######## No Correlation - Alt Spec - Testing Models Fit to 3 Extrapolation Bands
# - identical to specifying models with || 
modelName <- "e1_NC2_extrap_testVxBand"
e1_NC2_extrap_VxBMM <- brm(vx ~ condit * bandInt + (0 + bandInt|id) + (1|id),
                  data=test |> filter(expMode=="test-Nf"),file=paste0(here::here("data/model_cache",modelName)),
                  iter=5000,chains=4)



###########
# GAMS
############


library(mgcv)
mod_lm = gam(vx ~ 1 + bandInt, data = test)

mod_gam1 = gam(vx ~ condit+ s(bandInt, bs = "cr",k=4) + s(id, bandInt,bs='re'), data = test)
summary(mod_gam1)
plot(ggeffects::ggpredict(mod_gam1))
gratia::draw(mod_gam1)
gam.check(mod_gam1)


mod_gam2 = gam(vx ~ condit * s(bandInt, bs = "cr",k=4) + s(id, bandInt,bs='re'), data = test)
summary(mod_gam2)
plot(ggeffects::ggpredict(mod_gam1))
gratia::draw(mod_gam1)

brmGam <- brm(
  vx ~ condit + 
    s(bandInt,k=4) +
    s(bandInt, id, bs="re",k=4) + 
    (1 + bandInt || id), 
  data = test, 
  family = gaussian()
)

coef(brmGam)$id %>% as_tibble(rownames="id") 



brmGam2 <- brm(
  vx ~ condit + 
    s(bandInt,k=4) +
    s(bandInt, id, bs="re",k=4) + 
    s(id,bs="re") + 
    (1 + bandInt || id), 
  data = test, 
  family = gaussian()
)


coef(brmGam2)$id %>% as_tibble(rownames="id") 


brmGam2 <- brm(
  vx ~  
    s(bandInt,condit,k=2) +
    s(bandInt, id, bs="re",k=22) + 
    (1 + bandInt || id), 
  data = test, 
  family = gaussian()
)

gam1 <- brm(vx ~ s(bandInt, k=4) + (1 + bandInt | id), 
    data = test, 
    family = gaussian())

gam1

gam5 <- brm(vx ~ s(bandInt, k=5) + (1 + bandInt | id), 
            data = test, 
            family = gaussian())

gam5 <- brm(vx ~ s(bandInt, k=3) + (1 + bandInt | id), 
            data = test, 
            family = gaussian())

gam1


cr1 <- brm(vx ~ s(bandInt, bs="cr", k=5) + (1 + bandInt | id),
    data = test, 
    family = gaussian())



cr1 <- brm(vx ~ s(bandInt, bs="cr", k=1) + (1 + bandInt | id),
           data = test, 
           family = gaussian())




cr2 <- brm(vx ~ s(bandInt, bs="cr", k=2) + (1 + bandInt | id),
           data = test, 
           family = gaussian())




indv_coefs <- coef(gam1)$id |> 
  as_tibble(rownames="id") |> 
  select(id, starts_with("Est")) |>
  left_join(e1Sbjs, by=join_by(id) ) |> 
  group_by(condit) |> 
  mutate(rank = rank(desc(Estimate.bandInt))) 



