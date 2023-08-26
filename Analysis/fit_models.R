
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
mp1 <- GetModelStats(e1_distBMM) |> kable(escape=F,booktabs=T)
mp1

modelName <- "e1_testVxBand_RF_5k"
e1_vxBMM <- brm(vx ~ condit * bandInt + (1 + bandInt|id),
                        data=test,file=paste0(here::here("data/model_cache", modelName)),
                        iter=5000,chains=4,silent=0,, prior=prior, 
                        control=list(adapt_delta=0.94, max_treedepth=13))
mt2 <-GetModelStats(e1_vxBMM ) |> kable(escape=F,booktabs=T)
mt2


modelName <- "e2_testDistBand_RF_5K"
e2_distBMM <- brm(dist ~ condit * bandInt + (1 + bandInt|id),
                      data=testE2,file=paste0(here::here("data/model_cache",modelName)),
                      iter=5000,chains=4)
mp2 <- GetModelStats(e2_distBMM) |> kable(escape=F,booktabs=T)
mp2

e2_vxBMM <- brm(vx ~ condit * bandInt + (1 + bandInt|id),
                        data=testE2,file=paste0(here::here("data/model_cache", "e2_testVxBand_RF_5k")),
                        iter=5000,chains=4,silent=0,
                        control=list(adapt_delta=0.90, max_treedepth=11))
mt3<-GetModelStats(e2_vxBMM ) |> kable(escape=F,booktabs=T)
mt3



modelName <- "e3_testDistBand_RF_5K"
e3_distBMM <- brm(dist ~ condit * bandInt + (1 + bandInt|id),
                      data=testE3,file=paste0(here::here("data/model_cache",modelName)),
                      iter=5000,chains=4)
mp3 <- GetModelStats(e3_distBMM) |> kable(escape=F,booktabs=T)
mp3

e3_vxBMM <- brm(vx ~ condit * bandInt + (1 + bandInt|id),
                        data=testE3,file=paste0(here::here("data/model_cache", "e3_testVxBand_RF_5k")),
                        iter=5000,chains=4,silent=0,
                        control=list(adapt_delta=0.90, max_treedepth=11))
mt4<-GetModelStats(e3_vxBMM ) |> kable(escape=F,booktabs=T)
mt4





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






#### Grouped RF models
# - barely any difference from non-grouped models
#   - bf(vx ~ condit * bandInt + (1 + bandInt|gr(id, by = condit)))


modelName <- "e1_testVxBand_grRF"
b.eq <- bf(vx ~ condit * bandInt + (1 + bandInt|gr(id, by = condit)))
e1_gr_vxBMM <- brm(bf, data=test,file=paste0(here::here("data/model_cache", modelName)),
                   iter=5000,chains=4,silent=0,, prior=prior, 
                   control=list(adapt_delta=0.94, max_treedepth=13))


modelName <- "e1_testVxCondit_grRF2"
b.eq2 <- bf(vx ~ condit  + (bandInt|condit) + (1 + bandInt|gr(id, by = condit)))
e1_gr2_vxBMM <- brm(b.eq2, data=test,file=paste0(here::here("data/model_cache", modelName)),
                   iter=5000,chains=4,silent=0, prior=prior, 
                   control=list(adapt_delta=0.94, max_treedepth=13))





























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



