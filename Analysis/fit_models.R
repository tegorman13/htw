
source(here::here("Functions", "packages.R"))
options(brms.backend="cmdstanr",mc.cores=4)


test <- readRDS(here("data/e1_08-04-23.rds")) |> filter(expMode2 == "Test") 
testE2 <- readRDS(here("data/e2_08-04-23.rds")) |> filter(expMode2 == "Test") 
testE3 <- readRDS(here("data/e3_08-04-23.rds")) |> filter(expMode2 == "Test") 



# e1Sbjs <- test |> group_by(id,condit) |> summarise(n=n())
# testAvg <- test %>% group_by(id, condit, vb, bandInt,bandType,tOrder) %>%
#   summarise(nHits=sum(dist==0),vx=mean(vx),dist=mean(dist),sdist=mean(sdist),n=n(),Percent_Hit=nHits/n)


modelName <- "e1_testDistBand_RF_5K_Ml1"
e1_distBMM <- brm(dist ~ condit * bandInt + (1 + bandInt|id),
                      data=test,file=paste0(here::here("data/model_cache",modelName)),
                      iter=5000,chains=4)
mp1 <- GetModelStats(e1_distBMM) |> kable(escape=F,booktabs=T)
mp1

e1_vxBMM <- brm(vx ~ condit * bandInt + (1 + bandInt|id),
                        data=test,file=paste0(here::here("data/model_cache", "e1_testVxBand_RF_5k_Ml1")),
                        iter=5000,chains=4,silent=0,
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
