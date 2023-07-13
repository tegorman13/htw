---
title: "Testing Discrimination Analysis"
subtitle: "Fitting mixed effects models"
date: last-modified
categories: [Analysis, R]
code-fold: true
---


```{r}
#| output: false

lapply(c('tidyverse','data.table','lme4','lmerTest','knitr','kableExtra','cowplot','gghalves'),library,character.only = TRUE)
d <- readRDS('dPrune-01-19-23.rds')

dtest <- d %>% filter(expMode %in% c("test-Nf","test-train-nf")) %>% group_by(id,lowBound) %>% 
  mutate(nBand=n(),band=bandInt,id=factor(id)) %>% group_by(id) %>% mutate(nd=n_distinct(lowBound))
# unique(dtest[dtest$nd==4,]$sbjCode) # 7 in wrong condition
dtest <- dtest %>% group_by(id,lowBound) %>% filter(nBand>=5 & nd==6)
# for any id that has at least 1 nBand >=5, remove all rows with that id. 
dtest <- dtest %>% group_by(id) %>% filter(!id %in% unique(dtest$id[dtest$nBand<5]))

dtestAgg <- dtest %>% group_by(id,condit,catOrder,feedbackType,vb,band,lowBound,highBound,input) %>% mutate(vxCapped=ifelse(vx>1600,1600,vx)) %>%
  summarise(vxMean=mean(vx),devMean=mean(dist),vxMed=median(vx),devMed=median(dist),
            vxMeanCap=mean(vxCapped),.groups = "keep")
            
```

```{r fig.width=11, fig.height=9}
#| eval: false
#| include: false
dtest %>% group_by(id,lowBound,condit,catOrder,expMode) %>% summarise(n = n()) %>% ggplot(aes(x = n)) + geom_histogram(aes(fill = condit)) + facet_wrap(catOrder+expMode~lowBound,ncol=6)

dtest %>% group_by(id,lowBound,condit,catOrder,expMode) %>% summarise(n = n()) %>% mutate(nf=factor(n,levels=unique(n))) %>%
  group_by(condit,expMode,catOrder,nf,n) %>% summarise(c=n()) %>% arrange(n) %>% select(-nf) %>% DT::datatable()

dtest %>% group_by(id,lowBound,condit,catOrder) %>% summarise(n = n()) %>% group_by(lowBound,condit,catOrder) %>% summarise(mean = mean(n),sd = sd(n),n = n()) 

# print number of subjects in each condition combination with less than 5 trials in any lowBound. 
dtest %>% group_by(id,lowBound,condit,catOrder) %>% summarise(n = n()) %>% filter(n < 5) %>% DT::datatable()
dtest %>% group_by(id,expMode,catOrder,lowBound) %>% select(nTestF,nTestNf,nBand) %>% slice(1) %>% arrange(nBand)

```



## Quick Reminder of General Patterns

```{r fig.width=12, fig.height=10}
#| column: page-inset-right

fig1aCap=str_wrap("Figure 1a: Bands 100-300, 350-550 and 600-800 are novel extrapolations for both Original Order. Translucent rectangles indicate the correct band " ,width=170)

fig1bCap=str_wrap("Figure 1b: Bands 800-1000, 1000-1200,  and 1200-1400 are novel extrapolations for both Reverse Order. Translucent rectangles indicate the correct band " ,width=170)

plotDist <- function(df,title="",fcap=""){
  rectWidth=30
  df %>%ggplot()+aes(x = band, y = vxMeanCap, fill=vb) +
    # Set the color mapping in this layer so the points don't get a color
   geom_half_violin(color=NA)+ # remove border color
  geom_half_boxplot(position=position_nudge(x=-0.05),side="r",outlier.shape = NA,center=TRUE,
                    errorbar.draw = FALSE,width=20)+
  geom_half_point(transformation = position_jitter(width = 0.05, height = 0.05),size=.3,aes(color=vb))+
  facet_wrap(~condit,scale="free_x")+
    geom_rect(aes(xmin=band-rectWidth,xmax=band+rectWidth,ymin=band,ymax=highBound,fill=vb),alpha=.01)+
    geom_segment(aes(x=band-rectWidth,xend=band+rectWidth,y=highBound,yend=highBound),alpha=.8,linetype="dashed")+
    geom_segment(aes(x=band-rectWidth,xend=band+rectWidth,y=band,yend=band),alpha=.8,linetype="dashed")+
    labs(x = "Velocity Band", y = "vxMean",caption=fcap) +
    scale_y_continuous(expand=expansion(add=100),breaks=round(seq(0,2000,by=200),2))+
    scale_x_continuous(labels=sort(unique(df$band)),breaks=sort(unique(df$band)))+
    ggtitle(title) + theme(legend.position = "none")+theme_classic()+guides(fill="none",color="none")+
  theme(plot.caption=element_text(hjust=0,face="italic"))
}


#dtestAgg %>% plotDist()
dtestAgg %>% filter(catOrder=="orig") %>% plotDist(title="Empirical Vx - Original Order",fcap=fig1aCap)
dtestAgg %>% filter(catOrder=="rev") %>% plotDist(title="Empirical Vx - Reverse Order",fcap=fig1bCap)

```










## naive model that fits single slope and intercept to all subjects
```{r}
# Fit a model on all the data pooled together
m_pooled <- lm(vxMean ~ band, dtestAgg) 
# Repeat the intercept and slope terms for each participant
df_pooled <- tibble(
  Model = "Complete pooling",
  id = unique(dtestAgg$id),
  Intercept = coef(m_pooled)[1], 
  Slope_band = coef(m_pooled)[2]
)
#head(df_pooled)

# print the coefficents and residual of the model
summary(m_pooled)


```



## Fit no pooling model (individual fit for each subject)
```{r}

df_no_pooling <- lmList(vxMean ~ band | id, dtestAgg) %>% 
  coef() %>% rownames_to_column("id") %>% 
  rename(Intercept = `(Intercept)`, Slope_band = band) %>% 
  add_column(Model = "No pooling")

# print average coefficients and average residual for the model
summary(df_no_pooling)
# print average residual of no pooling model
summary(df_no_pooling$vxMean ~ df_no_pooling$band | df_no_pooling$id)

# sort the dataframe by the value of slope_band, highest to lowest
testSlopeIndv <- df_no_pooling %>% arrange(desc(Slope_band))

# Add a condit column to the dataframe, matching condition based on the value in dtestAgg for each sbjCode
testSlopeIndv <- testSlopeIndv %>% 
  left_join(dtestAgg %>% ungroup() %>% select(id, condit) %>% distinct(), by = "id") 

# Add a rank column to the dataframe, based on the value of slope_band. Smallest rank for highest value.
testSlopeIndv <- testSlopeIndv %>% group_by(condit) %>% 
  mutate(nGrp=n(),rank = nGrp -rank(Slope_band) +1,
         quantile = cut(rank, breaks = 4, labels = c("1st", "2nd", "3rd", "4th")),
         quintile=cut(rank,breaks=5,labels=c("1st", "2nd", "3rd", "4th","5th")),
         decile=cut(rank,breaks=10,labels=c(1:10))) %>% #select(-n)%>%
  arrange(rank)

# Reorder the sbjCode column so that the sbjCode with the highest slope_band is first
testSlopeIndv$id <- factor(testSlopeIndv$id, levels = testSlopeIndv$id)

#head(testSlopeIndv)
```


## Some individual plots showing the best fitting line against testing behavior (x velocity). 
-  Sample of high, and low discriminating subjects (i.e. highest and lowest slopes)
- Mean Vx for each band shown via dot. 
- correct bands shown with translucent rectangles

```{r fig.width=13, fig.height=11}
#| column: page-right
#| warning: false
# create plotting function that takes in a dataframe, and returns ggplot object
#rewrite plotSlope function to take line color as a function argument, and set the color of abline to that argument

plotSlope <- function(df,title="",colour=NULL){
  rectWidth=50
  df %>%ggplot()+aes(x = band, y = vxMean) +
    # Set the color mapping in this layer so the points don't get a color
    geom_abline(
      aes(intercept = Intercept, slope = Slope_band),
      size = .75,colour=colour,alpha=.2
    ) +geom_point(aes(color=vb)) +facet_wrap("id") +
    geom_rect(aes(xmin=band-rectWidth,xmax=band+rectWidth,ymin=band,ymax=highBound,fill=vb),alpha=.1)+
    geom_segment(aes(x=band-rectWidth,xend=band+rectWidth,y=highBound,yend=highBound),alpha=1,linetype="dashed")+
    geom_segment(aes(x=band-rectWidth,xend=band+rectWidth,y=band,yend=band),alpha=1,linetype="dashed")+
    labs(x = "Velocity Band", y = "vxMean") +
    scale_x_continuous(labels=sort(unique(df$band)),breaks=sort(unique(df$band)))+
    ggtitle(title) + theme(legend.position = "none")+theme_classic()+guides(fill="none",color="none")
}

tv<-testSlopeIndv %>% left_join(dtestAgg, by = c("id","condit")) %>% filter(condit=="Varied",rank<=6) %>% 
   plotSlope(.,colour="black",title="Largest Individually fit Varied Sbj. Slopes")
tc<-testSlopeIndv %>% left_join(dtestAgg, by = c("id","condit")) %>% filter(condit=="Constant",rank<=6) %>% 
   plotSlope(.,colour="black",title="Largest Individually fit Constant Sbj. Slopes")
bv<-testSlopeIndv %>% left_join(dtestAgg, by = c("id","condit")) %>% filter(condit=="Varied",rank>=nGrp-5) %>% 
   plotSlope(.,colour="black",title="Smallest Varied Sbj. Slopes")
bc<-testSlopeIndv %>% left_join(dtestAgg, by = c("id","condit")) %>% filter(condit=="Constant",rank>=nGrp-5) %>% 
   plotSlope(.,colour="black",title="Smallest Constant Sbj. Slopes.")
 
title = ggdraw()+draw_label("Highest and Lowest Slope Values",fontface = 'bold',x=0,hjust=0)+theme(plot.margin = margin(0, 0, 0, .5))
plot_grid(title,NULL,tv,tc,bv,bc,NULL,ncol=2,rel_heights = c(.1,1,1))

```





```{r}
#| eval: false
#| include: false

# m <- lmer(vxMed ~ 1 + band + (1 + band | id), dtestAgg)
# arm::display(m)
# 
# m2 <- lmer(vxMed ~ 1 + input + (1 + input | id), dtestAgg)
# arm::display(m2)

# fit model with more iterations and different optimizer

# m2 <- lmer(vxMed ~ 1 + band + (1 + band | id), dtestAgg, control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 3e5)))

# 
# bm1 <- lmer(vxMed ~ 1 + band + (1 + band | id), dtestAgg, control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 3e5)))
# arm::display(bm1)

# bm2 <- lmer(vxMed ~ 1 + input + (1 + input | id), dtestAgg, control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 3e5)))
# arm::display(bm2)
# 
# bm3 <- lmer(log(vxMean) ~ 1 + log(band) + (1 + log(band) | id), dtestAgg, control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 3e5)))
# arm::display(bm3)
```



## Fit partial pooling model (linear mixed model with random slope and intercept)
```{r}
#| warning: false
#| 
bm1 <- lmer(vxMed ~ 1 + band + (1 + band | id), dtestAgg, control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 3e5)))
arm::display(bm1)

df_partial_pooling <- coef(bm1)[["id"]] %>% 
  rownames_to_column("id") %>% 
  as_tibble() %>% 
  rename(Intercept = `(Intercept)`, Slope_band = band) %>% 
  add_column(Model = "Partial pooling")

head(df_partial_pooling)
summary(bm1)
```


```{r}
#| eval: false
#| include: false
#| echo: false
#| 
gbm <- lmer(vxMed ~ 1 + condit + band + (1 + band | id), dtestAgg, control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 3e5)))
arm::display(gbm)
anova(gbm)
summary(gbm)

gbd <- lmer(devMean ~ 1 + condit + band + (1 + band | id), dtestAgg, control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 3e5)))
arm::display(gbd)
anova(gbd)
summary(gbd)




gbmc <- lmer(vxMed ~ 1 + condit + band + catOrder + (1 + band | id), dtestAgg, control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 3e5)))
arm::display(gbmc)
anova(gbmc)
summary(gbmc)

gbdc <- lmer(devMean ~ 1 + condit + band + catOrder + (1 + band | id), dtestAgg, control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 3e5)))
arm::display(gbdc)
anova(gbdc)
summary(gbdc)


gbm2 <- dtestAgg %>% filter(catOrder=="orig") %>% lmer(vxMed ~ 1 + condit + band + (1 + band | id), data=., control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 3e5)))
arm::display(gbm2)
anova(gbm2)
summary(gbm2)

gbd2 <- dtestAgg %>% filter(catOrder=="orig") %>% lmer(devMean ~ 1 + condit + band + (1 + band | id), data=., control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 3e5)))
arm::display(gbd2)
anova(gbd2)
summary(gbd2)




gbm2.i <- dtestAgg %>% filter(catOrder=="orig") %>% lmer(vxMed ~ 1 + (condit * band) + (1 + band | id), data=., control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 3e5)))
arm::display(gbm2.i)
anova(gbm2.i)
summary(gbm2.i)

gbd2.i <- dtestAgg %>% filter(catOrder=="orig") %>% lmer(devMean ~ 1 + (condit * band) + (1 + band | id), data=., control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 3e5)))
arm::display(gbd2.i)
anova(gbd2.i)
summary(gbd2.i)

```



```{r fig.width=12, fig.height=9}
#| eval: false
#| include: false
#| 
#| column: page-inset-right
df_models <- bind_rows(df_pooled, df_no_pooling, df_partial_pooling) %>% 
  left_join(dtestAgg, by = c("id"))

# use base R to filter only testSlopeIndv with Rank < 3
testSlopeIndv[testSlopeIndv$rank<10,]$id

df_models %>% filter(id %in% testSlopeIndv[testSlopeIndv$rank<10,]$id) %>%
ggplot() + 
  aes(x = band, y = vxMed) + 
  # Set the color mapping in this layer so the points don't get a color
  geom_abline(
    aes(intercept = Intercept, slope = Slope_band, color = Model),
     linewidth= .75
  ) + 
  geom_point() +
  facet_wrap("id") +
  labs(x = "band", y = "vxMean") + 
  # Fix the color palette 
  scale_color_brewer(palette = "Dark2") + 
  theme(legend.position = "top", legend.justification = "left")
```



```{r fig.width=12,fig.height=13}
#| column: page-inset-right
#| message: false
# filter to only retain the no pooling and partial pooling models. 
# Compare the average slope and intercepts between constant and varied condits. Use barplots with standard error bars

df_models <- bind_rows(df_pooled, df_no_pooling, df_partial_pooling) %>% 
  left_join(dtestAgg, by = c("id"))

grpAvg<-  df_models %>% filter(Model %in% c("No pooling", "Partial pooling")) %>% group_by(id,Model) %>% slice(1) %>%
  group_by(Model, condit) %>% 
  summarise(
    n= n(),
    Intercept = mean(Intercept), 
    Slope_band = mean(Slope_band),
  ) %>% mutate( Intercept_se = sd(Intercept)/sqrt(n),
    Slope_band_se = sd(Slope_band)/sqrt(n), .groups="keep") 
#head(grpAvg)


 p1=grpAvg %>% ggplot() + 
  aes(x = Model, y = Slope_band, fill = condit) +
  geom_col(position = "dodge") + 
  geom_errorbar(aes(ymin = Slope_band - Slope_band_se, ymax = Slope_band + Slope_band_se), width = 0.2, position = position_dodge(0.9)) +
  labs(x = "Model", y = "Slope (band)") + 
  theme(legend.position = "top", legend.justification = "left")+ggtitle("Comparing Slopes between Conditions - Both pooling models")
 
 
 
 p2=grpAvg %>% ggplot() + 
  aes(x = Model, y = Intercept, fill = condit) +
  geom_col(position = "dodge") + 
  geom_errorbar(aes(ymin = Intercept - Intercept_se, ymax = Intercept + Intercept_se), width = 0.2, position = position_dodge(0.9)) +
  labs(x = "Model", y = "Intercept") + 
  theme(legend.position = "top", legend.justification = "left")+ggtitle("Comparing Intercepts between Conditions - Both pooling models")
  


# For the partial pooling model, visualize the correlation between the intercept and slope for each subject.
# Use geom_smooth to fit a linear model to the data, and plot the line of best fit.

p3=df_models %>% filter(Model == "Partial pooling") %>% 
  ggplot() + 
  aes(x = Intercept, y = Slope_band) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = "Intercept", y = "Slope (band)") + 
  theme(legend.position = "top", legend.justification = "left")+ggtitle("Correlation between Fit Slope and Intercept (Mixed Effects model)")


  # For the partial pooling model, visualize the correlation between slope and devMean for each subject.
# Use geom_smooth to fit a linear model to the data, and plot the line of best fit.

p4=df_models %>% filter(Model == "Partial pooling") %>% 
  ggplot() + 
  aes(x = Slope_band, y = devMean) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = "Slope (band)", y = "devMean") + 
  theme(legend.position = "top", legend.justification = "left")+ggtitle("Correlation between Fit Slope and testing performance (Mixed Effects model)")


# For the partial pooling model, visualize the correlation between Intercept and devMean for each subject.

p5=df_models %>% filter(Model == "Partial pooling") %>% 
  ggplot() + 
  aes(x = Intercept, y = devMean) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = "Intercept", y = "devMean") + 
  theme(legend.position = "top", legend.justification = "left")+ggtitle("Correlation between Fit Intercept and testing performance (Mixed Effects model)")



title = ggdraw()+draw_label("Examining the Fit Slopes and Intercepts",fontface = 'bold',x=0,hjust=0)+theme(plot.margin = margin(0, 0, 0, .5))
plot_grid(title,NULL,p1,p2,p3,p4,p5,ncol=2,rel_heights = c(.15,1,1,1))

```





## Correlation between fit parameters (Slope and Intercept) and testing Vx
* Noteworthy that The correlation between slope and Vx is strongest for the slowest bands (100-300 and 350-550), for both original and reverse ordered groups. 
The slow positions are extrapolation for the Original ordered group, and trained by the reverse ordered group. 
* Fairly similar patterns for Slope and Intercept


```{r fig.width=12,fig.height=10}
#| column: page-inset-right

# For the partial pooling model, visualize the correlation between slope and devMean for each subject. Facet by vb~catOrder. Group and color by condit. 
# Use geom_smooth to fit a linear model to the data, and plot the line of best fit.

df_models %>% filter(Model == "Partial pooling") %>% 
  ggplot() + 
  aes(x = Slope_band, y = vxMed, color = condit) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = "Slope (band)", y = "Median Vx") + 
  theme(legend.position = "top", legend.justification = "left") +
  facet_grid(catOrder~vb)+ ggtitle("Correlation between Slope and Median VX")



df_models %>% filter(Model == "Partial pooling") %>% 
  ggplot() + 
  aes(x = Intercept, y = vxMed, color = condit) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = "Intercept", y = "Median Vx") + 
  theme(legend.position = "top", legend.justification = "left") +
  facet_grid(catOrder~vb)+ ggtitle("Correlation between Intercept and Median Vx")

```



### Correlation between parameters and Mean Deviation. 
- Here  we see a powerful effect of slope for the slow bands (larger slopes tend to have smaller deviation)

```{r fig.width=12,fig.height=10}
#| column: page-inset-right

# For the partial pooling model, visualize the correlation between slope and devMean for each subject. Facet by vb~catOrder. Group and color by condit. 
# Use geom_smooth to fit a linear model to the data, and plot the line of best fit.

df_models %>% filter(Model == "Partial pooling") %>% 
  ggplot() + 
  aes(x = Slope_band, y = devMean, color = condit) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = "Slope (band)", y = "devMean") + 
  theme(legend.position = "top", legend.justification = "left") +
  facet_grid(catOrder~vb)+ ggtitle("Correlation between Slope and Mean Deviation")



df_models %>% filter(Model == "Partial pooling") %>% 
  ggplot() + 
  aes(x = Intercept, y = devMean, color = condit) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = "Intercept", y = "devMean") + 
  theme(legend.position = "top", legend.justification = "left") +
  facet_grid(catOrder~vb)+ ggtitle("Correlation between Intercept and Mean Deviation")

```





