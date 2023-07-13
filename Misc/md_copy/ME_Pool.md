---
title: "ME_Pool"
date: last-modified
categories: [Analysis, R]
page-layout: full
execute: 
  warning: false
  eval: false
---


https://www.tjmahr.com/plotting-partial-pooling-in-mixed-effects-models/



```{r}
source('HTW_Prep_Paper_Data.R')

dtestAgg <- dtest %>% group_by(sbjCode,vbLabel,condit,throwCategory) %>%
  summarise(vxMean=mean(vxCapped),lowBound=first(bandInt),highBound=first(highBound),
            vbLag=first(vbLag),vbLead=first(vbLead),.groups = 'keep')



# mixed effects model, with condit and lowBound as fixed effects, and sbjCode as random effects
m <- lmer(vxMean ~ condit + lowBound + (1 + lowBound | sbjCode), dtestAgg)




```



```{r}

df_no_pooling <- lmList(vxMean ~ lowBound | sbjCode, dtestAgg) %>% 
  coef() %>% rownames_to_column("sbjCode") %>% 
  rename(Intercept = `(Intercept)`, Slope_band = lowBound) %>% 
  add_column(Model = "No pooling")

# sort the dataframe by the value of slope_band, highest to lowest
testSlopeIndv <- df_no_pooling %>% arrange(desc(Slope_band))

# Add a condit column to the dataframe, matching condition based on the value in dtestAgg for each sbjCode
testSlopeIndv <- testSlopeIndv %>% 
  left_join(dtestAgg %>% ungroup() %>% select(sbjCode, condit) %>% distinct(), by = "sbjCode") 

# Add a rank column to the dataframe, based on the value of slope_band. Smallest rank for highest value.
testSlopeIndv <- testSlopeIndv %>% group_by(condit) %>% 
  mutate(n=n(),rank = n -rank(Slope_band) +1,
         quantile = cut(rank, breaks = 4, labels = c("1st", "2nd", "3rd", "4th")),
         quintile=cut(rank,breaks=5,labels=c("1st", "2nd", "3rd", "4th","5th"))) %>% select(-n) %>%
  arrange(rank)

# Reorder the sbjCode column so that the sbjCode with the highest slope_band is first
testSlopeIndv$sbjCode <- factor(testSlopeIndv$sbjCode, levels = testSlopeIndv$sbjCode)


head(testSlopeIndv)
```


# Create a dataframe with the slope_band and intercept for each sbjCode
df_no_pooling <- dtestAgg %>% 
  ungroup() %>% 
  as.data.table() %>% 
  setkey(sbjCode) %>% 
  by(sbjCode, .(Intercept = lmList(vxMean ~ lowBound | sbjCode, .)$coef[1],
                Slope_band = lmList(vxMean ~ lowBound | sbjCode, .)$coef[2]))

# Add a condit column to the dataframe, matching condition based on the value in dtestAgg for each sbjCode
df_no_pooling <- df_no_pooling[, c("Intercept", "Slope_band", condit = dtestAgg[.N, condit]), by = sbjCode]

# Add a rank column to the dataframe, based on the value of slope_band. Smallest rank for highest value.
df_no_pooling <- df_no_pooling[, c("Intercept", "Slope_band", condit, rank = rank(Slope_band)), by = condit]

# Reorder the sbjCode column so that the sbjCode with the highest slope_band is first
df_no_pooling$sbjCode <- factor(df_no_pooling$sbjCode, levels = df_no_pooling$sbjCode)

head(df_no_pooling)






### Varied Individual Subjects
should change to just show a subset of avg, high-discim,low-discim exemplar sbjs. 
```{r, fig.width=12, fig.height=11}

# dtestAgg %>% filter(condit=="Varied") %>% 
#   ggplot(aes(x=lowBound,y=vxMean)) + 
#   stat_smooth(method="lm",se=FALSE)+
#   geom_point()+facet_wrap("sbjCode")+labs(x="Band",y="Vx Mean")

 testSlopeIndv %>% left_join(dtestAgg, by = c("sbjCode","condit")) %>% filter(condit=="Varied",quintile=="1st" | quintile=="5th")%>%
   ggplot()+aes(x = lowBound, y = vxMean) +
  # Set the color mapping in this layer so the points don't get a color
  geom_abline(
    aes(intercept = Intercept, slope = Slope_band, color = quintile),
    size = .75
  ) +
  geom_point() +
  facet_wrap("sbjCode") +
  labs(x = "band", y = "vxMean") +
  # Fix the color palette
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "top", legend.justification = "left")

```

## Constant Subjects No-Pooling
```{r, fig.width=12, fig.height=11}

 testSlopeIndv %>% left_join(dtestAgg, by = c("sbjCode","condit")) %>% filter(condit=="Constant",quintile=="1st" | quintile=="5th")%>%
   ggplot()+aes(x = lowBound, y = vxMean) +
  # Set the color mapping in this layer so the points don't get a color
  geom_abline(
    aes(intercept = Intercept, slope = Slope_band, color = quintile),
    size = .75
  ) +
  geom_point() +
  facet_wrap("sbjCode") +
  labs(x = "band", y = "vxMean") +
  # Fix the color palette
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "top", legend.justification = "left")

```










```{r}
# Fit a model on all the data pooled together
m_pooled <- lm(vxMean ~ lowBound, dtestAgg) 
# Repeat the intercept and slope terms for each participant
df_pooled <- tibble(
  Model = "Complete pooling",
  sbjCode = unique(dtestAgg$sbjCode),
  Intercept = coef(m_pooled)[1], 
  Slope_band = coef(m_pooled)[2]
)
head(df_pooled)

```







```{r fig.width=14, fig.height=13}

# Join the raw data so we can use plot the points and the lines.
df_models <- bind_rows(df_pooled, df_no_pooling) %>% 
  left_join(dtestAgg, by = "sbjCode")

p_model_comparison <- ggplot(df_models) + 
  aes(x = lowBound, y = vxMean) + 
  # Set the color mapping in this layer so the points don't get a color
  geom_abline(
    aes(intercept = Intercept, slope = Slope_band, color = Model),
    size = .75
  ) + 
  geom_point() +
  facet_wrap("sbjCode") +
  labs(x = "band", y = "vxMean") + 
  # Fix the color palette 
  scale_color_brewer(palette = "Dark2") + 
  theme(legend.position = "top", legend.justification = "left")
p_model_comparison

```



```{r}
m <- lmer(vxMean ~ 1 + lowBound + (1 + lowBound | sbjCode), dtestAgg)
arm::display(m)


```


```{r }

# Make a dataframe with the fitted effects
df_partial_pooling <- coef(m)[["sbjCode"]] %>% 
  rownames_to_column("sbjCode") %>% 
  as_tibble() %>% 
  rename(Intercept = `(Intercept)`, Slope_band = lowBound) %>% 
  add_column(Model = "Partial pooling")

head(df_partial_pooling)
```



```{r fig.width=14, fig.height=13}
df_models <- bind_rows(df_pooled, df_no_pooling, df_partial_pooling) %>% 
  left_join(dtestAgg, by = c("sbjCode"))

# Replace the data-set of the last plot
#p_model_comparison %+replace% df_models

ggplot(df_models) + 
  aes(x = lowBound, y = vxMean) + 
  # Set the color mapping in this layer so the points don't get a color
  geom_abline(
    aes(intercept = Intercept, slope = Slope_band, color = Model),
    size = .75
  ) + 
  geom_point() +
  facet_wrap("sbjCode") +
  labs(x = "band", y = "vxMean") + 
  # Fix the color palette 
  scale_color_brewer(palette = "Dark2") + 
  theme(legend.position = "top", legend.justification = "left")

```

## Gravity
```{r}

## Summarize the intercept and slope_band of df_models by sbjCode, condit and model. 

model3 <- df_models %>% 
  group_by(sbjCode, condit, Model) %>% 
  summarise(
    Intercept = mean(Intercept),
    Slope_band = mean(Slope_band)
  )

# For each sbjCode, compute the difference in Slope_band between the 'No pooling' and 'Partial pooling' df_models
model3 <- model3 %>% 
  group_by(sbjCode) %>% 
  mutate(
    diff_Slope_band = Slope_band[Model == "No pooling"] - Slope_band[Model == "Partial pooling"]
  ) 




# Also visualize the point for the fixed effects
df_fixef <- tibble(
  Model = "Partial pooling (average)",
  Intercept = fixef(m)[1],
  Slope_band = fixef(m)[2]
)

# Complete pooling / fixed effects are center of gravity in the plot
df_gravity <- df_pooled %>% 
  distinct(Model, Intercept, Slope_band) %>% 
  bind_rows(df_fixef)
df_gravity



df_pulled <- bind_rows(df_no_pooling, df_partial_pooling) %>% 
  group_by(sbjCode) %>% 
  mutate(
    diff_Slope_band = Slope_band[Model == "No pooling"] - Slope_band[Model == "Partial pooling"]
  ) 

df_pulled %>% filter(abs(diff_Slope_band)>.16) %>% ggplot() + 
  aes(x = Intercept, y = Slope_band, color = Model, shape = Model) + 
  geom_point(size = 2) + 
  geom_point(
    data = df_gravity, 
    size = 5,
    # Prevent size-5 point from showing in legend keys
    show.legend = FALSE
  ) + 
  # Draw an arrow connecting the observations between models
  geom_path(
    aes(group = sbjCode, color = NULL), 
    arrow = arrow(length = unit(.01, "npc")),
    show.legend = FALSE
  ) + 
  # Use ggrepel to jitter the labels away from the points
  ggrepel::geom_text_repel(
    aes(label = sbjCode, color = NULL), 
   # data = df_no_pooling,
    show.legend = FALSE
  )+ 
  xlab("Intercept estimate") + 
  ylab("Slope estimate") + 
  scale_shape_manual(values = c(15:18)) +
  scale_color_brewer(palette = "Dark2") 


```


```{r}

# Extract the matrix
cov_mat <- VarCorr(m)[["sbjCode"]]

# Strip off some details so that just the useful part is printed
attr(cov_mat, "stddev") <- NULL
attr(cov_mat, "correlation") <- NULL
cov_mat


library(ellipse)
#> 
#> Attaching package: 'ellipse'
#> The following object is masked from 'package:graphics':
#> 
#>     pairs

# Helper function to make a data-frame of ellipse points that 
# includes the level as a column
make_ellipse <- function(cov_mat, center, level) {
  ellipse(cov_mat, centre = center, level = level) %>%
    as.data.frame() %>%
    add_column(level = level) %>% 
    as_tibble()
}

center <- fixef(m)
levels <- c(.1, .3, .5, .7, .9)

# Create an ellipse dataframe for each of the levels defined 
# above and combine them
df_ellipse <- levels %>%
  lapply(
    function(x) make_ellipse(cov_mat, center, level = x)
  ) %>% 
  bind_rows() %>% 
  rename(Intercept = `(Intercept)`, Slope_band = lowBound)

#df_ellipse

df_pulled %>% 
  #filter(abs(diff_Slope_band)<.1) %>% 
  ggplot() + 
  aes(x = Intercept, y = Slope_band, color = Model, shape = Model) + 
  # Draw contour lines from the distribution of effects
  geom_path(
    aes(group = level, color = NULL, shape = NULL), 
    data = df_ellipse, 
    linetype = "dashed", 
    color = "grey40"
  ) + 
  geom_point(
    aes(shape = Model),
    data = df_gravity, 
    size = 5,
    show.legend = FALSE
  ) + 
  geom_point(size = 2) + 
  geom_path(
    aes(group = sbjCode, color = NULL), 
    arrow = arrow(length = unit(.02, "npc")),
    show.legend = FALSE
  ) + 
  theme(
    legend.position = "bottom", 
    legend.justification = "right"
  ) + 
  ggtitle("Topographic map of regression parameters") + 
  xlab("Intercept estimate") + 
  ylab("Slope estimate") + 
  scale_color_brewer(palette = "Dark2") +
  scale_shape_manual(values = c(15:18))


last_plot() +
  coord_cartesian(
    xlim = range(df_pulled$Intercept), 
    ylim = range(df_pulled$Slope_band),
    expand = TRUE
  ) 

```

```{r}


# Euclidean distance
contour_dist <- function(xs, ys, center_x, center_y) {
  x_diff <- (center_x - xs) ^ 2
  y_diff <- (center_y - ys) ^ 2
  sqrt(x_diff + y_diff)
}

# Find the point to label in each ellipse.
df_label_locations <- df_ellipse %>% 
  group_by(level) %>%
  filter(
    Intercept < quantile(Intercept, .25), 
    Slope_band < quantile(Slope_band, .25)
  ) %>% 
  # Compute distance from center.
  mutate(
    dist = contour_dist(Intercept, Slope_band, fixef(m)[1], fixef(m)[2])
  ) %>% 
  # Keep smallest values.
  top_n(-1, wt = dist) %>% 
  ungroup()

# Tweak the last plot one more time!
last_plot() +
  geom_text(
    aes(label = level, color = NULL, shape = NULL), 
    data = df_label_locations, 
    nudge_x = .5, 
    nudge_y = .8, 
    size = 3.5, 
    color = "grey40"
  )

```



# Effect of condition on vxMean

```{r}

mod <- lme4::glmer(vxMean ~ lowBound * condit + (lowBound | sbjCode), data = dtestAgg, family = gaussian)




```








```{r}

# library(rstanarm)
# library(rstan)
# library(brms)

# # check the distribution of vxMean to see if it's normal. Use ggplot, color by throwCategory
# ggplot(dtestAgg, aes(x = vxMean, fill = throwCategory)) + 
#   geom_density() + 
#   theme_bw() + 
#   xlab("vxMean") + 
#   ylab("Density") + 
#   ggtitle("Distribution of vxMean by lowBound")


# # fit brms model, use gaussian family
# b <- brm(
#   vxMean ~ lowBound + (lowBound | sbjCode),
#   family = gaussian(),
#   data = dtestAgg)






# # fit stan model, use gaussian family
# b <- stan_glmer(
#   vxMean ~ lowBound + (lowBound | sbjCode),
#   family = gaussian(),
#   data = dtestAgg,
#   prior = normal(0, 2, autoscale = TRUE),
#   prior_intercept = normal(0, 5, autoscale = TRUE),
#   # reproducible blogging
# )




# b <- stan_glmer(
#   vxMean ~ lowBound + (lowBound | sbjCode),
#   family = gaussian(),
#   data = dtestAgg,
#   prior = normal(0, 2, autoscale = TRUE),
#   prior_intercept = normal(0, 5, autoscale = TRUE),
#   prior_covariance = decov(regularization = 2),
#   prior_aux = cauchy(0, 1, autoscale = TRUE), 
#   # reproducible blogging
#   seed = 20211116
# )

```






```{r}

# make a histogram, use ggplot2 
#ggplot(dtestAgg, aes(x = vxMean)) + geom_histogram(binwidth = 1) + facet_wrap("condit")




```


```{r}

```

