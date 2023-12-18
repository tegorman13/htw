---
title: "Combined Analysis of All 3 Experiments"
categories: [Analysis, R, Bayesian]
execute: 
  warning: false
  eval: false 
---

Questions
    - difficulty difference between bands
    - effect of ordinal feedback. 
    - effect of training with easier bands
    - do varied subjects who discriminate by end of training do better on test?



```{r}
source(here::here("Functions", "packages.R"))
dAll <- readRDS(here("data/dAll_08-21-23.rds"))


```

```{r}
dAll <- readRDS(here("data/dAll_08-21-23.rds"))

f <- (`Band` = vb)* expMode2 ~ dist * Exp* condit * ((`Avg.` = Mean)*Arguments(fmt='%.0f') )

datasummary(f,
            data = dAll,
            output = 'gt',
            sparse_header = TRUE) 



```


### Are there differences in difficulty between bands?

```{r}
mTv <- lmer(dist ~ 0+vb + Exp + (1|id), data = dAll |> filter(expMode2=="Train",condit=="Varied"))
summary(mTv)

mTvb <- lmer(dist ~ 0+bandInt + (1|id), data = dAll |> filter(expMode2=="Train",condit=="Varied"))
summary(mTvb)

summary(lmer(dist ~ 0+bandInt + (1|id), data = dAll |> filter(expMode2=="Train")))
summary(lmer(dist ~ 0+bandInt + (1|id), data = dAll |> filter(expMode=="test-feedback")))

summary(lmer(dist ~ 0+bandInt + (1|id), data = dAll |> filter(expMode=="test-train-nf")))

```



### Is there an effect of ordinal feedback?

```{r}
mfb <- lmer(dist ~ fb*vb + (1|id), data = dAll |> filter(expMode2=="Test"))
summary(mfb)

mfb2 <- lmer(dist ~ fb*vb*condit+ (1|id), data = dAll |> filter(expMode2=="Test"))
summary(mfb2)


mfb3 <- lmer(dist ~ fb*bandInt*condit+ (1|id), data = dAll |> filter(expMode2=="Test"))
summary(mfb3)

```




### Is there an effect of training with easier bands?

```{r}

mr <- lmer(dist ~ bandOrder*vb + (1|id), data = dAll |> filter(expMode2=="Test"))
summary(mr)

mr2 <- lmer(dist ~ bandOrder*vb*condit + (1|id), data = dAll |> filter(expMode2=="Test"))
summary(mr2)

mr3 <- lmer(dist ~ bandOrder*bandInt*condit + (1|id), data = dAll |> filter(expMode2=="Test"))
summary(mr3)

```




```{r}

dAll |> filter(condit=="Varied",expMode2=="Train") |> group_by(id,vb,bandInt,Exp) |> summarize(m=mean(dist)) |> ggplot(aes(x=vb,y=m,fill=vb)) + 
     stat_summary(geom = "bar", position=position_dodge(), fun = mean) +
    stat_summary(geom = "errorbar", position=position_dodge(.9), fun.data = mean_se, width = .4) +
  facet_wrap(~Exp) 


dAll |> filter(condit=="Varied",expMode2=="Test") |> group_by(id,vb,bandInt,Exp) |> summarize(m=mean(dist)) |> ggplot(aes(x=vb,y=m,fill=vb)) + 
     stat_summary(geom = "bar", position=position_dodge(), fun = mean) +
    stat_summary(geom = "errorbar", position=position_dodge(.9), fun.data = mean_se, width = .4) +
  facet_wrap(~Exp) 

dAll |> filter(condit=="Varied",expMode2=="Test") |> group_by(id,vb,bandInt,bandOrder) |> summarize(m=mean(dist)) |> ggplot(aes(x=vb,y=m,fill=vb)) + 
     stat_summary(geom = "bar", position=position_dodge(), fun = mean) +
    stat_summary(geom = "errorbar", position=position_dodge(.9), fun.data = mean_se, width = .4) +
  facet_wrap(~bandOrder) 

dAll |> filter(condit=="Varied") |> group_by(id,vb,bandInt,bandOrder) |> summarize(m=mean(dist)) |> ggplot(aes(x=vb,y=m,fill=vb)) + 
     stat_summary(geom = "bar", position=position_dodge(), fun = mean) +
    stat_summary(geom = "errorbar", position=position_dodge(.9), fun.data = mean_se, width = .4) +
  facet_wrap(~bandOrder) 



dAll  |> group_by(id,vb,bandInt,bandOrder,condit) |> summarize(m=mean(dist)) |> ggplot(aes(x=vb,y=m,fill=vb)) + 
     stat_summary(geom = "bar", position=position_dodge(), fun = mean) +
    stat_summary(geom = "errorbar", position=position_dodge(.9), fun.data = mean_se, width = .4) +
  facet_wrap(condit~bandOrder) 
```



```{r}


dAll  |> filter(expMode2=="Test") |> group_by(id,vb,bandInt,fb,condit) |> summarize(m=mean(dist)) |> ggplot(aes(x=vb,y=m,fill=vb)) + 
     stat_summary(geom = "bar", position=position_dodge(), fun = mean) +
    stat_summary(geom = "errorbar", position=position_dodge(.9), fun.data = mean_se, width = .4) +
  facet_wrap(condit~fb) 


dAll  |> filter(expMode2=="Test") |> group_by(id,vb,bandInt,fb,condit,Exp) |> summarize(m=mean(dist)) |> ggplot(aes(x=vb,y=m,fill=vb)) + 
     stat_summary(geom = "bar", position=position_dodge(), fun = mean) +
    stat_summary(geom = "errorbar", position=position_dodge(.9), fun.data = mean_se, width = .4) +
  facet_wrap(condit~Exp) 


dAll  |> filter(expMode2=="Test") |> group_by(id,vb,bandInt,fb,condit,Exp) |> summarize(m=mean(dist)) |> ggplot(aes(x=Exp,y=m,fill=Exp)) + 
     stat_summary(geom = "bar", position=position_dodge(), fun = mean) +
    stat_summary(geom = "errorbar", position=position_dodge(.9), fun.data = mean_se, width = .4) +
  facet_wrap(condit~vb) 

dAll  |> filter(expMode2=="Test") |> group_by(id,vb,bandInt,fb,condit,Exp) |> summarize(m=mean(dist)) |> ggplot(aes(x=Exp,y=m,fill=condit)) + 
     stat_summary(geom = "bar", position=position_dodge(), fun = mean) +
    stat_summary(geom = "errorbar", position=position_dodge(.9), fun.data = mean_se, width = .4) +
  facet_wrap(~vb)


dAll  |> filter(expMode2=="Test") |> group_by(id,vb,bandType,bandInt,fb,condit,Exp) |> summarize(m=mean(dist)) |> ggplot(aes(x=Exp,y=m,fill=condit)) + 
     stat_summary(geom = "bar", position=position_dodge(), fun = mean) +
    stat_summary(geom = "errorbar", position=position_dodge(.9), fun.data = mean_se, width = .4) +
  facet_wrap(~bandType)
```



```{r}

dAll|> group_by(id,vb,bandInt,fb,condit,expMode2) |> summarize(m=mean(dist)) |> filter(expMode2=="Test" | expMode2=="Train") |> group_by(id,fb,condit,expMode2) |> summarise(m=mean(m)) |> # pivot expMode into separate columns
pivot_wider( names_from = expMode2, values_from = m) |> ggplot(aes(x=Train,y=Test,fill=condit)) + geom_point() + geom_abline(intercept = 0, slope = 1) + facet_wrap(~fb)
  
# quantify correlation between train and test for varied and fixed
dAll|> group_by(id,vb,bandInt,fb,condit,expMode2) |> summarize(m=mean(dist)) |> filter(expMode2=="Test" | expMode2=="Train") |> group_by(id,fb,condit,expMode2) |> summarise(m=mean(m)) |> # pivot expMode into separate columns
pivot_wider( names_from = expMode2, values_from = m) |> filter(condit=="Varied") %>%cor.test(~Train+Test)

result <- dAll %>% 
  # Group and calculate mean
  group_by(id, vb, bandInt, fb, condit, expMode2) %>% 
  summarize(m = mean(dist), .groups = "drop") %>% 
  # Filter for Test and Train
  filter(expMode2 == "Test" | expMode2 == "Train") %>% 
  group_by(id, fb, condit, expMode2) %>% 
  summarise(m = mean(m), .groups = "drop") %>% 
  # Pivot expMode into separate columns
  pivot_wider(names_from = expMode2, values_from = m) %>% 
  # Filter for condition 'Varied'
  filter(condit == "Varied")

# Check if the resulting dataframe has 'Train' and 'Test' columns and they are numeric
if("Train" %in% colnames(result) & "Test" %in% colnames(result)){
  if(is.numeric(result$Train) & is.numeric(result$Test)){
    cor_result <- cor.test(result$Train, result$Test)
  } else {
    stop("Either 'Train' or 'Test' column is not numeric")
  }
} else {
  stop("The data frame does not have the expected 'Train' and 'Test' columns")
}

print(cor_result)


result <- dAll %>% 
  group_by(id, vb, bandInt, fb, condit, expMode2) %>% 
  summarize(m = mean(dist), .groups = "drop") %>% 
  filter(expMode2 == "Test" | expMode2 == "Train") %>% 
  group_by(id, fb, condit, expMode2) %>% 
  summarise(m = mean(m), .groups = "drop") %>% 
  pivot_wider(names_from = expMode2, values_from = m)

# Function to compute correlation for each condition
compute_correlation <- function(data) {
  if("Train" %in% colnames(data) & "Test" %in% colnames(data)){
    if(is.numeric(data$Train) & is.numeric(data$Test)){
      return(cor.test(data$Train, data$Test))
    } else {
      stop("Either 'Train' or 'Test' column is not numeric")
    }
  } else {
    stop("The data frame does not have the expected 'Train' and 'Test' columns")
  }
}



correlations <- result %>%
  group_by(condit) %>%
  group_map(~ compute_correlation(.x))
names(correlations) <- levels(result$condit)
print(correlations)


ggplot(result, aes(x=Train, y=Test)) + 
  geom_point(aes(color=condit), alpha=0.5) +  # Plot points with different colors based on condition
  geom_smooth(aes(color=condit), method='lm', se=FALSE) +  # Overlay a regression line
  facet_wrap(~condit) +  # Separate plot for each condition
  theme_minimal() + 
  labs(title="Correlation between Train and Test",
       subtitle="Faceted by Condition (Varied vs. Constant)",
       x="Train", y="Test", color="Condition") + 
  theme(legend.position="bottom")




correlations <- result %>%
  group_by(condit,fb) %>%
  group_map(~ compute_correlation(.x))
names(correlations) <- levels(result$condit)
print(correlations)


ggplot(result, aes(x=Train, y=Test)) + 
  geom_point(aes(color=condit), alpha=0.5) +  # Plot points with different colors based on condition
  geom_smooth(aes(color=condit), method='lm', se=FALSE) +  # Overlay a regression line
  facet_wrap(fb~condit) +  # Separate plot for each condition
  theme_minimal() + 
  labs(title="Correlation between Train and Test",
       subtitle="Faceted by Condition (Varied vs. Constant)",
       x="Train", y="Test", color="Condition") + 
  theme(legend.position="bottom")




```

