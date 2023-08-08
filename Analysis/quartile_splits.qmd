---
title: Comparing, High/Low performers
date: last-modified
categories: [Analysis, R]
page-layout: full
# fig-width: 15
# fig-height: 8
code-fold: true
code-tools: true
toc: false
execute: 
  warning: false
  eval: true
---


```{r setup}
pacman::p_load(tidyverse,lme4,emmeans,here,knitr,kableExtra,gt,gghalves)
e1 <- readRDS(here("data/e1_08-04-23.rds"))
source(here("Functions/Display_Functions.R"))

train <- e1 |> filter(expMode %in% c("train"))
trainAvg <- train %>% group_by(id, condit, vb, bandInt,trainStage) %>%
  summarise(vx=mean(vx),dist=mean(dist))
test <- e1 |> filter(expMode %in% c("test-Nf","test-train-nf"))
testAvg <- test %>% group_by(id, condit, vb, bandInt,bandType,tOrder) %>%
  summarise(vx=mean(vx),dist=mean(dist),sdist=mean(sdist))

```



```{r}


testAvg %>% 
  group_by(condit,vb) %>%  
  reframe(enframe(quantile(dist, c(0.0,0.25, 0.5, 0.75,1)), "quantile", "dist")) |> 
  pivot_wider(names_from=quantile,values_from=dist,names_prefix="Q_") |>
  group_by(vb,condit) |>
  summarise(across(starts_with("Q"), list(mean = mean)))



testAvg %>% 
  group_by(condit,vb) %>%  
  reframe(enframe(quantile(dist, seq(0,1,1/8)), "quantile", "dist")) |> 
  pivot_wider(names_from=quantile,values_from=dist,names_prefix="Q_") |>
  group_by(vb,condit) |>
  summarise(across(starts_with("Q"),round,1)) %>% kable(format="html",escape=F) %>% kable_styling() 
```



```{r}


raw_table <- testAvg %>% 
  group_by(condit,vb) %>%  
  reframe(enframe(quantile(dist, seq(0,1,1/8)), "quantile", "dist")) |> 
  pivot_wider(names_from=quantile,values_from=dist,names_prefix="Q_") |>
  group_by(vb,condit) |>
  summarise(across(starts_with("Q"), round,1))


long_data <- raw_table %>% 
  pivot_longer(
    cols = starts_with("Q"),
    names_to = "Quartile",
    values_to = "Value"
  ) %>%
  mutate(Quart = str_remove_all(Quartile, "Q_|%"), # Remove "Q_" and "%"
         Quart = as.numeric(Quart), # Convert Quart to numeric
         Quart = factor(Quart, levels = sort(unique(Quart)))) 






bold_lower <- function(data, by_group) {
  ifelse(data < by_group, cell_spec(data, "html", bold = T), as.character(data))
}

# Separate data by condition
constant_data <- raw_table %>% filter(condit == "Constant")
varied_data <- raw_table %>% filter(condit == "Varied")

# Apply function to varied data
varied_data <- varied_data %>%
  group_by(vb) %>%
  mutate(across(starts_with("Q"), function(.x) {
    col_name <- cur_column()
    by_group <- constant_data[[col_name]][constant_data$vb == first(vb)]
    bold_lower(.x, by_group)
  }, .names = "{.col}"))

# Format the constant_data to match the varied_data
constant_data <- constant_data %>%
  group_by(vb) %>%
  mutate(across(starts_with("Q"), ~as.character(.x), .names = "{.col}"))

# Join the constant and varied data frames back together
final_table <- bind_rows(constant_data, varied_data) %>%
  select(vb, condit, ends_with("%")) %>%
  arrange(vb, condit)

# Print the table
final_table %>% kable("html", escape = F) %>% kable_styling() %>%
  pack_rows(index = table(final_table$vb))


```


```{r}

ggplot(long_data, aes(x = Quart, y = Value, color = condit, group = condit)) +
  facet_wrap(~ vb) +
  geom_line() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(title = "Comparison of Varied and Constant Conditions across Octiles", 
       x = "Octile", 
       y = "Distance from target")
```


```{r}


model <- lmer(Value ~ condit * Quart + (1|vb), data = long_data)
summary(model)
library(multcomp)
comparison <- glht(model, linfct = mcp(condit = "Tukey"))
summary(comparison)


```


```{r}




raw_table <- test %>% 
  group_by(id,condit,vb) %>%  
  reframe(enframe(quantile(dist, seq(0,1,1/6)), "quantile", "dist")) |> 
  pivot_wider(names_from=quantile,values_from=dist,names_prefix="Q_") |>
  group_by(id,vb,condit) |>
  summarise(across(starts_with("Q"), round,1))


long_data <- raw_table %>% 
  pivot_longer(
    cols = starts_with("Q"),
    names_to = "Quartile",
    values_to = "Value"
  ) %>%
  mutate(Quart = str_remove_all(Quartile, "Q_|%"), # Remove "Q_" and "%"
         Quart = as.numeric(Quart), # Convert Quart to numeric
         Quart = factor(Quart, levels = sort(unique(Quart)))) 


ggplot(long_data, aes(x = Quart, y = Value, fill = condit, color=condit,group = condit)) +
  facet_wrap(~ vb) +
  stat_summary(geom="line",fun=mean)+
  stat_summary(geom="errorbar",fun.data=mean_se)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(title = "Comparison of Varied and Constant Conditions across Octiles", 
       x = "Octile", 
       y = "Distance from target")

model <- lmer(Value ~ condit * Quart + (1|vb)+(1|id), data = long_data)
summary(model)
comparison <- glht(model, linfct = mcp(condit = "Tukey"))
summary(comparison)
```