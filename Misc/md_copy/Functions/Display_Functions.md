library(tidyverse)



learn_curve_plot <- function(df, x_var, y_var, color_var, facet_var = NULL, groupVec, nbins, labels = FALSE) {
  df |> 
    group_by(pick({{ groupVec }})) |> 
    mutate(Trial_Bin = cut( {{x_var}} , breaks = seq(1,nbins+1), include.lowest = TRUE, labels = labels)) |>
    ggplot(aes(x = Trial_Bin, y = {{ y_var }}, col = {{ color_var }})) +
    stat_summary(aes(color = {{ color_var }}), geom = "line", fun = mean) +
    stat_summary(geom = "errorbar", fun.data = mean_se, width = .4, alpha = .7) +
    facet_wrap(vars({{facet_var}})) + 
    scale_x_continuous(breaks=seq(1,nbins+1))
}

learn_curve_table <- function(df, x_var, y_var,gw,groupVec, nbins, labels = FALSE,prefix="") {
  df |> 
    group_by(pick({{ groupVec }})) |> 
    mutate(Trial_Bin = cut( {{x_var}} , breaks = nbins, labels = labels)) |> 
    group_by(Trial_Bin,pick({{ groupVec }})) |>
    summarize(
      mean_y = round(mean({{ y_var }}, na.rm = TRUE),0),
      se_y = round(sd({{ y_var }}, na.rm = TRUE) / sqrt(n()),0),
      .groups = "drop"
    ) |>
    mutate(combined = paste0(mean_y, " (", se_y, ")")) |>
    select(-mean_y, -se_y) |>
    pivot_wider(names_from = {{ gw }}, values_from = combined,names_prefix=prefix) |>
    mutate_all(~as.character(.)) %>% replace(is.na(.), "..")
    # mutate(across(starts_with(c("combined")), as.character),
    #        across(starts_with(c("combined")), replace_na,"..")) 
    
}


plotWithTable <- function(p,t,arrange="H")
{
  require(patchwork); require(ggpmisc)
  ggp_table <- ggplot() +                            
    theme_void() +
    annotate(geom = "table",
             x = 1,
             y = 1,
             label = list(t))
  if (arrange=="H"){vp1+ggp_table}
  else {vp1/ggp_table}
   
  
}






pf <- function(df,groupVars=c("id","input","c","lr","inNodes","outNodes","noise_sd")){
  groupVars2=groupVars[!groupVars %in% c("id","input")]
  labVec = paste0(groupVars2, "=", "{", groupVars2, "}", collapse = ";")
  lc = df %>% group_by(across(all_of(groupVars)),trial) %>%
    summarize(almTrainDat = mean(almTrainDat), .groups = "drop") %>%
    mutate(simLab = as.factor(glue(labVec))) %>%
    ggplot(aes(x = trial, y = almTrainDat, color = input)) +
    geom_line() + ylim(c(0,1600))+facet_wrap(~simLab)
  wmp= df %>% filter(trial==1) %>% 
    mutate(k=map(weights,~ as.data.frame(matrix(unlist(.),nrow=7)) %>% 
                   rownames_to_column("Input") %>%
                   pivot_longer(-c(Input), names_to = "Output",names_prefix = "V", values_to = "weight") %>%
                   mutate(Output= fct_relevel(Output,unique(.$Output))))) %>%
    unnest(k) %>% ggplot(.,aes(x=Input, y=Output, fill=weight)) + 
    geom_raster() + 
    scale_fill_viridis_c()+facet_wrap(~id+c,scales="free")
  print(lc)
  
}


pfLc <- function(df,groupVars=c("id","input","c","lr","inNodes","outNodes","noise_sd")){
  groupVars2=groupVars[!groupVars %in% c("id","input")]
  labVec = paste0(groupVars2, "=", "{", groupVars2, "}", collapse = ";")
  lc = df %>% group_by(across(all_of(groupVars)),trial) %>%
    summarize(almTrainDat = mean(almTrainDat), .groups = "drop") %>%
    mutate(simLab = as.factor(glue(labVec))) %>%
    ggplot(aes(x = trial, y = almTrainDat, color = input)) +
    geom_line() + ylim(c(0,1600))+facet_wrap(~simLab)
  print(lc)
  
}


# function that takes training dataframe, bins trials into nBlocks blocks, and prints a table with a row for each model, and a column for each block
trainTab <- function(df,nBlocks=6,groupVars=c("id","c","lr","inNodes","outNodes","noise_sd")){
  #groupVars=c("id","c","lr","inNodes","outNodes","noise_sd")
  groupVars2=groupVars[!groupVars %in% "id"]
  labVec = paste0(groupVars2, "=", "{", groupVars2, "}", collapse = ";")
  vxTab=df %>%
    group_by(across(all_of(groupVars)),block = cut(trial, nBlocks, labels = FALSE)) %>%
    group_by(block, .add = TRUE) %>%
    summarize(almTrainDat = mean(almTrainDat, na.rm = TRUE)) %>%
    mutate(simLab = glue(labVec)) %>% 
    pivot_wider(names_from = block, values_from = almTrainDat)
  
  distTab = df %>% 
    group_by(across(all_of(groupVars))) %>%
    mutate(block = cut(trial, nBlocks, labels = FALSE),
           dist=vx-almTrainDat,
           absDist=abs(vx-almTrainDat)) %>%
    group_by(block, .add = TRUE) %>%
    summarize(meanAbsDist=round(mean(absDist,na.rm=TRUE),0)) %>%
    mutate(simParameters = glue(labVec)) %>% ungroup() %>% select(simParameters,block,meanAbsDist) %>% 
    pivot_wider(names_from = block, values_from = meanAbsDist)
  
  # colnames(distTab)[-1] <- paste0("Block ", colnames(distTab)[-1])
  disTab=distTab %>% # round all to 1 decimal place
    mutate(across(where(is.numeric), \(x) round(x, 0))) %>%
    kable() %>%
    add_header_above(c(" ", "Blocks" = nBlocks))  %>%
    column_spec(1, width = "5em") %>%  # adjust the width of the SimParameters column
    kable_styling(fixed_thead = TRUE)  
  
  print(distTab)
  return(distTab)
}
