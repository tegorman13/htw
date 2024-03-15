

#pacman::p_load(ggplot2,wesanderson,glue, ggdist,ggforce,patchwork,gghalves, ggh4x)
pacman::p_load(ggplot2,wesanderson,glue)
#pacman::p_load(ggplot2,ggpattern,viridis,scales,wesanderson,ggthemer)
options(dplyr.summarise.inform=FALSE)

select <- dplyr::select
mutate <- dplyr::mutate
filter <- dplyr::filter
map <- purrr::map

stat_bar <- list(stat_summary(fun=mean, geom="bar", position=position_dodge(), alpha=.75),
  stat_summary(fun.data=mean_se, geom="errorbar", position=position_dodge()))


stat_bar_sd <- list(stat_summary(fun=mean, geom="bar", position=position_dodge(), alpha=.75),
  stat_summary(fun.data=mean_sdl, geom="errorbar", ,fun.args = list(mult = 1), position=position_dodge(), alpha=.9))


#### Bayesian/Brms functions

rename_fm <- function(fm) {
  # Apply case_when to transform fm values
  transformed_fm <- dplyr::case_when(
    fm == "Test" ~ "Fit to Test Data",
    fm == "Train" ~ "Fit to Training Data",
    fm == "Test_Train" ~ "Fit to Test & Training Data",
    TRUE ~ "WRONG INPUT"
  )
  
  # Convert the transformed_fm into a factor with levels in the specific order
  factor(transformed_fm, levels = c("Fit to Test Data", "Fit to Test & Training Data","Fit to Training Data","WRONG INPUT"))
}


alm_plot <- function (){
  pacman::p_load(tidyverse,ggplot2,igraph,ggraph) 

inNodes <- c("exp(c * (100 - Stim)^2)", "exp(c * (350 - Stim)^2)", 
             "exp(c * (600 - Stim)^2)", "exp(c * (800 - Stim)^2)", 
             "exp(c * (1000 - Stim)^2)", "exp(c * (1200 - Stim)^2)")

outNodes <- c(100,350,600,800,1000,1200) %>% as.integer()
stim <- "Stim"
resp <- "Response"
inFlow <- tibble(expand.grid(from=stim,to=inNodes)) %>% mutate_all(as.character)
outFlow <- tibble(expand.grid(from=outNodes,to=resp)) %>% mutate_all(as.character)

gd <- tibble(expand.grid(from=inNodes,to=outNodes)) %>% mutate_all(as.character) %>%
  rbind(inFlow,.) %>% rbind(.,outFlow)

g = graph_from_data_frame(gd,directed=TRUE)
coords2=layout_as_tree(g)
colnames(coords2)=c("y","x")

odf <- as_tibble(coords2) %>% 
  mutate(label=vertex_attr(g,"name"),
         type=c("stim",rep("Input",length(inNodes)),rep("Output",length(outNodes)),"Resp"),
         x=x*-1) %>%
  mutate(y=ifelse(type=="Resp",0,y),
         # Adjust the width of input nodes
         xmin=ifelse(type=="Input", x-0.3, x-0.08),
         xmax=ifelse(type=="Input", x+0.3, x+0.08),
         ymin=y-.30, ymax=y+.30)

input_y <- odf %>% filter(type == "Input") %>% pull(y)
output_y <- odf %>% filter(type == "Output") %>% pull(y)
avg_input_y <- mean(input_y)
avg_output_y <- mean(output_y)
y_adjustment <- avg_input_y - avg_output_y
odf <- odf %>%
  mutate(y = ifelse(type == "Output", y + y_adjustment, y), 
         ymax = ifelse(type == "Output", ymax + y_adjustment, ymax),
         ymin = ifelse(type == "Output", ymin + y_adjustment, ymin))


plot_edges = gd %>% mutate(id=row_number()) %>%
  pivot_longer(cols=c("from","to"),names_to="s_e",values_to=("label")) %>%
  mutate(label=as.character(label)) %>% 
  group_by(id) %>%
  mutate(weight=sqrt(rnorm(1,mean=0,sd=10)^2)/10) %>%
  left_join(odf,by="label") %>%
  mutate(xmin=xmin+.02,xmax=xmax-.02)

lab <- inNodes
names(lab) <- inNodes
lab <- lapply(lab, function(x) paste0("bold(", x, ")")) # Bold the entire expression
lab <- sapply(lab, as.character)

g <- ggplot() + 
  geom_rect(data = odf, 
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = type), alpha = 0.2) +
  annotate("text", x = odf$x[odf$type == "Input"], y = odf$y[odf$type == "Input"], 
           label = lab, parse = TRUE, size = 3) +
  geom_text(data = odf %>% filter(type != "Input"), 
            aes(x = x, y = y, label = label), size = 3, fontface="bold") +
  geom_path(data = plot_edges, aes(x = x, y = y, group = id, alpha = weight), alpha=.2) +
  theme_void() +
  theme(legend.position = "none")
  print(g)
}




brms_posterior_checks <- function(brmModel, yvar, gvar, ndraws = 50) {
  
  dat <- brmModel$data
  
  yvar_name <- deparse(substitute(yvar))
  gvar_name <- deparse(substitute(gvar))
  
  y_data <- dat[[yvar_name]]
  g_data <- dat[[gvar_name]]
  
  msum <- summary(brmModel)
  pt <- posterior_table(brmModel)
  timeTab <- rstan::get_elapsed_time(brmModel[["fit"]]) |> 
    as_tibble(rownames = "chain") |>  mutate(total_seconds = warmup + sample) 
  
  print(msum)
  print(pt)
  print(timeTab)
  
  if (brmModel$family$family == "gaussian") {
   # print(conditional_effects(brmModel))
    print(pp_check(brmModel, type = "stat_grouped", group = "vb", ndraws = ndraws))
  #  print(bayes_R2(brmModel))
    
  }
  print(bayesplot::ppc_dens_overlay_grouped(y_data, posterior_predict(brmModel, ndraws = ndraws), g_data))
  #fixef(brmModel,summary=TRUE)
}



plot_subject_fits <- function(model, subject_code) {
  pattern <- glue("^r_id\\[{subject_code},.*\\]")
  plot <- mcmc_areas(model, prob = .5, regex_pars = c(pattern)) +
    ggtitle(glue("fit for subject #{subject_code}:"))
  return(plot)
}



posterior_table <- function(model){
  bayestestR::describe_posterior(model,
                                 verbose=FALSE,
                                 test=c("p_direction","p_significance"),
                                 centrality=c("median"),
                                 dispersion=FALSE)
  
}


indv_model_plot <- function(combined_df, indv_coefs, testAvg,slopeVar, rank_variable = "Estimate.Intercept", n_sbj = 5, order = "min") {
  slice_fn <- if (order == "min") slice_min else slice_max
  combined_df  |> 
    filter(id %in% (indv_coefs  |> 
                      slice_fn({{ rank_variable }}, n = n_sbj, by = condit) |> 
                      pull(id)))  |> 
    group_by(id, bandInt)  |> 
    sample_n(100)  |> 
    ggplot(aes(x = bandInt, y = estimate)) +
    geom_abline(aes(intercept = Intercept, slope = {{slopeVar}}), color = "grey50") +
    geom_abline(data = indv_coefs  |> 
                  slice_fn({{ rank_variable }}, n = n_sbj, by = condit),
                aes(intercept = Intercept, slope = {{slopeVar}}), color = "red") +
    stat_halfeye() +
    stat_halfeye(data = testAvg  |> 
                   filter(id %in% (indv_coefs |> slice_fn({{ rank_variable }}, n = n_sbj, by = condit) |> 
                                     pull(id))), 
                 aes(x = bandInt, y = vx), color = "blue") +
    scale_x_continuous(breaks = c(100, 350, 600, 800, 1000, 1200), 
                       labels = levels(testAvg$vb), 
                       limits = c(0, 1400)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5,size=8.5)) +
    ggh4x::facet_nested_wrap(vars(condit, id),ncol=3)
}


# bayes_R2(e1_testDistRF)
# 
# tidyMCMC(e1_testDistRF, conf.int = TRUE, conf.level = 0.95,
#          estimate.method = "median", conf.method = "HPDinterval")
# 
# (r_fit <- e1_testDistRF %>% 
#     tidy() %>% filter(effect=="fixed") |> select(-effect,-component, -group) |> 
#     mutate(term = janitor::make_clean_names(term)) |>
#     mutate(across(where(is.numeric), \(x) round(x, 0))) |> kbl() )

# brms_eq_tidy <-tidyMCMC(e1_testDistRF, conf.int = TRUE, conf.level = 0.95,
#                         estimate.method = "median", conf.method = "HPDinterval")

############



## Learning Plots
learn_curve_bins<- function(df, x_var, y_var,gw,groupVec, nbins, labels = FALSE,prefix="") {
  df |> 
    group_by(pick({{ groupVec }})) |> 
    mutate(Trial_Bin = cut( {{x_var}} , breaks = nbins, labels = labels)) |> 
    group_by(Trial_Bin,pick({{ groupVec }})) |>
    summarize(
      mean_y = round(mean({{ y_var }}, na.rm = TRUE),0),
      se_y = round(sd({{ y_var }}, na.rm = TRUE) / sqrt(n()),0),
      .groups = "drop"
     ) # |>
    # mutate(combined = paste0(mean_y, " (", se_y, ")")) |>
    # select(-mean_y, -se_y)
}


# learn_curve_plot <- function(df, x_var, y_var, color_var, facet_var = NULL, groupVec, nbins, labels = FALSE) {
#   df |> 
#     group_by(pick({{ groupVec }})) |> 
#     mutate(Trial_Bin = cut( {{x_var}} , breaks = seq(1,nbins+1), include.lowest = TRUE, labels = labels)) |>
#     ggplot(aes(x = Trial_Bin, y = {{ y_var }}, col = {{ color_var }})) +
#     stat_summary(aes(color = {{ color_var }}), geom = "line", fun = mean) +
#     stat_summary(geom = "errorbar", fun.data = mean_se, width = .4, alpha = .7) +
#     facet_wrap(vars({{facet_var}})) + 
#     scale_x_continuous(breaks=seq(1,nbins+1))
# }


learn_curve_plot <- function(df, x_var, y_var, color_var, facet_var = NULL, groupVec, nbins, labels = FALSE, y_label=NULL) {
  
  if (is.null(y_label)) {
    y_label <- deparse(substitute(y_var))
  }
  df |> 
    group_by({{ groupVec }}) |> 
    mutate(Trial_Bin = cut( {{x_var}}, breaks = seq(1, nbins + 1), include.lowest = TRUE, labels = labels)) |> 
    ggplot(aes(x = Trial_Bin, y = {{ y_var }}, color = {{ color_var }})) +
    stat_summary(aes(color = {{ color_var }}), geom = "line", fun = mean) +
    stat_summary(geom = "errorbar", fun.data = mean_se, width = .4, alpha = .7) +
    facet_wrap(vars({{facet_var}}), scales = 'free_y') + 
    labs(y = y_label) + # Set the y axis label dynamically
    scale_x_continuous(breaks = seq(1, nbins + 1)) 
}

learn_curve_plot2 <- function(df, x_var, y_var, color_var, facet_var = NULL, groupVec=NULL, labels = FALSE) {
  nbins= df |> ungroup() |> select({{x_var}}) %>% max()
  df |> 
    ggplot(aes(x = {{ x_var }}, y = {{ y_var }}, col = {{ color_var }})) +
    stat_summary(aes(color = {{ color_var }}), geom = "line", fun = mean) +
    stat_summary(geom = "errorbar", fun.data = mean_se, width = .4, alpha = .7) +
    facet_wrap(vars({{facet_var}})) + 
    scale_x_continuous(breaks=seq(1,nbins))
}


learn_curve_table <- function(df, x_var, y_var,gw,groupVec, nbins, nl=FALSE, labels = FALSE,prefix="") {
  separator = ifelse(nl, "<br>", " ")
  df |> 
    group_by(pick({{ groupVec }})) |> 
    mutate(Trial_Bin = cut( {{x_var}} , breaks = nbins, labels = labels)) |> 
    group_by(pick({{ groupVec }}),Trial_Bin) |> ungroup(id) |>
    summarize(
      mean_y = round(mean({{ y_var }}, na.rm = TRUE),0),
      se_y = round(sd({{ y_var }}, na.rm = TRUE) / sqrt(n()),0),
      .groups = "drop"
    ) |>
    mutate(combined = glue("{mean_y}", separator, "({se_y})")) |>
    select(-mean_y, -se_y) |>
    pivot_wider(names_from = {{ gw }}, values_from = combined,names_prefix=prefix) |>
    mutate_all(~as.character(.)) %>% replace(is.na(.), "..")
}


gt_temp <- function(gtable) {

if (!inherits(gtable, "gt_tbl")) {
  gtable <- gt::gt(gtable)
}

tmp <- tempfile(fileext = '.png') #generate path to temp .png file
gtsave(gtable, tmp) #save gt table as png
png::readPNG(tmp, native = TRUE) # read tmp png file
}



plotWithTable <- function(p,t,arrange="H")
{
  require(patchwork); require(ggpmisc)
  ggp_table <- ggplot() +                            
    theme_void() +
    annotate(geom = "table",
             x = 0,
             y = 0,
             label = list(t))
  if (arrange=="H"){vp1+ggp_table}
  else {vp1/ggp_table}
}




##### Distribution Plotting ########

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



plot_distByCondit <- function(df) {

  vbRect<- df |> group_by(vb) |>
    summarise(lowBound=first(bandInt),highBound=first(highBound)) |> 
    mutate(vbn=as.numeric(vb), rectWidth=.2,
           vbLag=vbn-rectWidth,vbLead=vbn+rectWidth)
  
  bandLines4 <- list(geom_segment(data=vbRect,aes(x=vbLag,xend=vbLead,y=highBound,yend=highBound),alpha=1,linetype="dashed"),
                   geom_segment(data=vbRect,aes(x=vbLag,xend=vbLead,y=lowBound,yend=lowBound),alpha=1,linetype="dashed"),
                   geom_text(data=vbRect,aes(x=vbLag-.03,y=lowBound+100,label=vb),angle=90,size=3.5,fontface="bold") )    

df %>% group_by(id,vb,condit,bandOrder) %>% 
  summarise(vxMean=mean(vx)) %>%
  ggplot(aes(x=vb,y=vxMean,fill=vb))+
  gghalves::geom_half_violin(color=NA)+ # remove border color
  gghalves::geom_half_boxplot(position=position_nudge(x=-0.05),side="r",outlier.shape = NA,center=TRUE, 
                              errorbar.draw = FALSE,width=.25)+
  gghalves::geom_half_point(transformation = position_jitter(width = 0.05, height = 0.05),size=.3,aes(color=vb))+
  facet_wrap(~condit,scale="free_x")+
  geom_rect(data=vbRect,aes(xmin=vbLag,xmax=vbLead,ymin=lowBound,ymax=highBound,fill=vb),alpha=.3,inherit.aes = FALSE)+
  bandLines4+
  #geom_text(data=sumStats,aes(x=vb,y=2100,label = groupMean),size=2, vjust = -0.5)+
  scale_y_continuous(expand=expansion(add=100),breaks=round(seq(0,2000,by=200),2))+
  theme(legend.position='none',
        #plot.title=element_text(face="bold"),
        #axis.title.x=element_text(face="bold"),
        # axis.title.y=element_text(face="bold"),
        axis.text.x = element_text(size = 8.5))+
  #ggtitle("Testing Performance (no-feedback) - X-Velocity Per Band")     
  ylab("Mean X Velocity")+xlab("Target  Band") 
  
 # geom_text(data=sumStats2,aes(y=2090,label = sumStatLab),size=1.9)

}


#####



# map(c("dist", "vx"), ~{
#   train %>%
#     learn_curve_bins(
#       gt.train,
#       !!sym(.x),
#       gw = Trial_Bin, 
#       groupVec = c(id, vb, condit),
#       nbins = 8
#     )
# }) %>%
#   walk(~ {
#     print(ggplot(data=.,aes(x=Trial_Bin, y=mean_y,fill=vb)) + 
#             ggdist::stat_halfeye(alpha=.5))
#   })
# 
# 
# titles = c("Title for dist", "Title for vx")
# map2(c("dist", "vx"), titles, ~{
#   train %>%
#     learn_curve_bins(
#       gt.train,
#       !!sym(.x),
#       gw = Trial_Bin, 
#       groupVec = c(id, vb, condit),
#       nbins = 8
#     ) %>%
#     ggplot(aes(x=Trial_Bin, y=mean_y,fill=vb)) + 
#     ggdist::stat_halfeye(alpha=.5) +
#     ggtitle(.y) 
# }) %>%
#   walk(~ print(.))



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





#################
#Themes
##################




# Custom theme for data visualizations
plot_theme <- function(title_size = NULL, 
                       xaxis_size = NULL, 
                       yaxis_size = NULL, 
                       strip_size = NULL, 
                       strip_face = NULL, 
                       caption.hjust = 1, 
                       caption.vjust = 0, 
                       x_axis_face = NULL, 
                       y_axis_face = NULL, 
                       transparent = FALSE, 
                       axis_text_size = NULL, 
                       legend_text_size = NULL,
                       subtitle_size = NULL,
                       caption_size = NULL,
                       base_size = 12,
                       ...) {
  .theme <- theme_minimal(base_size = base_size) + theme(
    # Specify the default settings for the plot title
    plot.title = element_text(
      size = title_size,
      #face = "bold",
      family = "serif"
    ),
     plot.tag = element_text(
      size = title_size,
      hjust=1,vjust=1,
      #face = "bold",
      family = "serif"
    ),
    # Specify the default settings for caption text
    plot.caption = element_text(
      size = caption_size,
      family = "serif",
      hjust = caption.hjust,
      vjust = caption.vjust
    ),
    # Specify the default settings for subtitle text
    plot.subtitle = element_text(
      size = subtitle_size,
      family = "serif"
    ),
    # Specify the default settings specific to the x axis title
    axis.title.y = element_text(
      size = yaxis_size, 
      face = y_axis_face, 
      family = "serif",
      margin = margin(r = 10, l = -10)
    ),
    # Specify the default settings specific to the y axis title
    axis.title.x = element_text(
      size = xaxis_size, 
      face = x_axis_face, 
      family = "serif",
      margin = margin(t = 10, b = -10)
    ),
    # Specify the default settings for x axis text
    axis.text.x = element_text(
      size = axis_text_size,
      family = "serif",
      face = x_axis_face
    ),
    # Specify the default settings for y axis text
    axis.text.y = element_text(
      size = axis_text_size,
      family = "serif",
      face = y_axis_face
    ),
    # Specify the default settings for legend titles
    legend.title = element_text(
      size = legend_text_size,
      face = "plain",
      family = "serif"
    ),
    # Specify the default settings for legend text
    legend.text = element_text(
      size = legend_text_size,
      family = "serif"
    ),
    # Set the strip background fill to blank
    strip.background = element_blank(),
    # Adjust the strip text size settings
    strip.text = element_text(
      family = "serif", 
      size = strip_size,
      face = strip_face
    ),
    # Additional Settings Passed to theme()
    ...
  )
  # Plot Transparency
  if (transparent == TRUE) {
    .theme <- .theme + theme(
      panel.background = element_rect(fill = "transparent", colour = NA),
      plot.background = element_rect(fill = "transparent", colour = NA),
      legend.background = element_rect(fill = "transparent", colour = NA),
      legend.key = element_rect(fill = "transparent", colour = NA)
    )
  }
  return(.theme)
}


#ggthemr::ggthemr('fresh')


## Set the base theme settings for ggplot2
theme_set(plot_theme(
  xaxis_size = 14,
  yaxis_size = 14,
  title_size = 16,
  caption_size = 12,
  axis_text_size = 10,
  strip_face = "plain",
  y_axis_face = "plain",
  x_axis_face = "plain",
  plot.margin = margin(2, 2, 4, 5, "mm"),
  plot.caption.position = "plot",
  plot.title.position = "plot",
  strip_size = 12,
  legend_text_size = 12,
  legend.position = "top",
  caption.hjust = 0, 
  caption.vjust = -1,
  transparent = TRUE
))

big_text <- function() 
{theme(
  xaxis_size = 18,
  yaxis_size = 18,
  title_size = 19,
  caption_size = 16,
  axis_text_size = 12,
  strip_face = "plain",
  y_axis_face = "plain",
  x_axis_face = "plain",
  strip_size = 12,
  legend_text_size = 16,
)
}

big_text <- function() {
  theme(panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "white", color = NA),
        plot.title = element_text(face = "bold",size=16),
        axis.title = element_text(face = "bold"),
        axis.title.x=element_text(face="bold",size=14),
        axis.title.y=element_text(face="bold",size=24),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        strip.text = element_text(face = "bold", size = rel(0.8), hjust = 0),
        strip.background = element_rect(fill = "grey80", color = NA),
        legend.title = element_text(face = "bold"))
}

big_text2 <- function() {
  theme(panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "white", color = NA),
        plot.title = element_text(face = "bold",size=16),
        axis.title = element_text(face = "bold"),
        axis.title.x=element_text(face="bold",size=14),
        axis.title.y=element_text(face="bold",size=24),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        strip.text = element_text(face = "bold", size = rel(0.8), hjust = 0),
        strip.background = element_rect(fill = "grey80", color = NA))
        #legend.title = element_text(face = "bold"))
}


# options(ggplot2.continuous.colour="viridis")
# options(ggplot2.continuous.fill = "viridis")

#ggokabeito::palette_okabe_ito()


# light_grey <- "#D3D3D3"  # This is a light grey color
# dark_grey <- "#000000"   # This is a dark grey color
# grey_palette <- colorRampPalette(colors = c("#D3D3D3", "#000000"))
# custom_grey_colors = grey_palette(6)



col_themes <- tibble::lst(darjeeling = c(wes_palette("Darjeeling1"),wes_palette("Darjeeling2"), ggokabeito::palette_okabe_ito(),
  wes_palette("AsteroidCity1"), wes_palette("AsteroidCity2")), 
                          wes2 = wes_palette("AsteroidCity1"), 
                          okabeito = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999", "#000000"))

# okabeito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999", "#000000")
# wes2 <- wes_palette("AsteroidCity1")


#### COLOR THEME ##### 
# darjeeling <- c(wes_palette("Darjeeling1"),wes_palette("Darjeeling2"))
#cat(darjeeling)


scale_colour_discrete <- function(...) {
  scale_colour_manual(..., values = col_themes$darjeeling)
}

scale_fill_discrete <- function(...) {
  scale_fill_manual(..., values = col_themes$darjeeling)
}

options(ggplot2.continuous.colour=col_themes$darjeeling)
options(ggplot2.continuous.fill = col_themes$darjeeling)

theme_clean <- function() {
  theme_minimal() +
    theme(panel.grid.minor = element_blank(),
          plot.background = element_rect(fill = "white", color = NA),
          plot.title = element_text(face = "bold"),
          axis.title = element_text(face = "bold"),
          strip.text = element_text(face = "bold", size = rel(0.8), hjust = 0),
          strip.background = element_rect(fill = "grey80", color = NA),
          legend.title = element_text(face = "bold"))
}

nested_settings <- ggh4x::strip_nested(
  background_x = list(element_rect(fill = "grey92"), NULL),
  by_layer_x = TRUE)



"%||%" <- function(a, b) {
  if (!is.null(a)) a else b
}

geom_flat_violin <- function(mapping = NULL, data = NULL, stat = "ydensity",
                             position = "dodge", trim = TRUE, scale = "area",
                             show.legend = NA, inherit.aes = TRUE, ...) {
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomFlatViolin,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      trim = trim,
      scale = scale,
      ...
    )
  )
}

#' @rdname ggplot2-ggproto
#' @format NULL
#' @usage NULL
#' @export
GeomFlatViolin <-
  ggproto("GeomFlatViolin", Geom,
          setup_data = function(data, params) {
            data$width <- data$width %||%
              params$width %||% (resolution(data$x, FALSE) * 0.9)
            
            # ymin, ymax, xmin, and xmax define the bounding rectangle for each group
            data %>%
              group_by(group) %>%
              mutate(
                ymin = min(y),
                ymax = max(y),
                xmin = x,
                xmax = x + width / 2
              )
          },
          
          draw_group = function(data, panel_scales, coord) {
            # Find the points for the line to go all the way around
            data <- transform(data,
                              xminv = x,
                              xmaxv = x + violinwidth * (xmax - x)
            )
            
            # Make sure it's sorted properly to draw the outline
            newdata <- rbind(
              plyr::arrange(transform(data, x = xminv), y),
              plyr::arrange(transform(data, x = xmaxv), -y)
            )
            
            # Close the polygon: set first and last point the same
            # Needed for coord_polar and such
            newdata <- rbind(newdata, newdata[1, ])
            
            ggplot2:::ggname("geom_flat_violin",
                             GeomPolygon$draw_panel(newdata, panel_scales, coord))
          },
          
          draw_key = draw_key_polygon,
          
          default_aes = aes(
            weight = 1, colour = "grey20", fill = "white", size = 0.5,
            alpha = NA, linetype = "solid"
          ),
          
          required_aes = c("x", "y")
  )
