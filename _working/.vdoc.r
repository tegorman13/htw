#
#
#
#
#
#
#
#
#
#
#
#
#
# load and view data
pacman::p_load(tidyverse,patchwork,here, pander, latex2exp)
purrr::walk(here::here(c("Functions/Display_Functions.R", "Functions/alm_core.R","Functions/misc_model_funs.R")),source)
select <- dplyr::select; mutate <- dplyr::mutate 

ds <- readRDS(here::here("data/e1_md_11-06-23.rds"))
dsAvg <- ds |> group_by(condit,expMode2,tr, x) |> 
  summarise(y=mean(y),.groups="keep") 

vAvg <- dsAvg |> filter(condit=="Varied")
cAvg <- dsAvg |> filter(condit=="Constant")

i1 <- ds |> filter(id=="3")

input.layer <- c(100,350,600,800,1000,1200)
output.layer <- c(100,350,600,800,1000,1200)

# dsAvg |> filter(condit=="Constant",expMode2=="Test") %>% nrow(.)

# number of distinct bands in the Test mode

dsAvg |> filter(condit=="Varied",expMode2=="Test") %>% nrow(.)
dsAvg |> filter(condit=="Constant",expMode2=="Test") %>% nrow(.)
dsAvg |> filter(condit=="Varied",expMode2=="Test") %>% pull(x) |> n_distinct()
dsAvg |> filter(condit=="Constant",expMode2=="Test") %>% pull(x) |> n_distinct()



purrr::walk(c("con_group_exam_fits", "var_group_exam_fits", "hybrid_group_exam_fits"), 
            ~ list2env(readRDS(here::here(paste0("data/model_cache/", .x, ".rds"))), 
            envir = .GlobalEnv))
#
#
#
#
#
#
#
#
#
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

ggplot() + 
  geom_rect(data = odf, 
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = type), alpha = 0.2) +
  annotate("text", x = odf$x[odf$type == "Input"], y = odf$y[odf$type == "Input"], 
           label = lab, parse = TRUE, size = 3) +
  geom_text(data = odf %>% filter(type != "Input"), 
            aes(x = x, y = y, label = label), size = 3, fontface="bold") +
  geom_path(data = plot_edges, aes(x = x, y = y, group = id, alpha = weight), alpha=.2) +
  theme_void() +
  theme(legend.position = "none")
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
almParamsV <- cbind(Model="ALM Test Only",pluck(a_te_v, "Fit"), pluck(a_te_v, "test") %>% summarise(Test_RMSE=RMSE(y,pred)) ) |>
  rbind(cbind(Model="ALM Test & Train", pluck(a_tetr_v,"Fit"), pluck(a_tetr_v, "test") %>% summarise(Test_RMSE=RMSE(y,pred)))) |>
  rbind(cbind(Model="ALM Train Only", pluck(a_tr_v, "Fit"), pluck(a_tr_v, "test") %>% summarise(Test_RMSE=RMSE(y,pred)))) |>
  mutate(across(where(is.numeric), \(x) round(x, 3)))

examParamsV <- cbind(Model="EXAM Test Only",pluck(ex_te_v, "Fit"), pluck(ex_te_v, "test") %>% summarise(Test_RMSE=RMSE(y,pred))) |>
  rbind(cbind(Model="EXAM Test & Train", pluck(ex_tetr_v,"Fit"), pluck(ex_tetr_v, "test") %>% summarise(Test_RMSE=RMSE(y,pred)))) |>
  rbind(cbind(Model="EXAM Train Only", pluck(ex_tr_v, "Fit"), pluck(ex_tr_v, "test") %>% summarise(Test_RMSE=RMSE(y,pred)))) |>
  mutate(across(where(is.numeric), \(x) round(x, 3)))

hybridParamsV <-cbind(Model="Hybrid Test Only",pluck(hybrid_te_v, "Fit"), pluck(hybrid_te_v, "test") %>% summarise(Test_RMSE=RMSE(y,pred))) |>
  rbind(cbind(Model="Hybrid Test & Train", pluck(hybrid_tetr_v,"Fit"), pluck(hybrid_tetr_v, "test") %>% summarise(Test_RMSE=RMSE(y,pred)))) |>
  rbind(cbind(Model="Hybrid Train Only", pluck(hybrid_tr_v, "Fit"), pluck(hybrid_tr_v, "test") %>% summarise(Test_RMSE=RMSE(y,pred)))) |>
  mutate(across(where(is.numeric), \(x) round(x, 3)))

almParamsC <- cbind(Model="ALM Test Only",pluck(a_te_c, "Fit"), pluck(a_te_c, "test") %>% summarise(Test_RMSE=RMSE(y,pred)) ) |>
  rbind(cbind(Model="ALM Test & Train", pluck(a_tetr_c,"Fit"), pluck(a_tetr_c, "test") %>% summarise(Test_RMSE=RMSE(y,pred)))) |>
  rbind(cbind(Model="ALM Train Only", pluck(a_tr_c, "Fit"), pluck(a_tr_c, "test") %>% summarise(Test_RMSE=RMSE(y,pred)))) |>
  mutate(across(where(is.numeric), \(x) round(x, 3)))

examParamsC <- cbind(Model="EXAM Test Only",pluck(ex0_te_c, "Fit"), pluck(ex0_te_c, "test") %>% summarise(Test_RMSE=RMSE(y,pred))) |>
  rbind(cbind(Model="EXAM Test & Train", pluck(ex0_tetr_c,"Fit"), pluck(ex0_tetr_c, "test") %>% summarise(Test_RMSE=RMSE(y,pred)))) |>
  rbind(cbind(Model="EXAM Train Only", pluck(ex0_tr_c, "Fit"), pluck(ex0_tr_c, "test") %>% summarise(Test_RMSE=RMSE(y,pred)))) |>
  mutate(across(where(is.numeric), \(x) round(x, 3)))

hybridParamsC <-cbind(Model="Hybrid Test Only",pluck(hybrid_te_c, "Fit"), pluck(hybrid_te_c, "test") %>% summarise(Test_RMSE=RMSE(y,pred))) |>
  rbind(cbind(Model="Hybrid Test & Train", pluck(hybrid_tetr_c,"Fit"), pluck(hybrid_tetr_c, "test") %>% summarise(Test_RMSE=RMSE(y,pred)))) |>
  rbind(cbind(Model="Hybrid Train Only", pluck(hybrid_tr_c, "Fit"), pluck(hybrid_tr_c, "test") %>% summarise(Test_RMSE=RMSE(y,pred)))) |>
  mutate(across(where(is.numeric), \(x) round(x, 3)))


calculate_AIC <- function(rmse, n, k) {
  return(n * log((rmse^2)) + 2 * k)
}

calculate_BIC <- function(rmse, n, k) {
  return(n * log((rmse^2)) + log(n) * k)
}

n <- dsAvg |> filter(expMode2=="Test") %>% pull(x) |> n_distinct()

# ALM AIC and BIC
almParamsV$AIC <- mapply(calculate_AIC, rmse = almParamsV$Test_RMSE, MoreArgs = list(n = n, k = 2))
almParamsV$BIC <- mapply(calculate_BIC, rmse = almParamsV$Test_RMSE, MoreArgs = list(n = n, k = 2))

# EXAM AIC and BIC
examParamsV$AIC <- mapply(calculate_AIC, rmse = examParamsV$Test_RMSE, MoreArgs = list(n = n, k = 2))
examParamsV$BIC <- mapply(calculate_BIC, rmse = examParamsV$Test_RMSE, MoreArgs = list(n = n, k = 2))

# Hybrid AIC and BIC
hybridParamsV$AIC <- mapply(calculate_AIC, rmse = hybridParamsV$Test_RMSE, MoreArgs = list(n = n, k = 3))
hybridParamsV$BIC <- mapply(calculate_BIC, rmse = hybridParamsV$Test_RMSE, MoreArgs = list(n = n, k = 3))


hybridParamsV


#
#
#
#

calculate_metrics <- function(df, n, k) {
  df %>% 
    mutate(
      AIC = calculate_AIC(Test_RMSE, n, k),
      BIC = calculate_BIC(Test_RMSE, n, k)
    ) %>% 
    mutate(across(where(is.numeric), round, 3))
}

# Number of distinct testing positions
n <- dsAvg |> filter(expMode2=="Test") %>% pull(x) |> n_distinct()

# Create a list of all model data frames
model_list <- list(
  almParamsV = list(data = almParamsV, k = 2),
  examParamsV = list(data = examParamsV, k = 2),
  hybridParamsV = list(data = hybridParamsV, k = 3),
  almParamsC = list(data = almParamsC, k = 2),
  examParamsC = list(data = examParamsC, k = 2),
  hybridParamsC = list(data = hybridParamsC, k = 3)
)

# Apply the calculate_metrics function to each model data frame
model_results <- map(model_list, ~calculate_metrics(.x$data, n, .x$k))

# Extract the results back into separate data frames
list2env(model_results, envir = .GlobalEnv)

# The final tables are now in the global environment with the same names as before
# For example, you can view the updated ALM params for the Varied group with:
almParamsV
examParamsV
#
#
#
#
#
pander(almParamsV, caption="ALM"); pander(examParamsV, caption="EXAM"); pander(hybridParamsV,caption="Hybrid") 
#
#
#
#
#

calculate_metrics <- function(df, n, k) {
  df %>%
    summarise(
      Test_RMSE = RMSE(y, pred),
      AIC = calculate_AIC(Test_RMSE, n, k),
      BIC = calculate_BIC(Test_RMSE, n, k)
    ) %>%
    mutate(across(where(is.numeric), round, 3))
}

n <- dsAvg |> filter(expMode2 == "Test") %>% pull(x) |> n_distinct()

# Bind the rows for each model and condition, then calculate AIC and BIC
almParams <- bind_rows(
  cbind(Model = "ALM Test Only", calculate_metrics(pluck(a_te_v, "test"), n, 2)),
  cbind(Model = "ALM Test & Train", calculate_metrics(pluck(a_tetr_v, "test"), n, 2)),
  cbind(Model = "ALM Train Only", calculate_metrics(pluck(a_tr_v, "test"), n, 2))
)

examParams <- bind_rows(
  cbind(Model = "EXAM Test Only", calculate_metrics(pluck(ex_te_v, "test"), n, 2)),
  cbind(Model = "EXAM Test & Train", calculate_metrics(pluck(ex_tetr_v, "test"), n, 2)),
  cbind(Model = "EXAM Train Only", calculate_metrics(pluck(ex_tr_v, "test"), n, 2))
)

hybridParams <- bind_rows(
  cbind(Model = "Hybrid Test Only", calculate_metrics(pluck(hybrid_te_v, "test"), n, 3)),
  cbind(Model = "Hybrid Test & Train", calculate_metrics(pluck(hybrid_tetr_v, "test"), n, 3)),
  cbind(Model = "Hybrid Train Only", calculate_metrics(pluck(hybrid_tr_v, "test"), n, 3))
)

# Example to show the result
head(hybridParams)
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
pander(almParamsC, caption="ALM"); pander(examParamsC, caption="EXAM"); pander(hybridParamsC,caption="Hybrid")
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
####

vte <-  pluck(a_te_v, "test") |> rename(ALM=pred,Observed=y) %>% 
  cbind(.,EXAM=pluck(ex_te_v, "test") |> pull(pred)) %>%
  cbind(., Hybrid=pluck(hybrid_te_v, "test") |> pull(pred)) |>  
  pivot_longer(Observed:Hybrid, names_to="Model", values_to = "vx") |> 
  ggplot(aes(x,vx,fill=Model, group=Model)) +geom_bar(position="dodge",stat="identity") +
  scale_fill_manual(values=col_themes$wes2)+
  scale_x_continuous(breaks=sort(unique(ds$x)), labels=sort(unique(ds$x)))+ylim(0,1500) +
  theme(legend.title = element_blank(), legend.position="top") +ggtitle("Fit to Test Only")

vtetr <-  pluck(a_tetr_v, "test") |> rename(ALM=pred,Observed=y) %>% 
  cbind(.,EXAM=pluck(ex_tetr_v, "test") |> pull(pred)) %>%
  cbind(., Hybrid=pluck(hybrid_tetr_v, "test") |> pull(pred)) |>  
  pivot_longer(Observed:Hybrid, names_to="Model", values_to = "vx") |> 
  ggplot(aes(x,vx,fill=Model, group=Model)) +geom_bar(position="dodge",stat="identity") + 
  scale_fill_manual(values=col_themes$wes2)+
  scale_x_continuous(breaks=sort(unique(ds$x)), labels=sort(unique(ds$x)))+ylim(0,1500) +
  theme(legend.title = element_blank(), legend.position="top") +ggtitle("Fit to Test and Train")

vtr <-  pluck(a_tr_v, "test") |> rename(ALM=pred,Observed=y) %>% 
  cbind(.,EXAM=pluck(ex_tr_v, "test") |> pull(pred)) %>%
  cbind(., Hybrid=pluck(hybrid_tr_v, "test") |> pull(pred)) |>  
  pivot_longer(Observed:Hybrid, names_to="Model", values_to = "vx") |> 
  ggplot(aes(x,vx,fill=Model, group=Model)) +geom_bar(position="dodge",stat="identity") +
  scale_fill_manual(values=col_themes$wes2)+
  scale_x_continuous(breaks=sort(unique(ds$x)), labels=sort(unique(ds$x)))+ylim(0,1500) +
  theme(legend.title = element_blank(), legend.position="top") +ggtitle("Fit to Train Only")
  

  


vte/vtetr/vtr
#
#
#
#
#
#
#
#
#
tvte<- pluck(a_te_v, "test") |> rename(ALM=pred,Observed=y) %>% 
  cbind(.,EXAM=pluck(ex_te_v, "test") |> pull(pred)) %>%
  cbind(., Hybrid=pluck(hybrid_te_v, "test") |> pull(pred))

tvtetr<-pluck(a_tetr_v, "test") |> rename(ALM=pred,Observed=y) %>% 
  cbind(.,EXAM=pluck(ex_tetr_v, "test") |> pull(pred)) %>%
  cbind(., Hybrid=pluck(hybrid_tetr_v, "test") |> pull(pred))

tvtr<- pluck(a_tr_v, "test") |> rename(ALM=pred,Observed=y) %>% 
  cbind(.,EXAM=pluck(ex_tr_v, "test") |> pull(pred)) %>%
  cbind(., Hybrid=pluck(hybrid_tr_v, "test") |> pull(pred))

pander(tvte, caption="Varied fit to test only")
pander(tvtetr,caption="Varied fit to train and test")
pander(tvtr,caption="Varied fit to train only")
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
####

cte <-  pluck(a_te_c, "test") |> rename(ALM=pred,Observed=y) %>% 
  cbind(.,EXAM=pluck(ex0_te_c, "test") |> pull(pred)) %>%
  cbind(., Hybrid=pluck(hybrid_te_c, "test") |> pull(pred)) |>  
  pivot_longer(Observed:Hybrid, names_to="Model", values_to = "vx") |> 
  ggplot(aes(x,vx,fill=Model, group=Model)) +geom_bar(position="dodge",stat="identity") +
  scale_fill_manual(values=col_themes$wes2)+
  scale_x_continuous(breaks=sort(unique(ds$x)), labels=sort(unique(ds$x)))+ylim(0,1500) +
  theme(legend.title = element_blank(), legend.position="top") +ggtitle("Fit to Test Only")

ctetr <-  pluck(a_tetr_c, "test") |> rename(ALM=pred,Observed=y) %>% 
  cbind(.,EXAM=pluck(ex0_tetr_c, "test") |> pull(pred)) %>%
  cbind(., Hybrid=pluck(hybrid_tetr_c, "test") |> pull(pred)) |>  
  pivot_longer(Observed:Hybrid, names_to="Model", values_to = "vx") |> 
  ggplot(aes(x,vx,fill=Model, group=Model)) +geom_bar(position="dodge",stat="identity") + 
  scale_fill_manual(values=col_themes$wes2)+
  scale_x_continuous(breaks=sort(unique(ds$x)), labels=sort(unique(ds$x)))+ylim(0,1500) +
  theme(legend.title = element_blank(), legend.position="top") +ggtitle("Fit to Test and Train")

ctr <-  pluck(a_tr_c, "test") |> rename(ALM=pred,Observed=y) %>% 
  cbind(.,EXAM=pluck(ex0_tr_c, "test") |> pull(pred)) %>%
  cbind(., Hybrid=pluck(hybrid_tr_c, "test") |> pull(pred)) |>  
  pivot_longer(Observed:Hybrid, names_to="Model", values_to = "vx") |> 
  ggplot(aes(x,vx,fill=Model, group=Model)) +geom_bar(position="dodge",stat="identity") +
  scale_fill_manual(values=col_themes$wes2)+
  scale_x_continuous(breaks=sort(unique(ds$x)), labels=sort(unique(ds$x)))+ylim(0,1500) +
  theme(legend.title = element_blank(), legend.position="top") +ggtitle("Fit to Train Only")
  
cte/ctetr/ctr
#
#
#
#
#
#
tcte<- pluck(a_te_c, "test") |> rename(ALM=pred,Observed=y) %>% 
  cbind(.,EXAM=pluck(ex0_te_c, "test") |> pull(pred)) %>%
  cbind(., Hybrid=pluck(hybrid_te_c, "test") |> pull(pred))

tctetr<-pluck(a_tetr_c, "test") |> rename(ALM=pred,Observed=y) %>% 
  cbind(.,EXAM=pluck(ex0_tetr_c, "test") |> pull(pred)) %>%
  cbind(., Hybrid=pluck(hybrid_tetr_c, "test") |> pull(pred))

tctr<- pluck(a_tr_c, "test") |> rename(ALM=pred,Observed=y) %>% 
  cbind(.,EXAM=pluck(ex0_tr_c, "test") |> pull(pred)) %>%
  cbind(., Hybrid=pluck(hybrid_tr_c, "test") |> pull(pred))

pander(tcte, caption="Constant fit to test only")
pander(tctetr,caption="Constant fit to train and test")
pander(tctr,caption="Constant fit to train only")
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
