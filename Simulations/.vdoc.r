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
# load and view data
pacman::p_load(tidyverse,patchwork,here, pander, latex2exp, flextable)
purrr::walk(here::here(c("Functions/Display_Functions.R", "Functions/alm_core.R","Functions/misc_model_funs.R")),source)
select <- dplyr::select; mutate <- dplyr::mutate 

ds <- readRDS(here::here("data/e1_md_11-06-23.rds"))
dsAvg <- ds |> group_by(condit,expMode2,tr, x) |> 
  summarise(y=mean(y),.groups="keep") 

vAvg <- dsAvg |> filter(condit=="Varied")
cAvg <- dsAvg |> filter(condit=="Constant")

#i1 <- ds |> filter(id=="3")

input.layer <- c(100,350,600,800,1000,1200)
output.layer <- c(100,350,600,800,1000,1200)


purrr::walk(c("con_group_exam_fits", "var_group_exam_fits", "hybrid_group_exam_fits"), 
            ~ list2env(readRDS(here::here(paste0("data/model_cache/", .x, ".rds"))), 
            envir = .GlobalEnv))

# pluck(ex_te_v, "Fit") |> mutate(w= ifelse(exists("w"), round(w,2),NA))
# pluck(hybrid_te_v, "Fit") |> mutate(w= ifelse(exists("w"), round(w,2), NA))
#
#
#
#
#| label: fig-alm-diagram
#| fig.cap: The basic structure of the ALM model. 
 alm_plot()
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
#| label: tbl-e1-cogmodel
#| tbl-cap: Fit Parameters and Model RMSE. The Test_RMSE column is the main performance indicator of interest, and represents the RMSE for just the testing data. The Fit_Method column indicates the data used to fit the model.
#| column: screen-inset-right
#| 
create_combined_df <- function(model_names, model_data, group) {
  combined_df <- do.call(rbind, Map(function(name, data) {
    model<- ifelse(grepl("Hybrid", name), "Hybrid", ifelse(grepl("EXAM", name), "EXAM", "ALM"))
    fit_method <- gsub(".*(Test Only|Test & Train|Train Only).*", "\\1", name)
    extract_params(model, data, group, fit_method)
  }, model_names, model_data))
  row.names(combined_df) <- NULL
  combined_df
}

# Adjust the extract_params function to accept and include the fit_method parameter
extract_params <- function(model_name, model_data, group, fit_method) {
  params_df <- cbind(Model = model_name, Group = group, Fit_Method = fit_method,
                     pluck(model_data, "Fit"), pluck(model_data, "test") |>
                     summarise(Test_RMSE = round(RMSE(y, pred),0))) |>
                     mutate(across(where(is.numeric), \(x) round(x, 3))) |>
                     mutate(w= ifelse(exists("w"), round(w,2), NA))
  return(params_df)
}

model_classes <- c("ALM", "EXAM", "Hybrid")

model_names <- c("ALM Test Only", "ALM Test & Train", "ALM Train Only", "EXAM Test Only", "EXAM Test & Train", "EXAM Train Only", "Hybrid Test Only", "Hybrid Test & Train", "Hybrid Train Only")

model_data_v <- list(a_te_v, a_tetr_v, a_tr_v, ex_te_v, ex_tetr_v, ex_tr_v, hybrid_te_v, hybrid_tetr_v, hybrid_tr_v)
model_data_c <- list(a_te_c, a_tetr_c, a_tr_c, ex0_te_c, ex0_tetr_c, ex0_tr_c, hybrid_te_c, hybrid_tetr_c, hybrid_tr_c)

combined_params_v <- create_combined_df(model_names, model_data_v, "Varied")
combined_params_c <- create_combined_df(model_names, model_data_c, "Constant")
all_combined_params <- rbind(combined_params_v, combined_params_c)



ft <- all_combined_params %>% select(-Value) |>
  rename("Fit Method" = Fit_Method, "Test RMSE" = Test_RMSE) |>
  pivot_longer(cols=c(c,lr,w,"Test RMSE"),names_to="Parameter") |>
    unite(Group, Group, Parameter) %>%
     pivot_wider(names_from = Group, values_from = value) %>%
     flextable::as_flextable() %>%
     ftExtra::span_header(sep = "\\_") |> align(i=1,align='center',part="header") |> 
  empty_blanks()  %>% 
  hline(part = "header", i = 1, j=3:6) %>% 
  hline(part = "header", i = 1, j=7:10)
ft
#
#
#
#
#| column: screen-inset-right
# Reshape all_combined_params
reshaped_df <- all_combined_params %>%
  select(-Value) |>
  rename("Fit Method" = Fit_Method, "Test RMSE" = Test_RMSE) |>
  pivot_longer(cols=c(c,lr,w,"Test RMSE"),names_to="Parameter") |>
  unite(Group, Group, Parameter) %>%
  pivot_wider(names_from = Group, values_from = value)

header_df <- data.frame(
  col_keys = c("Model", "Fit Method", "Varied_c", "Varied_lr", "Varied_w", "Varied_Test RMSE", 
               "Constant_c", "Constant_lr", "Constant_w", "Constant_Test RMSE"),
  line1 = c("", "", "Varied", "", "", "", "Constant", "", "", ""),
  line2 = c("Model", "Fit Method", "c", "lr", "w", "Test RMSE", "c", "lr", "w", "Test RMSE")
)

# Apply the custom header to the flextable
ft2 <- flextable(reshaped_df) %>% 
  set_header_df(
    mapping = header_df,
    key = "col_keys"
  ) %>% add_header_lines(values = " ") %>%
  theme_booktabs() %>% 
  merge_v(part = "header") %>% 
  merge_h(part = "header") %>%
  merge_h(part = "header") %>%
  align(align = "center", part = "all") %>% 
  autofit() %>% 
  empty_blanks() %>% 
  fix_border_issues() %>% 
  hline(part = "header", i = 2, j=3:6) %>% 
  hline(part = "header", i = 2, j=7:10)



# Display the flextable
ft2


#
#
#
#
#
#| column: screen-inset-right
#| 
flextable_to_markdown <- function(ft) {
  header <- ft$header$dataset %>% `rownames<-`( NULL )  %>% `colnames<-`( NULL )
  body <- ft$body$dataset %>% `rownames<-`( NULL )  %>% `colnames<-`( NULL )
  col_names <- paste0("V", 1:ncol(header))
  colnames(header) <- col_names
  colnames(body) <- col_names
  combined_df <- rbind(header, body) %>% `rownames<-`( NULL )  %>% `colnames<-`( NULL )
  markdown_table_k <- knitr::kable(combined_df, format="markdown")
  markdown_table_p <- pander::pandoc.table(combined_df,style="rmarkdown", split.table = Inf)
  return(tibble::lst(markdown_table_k,markdown_table_p))
}
flextable_to_markdown(ft)

#
#
#
#
#
#
#| label: tbl-e1-hybridParam
#| tbl-cap: Hybrid Model - w parameter. The $w$ parameter determines the balance between the ALM and EXAM response generation processes, and is only included for the hybrid model. A weight of .5 would indicate equal contribution from both models. $w$ values approaching 1 indicate stronger weight for EXAM. 
#| eval: false
all_combined_params |> filter(Model =="Hybrid") |> 
  group_by(Group, Fit_Method) |>
  summarise(w=fmt_iden(first(w)), Test_RMSE=(fmt_iden(Test_RMSE))) |>
  group_by(Group,Fit_Method) |> flextable()

#
#
#
#
#
#| label: fig-model-preds-varied
#| fig-cap: Varied Group - Mean Model predictions vs. observations
#| fig-height: 12
#| fig-width: 14
#| column: screen-inset-right

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
#| label: tbl-e1-predsV
#| tbl-cap: Varied group - mean model predictions vs. observations
#| 
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
#| label: fig-model-preds-constant
#| fig-cap: Constant Group - Mean Model predictions vs. observations
#| fig-height: 12
#| fig-width: 14
#| column: screen-inset-right

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
#| label: tbl-e1-predsC
#| tbl-cap: Constant group - mean model predictions vs. observations
#| 
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
#| eval: false
pluck(a_te_v, "train") |> pivot_longer(y:almResp, names_to="Resp", values_to = "vx") |> 
  mutate(dev=x-vx,abs_dev=abs(x-vx)) |>
  learn_curve_plot(tr, vx, Resp,facet_var=x, groupVec=Resp,nbins=8)

pluck(a_te_v, "train") |> pivot_longer(y:almResp, names_to="Resp", values_to = "vx") |> 
  mutate(dev=x-vx,abs_dev=abs(x-vx)) |>
  ungroup() %>%
  gather(key = "variable", value = "y_value", dev, abs_dev, vx) %>%
  group_by(variable) %>%
  group_map(~ learn_curve_plot(.x, x_var = tr, y_var = y_value, color_var = Resp, facet_var = x, groupVec = Resp, nbins = 8, y_label = .y$variable), .keep = TRUE)

list(a_tr_v, a_te_v,a_tetr_v) |> map( ~{pluck(.x, "train") |> pivot_longer(y:almResp, names_to="Resp", values_to = "vx") |> mutate(dev=x-vx,abs_dev=abs(x-vx)) |>
  ungroup() %>%
  gather(key = "variable", value = "y_value", dev, abs_dev, vx) %>%
  group_by(variable) %>%
  group_map(~ learn_curve_plot(.x, x_var = tr, y_var = y_value, color_var = Resp, facet_var = x, groupVec = Resp, nbins = 8, y_label = .y$variable), .keep = TRUE) })

#
#
#
#
#
#
#| eval: false
pluck(ex_te_v, "train") |> pivot_longer(y:almResp, names_to="Resp", values_to = "vx") |> 
  mutate(dev=x-vx,abs_dev=abs(x-vx)) |>
  learn_curve_plot(tr, vx, Resp,facet_var=x, groupVec=Resp,nbins=8)

pluck(ex_te_v, "train") |> pivot_longer(y:almResp, names_to="Resp", values_to = "vx") |> 
  mutate(dev=x-vx,abs_dev=abs(x-vx)) |>
  ungroup() %>%
  gather(key = "variable", value = "y_value", dev, abs_dev, vx) %>%
  group_by(variable) %>%
  group_map(~ learn_curve_plot(.x, x_var = tr, y_var = y_value, color_var = Resp, facet_var = x, groupVec = Resp, nbins = 8, y_label = .y$variable), .keep = TRUE)

#
#
#
#
#
#| eval: false


optimize_params_weighted_individual <- function(ds, c_values, lr_values, weight_exam_values, input.layer, output.layer) {
    all_results <- list()
    
    # Loop through each unique id
    for (individual in unique(ds$id)) {
        indiv_data <- ds[ds$id == individual, ]
        
        # Run the optimization function for the individual's data
        result <- optimize_params_weighted(indiv_data, c_values, lr_values, weight_exam_values, input.layer, output.layer)
        
        all_results[[as.character(individual)]] <- result
    }
    
    all_results
}

dss <- ds |> filter(id %in% c(1,2))

all_results_weighted_hybrid <- readRDS(here::here('data/model_cache/indv_hybrid_fits.rds'))
ma <- map(all_results_weighted_hybrid, "best_params") |> map("c")


data.frame(id=names(ma),c=as.numeric(ma))

ma = cbind(id=names(all_results_weighted_hybrid),map(all_results_weighted_hybrid, "best_params") |> map_dfr(magrittr::extract,c("c","lr","weight_exam")))

ds |> group_by(id,condit) |> distinct(id,condit) |> left_join(ma,by=join_by(id))


map(all_results_weighted_hybrid,"best_params") |> pluck("c")
all_results_weighted_hybrid[["1"]]$b

 map(~ map(.x$best_params, pluck, "c"))

 map_df(~ map_df(.x$train, pluck, "d"), .id = "density")
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
