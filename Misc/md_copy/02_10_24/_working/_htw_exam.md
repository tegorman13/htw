---
title: EXAM Fits and Predictions
date: last-modified
categories: [Simulation, ALM, EXAM, R]
code-fold: true
code-tools: true
execute: 
  warning: false
  eval: true
---

<link href="../site_libs/tabwid-1.1.3/tabwid.css" rel="stylesheet" />
<script src="../site_libs/tabwid-1.1.3/tabwid.js"></script>


<details class="code-fold">
<summary>Code</summary>

```{r}
# load and view data
pacman::p_load(tidyverse,patchwork,here, pander, latex2exp)
purrr::walk(here::here(c("Functions/Display_Functions.R", "Functions/alm_core.R","Functions/fun_model.R")),source)
select <- dplyr::select; mutate <- dplyr::mutate 

ds <- readRDS(here::here("data/e1_md_11-06-23.rds"))
dsAvg <- ds |> group_by(condit,expMode2,tr, x) |> 
  summarise(y=mean(y),.groups="keep") 

vAvg <- dsAvg |> filter(condit=="Varied")
cAvg <- dsAvg |> filter(condit=="Constant")

i1 <- ds |> filter(id=="3")

input.layer <- c(100,350,600,800,1000,1200)
output.layer <- c(100,350,600,800,1000,1200)


purrr::walk(c("con_group_exam_fits", "var_group_exam_fits", "hybrid_group_exam_fits"), 
            ~ list2env(readRDS(here::here(paste0("data/model_cache/", .x, ".rds"))), 
            envir = .GlobalEnv))
```

</details>
<details class="code-fold">
<summary>Code</summary>

```{r}
 alm_plot()
```

</details>
![](htw_exam.markdown_strict_files/figure-markdown_strict/fig-alm-diagram-1.jpeg)

# Modeling

In project 1, we applied model-based techniques to quantify and control for the similarity between training and testing experience, which in turn enabled us to account for the difference between varied and constant training via an extended version of a similarity based generalization model. In project 2, we will go a step further, implementing a full process model capable of both 1) producing novel responses and 2) modeling behavior in both the learning and testing stages of the experiment. For this purpose, we will apply the associative learning model (ALM) and the EXAM model of function learning (DeLosh 1997). ALM is a simple connectionist learning model which closely resembles Kruschke's ALCOVE model (Kruscke 1992), with modifications to allow for the generation of continuous responses.

## ALM & Exam Description

DeLosh et al. (1997) introduced the associative learning model (ALM), a connectionist model within the popular class of radial-basis networks. ALM was inspired by, and closely resembles Kruschke's influential ALCOVE model of categorization (Kruschke, 1992).

ALM is a localist neural network model, with each input node corresponding to a particular stimulus, and each output node corresponding to a particular response value. The units in the input layer activate as a function of their Gaussian similarity to the input stimulus. So, for example, an input stimulus of value 55 would induce maximal activation of the input unit tuned to 55. Depending on thevalue of the generalization parameter, the nearby units (e.g. 54 and 56; 53 and 57) may also activate to some degree. ALM is structured with input and output nodes that correspond to regions of the stimulus space, and response space, respectively. The units in the input layer activate as a function of their similarity to a presented stimulus. As was the case with the exemplar-based models, similarity in ALM is exponentially decaying function of distance. The input layer is fully connected to the output layer, and the activation for any particular output node is simply the weighted sum of the connection weights between that node and the input activations. The network then produces a response by taking the weighted average of the output units (recall that each output unit has a value corresponding to a particular response). During training, the network receives feedback which activates each output unit as a function of its distance from the ideal level of activation necessary to produce the correct response. The connection weights between input and output units are then updated via the standard delta learning rule, where the magnitude of weight changes are controlled by a learning rate parameter.

See <a href="#tbl-alm-exam" class="quarto-xref">Table 5</a> for a full specification of the equations that define ALM and EXAM.

## Model Fitting and Comparison

Following the procedure used by Mcdaniel et al. (2009), we will assess the ability of both ALM and EXAM to account for the empirical data when fitting the models to 1) only the training data, and 2) both training and testing data. Models were fit to the aggregated participant data by minimizing the root-mean squared deviation (RMSE). Because ALM has been shown to do poorly at accounting for human patterns extrapolation (DeLosh et al., 1997), we will also generate predictions from the EXAM model for the testing stage. EXAM which operates identically to ALM during training, but includes a linear extrapolation mechanism for generating novel responses during testing.

For the hybrid model, predictions are computed by first generating separate predictions from ALM and EXAM, and then combining them using the following equation: $\hat{y} = (1 - w) \cdot alm.pred + w \cdot exam.pred$. For the grid search, the weight parameter is varied from 0 to 1, and the resulting RMSE is recorded.

Each model was fit to the data in 3 different ways. 1) To just the testing data, 2) Both the training and testing data, 3) Only the training data. In all cases, the model only updates its weights during the training phase, and the weights are frozen during the testing phase. In all cases, only the ALM model generates predictions during the training phase. For the testing phase, all 3 models are used to generate predictions.



<details class="code-fold">
<summary>Code</summary>

```{r}
create_combined_df <- function(model_names, model_data, group) {
  do.call(rbind, Map(function(name, data) {
    model<- ifelse(grepl("Hybrid", name), "Hybrid", ifelse(grepl("EXAM", name), "EXAM", "ALM"))
    fit_method <- gsub(".*(Test Only|Test & Train|Train Only).*", "\\1", name)
    extract_params(model, data, group, fit_method)
  }, model_names, model_data))
}

# Adjust the extract_params function to accept and include the fit_method parameter
extract_params <- function(model_name, model_data, group, fit_method) {
  params_df <- cbind(Model = model_name, Group = group, Fit_Method = fit_method,
                     pluck(model_data, "Fit"), pluck(model_data, "test") %>% summarise(Test_RMSE = round(RMSE(y, pred)),1)) %>% 
                     mutate(across(where(is.numeric), round, 4))
  if (!"w" %in% names(params_df)) {
    params_df$w <- NA
  }
  return(params_df)
}


model_classes <- c("ALM", "EXAM", "Hybrid")
groups <- c("V", "C")

model_names_v <- c("ALM Test Only", "ALM Test & Train", "ALM Train Only", "EXAM Test Only", "EXAM Test & Train", "EXAM Train Only", "Hybrid Test Only", "Hybrid Test & Train", "Hybrid Train Only")

#model_names_v <- c("Test Only", "Test & Train", "Train Only")
model_names_c <- model_names_v

model_data_v <- list(a_te_v, a_tetr_v, a_tr_v, ex_te_v, ex_tetr_v, ex_tr_v, hybrid_te_v, hybrid_tetr_v, hybrid_tr_v)
model_data_c <- list(a_te_c, a_tetr_c, a_tr_c, ex0_te_c, ex0_tetr_c, ex0_tr_c, hybrid_te_c, hybrid_tetr_c, hybrid_tr_c)

combined_params_v <- create_combined_df(model_names_v, model_data_v, "Varied")
combined_params_c <- create_combined_df(model_names_c, model_data_c, "Constant")

all_combined_params <- rbind(combined_params_v, combined_params_c)


library(flextable)

identity <- function (x){x}
fmt_iden <- function(x, digits = 2) {
  paste0(x)
}

# all_combined_params %>% 
#   group_by(Group, Model) |> summarise(across(c(c,lr,Test_RMSE), ~max(.x,na.rm=TRUE)))


grouped_params <- all_combined_params %>% 
  group_by(Group, Model) %>%
  summarise(across(c(c, lr, w, Value, Test_RMSE), ~identity(.x))) %>%
  ungroup()

ft <- flextable(all_combined_params) %>% 
  theme_vanilla() %>% 
  align(align = "left", part = "all") %>% 
  colformat_double(digits = 2) %>% 
  autofit()


# Pre-calculate the averages and standard deviations
all_combined_params2 <- all_combined_params %>%
  group_by(Fit_Method, Group, Model) %>%
  summarise(
    c_stats = fmt_iden(c),
    lr_stats = fmt_iden(lr),
     rmse_stats = fmt_iden(Test_RMSE),
    .groups = "drop"
  )

# Use the pre-calculated stats in tabulator
tab <- tabulator(
  x = all_combined_params2, rows = c("Fit_Method", "Model"),
  columns = "Group",
  `c` = as_paragraph(c_stats),
  `lr` = as_paragraph(lr_stats),
  `Test_RMSE` = as_paragraph(rmse_stats)
) %>%
  as_flextable() |> align(align="left",part="all")


tab
# a=aggregate(c ~ Group+Model+Fit_Method, data=all_combined_params,identity)
# 
# all_combined_params |> group_by(Model, Group,Fit_Method) %>% 
#   summarise(
#     across(
#       all_of(c("c", "lr")),
#       list(
#         c = ~ identity(.x),
#         lr = ~ identity(.x)
#       )
#     )
#   ) %>% flextable() %>% separate_header()
# 
# 
# ft_1 <-  tabulator(
#   x = all_combined_params,
#   rows = "Group",
#   columns = c("Fit_Method","Model"),
#   `c1` = c,
#   `$lr$` = as_paragraph(lr) ) %>% as_flextable()
# 
# ft_1
```

</details>

<div class="tabwid"><style>.cl-d4b05f08{}.cl-d4ac1f4c{font-family:'Helvetica';font-size:11pt;font-weight:normal;font-style:normal;text-decoration:none;color:rgba(0, 0, 0, 1.00);background-color:transparent;}.cl-d4ae178e{margin:0;text-align:left;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);padding-bottom:5pt;padding-top:5pt;padding-left:5pt;padding-right:5pt;line-height: 1;background-color:transparent;}.cl-d4ae178f{margin:0;text-align:left;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);padding-bottom:0;padding-top:0;padding-left:0;padding-right:0;line-height: 1;background-color:transparent;}.cl-d4ae1790{margin:0;text-align:left;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);padding-bottom:5pt;padding-top:5pt;padding-left:5pt;padding-right:5pt;line-height: 1;background-color:transparent;}.cl-d4ae1798{margin:0;text-align:left;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);padding-bottom:5pt;padding-top:5pt;padding-left:5pt;padding-right:5pt;line-height: 1;background-color:transparent;}.cl-d4ae1799{margin:0;text-align:left;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);padding-bottom:0;padding-top:0;padding-left:0;padding-right:0;line-height: 1;background-color:transparent;}.cl-d4ae2210{width:1.096in;background-color:transparent;vertical-align: bottom;border-bottom: 1.5pt solid rgba(102, 102, 102, 1.00);border-top: 1.5pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-d4ae221a{width:0.705in;background-color:transparent;vertical-align: bottom;border-bottom: 1.5pt solid rgba(102, 102, 102, 1.00);border-top: 1.5pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-d4ae221b{width:0.079in;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-d4ae221c{width:0.731in;background-color:transparent;vertical-align: middle;border-bottom: 0.75pt solid rgba(102, 102, 102, 1.00);border-top: 1.5pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-d4ae221d{width:1.087in;background-color:transparent;vertical-align: middle;border-bottom: 0.75pt solid rgba(102, 102, 102, 1.00);border-top: 1.5pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-d4ae2224{width:1.096in;background-color:transparent;vertical-align: bottom;border-bottom: 1.5pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-d4ae2225{width:0.705in;background-color:transparent;vertical-align: bottom;border-bottom: 1.5pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-d4ae2226{width:0.079in;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-d4ae2227{width:0.731in;background-color:transparent;vertical-align: middle;border-bottom: 1.5pt solid rgba(102, 102, 102, 1.00);border-top: 0.75pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-d4ae2228{width:1.087in;background-color:transparent;vertical-align: middle;border-bottom: 1.5pt solid rgba(102, 102, 102, 1.00);border-top: 0.75pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-d4ae222e{width:1.096in;background-color:transparent;vertical-align: top;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-d4ae222f{width:0.705in;background-color:transparent;vertical-align: top;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-d4ae2230{width:0.079in;background-color:transparent;vertical-align: top;border-bottom: 0 solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-d4ae2231{width:0.731in;background-color:transparent;vertical-align: top;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-d4ae2232{width:1.087in;background-color:transparent;vertical-align: top;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-d4ae2233{width:1.096in;background-color:transparent;vertical-align: top;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-d4ae2234{width:0.705in;background-color:transparent;vertical-align: top;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-d4ae2238{width:0.079in;background-color:transparent;vertical-align: top;border-bottom: 0 solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-d4ae2239{width:0.731in;background-color:transparent;vertical-align: top;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-d4ae223a{width:1.087in;background-color:transparent;vertical-align: top;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-d4ae223b{width:1.096in;background-color:transparent;vertical-align: top;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-d4ae223c{width:0.705in;background-color:transparent;vertical-align: top;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-d4ae223d{width:0.079in;background-color:transparent;vertical-align: top;border-bottom: 0 solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-d4ae2242{width:0.731in;background-color:transparent;vertical-align: top;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-d4ae2243{width:1.087in;background-color:transparent;vertical-align: top;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-d4ae2244{width:1.096in;background-color:transparent;vertical-align: top;border-bottom: 1.5pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-d4ae2245{width:0.705in;background-color:transparent;vertical-align: top;border-bottom: 1.5pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-d4ae2246{width:0.731in;background-color:transparent;vertical-align: top;border-bottom: 1.5pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-d4ae224c{width:1.087in;background-color:transparent;vertical-align: top;border-bottom: 1.5pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}</style><table data-quarto-disable-processing='true' class='cl-d4b05f08'><thead><tr style="overflow-wrap:break-word;"><th  rowspan="2"class="cl-d4ae2210"><p class="cl-d4ae178e"><span class="cl-d4ac1f4c">Fit_Method</span></p></th><th  rowspan="2"class="cl-d4ae221a"><p class="cl-d4ae178e"><span class="cl-d4ac1f4c">Model</span></p></th><th class="cl-d4ae221b"><p class="cl-d4ae178f"><span class="cl-d4ac1f4c"></span></p></th><th  colspan="3"class="cl-d4ae221c"><p class="cl-d4ae1790"><span class="cl-d4ac1f4c">Constant</span></p></th><th class="cl-d4ae221b"><p class="cl-d4ae178f"><span class="cl-d4ac1f4c"></span></p></th><th  colspan="3"class="cl-d4ae221c"><p class="cl-d4ae1790"><span class="cl-d4ac1f4c">Varied</span></p></th></tr><tr style="overflow-wrap:break-word;"><th class="cl-d4ae2226"><p class="cl-d4ae178f"><span class="cl-d4ac1f4c"></span></p></th><th class="cl-d4ae2227"><p class="cl-d4ae1790"><span class="cl-d4ac1f4c">c</span></p></th><th class="cl-d4ae2227"><p class="cl-d4ae1790"><span class="cl-d4ac1f4c">lr</span></p></th><th class="cl-d4ae2228"><p class="cl-d4ae1790"><span class="cl-d4ac1f4c">Test_RMSE</span></p></th><th class="cl-d4ae2226"><p class="cl-d4ae178f"><span class="cl-d4ac1f4c"></span></p></th><th class="cl-d4ae2227"><p class="cl-d4ae1790"><span class="cl-d4ac1f4c">c</span></p></th><th class="cl-d4ae2227"><p class="cl-d4ae1790"><span class="cl-d4ac1f4c">lr</span></p></th><th class="cl-d4ae2228"><p class="cl-d4ae1790"><span class="cl-d4ac1f4c">Test_RMSE</span></p></th></tr></thead><tbody><tr style="overflow-wrap:break-word;"><td  rowspan="3"class="cl-d4ae222e"><p class="cl-d4ae1798"><span class="cl-d4ac1f4c">Test &amp; Train</span></p></td><td class="cl-d4ae222f"><p class="cl-d4ae1798"><span class="cl-d4ac1f4c">ALM</span></p></td><td class="cl-d4ae2230"><p class="cl-d4ae1799"><span class="cl-d4ac1f4c"></span></p></td><td class="cl-d4ae2231"><p class="cl-d4ae1798"><span class="cl-d4ac1f4c">0.047</span></p></td><td class="cl-d4ae2231"><p class="cl-d4ae1798"><span class="cl-d4ac1f4c">0.0804</span></p></td><td class="cl-d4ae2232"><p class="cl-d4ae1798"><span class="cl-d4ac1f4c">329</span></p></td><td class="cl-d4ae2230"><p class="cl-d4ae1799"><span class="cl-d4ac1f4c"></span></p></td><td class="cl-d4ae2231"><p class="cl-d4ae1798"><span class="cl-d4ac1f4c">0.0671</span></p></td><td class="cl-d4ae2231"><p class="cl-d4ae1798"><span class="cl-d4ac1f4c">0.1005</span></p></td><td class="cl-d4ae2232"><p class="cl-d4ae1798"><span class="cl-d4ac1f4c">107</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-d4ae222f"><p class="cl-d4ae1798"><span class="cl-d4ac1f4c">EXAM</span></p></td><td class="cl-d4ae2230"><p class="cl-d4ae1799"><span class="cl-d4ac1f4c"></span></p></td><td class="cl-d4ae2231"><p class="cl-d4ae1798"><span class="cl-d4ac1f4c">0.0805</span></p></td><td class="cl-d4ae2231"><p class="cl-d4ae1798"><span class="cl-d4ac1f4c">0.1608</span></p></td><td class="cl-d4ae2232"><p class="cl-d4ae1798"><span class="cl-d4ac1f4c">132</span></p></td><td class="cl-d4ae2230"><p class="cl-d4ae1799"><span class="cl-d4ac1f4c"></span></p></td><td class="cl-d4ae2231"><p class="cl-d4ae1798"><span class="cl-d4ac1f4c">0.0738</span></p></td><td class="cl-d4ae2231"><p class="cl-d4ae1798"><span class="cl-d4ac1f4c">0.1005</span></p></td><td class="cl-d4ae2232"><p class="cl-d4ae1798"><span class="cl-d4ac1f4c">60</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-d4ae2234"><p class="cl-d4ae1798"><span class="cl-d4ac1f4c">Hybrid</span></p></td><td class="cl-d4ae2238"><p class="cl-d4ae1799"><span class="cl-d4ac1f4c"></span></p></td><td class="cl-d4ae2239"><p class="cl-d4ae1798"><span class="cl-d4ac1f4c">0.0672</span></p></td><td class="cl-d4ae2239"><p class="cl-d4ae1798"><span class="cl-d4ac1f4c">0.1345</span></p></td><td class="cl-d4ae223a"><p class="cl-d4ae1798"><span class="cl-d4ac1f4c">136</span></p></td><td class="cl-d4ae2238"><p class="cl-d4ae1799"><span class="cl-d4ac1f4c"></span></p></td><td class="cl-d4ae2239"><p class="cl-d4ae1798"><span class="cl-d4ac1f4c">0.1345</span></p></td><td class="cl-d4ae2239"><p class="cl-d4ae1798"><span class="cl-d4ac1f4c">2.0168</span></p></td><td class="cl-d4ae223a"><p class="cl-d4ae1798"><span class="cl-d4ac1f4c">47</span></p></td></tr><tr style="overflow-wrap:break-word;"><td  rowspan="3"class="cl-d4ae223b"><p class="cl-d4ae1798"><span class="cl-d4ac1f4c">Test Only</span></p></td><td class="cl-d4ae223c"><p class="cl-d4ae1798"><span class="cl-d4ac1f4c">ALM</span></p></td><td class="cl-d4ae223d"><p class="cl-d4ae1799"><span class="cl-d4ac1f4c"></span></p></td><td class="cl-d4ae2242"><p class="cl-d4ae1798"><span class="cl-d4ac1f4c">0</span></p></td><td class="cl-d4ae2242"><p class="cl-d4ae1798"><span class="cl-d4ac1f4c">0.1005</span></p></td><td class="cl-d4ae2243"><p class="cl-d4ae1798"><span class="cl-d4ac1f4c">348</span></p></td><td class="cl-d4ae223d"><p class="cl-d4ae1799"><span class="cl-d4ac1f4c"></span></p></td><td class="cl-d4ae2242"><p class="cl-d4ae1798"><span class="cl-d4ac1f4c">0.1342</span></p></td><td class="cl-d4ae2242"><p class="cl-d4ae1798"><span class="cl-d4ac1f4c">2.0302</span></p></td><td class="cl-d4ae2243"><p class="cl-d4ae1798"><span class="cl-d4ac1f4c">95</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-d4ae223c"><p class="cl-d4ae1798"><span class="cl-d4ac1f4c">EXAM</span></p></td><td class="cl-d4ae223d"><p class="cl-d4ae1799"><span class="cl-d4ac1f4c"></span></p></td><td class="cl-d4ae2242"><p class="cl-d4ae1798"><span class="cl-d4ac1f4c">0.0067</span></p></td><td class="cl-d4ae2242"><p class="cl-d4ae1798"><span class="cl-d4ac1f4c">1.3266</span></p></td><td class="cl-d4ae2243"><p class="cl-d4ae1798"><span class="cl-d4ac1f4c">127</span></p></td><td class="cl-d4ae223d"><p class="cl-d4ae1799"><span class="cl-d4ac1f4c"></span></p></td><td class="cl-d4ae2242"><p class="cl-d4ae1798"><span class="cl-d4ac1f4c">0.4094</span></p></td><td class="cl-d4ae2242"><p class="cl-d4ae1798"><span class="cl-d4ac1f4c">1.9096</span></p></td><td class="cl-d4ae2243"><p class="cl-d4ae1798"><span class="cl-d4ac1f4c">46</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-d4ae223c"><p class="cl-d4ae1798"><span class="cl-d4ac1f4c">Hybrid</span></p></td><td class="cl-d4ae223d"><p class="cl-d4ae1799"><span class="cl-d4ac1f4c"></span></p></td><td class="cl-d4ae2242"><p class="cl-d4ae1798"><span class="cl-d4ac1f4c">0.0084</span></p></td><td class="cl-d4ae2242"><p class="cl-d4ae1798"><span class="cl-d4ac1f4c">1.5798</span></p></td><td class="cl-d4ae2243"><p class="cl-d4ae1798"><span class="cl-d4ac1f4c">127</span></p></td><td class="cl-d4ae223d"><p class="cl-d4ae1799"><span class="cl-d4ac1f4c"></span></p></td><td class="cl-d4ae2242"><p class="cl-d4ae1798"><span class="cl-d4ac1f4c">0.395</span></p></td><td class="cl-d4ae2242"><p class="cl-d4ae1798"><span class="cl-d4ac1f4c">2.0168</span></p></td><td class="cl-d4ae2243"><p class="cl-d4ae1798"><span class="cl-d4ac1f4c">34</span></p></td></tr><tr style="overflow-wrap:break-word;"><td  rowspan="3"class="cl-d4ae2244"><p class="cl-d4ae1798"><span class="cl-d4ac1f4c">Train Only</span></p></td><td class="cl-d4ae223c"><p class="cl-d4ae1798"><span class="cl-d4ac1f4c">ALM</span></p></td><td class="cl-d4ae223d"><p class="cl-d4ae1799"><span class="cl-d4ac1f4c"></span></p></td><td class="cl-d4ae2242"><p class="cl-d4ae1798"><span class="cl-d4ac1f4c">0.0604</span></p></td><td class="cl-d4ae2242"><p class="cl-d4ae1798"><span class="cl-d4ac1f4c">0.1005</span></p></td><td class="cl-d4ae2243"><p class="cl-d4ae1798"><span class="cl-d4ac1f4c">329</span></p></td><td class="cl-d4ae223d"><p class="cl-d4ae1799"><span class="cl-d4ac1f4c"></span></p></td><td class="cl-d4ae2242"><p class="cl-d4ae1798"><span class="cl-d4ac1f4c">0.047</span></p></td><td class="cl-d4ae2242"><p class="cl-d4ae1798"><span class="cl-d4ac1f4c">0.0804</span></p></td><td class="cl-d4ae2243"><p class="cl-d4ae1798"><span class="cl-d4ac1f4c">109</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-d4ae223c"><p class="cl-d4ae1798"><span class="cl-d4ac1f4c">EXAM</span></p></td><td class="cl-d4ae223d"><p class="cl-d4ae1799"><span class="cl-d4ac1f4c"></span></p></td><td class="cl-d4ae2242"><p class="cl-d4ae1798"><span class="cl-d4ac1f4c">0.0604</span></p></td><td class="cl-d4ae2242"><p class="cl-d4ae1798"><span class="cl-d4ac1f4c">0.1005</span></p></td><td class="cl-d4ae2243"><p class="cl-d4ae1798"><span class="cl-d4ac1f4c">200</span></p></td><td class="cl-d4ae223d"><p class="cl-d4ae1799"><span class="cl-d4ac1f4c"></span></p></td><td class="cl-d4ae2242"><p class="cl-d4ae1798"><span class="cl-d4ac1f4c">0.047</span></p></td><td class="cl-d4ae2242"><p class="cl-d4ae1798"><span class="cl-d4ac1f4c">0.0804</span></p></td><td class="cl-d4ae2243"><p class="cl-d4ae1798"><span class="cl-d4ac1f4c">65</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-d4ae2245"><p class="cl-d4ae1798"><span class="cl-d4ac1f4c">Hybrid</span></p></td><td class="cl-d4ae223d"><p class="cl-d4ae1799"><span class="cl-d4ac1f4c"></span></p></td><td class="cl-d4ae2246"><p class="cl-d4ae1798"><span class="cl-d4ac1f4c">0.042</span></p></td><td class="cl-d4ae2246"><p class="cl-d4ae1798"><span class="cl-d4ac1f4c">0.0672</span></p></td><td class="cl-d4ae224c"><p class="cl-d4ae1798"><span class="cl-d4ac1f4c">330</span></p></td><td class="cl-d4ae223d"><p class="cl-d4ae1799"><span class="cl-d4ac1f4c"></span></p></td><td class="cl-d4ae2246"><p class="cl-d4ae1798"><span class="cl-d4ac1f4c">0.042</span></p></td><td class="cl-d4ae2246"><p class="cl-d4ae1798"><span class="cl-d4ac1f4c">0.0672</span></p></td><td class="cl-d4ae224c"><p class="cl-d4ae1798"><span class="cl-d4ac1f4c">110</span></p></td></tr></tbody></table></div>

<details class="code-fold">
<summary>Code</summary>

```{r}
extract_params <- function(model_name, model_data) {
  cbind(Model = model_name,
        pluck(model_data, "Fit"),
        pluck(model_data, "test") %>% summarise(Test_RMSE = RMSE(y, pred))
       ) %>%
  mutate(across(where(is.numeric), \(x) round(x, 3)))
}

# Define model names and data
# model_names <- c("ALM Test Only", "ALM Test & Train", "ALM Train Only", "EXAM Test Only", "EXAM Test & Train", "EXAM Train Only", "Hybrid Test Only", "Hybrid Test & Train", "Hybrid Train Only")
# model_data_v <- list(a_te_v, a_tetr_v, a_tr_v, ex_te_v, ex_tetr_v, ex_tr_v, hybrid_te_v, hybrid_tetr_v, hybrid_tr_v)
# model_data_c <- list(a_te_c, a_tetr_c, a_tr_c, ex0_te_c, ex0_tetr_c, ex0_tr_c, hybrid_te_c, hybrid_tetr_c, hybrid_tr_c)


# params_list_v <- Map(extract_params, model_names, model_data_v)
# params_list_c <- Map(extract_params, model_names, model_data_c)

# almParamsV <- do.call(rbind, params_list_v[1:3])
# examParamsV <- do.call(rbind, params_list_v[4:6])
# hybridParamsV <- do.call(rbind, params_list_v[7:9])

# almParamsC <- do.call(rbind, params_list_c[1:3])
# examParamsC <- do.call(rbind, params_list_c[4:6])
# hybridParamsC <- do.call(rbind, params_list_c[7:9])



all_combined_params |> filter(Model =="Hybrid") |> 
  group_by(Group, Fit_Method) |>
  summarise(w=fmt_iden(first(w)), Test_RMSE=(fmt_iden(Test_RMSE))) |>
  group_by(Group,Fit_Method) |> flextable()
```

</details>

<div class="tabwid"><style>.cl-d4bdbe14{}.cl-d4ba2b5a{font-family:'Helvetica';font-size:11pt;font-weight:normal;font-style:normal;text-decoration:none;color:rgba(0, 0, 0, 1.00);background-color:transparent;}.cl-d4bb7b36{margin:0;text-align:left;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);padding-bottom:5pt;padding-top:5pt;padding-left:5pt;padding-right:5pt;line-height: 1;background-color:transparent;}.cl-d4bb83b0{width:0.75in;background-color:transparent;vertical-align: middle;border-bottom: 1.5pt solid rgba(102, 102, 102, 1.00);border-top: 1.5pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-d4bb83ba{width:0.75in;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-d4bb83bb{width:0.75in;background-color:transparent;vertical-align: middle;border-bottom: 1.5pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}</style><table data-quarto-disable-processing='true' class='cl-d4bdbe14'><thead><tr style="overflow-wrap:break-word;"><th class="cl-d4bb83b0"><p class="cl-d4bb7b36"><span class="cl-d4ba2b5a">Group</span></p></th><th class="cl-d4bb83b0"><p class="cl-d4bb7b36"><span class="cl-d4ba2b5a">Fit_Method</span></p></th><th class="cl-d4bb83b0"><p class="cl-d4bb7b36"><span class="cl-d4ba2b5a">w</span></p></th><th class="cl-d4bb83b0"><p class="cl-d4bb7b36"><span class="cl-d4ba2b5a">Test_RMSE</span></p></th></tr></thead><tbody><tr style="overflow-wrap:break-word;"><td class="cl-d4bb83ba"><p class="cl-d4bb7b36"><span class="cl-d4ba2b5a">Constant</span></p></td><td class="cl-d4bb83ba"><p class="cl-d4bb7b36"><span class="cl-d4ba2b5a">Test &amp; Train</span></p></td><td class="cl-d4bb83ba"><p class="cl-d4bb7b36"><span class="cl-d4ba2b5a">1</span></p></td><td class="cl-d4bb83ba"><p class="cl-d4bb7b36"><span class="cl-d4ba2b5a">136</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-d4bb83ba"><p class="cl-d4bb7b36"><span class="cl-d4ba2b5a">Constant</span></p></td><td class="cl-d4bb83ba"><p class="cl-d4bb7b36"><span class="cl-d4ba2b5a">Test Only</span></p></td><td class="cl-d4bb83ba"><p class="cl-d4bb7b36"><span class="cl-d4ba2b5a">1</span></p></td><td class="cl-d4bb83ba"><p class="cl-d4bb7b36"><span class="cl-d4ba2b5a">127</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-d4bb83ba"><p class="cl-d4bb7b36"><span class="cl-d4ba2b5a">Constant</span></p></td><td class="cl-d4bb83ba"><p class="cl-d4bb7b36"><span class="cl-d4ba2b5a">Train Only</span></p></td><td class="cl-d4bb83ba"><p class="cl-d4bb7b36"><span class="cl-d4ba2b5a">0</span></p></td><td class="cl-d4bb83ba"><p class="cl-d4bb7b36"><span class="cl-d4ba2b5a">330</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-d4bb83ba"><p class="cl-d4bb7b36"><span class="cl-d4ba2b5a">Varied</span></p></td><td class="cl-d4bb83ba"><p class="cl-d4bb7b36"><span class="cl-d4ba2b5a">Test &amp; Train</span></p></td><td class="cl-d4bb83ba"><p class="cl-d4bb7b36"><span class="cl-d4ba2b5a">0.7857</span></p></td><td class="cl-d4bb83ba"><p class="cl-d4bb7b36"><span class="cl-d4ba2b5a">47</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-d4bb83ba"><p class="cl-d4bb7b36"><span class="cl-d4ba2b5a">Varied</span></p></td><td class="cl-d4bb83ba"><p class="cl-d4bb7b36"><span class="cl-d4ba2b5a">Test Only</span></p></td><td class="cl-d4bb83ba"><p class="cl-d4bb7b36"><span class="cl-d4ba2b5a">0.6429</span></p></td><td class="cl-d4bb83ba"><p class="cl-d4bb7b36"><span class="cl-d4ba2b5a">34</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-d4bb83bb"><p class="cl-d4bb7b36"><span class="cl-d4ba2b5a">Varied</span></p></td><td class="cl-d4bb83bb"><p class="cl-d4bb7b36"><span class="cl-d4ba2b5a">Train Only</span></p></td><td class="cl-d4bb83bb"><p class="cl-d4bb7b36"><span class="cl-d4ba2b5a">0</span></p></td><td class="cl-d4bb83bb"><p class="cl-d4bb7b36"><span class="cl-d4ba2b5a">110</span></p></td></tr></tbody></table></div>

## Varied Testing Predictions

<details class="code-fold">
<summary>Code</summary>

```{r}
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
```

</details>
![](htw_exam.markdown_strict_files/figure-markdown_strict/fig-model-preds-varied-1.jpeg)

## Varied Testing

<details class="code-fold">
<summary>Code</summary>

```{r}
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
```

</details>

|  x   | Observed | ALM  | EXAM  | Hybrid |
|:----:|:--------:|:----:|:-----:|:------:|
| 100  |   663    | 675  | 715.6 | 708.5  |
| 350  |  764.2   | 675  | 817.2 | 792.1  |
| 600  |  883.9   | 675  | 895.1 | 874.7  |
| 800  |   1083   | 1078 | 1000  |  1091  |
| 1000 |   1196   | 1202 | 1199  |  1204  |
| 1200 |   1283   | 1230 | 1282  |  1221  |

Varied fit to test only

|  x   | Observed | ALM  | EXAM  | Hybrid |
|:----:|:--------:|:----:|:-----:|:------:|
| 100  |   663    | 675  | 715.6 | 707.3  |
| 350  |  764.2   | 675  | 817.2 |  788   |
| 600  |  883.9   | 675  |  902  | 851.5  |
| 800  |   1083   | 1000 | 1000  |  1004  |
| 1000 |   1196   | 1163 | 1165  |  1196  |
| 1200 |   1283   | 1191 | 1194  |  1227  |

Varied fit to train and test

|  x   | Observed |  ALM  | EXAM  | Hybrid |
|:----:|:--------:|:-----:|:-----:|:------:|
| 100  |   663    |  675  | 715.6 |  675   |
| 350  |  764.2   |  675  | 817.1 |  675   |
| 600  |  883.9   |  675  | 904.8 |  675   |
| 800  |   1083   | 999.8 | 999.8 | 999.3  |
| 1000 |   1196   | 1150  | 1150  |  1143  |
| 1200 |   1283   | 1180  | 1180  |  1176  |

Varied fit to train only

## Constant Testing Predictions

<details class="code-fold">
<summary>Code</summary>

```{r}
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
```

</details>
![](htw_exam.markdown_strict_files/figure-markdown_strict/fig-model-preds-constant-1.jpeg)

<details class="code-fold">
<summary>Code</summary>

```{r}
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
```

</details>

|  x   | Observed | ALM | EXAM  | Hybrid |
|:----:|:--------:|:---:|:-----:|:------:|
| 100  |  526.7   | 675 | 716.9 | 716.8  |
| 350  |  666.3   | 675 | 821.7 | 821.3  |
| 600  |  779.6   | 675 | 926.6 | 925.7  |
| 800  |   980    | 675 | 1010  |  1009  |
| 1000 |   1163   | 675 | 1094  |  1093  |
| 1200 |   1277   | 675 | 1178  |  1176  |

Constant fit to test only

|  x   | Observed |  ALM  | EXAM  | Hybrid |
|:----:|:--------:|:-----:|:-----:|:------:|
| 100  |  526.7   |  675  | 712.4 | 710.6  |
| 350  |  666.3   |  675  | 806.1 | 799.6  |
| 600  |  779.6   |  675  | 899.7 | 888.6  |
| 800  |   980    | 858.9 | 974.6 | 959.8  |
| 1000 |   1163   |  675  | 1049  |  1031  |
| 1200 |   1277   |  675  | 1124  |  1102  |

Constant fit to train and test

|  x   | Observed | ALM | EXAM | Hybrid |
|:----:|:--------:|:---:|:----:|:------:|
| 100  |  526.7   | 675 | 697  |  675   |
| 350  |  666.3   | 675 | 752  |  675   |
| 600  |  779.6   | 675 | 807  |  675   |
| 800  |   980    | 851 | 851  | 832.7  |
| 1000 |   1163   | 675 | 895  |  675   |
| 1200 |   1277   | 675 | 939  |  675   |

Constant fit to train only

<details class="code-fold">
<summary>Code</summary>

```{r}
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
```

</details>

# EXAM fit learning curves

<details class="code-fold">
<summary>Code</summary>

```{r}
pluck(ex_te_v, "train") |> pivot_longer(y:almResp, names_to="Resp", values_to = "vx") |> 
  mutate(dev=x-vx,abs_dev=abs(x-vx)) |>
  learn_curve_plot(tr, vx, Resp,facet_var=x, groupVec=Resp,nbins=8)

pluck(ex_te_v, "train") |> pivot_longer(y:almResp, names_to="Resp", values_to = "vx") |> 
  mutate(dev=x-vx,abs_dev=abs(x-vx)) |>
  ungroup() %>%
  gather(key = "variable", value = "y_value", dev, abs_dev, vx) %>%
  group_by(variable) %>%
  group_map(~ learn_curve_plot(.x, x_var = tr, y_var = y_value, color_var = Resp, facet_var = x, groupVec = Resp, nbins = 8, y_label = .y$variable), .keep = TRUE)
```

</details>
<details class="code-fold">
<summary>Code</summary>

```{r}
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
```

</details>



## Model Table

### ALM Activation & Response

| Step                          | Equation                                                                            | Description                                                                                                                            |
|----------------|-------------------------|--------------------------------|
| **ALM Activation & Response** |                                                                                     |                                                                                                                                        |
| Input Activation              | $a_i(X) = \frac{e^{-c(X-X_i)^2}}{\sum_{k=1}^M e^{-c(X-X_k)^2}}$                     | Activation of each input node $X_i$, is a function of the Gaussian similarity between the node value and stimulus X.                   |
| Output Activation             | $O_j(X) = \sum_{k=1}^M w_{ji} \cdot a_i(X)$                                         | Activation of each Output unit $O_j$ is the weighted sum of the input activations and association weights.                             |
| Output Probability            | $P[Y_j|X] = \frac{O_j(X)}{\sum_{k=1}^M O_k(X)}$                                     | Each output node has associated response, $Y_j$. The probability of response $Y_j$ is determined by the ratio of output activations.   |
| Mean Output                   | $m(x) = \sum_{j=1}^L Y_j \cdot \frac{O_j(x)}{\sum_{k=1}^M O_k(X)}$                  | The response to stimulus x is the weighted average of the response probabilities.                                                      |
| **ALM Learning**              |                                                                                     |                                                                                                                                        |
| Feedback Activation           | $f_j(Z) = e^{-c(Z-Y_j)^2}$                                                          | After responding, feedback signal Z is presented, activating each output node via the Gaussian similarity to the ideal response.       |
| Update Weights                | $w_{ji}^{new}=w_{ji}+\eta\alpha_{ji}$                                               | Delta rule to update weights. Magnitude of weight changes controlled by learning rate parameter $\alpha$.                              |
| **EXAM**                      |                                                                                     |                                                                                                                                        |
| Extrapolation                 | $P[X_i|X] = \frac{a_i(X)}{\sum_{k=1}^M a_k(X)}$                                     | Novel test stimulus X activates input nodes associated with trained stimuli.                                                           |
|                               | $E[Y|X_i] = m(X_i) + \frac{m(X_{i+1})-m(X_{i-1})}{X_{i+1}-X_{i-1}} \cdot [X - X_i]$ | Slope value computed from nearest training instances and then added to the response associated with the nearest training instance,m(x) |



DeLosh, E. L., McDaniel, M. A., & Busemeyer, J. R. (1997). Extrapolation: The Sine Qua Non for Abstraction in Function Learning. *Journal of Experimental Psychology: Learning, Memory, and Cognition*, *23*(4), 19. <https://doi.org/10.1037/0278-7393.23.4.968>

Kruschke, J. K. (1992). ALCOVE: An exemplar-based connectionist model of Category Learning. *Psychological Review*, *99*(1). <https://doi.org/10.1037/0033-295X.99.1.22>

Mcdaniel, M. A., Dimperio, E., Griego, J. A., & Busemeyer, J. R. (2009). Predicting transfer performance: A comparison of competing function learning models. *Journal of Experimental Psychology. Learning, Memory, and Cognition*, *35*, 173--195. <https://doi.org/10.1037/a0013982>
