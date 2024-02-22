# HTW Model
Thomas Gorman
2024-02-22

<details class="code-fold">
<summary>Code</summary>

``` r
pacman::p_load(dplyr,purrr,tidyr,ggplot2, data.table, here, patchwork, conflicted, 
               stringr,future,furrr, knitr, reactable, flextable,ggstance, htmltools)
#conflict_prefer_all("dplyr", quiet = TRUE)
options(scipen = 999)
walk(c("Display_Functions","fun_alm","fun_indv_fit","fun_model"), ~ source(here::here(paste0("Functions/", .x, ".R"))))
```

</details>
<details class="code-fold">
<summary>Code</summary>

``` r
# alm_plot()
```

</details>



# Modeling Approach

In project 1, I applied model-based techniques to quantify and control
for the similarity between training and testing experience, which in
turn enabled us to account for the difference between varied and
constant training via an extended version of a similarity based
generalization model. In project 2, I will go a step further,
implementing a full process model capable of both 1) producing novel
responses and 2) modeling behavior in both the learning and testing
stages of the experiment. For this purpose, we will apply the
associative learning model (ALM) and the EXAM model of function learning
(DeLosh et al., 1997). ALM is a simple connectionist learning model
which closely resembles Kruschke’s ALCOVE model (Kruschke, 1992), with
modifications to allow for the generation of continuous responses.

## ALM & Exam Description

ALM is a localist neural network model (Page, 2000), with each input
node corresponding to a particular stimulus, and each output node
corresponding to a particular response value. The units in the input
layer activate as a function of their Gaussian similarity to the input
stimulus. So, for example, an input stimulus of value 55 would induce
maximal activation of the input unit tuned to 55. Depending on the value
of the generalization parameter, the nearby units (e.g. 54 and 56; 53
and 57) may also activate to some degree. ALM is structured with input
and output nodes that correspond to regions of the stimulus space, and
response space, respectively. The units in the input layer activate as a
function of their similarity to a presented stimulus. As was the case
with the exemplar-based models, similarity in ALM is exponentially
decaying function of distance. The input layer is fully connected to the
output layer, and the activation for any particular output node is
simply the weighted sum of the connection weights between that node and
the input activations. The network then produces a response by taking
the weighted average of the output units (recall that each output unit
has a value corresponding to a particular response). During training,
the network receives feedback which activates each output unit as a
function of its distance from the ideal level of activation necessary to
produce the correct response. The connection weights between input and
output units are then updated via the standard delta learning rule,
where the magnitude of weight changes are controlled by a learning rate
parameter. The EXAM model is an extension of ALM, with the same learning
rule and representational scheme for input and output units. The primary
difference is that EXAM includes a linear extrapolation mechanism for
generating novel responses during testing, a modification necessary to
account for human extrapolation patterns in past research Brown &
Lacroix (2017). Although this extrapolation rule departs from a strictly
similarity-based generalization mechanism, EXAM is distinct from pure
rule-based models in that it remains constrained by the weights learned
during training.

See <a href="#tbl-alm-exam" class="quarto-xref">Table 1</a> for a full
specification of the equations that define ALM and EXAM.

<div class="column-page-inset-right">

|                    | **ALM Response Generation**                                        |                                                                                               |
|--------------------|--------------------------------------------------------------------|-----------------------------------------------------------------------------------------------|
| Input Activation   | $a_i(X) = \frac{e^{-c(X-X_i)^2}}{\sum_{k=1}^M e^{-c(X-X_k)^2}}$    | Input nodes activate as a function of Gaussian similarity to stimulus                         |
| Output Activation  | $O_j(X) = \sum_{k=1}^M w_{ji} \cdot a_i(X)$                        | Output unit $O_j$ activation is the weighted sum of input activations and association weights |
| Output Probability | $P[Y_j|X] = \frac{O_j(X)}{\sum_{k=1}^M O_k(X)}$                    | The response, $Y_j$ probabilites computed via Luce’s choice rule                              |
| Mean Output        | $m(X) = \sum_{j=1}^L Y_j \cdot \frac{O_j(x)}{\sum_{k=1}^M O_k(X)}$ | Weighted average of probabilities determines response to X                                    |
|                    | **ALM Learning**                                                   |                                                                                               |
| Feedback           | $f_j(Z) = e^{-c(Z-Y_j)^2}$                                         | feedback signal Z computed as similarity between ideal response and observed response         |
| magnitude of error | $\Delta_{ji}=(f_{j}(Z)-o_{j}(X))a_{i}(X)$                          | Delta rule to update weights.                                                                 |
| Update Weights     | $w_{ji}^{new}=w_{ji}+\eta\Delta_{ji}$                              | Updates scaled by learning rate parameter $\eta$.                                             |
|                    | **EXAM Extrapolation**                                             |                                                                                               |
| Instance Retrieval | $P[X_i|X] = \frac{a_i(X)}{\sum_{k=1}^M a_k(X)}$                    | Novel test stimulus $X$ activates input nodes $X_i$                                           |
| Slope Computation  | $S =$ $\frac{m(X_{1})-m(X_{2})}{X_{1}-X_{2}}$                      | Slope value, $S$ computed from nearest training instances                                     |
| Response           | $E[Y|X_i] = m(X_i) + S \cdot [X - X_i]$                            | ALM response $m(X_i)$ adjusted by slope.                                                      |

</div>

## Model Fitting Strategy

To fit ALM and EXAM to our participant data, we employ a similar method
to Mcdaniel et al. (2009), wherein we examine the performance of each
model after being fit to various subsets of the data. Each model was fit
to the data in with separate procedures: 1) fit to maximize predictions
of the testing data, 2) fit to maximize predictions of both the training
and testing data, 3) fit to maximize predictions of the just the
training data. We refer to this fitting manipulations as “Fit Method” in
the tables and figures below. It should be emphasized that for all three
fit methods, the ALM and EXAM models behave identically - with weights
updating only during the training phase.Models to were fit separately to
the data of each individual participant. The free parameters for both
models are the generalization ($c$) and learning rate ($lr$) parameters.
Parameter estimation was performed using approximate bayesian
computation (ABC), which we describe in detail below.

### Approximate Bayesian Computation

To estimate parameters, we used approximate bayesian computation (ABC),
enabling us to obtain an estimate of the posterior distribution of the
generalization and learning rate parameters for each individual. ABC
belongs to the class of simulation-based inference methods (Cranmer et
al., 2020), which have begun being used for parameter estimation in
cognitive modeling relatively recently (Kangasrääsiö et al., 2019;
Turner et al., 2016; Turner & Van Zandt, 2012). Although they can be
applied to any model from which data can be simulated, ABC methods are
most useful for complex models that lack an explicit likelihood function
(e.g. many neural network and evidence accumulation models).

The general ABC procedure is to 1) define a prior distribution over
model parameters. 2) sample candidate parameter values, $\theta^*$, from
the prior. 3) Use $\theta^*$ to generate a simulated dataset,
$Data_{sim}$. 4) Compute a measure of discrepancy between the simulated
and observed datasets, $discrep$($Data_{sim}$, $Data_{obs}$). 5) Accept
$\theta^*$ if the discrepancy is less than the tolerance threshold,
$\epsilon$, otherwise reject $\theta^*$. 6) Repeat until desired number
of posterior samples are obtained.

Although simple in the abstract, implementations of ABC require
researchers to make a number of non-trivial decisions as to i) the
discrepancy function between observed and simulated data, ii) whether to
compute the discrepancy between trial level data, or a summary statistic
of the datasets, iii) the value of the minimum tolerance $\epsilon$
between simulated and observed data. For the present work, we follow the
guidelines from previously published ABC tutorials (Farrell &
Lewandowsky, 2018; Turner & Van Zandt, 2012). For the test stage, we
summarized datasets with mean velocity of each band in the observed
dataset as $V_{obs}^{(k)}$ and in the simulated dataset as
$V_{sim}^{(k)}$, where $k$ represents each of the six velocity bands.
For computing the discrepancy between datasets in the training stage, we
aggregated training trials into three equally sized blocks (separately
for each velocity band in the case of the varied group). After obtaining
the summary statistics of the simulated and observed datasets, the
discrepancy was computed as the mean of the absolute difference between
simulated and observed datasets
(<a href="#eq-discrep-test" class="quarto-xref">Equation 1</a> and
<a href="#eq-discrep-train" class="quarto-xref">Equation 2</a>). For the
models fit to both training and testing data, discrepancies were
computed for both stages, and then averaged together.

<div class="column-page-inset-left">

<span id="eq-discrep-test">$$
discrep_{Test}(Data_{sim}, Data_{obs}) = \frac{1}{6} \sum_{k=1}^{6} |V_{obs}^{(k)} - V_{sim}^{(k)}|
 \qquad(1)$$</span>

<span id="eq-discrep-train">$$
\begin{aligned} \\
discrep_{Train,constant}(Data_{sim}, Data_{obs}) = \frac{1}{N_{blocks}} \sum_{j=1}^{N_{blocks}} |V_{obs,constant}^{(j)} - V_{sim,constant}^{(j)}| \\ \\
discrep_{Train,varied}(Data_{sim}, Data_{obs}) = \frac{1}{N_{blocks} \times 3} \sum_{j=1}^{N_{blocks}} \sum_{k=1}^{3} |V_{obs,varied}^{(j,k)} - V_{sim,varied}^{(j,k)}|
\end{aligned}
 \qquad(2)$$</span>

</div>

The final component of our ABC implementation is the determination of an
appropriate value of $\epsilon$. The setting of $\epsilon$ exerts strong
influence on the approximated posterior distribution. Smaller values of
$\epsilon$ increase the rejection rate, and improve the fidelity of the
approximated posterior, while larger values result in an ABC sampler
that simply reproduces the prior distribution. Because the individual
participants in our dataset differed substantially in terms of the
noisiness of their data, we employed an adaptive tolerance setting
strategy to tailor $\epsilon$ to each individual. The initial value of
$\epsilon$ was set to the overall standard deviation of each individuals
velocity values. Thus, sampled parameter values that generated simulated
data within a standard deviation of the observed data were accepted,
while worse performing parameters were rejected. After every 300 samples
the tolerance was allowed to increase only if the current acceptance
rate of the algorithm was less than 1%. In such cases, the tolerance was
shifted towards the average discrepancy of the 5 best samples obtained
thus far. To ensure the acceptance rate did not become overly
permissive, $\epsilon$ was also allowed to decrease every time a sample
was accepted into the posterior.

For each of the 156 participants from Experiment 1, the ABC algorithm
was run until 200 samples of parameters were accepted into the posterior
distribution. Obtaining this number of posterior samples required an
average of 205,000 simulation runs per participant. Fitting each
combination of participant, Model (EXAM & ALM), and fitting method (Test
only, Train only, Test & Train) required a total of 192 million
simulation runs. To facilitate these intensive computational demands, we
used the Future Package in R (Bengtsson, 2021), allowing us to
parallelize computations across a cluster of ten M1 iMacs, each with 8
cores.

### Modelling Results

<details class="code-fold">
<summary>Code</summary>

``` r
ds <- readRDS(here::here("data/e1_md_11-06-23.rds"))  |> as.data.table()
nbins <- 3

fd <- readRDS(here("data/e1_08-21-23.rds"))
test <- fd |> filter(expMode2 == "Test") 
testAvg <- test %>% group_by(id, condit, vb, bandInt,bandType,tOrder) %>%
  summarise(nHits=sum(dist==0),vx=mean(vx),dist=mean(dist),sdist=mean(sdist),n=n(),Percent_Hit=nHits/n)

trainAvg <- fd |> filter(expMode2 == "Train") |> group_by(id) |> 
  mutate(tr=trial,x=vb,Block=case_when(expMode2=="Train" ~ cut(tr,breaks=seq(1,max(tr), length.out=nbins+1),include.lowest=TRUE,labels=FALSE),
                                         expMode2=="Test" ~ 4)) |> 
  group_by(id,condit,vb,x,Block) |> 
  summarise(dist=mean(dist),y=mean(vx))

input_layer <<- output_layer <<-  c(100,350,600,800,1000,1200)
ids2 <- c(1,66,36)
file_name <- "n_iter_200_ntry_300_5354"
#file_name <- "n_iter_400_ntry_100_2944"


ind_fits <- map(list.files(here(paste0('data/abc_reject/'),file_name),full.names=TRUE), readRDS)
ind_fits_df <- ind_fits |> map(~list(dat=.x[[1]], Model = .x[["Model"]], Fit_Method=.x[["Fit_Method"]]))
ind_fits_df <- ind_fits_df |> map(~rbindlist(.x$dat) |> mutate(Model = .x$Model, Fit_Method = .x$Fit_Method)) |> rbindlist() 
```

</details>
<details class="code-fold">
<summary>Code</summary>

``` r
generate_data <- function(Model, post_samples, data, num_samples = 1, return_dat = "train_data, test_data") {
  # Filter data for the specific id without invalidating selfref
  sbj_data <- copy(data[id == post_samples$id[1]])
  simulation_function <- ifelse(Model == "EXAM", full_sim_exam, full_sim_alm)

  target_data <- switch(return_dat,
                        "test_data" = copy(sbj_data[expMode2 == "Test"]),
                        "train_data" = copy(sbj_data[expMode2 == "Train"]),
                        "train_data, test_data" = copy(sbj_data[expMode2 %in% c("Test", "Train")]))
  
  post_samples <- post_samples[order(mean_error)][1:num_samples, .(c, lr, mean_error, rank = .I)]

  simulated_data_list <- lapply(1:nrow(post_samples), function(i) {
    params <- post_samples[i]
    sim_data <- simulation_function(sbj_data, params$c, params$lr, input_layer = input_layer, 
                                    output_layer = output_layer, return_dat = return_dat)
    sim_data_dt <- data.table(id = sbj_data$id[1], condit = sbj_data$condit[1], 
                              expMode2 = target_data$expMode2, Model = Model,tr=target_data$tr,
                              y = target_data$y, x = target_data$x, c = params$c, 
                              lr = params$lr, mean_error = params$mean_error, rank = i,
                              pred = sim_data)
    return(sim_data_dt)
  })
  
  result_dt <- rbindlist(simulated_data_list)
  setcolorder(result_dt, c("id", "condit", "expMode2","tr", "c", "lr", "x", "y", "pred"))
  return(result_dt)
}

#future::plan(multisession)

nestSbjModelFit <- ind_fits_df %>% nest(.by=c(id,Model,Fit_Method))

# organize test data predictions
# post_dat <- nestSbjModelFit |> mutate(pp=furrr::future_pmap(list(id,Model,Fit_Method,data), ~{
#    generate_data(..2, ..4 |> mutate(id=..1), ds, num_samples = 50, return_dat="test_data")
#    })) |> 
#   select(Fit_Method,pp,-data) |>  
#   unnest(pp) |>  filter(expMode2=="Test") |> as.data.table()
# 
# saveRDS(post_dat, here("data/model_cache/post_dat.rds"))

post_dat <- readRDS(here("data/model_cache/post_dat.rds"))

post_dat_avg <- post_dat |> group_by(id, condit, Model, Fit_Method, x, c, lr, rank) |> 
  mutate(error2 = y - pred) |>
  summarise(y = mean(y), pred = mean(pred), error = y - pred, error2=mean(error2)) |> as.data.table()

setorder(post_dat_avg, id, x, rank)
post_dat_l <- melt(post_dat_avg, id.vars = c("id", "condit", "Model", "Fit_Method", "x", "c", "lr", "rank","error"),
                   measure.vars = c("pred", "y"), variable.name = "Resp", value.name = "val")
post_dat_l[, Resp := fifelse(Resp == "y", "Observed",
                             fifelse(Model == "ALM", "ALM", "EXAM"))]
setorder(post_dat_l, id, Resp)
#rm(post_dat_avg)


# organize training data predictions
# pd_train <- nestSbjModelFit |> mutate(pp=furrr::future_pmap(list(id,Model,Fit_Method,data), ~{
#    generate_data(..2, ..4 |> mutate(id=..1), ds, num_samples = 20, return_dat="train_data")
#    })) |>
#   select(Fit_Method,pp,-data) |>
#   unnest(pp) |> as.data.table() |> filter(expMode2=="Train")

#saveRDS(pd_train, here("data/model_cache/pd_train.rds"))

pd_train <- readRDS(here("data/model_cache/pd_train.rds"))






nbins <- 3
pd_train <- pd_train |> group_by(id,condit,Model,Fit_Method) |>
  mutate(Block=cut(tr,breaks=seq(1,max(tr), length.out=nbins+1),include.lowest=TRUE,labels=FALSE))
setorder(pd_train, id, x,Block, rank)

pd_train_l <- melt(pd_train, id.vars = c("id", "condit", "Model","Block", "Fit_Method", "x", "c", "lr", "rank"),
                   measure.vars = c("pred", "y"), variable.name = "Resp", value.name = "val") |> as.data.table()
pd_train_l[, Resp := fifelse(Resp == "y", "Observed",
                             fifelse(Model == "ALM", "ALM", "EXAM"))] 
setorder(pd_train_l, id,Block, Resp) 

pd_train_l <- pd_train_l  |>
  mutate(dist = case_when(
    val >= x & val <= x + 200 ~ 0,                 
    val < x ~ abs(x - val),                       
    val > x + 200 ~ abs(val - (x + 200)),           
    TRUE ~ NA_real_                                 
  ))

#plan(sequential)
```

</details>

### Group level aggregations

<details class="code-fold">
<summary>Code</summary>

``` r
post_tabs <- abc_tables(post_dat,post_dat_l)

# post_tabs$agg_pred_full |> flextable::tabulator(rows=c("Fit_Method","Model"), columns=c("condit"), 
#                        `ME` = as_paragraph(mean_error)) |> as_flextable() 

post_tabs$agg_pred_full |> pander::pandoc.table()
```

</details>

<div class="cell-output cell-output-stdout">


    --------------------------------------------
      condit    Model   Fit_Method   mean_error 
    ---------- ------- ------------ ------------
     Constant    ALM       Test        276.7    

     Constant    ALM    Test_Train     288.2    

     Constant    ALM      Train        528.1    

     Constant   EXAM       Test        215.9    

     Constant   EXAM    Test_Train     228.6    

     Constant   EXAM      Train        340.3    

      Varied     ALM       Test        231.2    

      Varied     ALM    Test_Train     268.3    

      Varied     ALM      Train        368.7    

      Varied    EXAM       Test         215     

      Varied    EXAM    Test_Train     250.7    

      Varied    EXAM      Train        370.9    
    --------------------------------------------

</div>

The posterior distributions of the $c$ and $lr$ parameters are shown
<a href="#fig-htw-post-dist" class="quarto-xref">Figure 2</a>
(i.e. these plots combine all the posterior samples from all of the
subjects). There were substantial individual differences in the
posteriors of both parameters, with the within-group individual
differences generally swamped any between-group or between-model
differences. The magnitude of these individual differences remains even
if we consider only the single best parameter set for each subject.

We used the posterior distribution of $c$ and $lr$ parameters to
generate a posterior predictive distribution of the observed data for
each participant, which then allows us to compare the empirical data to
the full range of predictions from each model. Model residuals are shown
in the upper panels of
<a href="#fig-htw-resid-pred" class="quarto-xref">Figure 1</a>. The
pattern of training stage residual errors are unsurprising across the
combinations of models and fitting method . Differences between ALM and
EXAM are generally minor (the two models have identical learning
mechanisms). The differences in the magnitude of residuals across the
three fitting methods are also straightforward, with massive errors for
the ‘fit to Test Only’ model, and the smallest errors for the ‘fit to
train only’ models. It is also noteworthy that the residual errors are
generally larger for the first block of training, which is likely due to
the initial values of the ALM weights being unconstrained by whatever
initial biases participants tend to bring to the task. Future work may
explore the ability of the models to capture more fine grained aspects
of the learning trajectories. However for the present purposes, our
primary interest is in the ability of ALM and EXAM to account for the
testing patterns while being constrained, or not constrained, by the
training data. All subsequent analyses and discussion will thus focus on
the testing stage.

The residuals of the model predictions for the testing stage
(<a href="#fig-htw-resid-pred" class="quarto-xref">Figure 1</a>) also
show a sensible pattern across fitting methods - with models fit only to
the test data showing the best performance, followed by models fit to
both training and test data, and with models fit only to the training
data showing the worst performance (note that y axes are scaled
different between plots). Unsurprisingly, the advantage of EXAM is
strongest for extrapolation positions (the three smallest bands for both
groups - as well as the two highest bands for the Constant group).
Although EXAM tends to perform better for both Constant and Varied
participants (see also
<a href="#tbl-htw-modelError" class="quarto-xref">Table 2</a>), the
relative advantage of EXAM is generally larger for the Constant group -
a pattern consistent across all three fitting methods.

Panel B of
<a href="#fig-htw-resid-pred" class="quarto-xref">Figure 1</a> directly
compares the aggregated observed data to the posterior predictive
distributions for the testing stage. Of interest are a) the extent to
which the median estimates of the ALM and EXAM posteriors deviate from
the observed medians for each velocity band; b) the ability of ALM and
EXAM to discriminate between velocity bands; c) the relative performance
of models that are constrained by the training data (i.e. the ‘fit to
train only’ and ‘fit to both’ models) compared to the ‘fit to test only’
models; and d) the extent to which the variance of the posterior
predictive distributions mimics the variance of the observed data.

Considering first the models fit to only the testing data, which reflect
the best possible performance of ALM and EXAM at capturing the
group-aggregated testing patterns. For the varied group, both ALM and
EXAM are able to capture the median values of the observed data within
the 66% credible intervals, and the spread of model predictions
generally matches that of the observed data. For the constant group,
only EXAM is able to capture the median range of values across the
velocity bands, with ALM generally underestimating human velocoties in
the upper bands, and overestimating in the lower bands. In the case of
band 100, the median ALM prediction appears to match that of our
participants - however this is due to a large subset of participants
have ALM predictions near 0 for band 100, a pattern we will explore
further in our considertation of individual patterns below. Models fit
to both training and testing data show a similar pattern to only the
testing data display the same basic pattern as those fit to only the
testing data, albeit with slightly larger residuals. However models fit
to only the training data display markedly worse performance at
accounting for the key testing patterns.

- \*\* explain how the constant group ALM predictions for band 100 look
  deceptively good due to aggregation of a large subset of subjects
  having ALM predictions of 0 for vb100, and a large subset with ALM
  predictions close to their position 800 value. This is relected by
  much greater variance of the ALM esimates in the posterior predictive
  plot

- \*\* comment on how much constrained by the training data has a worse
  impact on the EXAM predictions for varied than for constant - perhaps
  due to the varied training data being much noisier than the constant
  training data.

- \*\* comment on EXAM doing a better job mimicing the within-condition
  variance of the observed data

- \*\* comment on the % of Constant subjects being best accounted for by
  EXAM being higher.

- \*\* does EXAM do better for the Constant group because the constant
  group performs better? Or does training with a single example
  encourage an exam sort of strategy?

<details class="code-fold">
<summary>Code</summary>

``` r
##| layout: [[45,-5, 45], [100]]
##| fig-subcap: ["Model Residuals - training data", "Model Residuals - testing data","Full posterior predictive distributions vs. observed data from participants."]
train_resid <- pd_train |> group_by(id,condit,Model,Fit_Method, Block) |> 
  summarise(y = mean(y), pred = mean(pred), error = y - pred) |>
  ggplot(aes(x = Block, y = abs(error), fill=Model)) + 
  stat_bar + 
  ggh4x::facet_nested_wrap(rename_fm(Fit_Method)~condit, scales="free",ncol=2) +
  scale_fill_manual(values=wes_palette("AsteroidCity2"))+
  labs(title="Model Residual Errors - Training Stage", y="RMSE", x= "Training Block") +
  theme(legend.title = element_blank(), legend.position="top")

test_resid <-  post_dat |> 
   group_by(id,condit,x,Model,Fit_Method,rank) |>
   summarize(error=mean(abs(y-pred)),n=n()) |>
   group_by(id,condit,x,Model,Fit_Method) |>
   summarize(error=mean(error)) |>
  mutate(vbLab = factor(paste0(x,"-",x+200))) |>
  ggplot(aes(x = vbLab, y = abs(error), fill=Model)) + 
  stat_bar + 
  scale_fill_manual(values=wes_palette("AsteroidCity2"))+
  ggh4x::facet_nested_wrap(rename_fm(Fit_Method)~condit, axes = "all",ncol=2,scale="free") +
  labs(title="Model Residual Errors - Testing Stage",y="RMSE", x="Velocity Band") +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5)) 

group_pred <- post_dat_l |> 
  mutate(vbLab = factor(paste0(x,"-",x+200),levels=levels(testAvg$vb))) |>
  ggplot(aes(x=val,y=vbLab,col=Resp)) + 
  stat_pointinterval(position=position_dodge(.5), alpha=.9) + 
  scale_color_manual(values=wes_palette("AsteroidCity2"))+
  ggh4x::facet_nested_wrap(rename_fm(Fit_Method)~condit, axes = "all",ncol=2,scale="free") +
  labs(title="Posterior Predictions - Testing Stage",y="Velocity Band (lower bound)", x="X Velocity") +
theme(legend.title=element_blank(),axis.text.y = element_text(angle = 45, hjust = 0.5, vjust = 0.5))


((train_resid | test_resid) / group_pred) +
  plot_layout(heights=c(1,1.5)) & 
  plot_annotation(tag_levels = list(c('A1','A2','B')),tag_suffix = ') ') & 
  theme(plot.tag.position = c(0, 1))
```

</details>
![](htw_model_files/figure-commonmark/fig-htw-resid-pred-1.jpeg)

<details class="code-fold">
<summary>Code</summary>

``` r
c_post <- post_dat_avg %>%
    group_by(id, condit, Model, Fit_Method, rank) %>%
    slice_head(n = 1) |>
    ggplot(aes(y=log(c), x = Fit_Method,col=condit)) + stat_pointinterval(position=position_dodge(.2)) +
    ggh4x::facet_nested_wrap(~Model) + labs(title="c parameter") +
  theme(legend.title = element_blank(), legend.position="right",plot.title=element_text(hjust=.4))

lr_post <- post_dat_avg %>%
    group_by(id, condit, Model, Fit_Method, rank) %>%
    slice_head(n = 1) |>
    ggplot(aes(y=lr, x = Fit_Method,col=condit)) + stat_pointinterval(position=position_dodge(.2)) +
    ggh4x::facet_nested_wrap(~Model) + labs(title="learning rate parameter") +
  theme(legend.title = element_blank(), legend.position = "none",plot.title=element_text(hjust=.4))
c_post + lr_post
```

</details>
![](htw_model_files/figure-commonmark/fig-htw-post-dist-1.jpeg)

### Accounting for individual patterns

To more accurately assess the relative abilities of ALM and EXAM to
capture important empirical patterns - we will now examine the
predictions of both models for the subset of individual participants
shown in <a href="#fig-htw-indv-pred" class="quarto-xref">Figure 3</a>.
Panel A presents three varied and constant participants who demonstrated
a reasonable degree of discrimination between the 6 velocity bands
during testing.

- \*\* comment on the different ways ALM can completely fail to mimic
  discrimination patterns (sbj. 35; sbj. 137),and on how it can
  sometimes partially succeed (sbj. 11; 14,74)

- \*\* comment on how EXAM can somtimes mimic non-monotonic spacing
  between bands due to associative stregth from training (i.e. subject
  47)

- \*\* compare c values to slope parameters from the statistical models
  earlier in paper

<details class="code-fold">
<summary>Code</summary>

``` r
cId_tr <- c(137, 181, 11)
vId_tr <- c(14, 193, 47)
cId_tt <- c(11, 93, 35)
vId_tt <- c(1,14,74)
# filter(id %in% (filter(bestTestEXAM,group_rank<=9, Fit_Method=="Test")

testIndv <- post_dat_l |> filter(id %in% c(cId_tt,vId_tt), Fit_Method=="Test_Train") |> 
   mutate(x=as.factor(x), Resp=as.factor(Resp)) |>
  group_by(id,condit,Fit_Method,Model,Resp) |>
   mutate(flab=paste0("Subject: ",id)) |>
  ggplot(aes(x = Resp, y = val, fill=x)) + 
  stat_bar_sd + ggh4x::facet_nested_wrap(condit~flab, axes = "all",ncol=3) +
  labs(title="Individual Participant fits from Test & Train Fitting Method",
       y="X Velocity",fill="Target Velocity") +
   guides(fill = guide_legend(nrow = 1)) + 
  theme(legend.position = "bottom",axis.title.x = element_blank())


trainIndv <- post_dat_l |> filter(id %in% c(cId_tr,vId_tr), Fit_Method=="Train") |> 
   mutate(x=as.factor(x), Resp=as.factor(Resp), flab=paste0("Subject: ",id)) |>
  group_by(id,condit,Fit_Method,Model,Resp) |>
  ggplot(aes(x = Resp, y = val, fill=x)) + 
  stat_bar + 
  ggh4x::facet_nested_wrap(condit~flab, axes = "all",ncol=3) +
  labs(title="Individual Participant fits from Train Only Fitting Method", y="X Velocity",
       fill="Target Velocity") +
     guides(fill = guide_legend(nrow = 1)) + 
  theme(legend.position = "bottom",axis.title.x = element_blank())


(testIndv  / trainIndv) +
  plot_annotation(tag_levels = list(c('A','B')),tag_suffix = ') ') & 
  theme(plot.tag.position = c(0, 1))
```

</details>
![](htw_model_files/figure-commonmark/fig-htw-indv-pred-1.jpeg)

<details class="code-fold">
<summary>Code</summary>

``` r
# 
# 
# full_comparison |> filter(Fit_Method=="Test_Train") |>
#   ungroup() |>
#   mutate(id = reorder(id, Error_difference)) %>%
#   ggplot(aes(y=id,x=Error_difference,fill=Best_Model))+
#   geom_col()+
#   ggh4x::facet_grid2(~condit,axes="all",scales="free_y", independent = "y")
# 

# d <- testAvg |> left_join(full_comparison, by=c("id","condit")) |> filter(Fit_Method=="Test_Train")
# 
# 
# 
# d |> ggplot(aes(x=vb,y=vx,fill=condit)) + stat_bar + facet_wrap(Fit_Method~Best_Model2,ncol=2)
# 
# d |> 
#   group_by(condit,Fit_Method,Best_Model) |> 
#   mutate(nGrp2=n()) |>
#   filter(abs(Error_difference)>15) |> 
#   ggplot(aes(x=vb,y=dist,fill=condit)) + 
#   stat_bar + facet_wrap(Fit_Method~Best_Model,ncol=2)
# 
# d |> group_by(condit,Fit_Method,Best_Model) %>% tally() |> mutate(n=n/6)
# 
# d |> group_by(condit,Fit_Method,Best_Model) |> filter(abs(Error_difference)>15) |> tally() |> mutate(n=n/6)
# 
# d |> group_by(condit,Fit_Method) |> mutate(m=mean(Error_difference), 
#                                               sd=sd(Error_difference), 
#                                               n=n()/6,se=sd/sqrt(n)) |>
#    group_by(condit,Fit_Method,Best_Model) |> 
#  # filter(abs(Error_difference)>(2.5*se)) |> 
#   ggplot(aes(x=vb,y=dist,fill=condit)) + 
#   stat_bar + facet_wrap(Fit_Method~Best_Model,ncol=2)
#   
```

</details>

### To add to appendix

<details class="code-fold">
<summary>Code</summary>

``` r
# post_tabs$agg_x_full |> flextable::tabulator(rows=c("Fit_Method","x"), columns=c("condit","Model"), 
#                        `X` = as_paragraph(mean_error)) |> as_flextable() 


post_dat  |> group_by(id,condit,Model,Fit_Method,x) |> 
  mutate(e2=abs(y-pred)) |> 
  summarise(y1=mean(y), pred1=mean(pred)) |>
  group_by(condit,Model,Fit_Method,x) |> 
  summarise(y=mean(y1), pred=mean(pred1),mean_error=abs(y-pred)) |> 
  round_tibble(1) |> pander::pandoc.table()
```

</details>


    -------------------------------------------------------------------
      condit    Model   Fit_Method    x       y     pred    mean_error 
    ---------- ------- ------------ ------ ------- ------- ------------
     Constant    ALM       Test      100    522.8   447.1      75.7    

     Constant    ALM       Test      350    660.6   797.6     136.9    

     Constant    ALM       Test      600    771.7    944      172.3    

     Constant    ALM       Test      800    1009    944.3      64.4    

     Constant    ALM       Test      1000   1173     944      229.4    

     Constant    ALM       Test      1200   1307    962.7     344.3    

     Constant    ALM    Test_Train   100    522.8   454.5      68.4    

     Constant    ALM    Test_Train   350    660.6   774.9     114.3    

     Constant    ALM    Test_Train   600    771.7   909.1     137.4    

     Constant    ALM    Test_Train   800    1009    909.5      99.3    

     Constant    ALM    Test_Train   1000   1173    909.1     264.3    

     Constant    ALM    Test_Train   1200   1307    926.7     380.3    

     Constant    ALM      Train      100    522.8   397.3     125.5    

     Constant    ALM      Train      350    660.6   496.7     163.9    

     Constant    ALM      Train      600    771.7   573.1     198.6    

     Constant    ALM      Train      800    1009    895.9     112.9    

     Constant    ALM      Train      1000   1173    573.1     600.3    

     Constant    ALM      Train      1200   1307    515.3     791.7    

     Constant   EXAM       Test      100    522.8   439.9      82.9    

     Constant   EXAM       Test      350    660.6   630.6      30.1    

     Constant   EXAM       Test      600    771.7   821.2      49.5    

     Constant   EXAM       Test      800    1009    973.7       35     

     Constant   EXAM       Test      1000   1173    1126       47.2    

     Constant   EXAM       Test      1200   1307    1279       28.2    

     Constant   EXAM    Test_Train   100    522.8   406.9     115.9    

     Constant   EXAM    Test_Train   350    660.6   596.1      64.5    

     Constant   EXAM    Test_Train   600    771.7   785.3      13.6    

     Constant   EXAM    Test_Train   800    1009    936.7      72.1    

     Constant   EXAM    Test_Train   1000   1173    1088       85.4    

     Constant   EXAM    Test_Train   1200   1307    1239       67.6    

     Constant   EXAM      Train      100    522.8   426.2      96.6    

     Constant   EXAM      Train      350    660.6   593.8      66.8    

     Constant   EXAM      Train      600    771.7   761.4      10.3    

     Constant   EXAM      Train      800    1009    895.5     113.3    

     Constant   EXAM      Train      1000   1173    1030      143.8    

     Constant   EXAM      Train      1200   1307    1164      143.3    

      Varied     ALM       Test      100    665.5   676.8      11.3    

      Varied     ALM       Test      350    772.5   832.3      59.9    

      Varied     ALM       Test      600    879.5   902.8      23.2    

      Varied     ALM       Test      800    1070     981       89.2    

      Varied     ALM       Test      1000   1181    1133       47.6    

      Varied     ALM       Test      1200   1269    1202        67     

      Varied     ALM    Test_Train   100    665.5   634.1      31.4    

      Varied     ALM    Test_Train   350    772.5   844.3      71.8    

      Varied     ALM    Test_Train   600    879.5   939.9      60.4    

      Varied     ALM    Test_Train   800    1070    964.7     105.6    

      Varied     ALM    Test_Train   1000   1181    1035      145.8    

      Varied     ALM    Test_Train   1200   1269    1082       187     

      Varied     ALM      Train      100    665.5   842.8     177.3    

      Varied     ALM      Train      350    772.5   956.3     183.8    

      Varied     ALM      Train      600    879.5   968.1      88.6    

      Varied     ALM      Train      800    1070    986.5      83.8    

      Varied     ALM      Train      1000   1181    1019      161.4    

      Varied     ALM      Train      1200   1269    1049      220.4    

      Varied    EXAM       Test      100    665.5    642       23.5    

      Varied    EXAM       Test      350    772.5   777.5      5.1     

      Varied    EXAM       Test      600    879.5    913       33.5    

      Varied    EXAM       Test      800    1070    1022       48.8    

      Varied    EXAM       Test      1000   1181    1130       50.8    

      Varied    EXAM       Test      1200   1269    1222       47.7    

      Varied    EXAM    Test_Train   100    665.5   684.1      18.6    

      Varied    EXAM    Test_Train   350    772.5   785.2      12.8    

      Varied    EXAM    Test_Train   600    879.5   886.4      6.8     

      Varied    EXAM    Test_Train   800    1070    967.3     102.9    

      Varied    EXAM    Test_Train   1000   1181    1048      132.5    

      Varied    EXAM    Test_Train   1200   1269    1075      194.3    

      Varied    EXAM      Train      100    665.5   859.6     194.1    

      Varied    EXAM      Train      350    772.5   904.1     131.6    

      Varied    EXAM      Train      600    879.5   948.6      69.1    

      Varied    EXAM      Train      800    1070    984.2       86     

      Varied    EXAM      Train      1000   1181    1020      160.9    

      Varied    EXAM      Train      1200   1269    1051      218.6    
    -------------------------------------------------------------------

<details class="code-fold">
<summary>Code</summary>

``` r
post_dat  |> group_by(id,condit,Model,Fit_Method,x) |> 
  mutate(e2=abs(y-pred)) |> 
  summarise(y1=mean(y), pred1=mean(pred)) |>
  group_by(condit,Model,Fit_Method) |> 
  summarise(y=mean(y1), pred=mean(pred1),mean_error=abs(y-pred)) |> 
  round_tibble(1) |> pander::pandoc.table()
```

</details>


    ------------------------------------------------------------
      condit    Model   Fit_Method     y     pred    mean_error 
    ---------- ------- ------------ ------- ------- ------------
     Constant    ALM       Test      907.4   839.9      67.4    

     Constant    ALM    Test_Train   907.4    814       93.4    

     Constant    ALM      Train      907.4   575.2     332.1    

     Constant   EXAM       Test      907.4   878.4       29     

     Constant   EXAM    Test_Train   907.4   842.1      65.3    

     Constant   EXAM      Train      907.4   811.7      95.7    

      Varied     ALM       Test       973    954.7      18.2    

      Varied     ALM    Test_Train    973    916.7      56.3    

      Varied     ALM      Train       973    970.3      2.6     

      Varied    EXAM       Test       973    950.9       22     

      Varied    EXAM    Test_Train    973    907.7      65.2    

      Varied    EXAM      Train       973    961.2      11.8    
    ------------------------------------------------------------

<details class="code-fold">
<summary>Code</summary>

``` r
post_dat  |> group_by(id,condit,Model,Fit_Method,x) |> 
  mutate(e2=abs(y-pred),Fit_Method=rename_fm(Fit_Method)) |> 
  summarise(y1=mean(y), pred1=mean(pred)) |>
  group_by(condit,Model,Fit_Method) |> 
  summarise(y=mean(y1), pred=mean(pred1),mean_error=abs(y-pred)) |> 
  select(-y,-pred) |>
  arrange(condit,Fit_Method) |>
  round_tibble(1) |> pander::pandoc.table()
```

</details>


    ---------------------------------------------------------------
      condit    Model            Fit_Method             mean_error 
    ---------- ------- ------------------------------- ------------
     Constant    ALM          Fit to Test Data             67.4    

     Constant   EXAM          Fit to Test Data              29     

     Constant    ALM    Fit to Test and Training Data      93.4    

     Constant   EXAM    Fit to Test and Training Data      65.3    

     Constant    ALM        Fit to Training Data          332.1    

     Constant   EXAM        Fit to Training Data           95.7    

      Varied     ALM          Fit to Test Data             18.2    

      Varied    EXAM          Fit to Test Data              22     

      Varied     ALM    Fit to Test and Training Data      56.3    

      Varied    EXAM    Fit to Test and Training Data      65.2    

      Varied     ALM        Fit to Training Data           2.6     

      Varied    EXAM        Fit to Training Data           11.8    
    ---------------------------------------------------------------

<details class="code-fold">
<summary>Code</summary>

``` r
best_id_x <- post_dat_avg |> group_by(id,condit,Model,Fit_Method,x) |>
  mutate(e2=(y-pred)) |>
  summarise(y=mean(y), pred=mean(pred), mean_error=mean(e2),abs_me=abs(mean_error)) |>
  group_by(id,condit,Fit_Method,x) |> mutate(best=ifelse(mean_error==min(mean_error),1,0)) |>
  group_by(id,condit,Fit_Method,Model) |> mutate(n_best=sum(best)) 

best_id <- best_id_x |> group_by(id,condit,Fit_Method,Model) |>
  summarise(mean_error=mean(mean_error), n_best=dplyr::first(n_best),abs_me=mean(abs_me))
  
lowest_error_model <- best_id %>%
  group_by(id, condit,Fit_Method) %>%
  summarise(Best_Model = Model[which.min(mean_error)],
            Best_Model2 = Model[which.min(abs_me)],
            n_best = n_best[which.min(mean_error)],
            Lowest_error = min(mean_error),
            differential = min(mean_error) - max(mean_error)) %>%
  ungroup()

error_difference <- best_id %>%
  select(id, condit, Model,Fit_Method, mean_error,abs_me) %>%
  pivot_wider(names_from = Model, values_from = c(abs_me,mean_error)) %>%
  mutate(Error_difference = (mean_error_ALM - mean_error_EXAM), abs_error_dif = (abs_me_ALM - abs_me_EXAM))

full_comparison <- lowest_error_model |> left_join(error_difference, by=c("id","condit","Fit_Method"))  |> 
  group_by(condit,Fit_Method,Best_Model) |> mutate(nGrp=n(), model_rank = nGrp - rank(Error_difference) ) |> 
  arrange(Fit_Method,-Error_difference)

full_comparison |> filter(Fit_Method=="Test_Train") |> 
  ungroup() |>
  mutate(id = reorder(id, abs_error_dif)) %>%
  ggplot(aes(y=id,x=abs_error_dif,fill=Best_Model2))+
  geom_col()+
  ggh4x::facet_grid2(~condit,axes="all",scales="free_y", independent = "y")+
  labs(fill="Best Model",x="Mean Model Error Difference (ALM - EXAM)",y="Participant")
```

</details>
![](htw_model_files/figure-commonmark/fig-htw-best-model-1.jpeg)

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0" line-spacing="2">

<div id="ref-bengtssonUnifyingFrameworkParallel2021" class="csl-entry">

Bengtsson, H. (2021). A Unifying Framework for Parallel and Distributed
Processing in R using Futures. *The R Journal*, *13*(2), 208.
<https://doi.org/10.32614/RJ-2021-048>

</div>

<div id="ref-brownUnderestimationLinearFunction2017" class="csl-entry">

Brown, M. A., & Lacroix, G. (2017). Underestimation in linear function
learning: Anchoring to zero or x-y similarity? *Canadian Journal of
Experimental Psychology/Revue Canadienne de Psychologie Expérimentale*,
*71*(4), 274–282. <https://doi.org/10.1037/cep0000129>

</div>

<div id="ref-cranmerFrontierSimulationbasedInference2020"
class="csl-entry">

Cranmer, K., Brehmer, J., & Louppe, G. (2020). The frontier of
simulation-based inference. *Proceedings of the National Academy of
Sciences*, *117*(48), 30055–30062.
<https://doi.org/10.1073/pnas.1912789117>

</div>

<div id="ref-deloshExtrapolationSineQua1997" class="csl-entry">

DeLosh, E. L., McDaniel, M. A., & Busemeyer, J. R. (1997).
Extrapolation: The Sine Qua Non for Abstraction in Function Learning.
*Journal of Experimental Psychology: Learning, Memory, and Cognition*,
*23*(4), 19. <https://doi.org/10.1037/0278-7393.23.4.968>

</div>

<div id="ref-farrellComputationalModelingCognition2018"
class="csl-entry">

Farrell, S., & Lewandowsky, S. (2018). *Computational Modeling of
Cognition and Behavior:* (1st ed.). Cambridge University Press.
<https://doi.org/10.1017/CBO9781316272503>

</div>

<div id="ref-kangasraasioParameterInferenceComputational2019"
class="csl-entry">

Kangasrääsiö, A., Jokinen, J. P. P., Oulasvirta, A., Howes, A., & Kaski,
S. (2019). Parameter Inference for Computational Cognitive Models with
Approximate Bayesian Computation. *Cognitive Science*, *43*(6), e12738.
<https://doi.org/10.1111/cogs.12738>

</div>

<div id="ref-kruschkeALCOVEExemplarbasedConnectionist1992"
class="csl-entry">

Kruschke, J. K. (1992). ALCOVE: An exemplar-based connectionist model of
Category Learning. *Psychological Review*, *99*(1).
<https://doi.org/10.1037/0033-295X.99.1.22>

</div>

<div id="ref-mcdanielPredictingTransferPerformance2009"
class="csl-entry">

Mcdaniel, M. A., Dimperio, E., Griego, J. A., & Busemeyer, J. R. (2009).
Predicting transfer performance: A comparison of competing function
learning models. *Journal of Experimental Psychology. Learning, Memory,
and Cognition*, *35*, 173–195. <https://doi.org/10.1037/a0013982>

</div>

<div id="ref-pageConnectionistModellingPsychology2000a"
class="csl-entry">

Page, M. (2000). Connectionist modelling in psychology: A localist
manifesto. *Behavioral and Brain Sciences*, *23*(4), 443–467.
<https://doi.org/10.1017/S0140525X00003356>

</div>

<div id="ref-turnerBayesianAnalysisSimulationbased2016"
class="csl-entry">

Turner, B. M., Sederberg, P. B., & McClelland, J. L. (2016). Bayesian
analysis of simulation-based models. *Journal of Mathematical
Psychology*, *72*, 191–199. <https://doi.org/10.1016/j.jmp.2014.10.001>

</div>

<div id="ref-turnerTutorialApproximateBayesian2012" class="csl-entry">

Turner, B. M., & Van Zandt, T. (2012). A tutorial on approximate
Bayesian computation. *Journal of Mathematical Psychology*, *56*(2),
69–85. <https://doi.org/10.1016/j.jmp.2012.02.005>

</div>

</div>
