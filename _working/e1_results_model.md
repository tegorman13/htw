
# HTW Project - HTW E1 Testing

We compared varied and constant performance across two measures, deviation and discrimination. Deviation was quantified as the absolute deviation from the nearest boundary of the velocity band, or set to 0 if the throw velocity fell anywhere inside the target band. Thus, when the target band was 600-800, throws of 400, 650, and 1100 would result in deviation values of 200, 0, and 300, respectively. Discrimination was measured by fitting a linear model to the testing throws of each subjects, with the lower end of the target velocity band as the predicted variable, and the x velocity produced by the participants as the predictor variable. Participants who reliably discriminated between velocity bands tended to have positive slopes with values ~1, while participants who made throws irrespective of the current target band would have slopes ~0.


### Testing Phase - No feedback.

In the first part of the testing phase, participants are tested from each of the velocity bands, and receive no feedback after each throw.

#### Deviation From Target Band

Descriptive summaries testing deviation data are provided in [Table 2](http://localhost:4200/Analysis/e1_test.html#tbl-e1-test-nf-deviation) and [Figure 1](http://localhost:4200/Analysis/e1_test.html#fig-e1-test-dev). To model differences in accuracy between groups, we used Bayesian mixed effects regression models to the trial level data from the testing phase. The primary model predicted the absolute deviation from the target velocity band (dist) as a function of training condition (condit), target velocity band (band), and their interaction, with random intercepts and slopes for each participant (id).

\\\[\\begin{equation} dist\_{ij} = \\beta\_0 + \\beta\_1 \\cdot condit\_{ij} + \\beta\_2 \\cdot band\_{ij} + \\beta\_3 \\cdot condit\_{ij} \\cdot band\_{ij} + b\_{0i} + b\_{1i} \\cdot band\_{ij} + \\epsilon\_{ij} \\end{equation}\\\] Code

````
```{r}
#| label: tbl-e1-test-nf-deviation
#| tbl-cap: "Testing Deviation - Empirical Summary"
#| tbl-subcap: ["Constant Testing - Deviation", "Varied Testing - Deviation"]
#| layout-ncol: 2

result <- test_summary_table(test, "dist","Deviation", mfun = list(mean = mean, median = median, sd = sd))
result$constant
result$varied
```
````

Table 2: Testing Deviation - Empirical Summary

Summary of Deviation - Constant
| Band | Band Type | Mean | Median | Sd |
| --- | --- | --- | --- | --- |
| 100-300 | Extrapolation | 254 | 148 | 298 |
| 350-550 | Extrapolation | 191 | 110 | 229 |
| 600-800 | Extrapolation | 150 | 84 | 184 |
| 800-1000 | Trained | 184 | 106 | 242 |
| 1000-1200 | Extrapolation | 233 | 157 | 282 |
| 1200-1400 | Extrapolation | 287 | 214 | 290 |

Summary of Deviation - Varied
| Band | Band Type | Mean | Median | Sd |
| --- | --- | --- | --- | --- |
| 100-300 | Extrapolation | 386 | 233 | 426 |
| 350-550 | Extrapolation | 285 | 149 | 340 |
| 600-800 | Extrapolation | 234 | 144 | 270 |
| 800-1000 | Trained | 221 | 149 | 248 |
| 1000-1200 | Trained | 208 | 142 | 226 |
| 1200-1400 | Trained | 242 | 182 | 235 |



Table 3: Experiment 1. Bayesian Mixed Model predicting absolute deviation as a function of condition (Constant vs. Varied) and Velocity Band

Model Coefficients
| Term | Estimate | 95% CrI Lower | 95% CrI Upper | pd |
| --- | --- | --- | --- | --- |
| Intercept | 205.09 | 136.86 | 274.06 | 1.00 |
| conditVaried | 157.44 | 60.53 | 254.90 | 1.00 |
| Band | 0.01 | \-0.07 | 0.08 | 0.57 |
| condit\*Band | \-0.16 | \-0.26 | \-0.06 | 1.00 |



```r
e1_distBMM |> 
  emmeans("condit",by="bandInt",at=list(bandInt=c(100,350,600,800,1000,1200)),
          epred = TRUE, re_formula = NA) |> 
  pairs() |> gather_emmeans_draws()  |> 
   summarize(median_qi(.value),pd=sum(.value>0)/n()) |>
   select(contrast,Band=bandInt,value=y,lower=ymin,upper=ymax,pd) |> 
   mutate(across(where(is.numeric), \(x) round(x, 2)),
          pd=ifelse(value<0,1-pd,pd)) |>
   kbl(caption="Contrasts")
```

Contrasts
| contrast | Band | value | lower | upper | pd |
| --- | --- | --- | --- | --- | --- |
| Constant - Varied | 100 | \-141.49 | \-229.19 | \-53.83 | 1.00 |
| Constant - Varied | 350 | \-101.79 | \-165.62 | \-36.32 | 1.00 |
| Constant - Varied | 600 | \-62.02 | \-106.21 | \-14.77 | 1.00 |
| Constant - Varied | 800 | \-30.11 | \-65.08 | 6.98 | 0.94 |
| Constant - Varied | 1000 | 2.05 | \-33.46 | 38.41 | 0.54 |
| Constant - Varied | 1200 | 33.96 | \-11.94 | 81.01 | 0.92 |


The model predicting absolute deviation (dist) showed clear effects of both training condition and target velocity band (Table X). Overall, the varied training group showed a larger deviation relative to the constant training group (β = 157.44, 95% CI \[60.53, 254.9\]). Deviation also depended on target velocity band, with lower bands showing less deviation. See [Table 3](http://localhost:4200/Analysis/e1_test.html#tbl-e1-bmm-dist) for full model output.


### Discrimination between bands

In addition to accuracy/deviation, we also assessed the ability of participants to reliably discriminate between the velocity bands (i.e. responding differently when prompted for band 600-800 than when prompted for band 150-350). [Table 4](http://localhost:4200/Analysis/e1_test.html#tbl-e1-test-nf-vx) shows descriptive statistics of this measure, and Figure 1 visualizes the full distributions of throws for each combination of condition and velocity band. To quantify discrimination, we again fit Bayesian Mixed Models as above, but this time the dependent variable was the raw x velocity generated by participants on each testing trial.

\\\[\\begin{equation} vx\_{ij} = \\beta\_0 + \\beta\_1 \\cdot condit\_{ij} + \\beta\_2 \\cdot bandInt\_{ij} + \\beta\_3 \\cdot condit\_{ij} \\cdot bandInt\_{ij} + b\_{0i} + b\_{1i} \\cdot bandInt\_{ij} + \\epsilon\_{ij} \\end{equation}\\\]


````
```{r}
#| label: tbl-e1-test-nf-vx
#| tbl-cap: "Testing vx - Empirical Summary"
#| tbl-subcap: ["Constant Testing - vx", "Varied Testing - vx"]
#| layout-ncol: 2

result <- test_summary_table(test, "vx","X Velocity", mfun = list(mean = mean, median = median, sd = sd))
result$constant 
result$varied 
```
````

Table 4: Testing vx - Empirical Summary

Summary of X Velocity - Constant
| Band | Band Type | Mean | Median | Sd |
| --- | --- | --- | --- | --- |
| 100-300 | Extrapolation | 524 | 448 | 327 |
| 350-550 | Extrapolation | 659 | 624 | 303 |
| 600-800 | Extrapolation | 770 | 724 | 300 |
| 800-1000 | Trained | 1001 | 940 | 357 |
| 1000-1200 | Extrapolation | 1167 | 1104 | 430 |
| 1200-1400 | Extrapolation | 1283 | 1225 | 483 |

Summary of X Velocity - Varied
| Band | Band Type | Mean | Median | Sd |
| --- | --- | --- | --- | --- |
| 100-300 | Extrapolation | 664 | 533 | 448 |
| 350-550 | Extrapolation | 768 | 677 | 402 |
| 600-800 | Extrapolation | 876 | 813 | 390 |
| 800-1000 | Trained | 1064 | 1029 | 370 |
| 1000-1200 | Trained | 1180 | 1179 | 372 |
| 1200-1400 | Trained | 1265 | 1249 | 412 |



Table 5: Experiment 1. Bayesian Mixed Model Predicting Vx as a function of condition (Constant vs. Varied) and Velocity Band

Fit to all 6 bands
| Term | Estimate | 95% CrI Lower | 95% CrI Upper | pd |
| --- | --- | --- | --- | --- |
| Intercept | 408.55 | 327.00 | 490.61 | 1.00 |
| conditVaried | 164.05 | 45.50 | 278.85 | 1.00 |
| Band | 0.71 | 0.62 | 0.80 | 1.00 |
| condit\*Band | \-0.14 | \-0.26 | \-0.01 | 0.98 |

Fit to 3 extrapolation bands
| Term | Estimate | 95% CrI Lower | 95% CrI Upper | pd |
| --- | --- | --- | --- | --- |
| Intercept | 478.47 | 404.00 | 551.45 | 1.00 |
| conditVaried | 142.04 | 37.17 | 247.59 | 1.00 |
| Band | 0.50 | 0.42 | 0.57 | 1.00 |
| condit\*Band | \-0.07 | \-0.17 | 0.04 | 0.89 |

See [Table 5](http://localhost:4200/Analysis/e1_test.html#tbl-e1-bmm-vx) for the full model results. The estimated coefficient for training condition (β = 164.05, 95% CrI \[45.5, 278.85\]) suggests that the varied group tends to produce harder throws than the constant group, but is not in and of itself useful for assessing discrimination. Most relevant to the issue of discrimination is the slope on Velocity Band (β = 0.71, 95% CrI \[0.62, 0.8\]). Although the median slope does fall underneath the ideal of value of 1, the fact that the 95% credible interval does not contain 0 provides strong evidence that participants exhibited some discrimination between bands. The estimate for the interaction between slope and condition (β = -0.14, 95% CrI \[-0.26, -0.01\]), suggests that the discrimination was somewhat modulated by training condition, with the varied participants showing less senitivity between vands than the constant condition. This difference is depicted visually in [Figure 4](http://localhost:4200/Analysis/e1_test.html#fig-e1-bmm-vx).@tbl-e1-slope-quartile shows the average slope coefficients for varied and constant participants separately for each quartile. The constant participant participants appear to have larger slopes across quartiles, but the difference between conditions may be less pronounced for the top quartiles of subjects who show the strongest discrimination. Figure [Figure 5](http://localhost:4200/Analysis/e1_test.html#fig-e1-bmm-bx2) shows the distributions of slope values for each participant, and the compares the probability density of slope coefficients between training conditions. [Figure 6](http://localhost:4200/Analysis/e1_test.html#fig-e1-indv-slopes)

The second model, which focused solely on extrapolation bands, revealed similar patterns. The Velocity Band term (β = 0.5, 95% CrI \[0.42, 0.57\]) still demonstrates a high degree of discrimination ability. However, the posterior distribution for interaction term (β = -0.07, 95% CrI \[-0.17, 0.04\] ) does across over 0, suggesting that the evidence for decreased discrimination ability for the varied participants is not as strong when considering only the three extrapolation bands.


Table 6: Slope coefficients by quartile, per condition

| Condition | Q\_0%\_mean | Q\_25%\_mean | Q\_50%\_mean | Q\_75%\_mean | Q\_100%\_mean |
| --- | --- | --- | --- | --- | --- |
| Constant | \-0.1025710 | 0.4796476 | 0.6891841 | 0.9284214 | 1.394802 |
| Varied | \-0.2015313 | 0.2673139 | 0.5888388 | 0.8980842 | 1.295650 |

[Figure 5](http://localhost:4200/Analysis/e1_test.html#fig-e1-bmm-bx2) shows the distributions of estimated slopes relating velocity band to x velocity for each participant, ordered from lowest to highest within condition. Slope values are lower overall for varied training compared to constant training. Figure Xb plots the density of these slopes for each condition. The distribution for varied training has more mass at lower values than the constant training distribution. Both figures illustrate the model’s estimate that varied training resulted in less discrimination between velocity bands, evidenced by lower slopes on average.

# Computational Modeling

In project 1, we applied model-based techniques to quantify and control for the similarity between training and testing experience, which in turn enabled us to account for the difference between varied and constant training via an extended version of a similarity based generalization model. In project 2, we will go a step further, implementing a full process model capable of both 1) producing novel responses and 2) modeling behavior in both the learning and testing stages of the experiment. For this purpose, we will apply the associative learning model (ALM) and the EXAM model of function learning (DeLosh 1997). ALM is a simple connectionist learning model which closely resembles Kruschke's ALCOVE model (Kruscke 1992), with modifications to allow for the generation of continuous responses.

## ALM & Exam Description

DeLosh et al. (1997) introduced the associative learning model (ALM), a connectionist model within the popular class of radial-basis networks. ALM was inspired by, and closely resembles Kruschke's influential ALCOVE model of categorization (Kruschke, 1992).

ALM is a localist neural network model, with each input node corresponding to a particular stimulus, and each output node corresponding to a particular response value. The units in the input layer activate as a function of their Gaussian similarity to the input stimulus. So, for example, an input stimulus of value 55 would induce maximal activation of the input unit tuned to 55. Depending on thevalue of the generalization parameter, the nearby units (e.g. 54 and 56; 53 and 57) may also activate to some degree. ALM is structured with input and output nodes that correspond to regions of the stimulus space, and response space, respectively. The units in the input layer activate as a function of their similarity to a presented stimulus. As was the case with the exemplar-based models, similarity in ALM is exponentially decaying function of distance. The input layer is fully connected to the output layer, and the activation for any particular output node is simply the weighted sum of the connection weights between that node and the input activations. The network then produces a response by taking the weighted average of the output units (recall that each output unit has a value corresponding to a particular response). During training, the network receives feedback which activates each output unit as a function of its distance from the ideal level of activation necessary to produce the correct response. The connection weights between input and output units are then updated via the standard delta learning rule, where the magnitude of weight changes are controlled by a learning rate parameter.

See for a full specification of the equations that define ALM and EXAM.



## Model Table

### ALM Activation & Response

| Step                          | Equation                                                                            | Description                                                                                                                            |
|------------------|------------------------|------------------------------|
| **ALM Activation & Response** |                                                                                     |                                                                                                                                        |
| Input Activation              | $a_i(X) = \frac{e^{-c(X-X_i)^2}}{\sum_{k=1}^M e^{-c(X-X_k)^2}}$                     | Activation of each input node $X_i$, is a function of the Gaussian similarity between the node value and stimulus X.                   |
| Output Activation             | $O_j(X) = \sum_{k=1}^M w_{ji} \cdot a_i(X)$                                         | Activation of each Output unit $O_j$ is the weighted sum of the input activations and association weights.                             |
| Output Probability            | $P[Y_j|X] = \frac{O_j(X)}{\sum_{k=1}^M O_k(X)}$                                     | Each output node has associated response, $Y_j$. The probability of response $Y_j$ is determined by the ratio of output activations.   |
| Mean Output                   | $m(x) = \sum_{j=1}^L Y_j \cdot \frac{O_j(x)}{\sum_{k=1}^M O_k(X)}$                  | The response to stimulus x is the weighted average of the response probabilities.                                                      |
| **ALM Learning**              |                                                                                     |                                                                                                                                        |
| Feedback Activation           | $f_j(Z) = e^{-c(Z-Y_j)^2}$                                                          | After responding, feedback signal Z is presented, activating each output node via the Gaussian similarity to the ideal response.       |
| Update Weights                | $w_{ji}(t + 1) = w_{ji}(t) + \alpha \cdot (f_j(Z(t)) - O_j(X(t)) \cdot a_i(X(t))$   | Delta rule to update weights. Magnitude of weight changes controlled by learning rate parameter alpha.                                 |
| **EXAM**                      |                                                                                     |                                                                                                                                        |
| Extrapolation                 | $P[X_i|X] = \frac{a_i(X)}{\sum_{k=1}^M a_k(X)}$                                     | Novel test stimulus X activates input nodes associated with trained stimuli.                                                           |
|                               | $E[Y|X_i] = m(X_i) + \frac{m(X_{i+1})-m(X_{i-1})}{X_{i+1}-X_{i-1}} \cdot [X - X_i]$ | Slope value computed from nearest training instances and then added to the response associated with the nearest training instance,m(x) |


## Model Fitting and Comparison

Following the procedure used by Mcdaniel et al. (2009), we will assess the ability of both ALM and EXAM to account for the empirical data when fitting the models to 1) only the training data, and 2) both training and testing data. Models will be fit directly to the trial by trial data of each individual participants, both by minimizing the root-mean squared deviation (RMSE). Because ALM has been shown to do poorly at accounting for human patterns extrapolation (DeLosh et al., 1997), we will also fit the extended EXAM version of the model, which operates identically to ALM during training, but includes a linear extrapolation mechanism for generating novel responses during testing.

We also fit a 3rd, hybrid model, which generates responses as a weighted sum of ALM and EXAM predictions. For the hybrid model, predictions are computed by first generating separate predictions from ALM and EXAM, and then combining them using the following equation: $\hat{y} = (1 - w) \cdot alm_pred + w \cdot exam_pred$, where $w$ is a third fit parameters that sets the relative contribution between the two models. For the grid search, the weight parameter is varied from 0 to 1, and the resulting RMSE is recorded. 

Each model was fit to the data in 3 different ways. 1) To just the testing data, 2) Both the training and testing data, 3) Only the training data. In all cases, the model only updates its weights during the training phase, and the weights are frozen during the testing phase. In all cases, only the ALM model generates predictions during the training phase. For the testing phase, all 3 models are used to generate predictions. 


### Varied Parameter fits and RMSE

|      Model       |   c   |  lr  | Value | Test_RMSE |
|:----------------:|:-----:|:----:|:-----:|:---------:|
|  ALM Test Only   | 0.134 | 2.03 | 95.46 |   95.46   |
| ALM Test & Train | 0.067 | 0.1  | 247.3 |   106.5   |
|  ALM Train Only  | 0.047 | 0.08 | 139.2 |    109    |

ALM

|       Model       |   c   |  lr  | Value | Test_RMSE |
|:-----------------:|:-----:|:----:|:-----:|:---------:|
|  EXAM Test Only   | 0.409 | 1.91 | 45.84 |   45.84   |
| EXAM Test & Train | 0.074 | 0.1  | 201.4 |   60.18   |
|  EXAM Train Only  | 0.047 | 0.08 | 139.2 |   65.31   |

EXAM

|        Model        |   c   |  lr   |   w   | Value | Test_RMSE |
|:-------------------:|:-----:|:-----:|:-----:|:-----:|:---------:|
|  Hybrid Test Only   | 0.395 | 2.017 | 0.643 | 33.88 |   33.88   |
| Hybrid Test & Train | 0.134 | 2.017 | 0.786 | 197.2 |   46.51   |
|  Hybrid Train Only  | 0.042 | 0.067 |   0   | 139.2 |   110.3   |

Hybrid



```{r}
pander(almParamsC, caption="ALM"); pander(examParamsC, caption="EXAM"); pander(hybridParamsC,caption="Hybrid")
```

### Constant Parameter fits and RMSE

|      Model       |   c   |  lr  | Value | Test_RMSE |
|:----------------:|:-----:|:----:|:-----:|:---------:|
|  ALM Test Only   |   0   | 0.1  | 309.5 |   347.8   |
| ALM Test & Train | 0.047 | 0.08 |  361  |   328.5   |
|  ALM Train Only  | 0.06  | 0.1  | 32.44 |    329    |

ALM

|       Model       |   c   |  lr   | Value | Test_RMSE |
|:-----------------:|:-----:|:-----:|:-----:|:---------:|
|  EXAM Test Only   | 0.007 | 1.327 | 127.3 |   127.3   |
| EXAM Test & Train | 0.081 | 0.161 | 194.6 |    132    |
|  EXAM Train Only  | 0.06  |  0.1  | 32.44 |   199.8   |

EXAM

|        Model        |   c   |  lr   |  w  | Value | Test_RMSE |
|:-------------------:|:-----:|:-----:|:---:|:-----:|:---------:|
|  Hybrid Test Only   | 0.008 | 1.58  |  1  | 127.3 |   127.3   |
| Hybrid Test & Train | 0.067 | 0.134 |  1  | 194.5 |   136.4   |
|  Hybrid Train Only  | 0.042 | 0.067 |  0  | 31.5  |   330.3   |

Hybrid


## Varied Testing Predictions

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


### References 

DeLosh, E. L., McDaniel, M. A., & Busemeyer, J. R. (1997). Extrapolation: The Sine Qua Non for Abstraction in Function Learning. *Journal of Experimental Psychology: Learning, Memory, and Cognition*, *23*(4), 19. <https://doi.org/10.1037/0278-7393.23.4.968>

Kruschke, J. K. (1992). ALCOVE: An exemplar-based connectionist model of Category Learning. *Psychological Review*, *99*(1). <https://doi.org/10.1037/0033-295X.99.1.22>

Mcdaniel, M. A., Dimperio, E., Griego, J. A., & Busemeyer, J. R. (2009). Predicting transfer performance: A comparison of competing function learning models. *Journal of Experimental Psychology. Learning, Memory, and Cognition*, *35*, 173--195. <https://doi.org/10.1037/a0013982>


