---
title: "HTW"
short-title: "Variability and Extrapolation"
date: last-modified
author:
- name: Thomas Gorman
  affiliation: Indiana University
  url: www.tegorman13.github.io
  email: tegorman@iu.edu
  orcid: 0000-0001-5366-5442
- name: Rob Goldstone
  affiliation: Indiana University
  url: www.pc.cogs.indiana.edu 
  email: rgoldsto@indiana.edu
  orcid: 0000-0001-8357-8358
abstract: |
  In project 1, we applied model-based techniques to quantify and control for the similarity between training and testing experience, which in turn enabled us to account for the difference between varied and constant training via an extended version of a similarity based generalization model. In project 2, we will go a step further, implementing a full process model capable of both 1) producing novel responses and 2) modeling behavior in both the learning and testing stages of the experiment. Project 2 also places a greater emphasis on extrapolation performance following training - as varied training has often been purported to be particularly beneficial in such situations. 
keywords:
  - Learning Generalization
  - Function Learning
  - Visuomotor learning
  - Training Variability
code-repo: "Access the code, data, and analysis at <https://github.com/tegorman13/htw>"
bibliography: ../Assets/Bib/HTW.bib
link-citations: true
keep-md: true
---





# Introduction

In project 1, we applied model-based techniques to quantify and control for the similarity between training and testing experience, which in turn enabled us to account for the difference between varied and constant training via an extended version of a similarity based generalization model. In project 2, we will go a step further, implementing a full process model capable of both 1) producing novel responses and 2) modeling behavior in both the learning and testing stages of the experiment. Project 2 also places a greater emphasis on extrapolation performance following training - as varied training has often been purported to be particularly beneficial in such situations. Extrapolation has long been a focus of the literature on function learning [@brehmerHypothesesRelationsScaled1974; @carrollFunctionalLearningLearning1963]. Central questions of the function learning literature have included the relative difficulties of learning various functional forms (e.g. linear vs.bilinear vs. quadratic), and the relative effectiveness of rule-based vs. association-based exemplar models vs. various hybrid models [@bottNonmonotonicExtrapolationFunction2004; @deloshExtrapolationSineQua1997; @jonesActiveFunctionLearning2018; @kalishPopulationLinearExperts2004; @mcdanielConceptualBasisFunction2005; @mcdanielPredictingTransferPerformance2009]. However the issue of training variation has received surprisingly little attention in this area.


# Methods

## Participants

Data was collected from 647 participants (after exclusions). The results shown below consider data from subjects in our initial experiment, which consisted of 196 participants (106 constant, 90 varied). The follow-up experiments entailed minor manipulations: 1) reversing the velocity bands that were trained on vs. novel during testing; 2) providing ordinal rather than numerical feedback during training (e.g. correct, too low, too high). The data from these subsequent experiments are largely consistently with our initial results shown below.

## Task

We developed a novel visuomotor extrapolation task, termed the Hit The Wall task, wherein participants learned to launch a projectile such that it hit a rectangle at the far end of the screen with an appropriate amount of force. Although the projectile had both x and y velocity components, only the x-dimension was relevant for the task.  [Link to task demo](https://pcl.sitehost.iu.edu/tg/HTW/HTW_Index.html?sonaid=){target="_blank"}

## Procedure
Upon arrival at the laboratory, participants were provided with a description of the experiment and signed informed consent forms. They were then seated in front of a computer equipped with a mouse and were given instructions on how to perform the "Hit The Wall" (HTW) visuomotor extrapolation task.

The HTW task involved launching projectiles to hit a target displayed on the computer screen. Participants completed a total of 90 trials during the training stage. In the varied training condition, participants encountered three velocity bands (800-1000, 1000-1200, and 1200-1400). In contrast, participants in the constant training condition encountered only one velocity band (800-1000).

During the training stage, participants in both conditions also completed "no feedback" trials, where they received no information about their performance. These trials were randomly interleaved with the regular training trials.

Following the training stage, participants proceeded to the testing stage, which consisted of three phases. In the first phase, participants completed "no-feedback" testing from three novel extrapolation bands (100-300, 350-550, and 600-800), with each band consisting of 15 trials.

In the second phase of testing, participants completed "no-feedback" testing from the three velocity bands used during the training stage (800-1000, 1000-1200, and 1200-1400). In the constant training condition, two of these bands were novel, while in the varied training condition, all three bands were encountered during training.

The third and final phase of testing involved "feedback" testing for each of the three extrapolation bands (100-300, 350-550, and 600-800), with each band consisting of 10 trials. Participants received feedback on their performance during this phase.

Throughout the experiment, participants' performance was measured by calculating the distance between the produced x-velocity of the projectiles and the closest edge of the current velocity band. Lower distances indicated better performance.

After completing the experiment, participants were debriefed and provided with an opportunity to ask questions about the study.




:::{.cell .column-screen-inset-right fig-width="6" fig-height="2.5" fig-responsive=false}

:::{.cell-output-display}

:::{#fig-design-e1}

:::{}

![](manuscript_files/figure-odt/dot-figure-3.png){width="6in" height="2.5in"}
:::


Experiment 1 Design. Constant and Varied participants complete different training conditions.
:::
:::
:::






## Analyses Strategy

All data processing and statistical analyses were performed in R version 4.31 @rcoreteamLanguageEnvironmentStatistical2020. To assess differences between groups, we used Bayesian Mixed Effects Regression. Model fitting was performed with the brms package in R @burknerBrmsPackageBayesian2017, and descriptive stats and tables were extracted with the BayestestR package @makowskiBayestestRDescribingEffects2019a. Mixed effects regression enables us to take advantage of partial pooling, simultaneously estimating parameters at the individual and group level. Our use of Bayesian, rather than frequentist methods allows us to directly quantify the uncertainty in our parameter estimates, as well as circumventing convergence issues common to the frequentist analogues of our mixed models. For each model, we report the median values of the posterior distribution, and 95% credible intervals.

Each model was set to run with 4 chains, 5000 iterations per chain, with the first 2500 of which were discarded as warmup chains. Rhat values were generally within an acceptable range, with values \<=1.02 (see appendix for diagnostic plots). We used uninformative priors for the fixed effects of the model (condition and velocity band), and weakly informative Student T distributions for for the random effects.

We compared varied and constant performance across two measures, deviation and discrimination. Deviation was quantified as the absolute deviation from the nearest boundary of the velocity band, or set to 0 if the throw velocity fell anywhere inside the target band. Thus, when the target band was 600-800, throws of 400, 650, and 1100 would result in deviation values of 200, 0, and 300, respectively. Discrimination was measured by fitting a linear model to the testing throws of each subjects, with the lower end of the target velocity band as the predicted variable, and the x velocity produced by the participants as the predictor variable. Participants who reliably discriminated between velocity bands tended to have positive slopes with values \~1, while participants who made throws irrespective of the current target band would have slopes \~0.








::: {#tbl-e1-test-nf-deviation .cell layout-ncol="1" layout-align="center" tbl-cap='Testing Deviation - Empirical Summary' tbl-subcap='["Full datasets","Intersection of samples with all labels available"]'}
::: {.cell-output-display}

`````{=html}
<table>
<caption>Summary of Deviation- Constant</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> Band </th>
   <th style="text-align:left;"> Band Type </th>
   <th style="text-align:right;"> Mean </th>
   <th style="text-align:right;"> Median </th>
   <th style="text-align:right;"> Sd </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> 100-300 </td>
   <td style="text-align:left;"> Extrapolation </td>
   <td style="text-align:right;"> 254 </td>
   <td style="text-align:right;"> 148 </td>
   <td style="text-align:right;"> 298 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 350-550 </td>
   <td style="text-align:left;"> Extrapolation </td>
   <td style="text-align:right;"> 191 </td>
   <td style="text-align:right;"> 110 </td>
   <td style="text-align:right;"> 229 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 600-800 </td>
   <td style="text-align:left;"> Extrapolation </td>
   <td style="text-align:right;"> 150 </td>
   <td style="text-align:right;"> 84 </td>
   <td style="text-align:right;"> 184 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 800-1000 </td>
   <td style="text-align:left;"> Trained </td>
   <td style="text-align:right;"> 184 </td>
   <td style="text-align:right;"> 106 </td>
   <td style="text-align:right;"> 242 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 1000-1200 </td>
   <td style="text-align:left;"> Extrapolation </td>
   <td style="text-align:right;"> 233 </td>
   <td style="text-align:right;"> 157 </td>
   <td style="text-align:right;"> 282 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 1200-1400 </td>
   <td style="text-align:left;"> Extrapolation </td>
   <td style="text-align:right;"> 287 </td>
   <td style="text-align:right;"> 214 </td>
   <td style="text-align:right;"> 290 </td>
  </tr>
</tbody>
</table>

`````

:::

::: {.cell-output-display}

`````{=html}
<table>
<caption>Summary of Deviation- Varied</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> Band </th>
   <th style="text-align:left;"> Band Type </th>
   <th style="text-align:right;"> Mean </th>
   <th style="text-align:right;"> Median </th>
   <th style="text-align:right;"> Sd </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> 100-300 </td>
   <td style="text-align:left;"> Extrapolation </td>
   <td style="text-align:right;"> 386 </td>
   <td style="text-align:right;"> 233 </td>
   <td style="text-align:right;"> 426 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 350-550 </td>
   <td style="text-align:left;"> Extrapolation </td>
   <td style="text-align:right;"> 285 </td>
   <td style="text-align:right;"> 149 </td>
   <td style="text-align:right;"> 340 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 600-800 </td>
   <td style="text-align:left;"> Extrapolation </td>
   <td style="text-align:right;"> 234 </td>
   <td style="text-align:right;"> 144 </td>
   <td style="text-align:right;"> 270 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 800-1000 </td>
   <td style="text-align:left;"> Trained </td>
   <td style="text-align:right;"> 221 </td>
   <td style="text-align:right;"> 149 </td>
   <td style="text-align:right;"> 248 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 1000-1200 </td>
   <td style="text-align:left;"> Trained </td>
   <td style="text-align:right;"> 208 </td>
   <td style="text-align:right;"> 142 </td>
   <td style="text-align:right;"> 226 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 1200-1400 </td>
   <td style="text-align:left;"> Trained </td>
   <td style="text-align:right;"> 242 </td>
   <td style="text-align:right;"> 182 </td>
   <td style="text-align:right;"> 235 </td>
  </tr>
</tbody>
</table>

`````

:::
:::



## Results

### Testing Phase - No feedback.

In the first part of the testing phase, participants are tested from each of the velocity bands, and receive no feedback after each throw.

#### Deviation From Target Band

Descriptive summaries testing deviation data are provided in @tbl-e1-test-nf-deviation and @fig-e1-test-dev. To model differences in accuracy between groups, we used Bayesian mixed effects regression models to the trial level data from the testing phase. The primary model predicted the absolute deviation from the target velocity band (dist) as a function of training condition (condit), target velocity band (band), and their interaction, with random intercepts and slopes for each participant (id).




```{=tex}
\begin{equation}
dist_{ij} = \beta_0 + \beta_1 \cdot condit_{ij} + \beta_2 \cdot band_{ij} + \beta_3 \cdot condit_{ij} \cdot band_{ij} + b_{0i} + b_{1i} \cdot band_{ij} + \epsilon_{ij}
\end{equation}
```


::: {.cell layout-align="center"}
::: {.cell-output-display}
![E1. Deviations from target band during testing without feedback stage.](manuscript_files/figure-odt/fig-e1-test-dev-1.png){#fig-e1-test-dev fig-align='center' width=100%}
:::
:::

::: {#tbl-e1-bmm-dist .cell layout-align="center" tbl-cap='Experiment 1. Bayesian Mixed Model predicting absolute deviation as a function of condition (Constant vs. Varied) and Velocity Band' tbl-subcap='["Constant Testing1 - Deviation","Varied Testing - Deviation"]'}
::: {.cell-output-display}

`````{=html}
<table>
<caption>Coefficients</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> Term </th>
   <th style="text-align:right;"> Estimate </th>
   <th style="text-align:right;"> 95% CrI Lower </th>
   <th style="text-align:right;"> 95% CrI Upper </th>
   <th style="text-align:right;"> pd </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Intercept </td>
   <td style="text-align:right;"> 205.09 </td>
   <td style="text-align:right;"> 136.86 </td>
   <td style="text-align:right;"> 274.06 </td>
   <td style="text-align:right;"> 1.00 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> conditVaried </td>
   <td style="text-align:right;"> 157.44 </td>
   <td style="text-align:right;"> 60.53 </td>
   <td style="text-align:right;"> 254.90 </td>
   <td style="text-align:right;"> 1.00 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Band </td>
   <td style="text-align:right;"> 0.01 </td>
   <td style="text-align:right;"> -0.07 </td>
   <td style="text-align:right;"> 0.08 </td>
   <td style="text-align:right;"> 0.57 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> condit*Band </td>
   <td style="text-align:right;"> -0.16 </td>
   <td style="text-align:right;"> -0.26 </td>
   <td style="text-align:right;"> -0.06 </td>
   <td style="text-align:right;"> 1.00 </td>
  </tr>
</tbody>
</table>

`````

:::

::: {.cell-output-display}

`````{=html}
<table>
 <thead>
  <tr>
   <th style="text-align:left;"> contrast </th>
   <th style="text-align:right;"> Band </th>
   <th style="text-align:right;"> value </th>
   <th style="text-align:right;"> lower </th>
   <th style="text-align:right;"> upper </th>
   <th style="text-align:right;"> pd </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Constant - Varied </td>
   <td style="text-align:right;"> 100 </td>
   <td style="text-align:right;"> -141.49 </td>
   <td style="text-align:right;"> -229.2 </td>
   <td style="text-align:right;"> -53.83 </td>
   <td style="text-align:right;"> 1.00 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Constant - Varied </td>
   <td style="text-align:right;"> 350 </td>
   <td style="text-align:right;"> -101.79 </td>
   <td style="text-align:right;"> -165.6 </td>
   <td style="text-align:right;"> -36.32 </td>
   <td style="text-align:right;"> 1.00 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Constant - Varied </td>
   <td style="text-align:right;"> 600 </td>
   <td style="text-align:right;"> -62.02 </td>
   <td style="text-align:right;"> -106.2 </td>
   <td style="text-align:right;"> -14.77 </td>
   <td style="text-align:right;"> 1.00 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Constant - Varied </td>
   <td style="text-align:right;"> 800 </td>
   <td style="text-align:right;"> -30.11 </td>
   <td style="text-align:right;"> -65.1 </td>
   <td style="text-align:right;"> 6.98 </td>
   <td style="text-align:right;"> 0.94 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Constant - Varied </td>
   <td style="text-align:right;"> 1000 </td>
   <td style="text-align:right;"> 2.05 </td>
   <td style="text-align:right;"> -33.5 </td>
   <td style="text-align:right;"> 38.41 </td>
   <td style="text-align:right;"> 0.54 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Constant - Varied </td>
   <td style="text-align:right;"> 1200 </td>
   <td style="text-align:right;"> 33.96 </td>
   <td style="text-align:right;"> -11.9 </td>
   <td style="text-align:right;"> 81.01 </td>
   <td style="text-align:right;"> 0.92 </td>
  </tr>
</tbody>
</table>

`````

:::
:::




The model predicting absolute deviation (dist) showed clear effects of both training condition and target velocity band (Table X). Overall, the varied training group showed a larger deviation relative to the constant training group (β = 157.44, 95% CI \[60.53, 254.9\]). Deviation also depended on target velocity band, with lower bands showing less deviation. See @tbl-e1-bmm-dist for full model output.


#### Discrimination between bands

In addition to accuracy/deviation, we also assessed the ability of participants to reliably discriminate between the velocity bands (i.e. responding differently when prompted for band 600-800 than when prompted for band 150-350). @tbl-e1-test-nf-vx shows descriptive statistics of this measure, and Figure 1 visualizes the full distributions of throws for each combination of condition and velocity band. To quantify discrimination, we again fit Bayesian Mixed Models as above, but this time the dependent variable was the raw x velocity generated by participants on each testing trial.




```{=tex}
\begin{equation}
vx_{ij} = \beta_0 + \beta_1 \cdot condit_{ij} + \beta_2 \cdot bandInt_{ij} + \beta_3 \cdot condit_{ij} \cdot bandInt_{ij} + b_{0i} + b_{1i} \cdot bandInt_{ij} + \epsilon_{ij}
\end{equation}
```


::: {.cell layout-align="center"}
::: {.cell-output-display}
![E1 testing x velocities. Translucent bands with dash lines indicate the correct range for each velocity band.](manuscript_files/figure-odt/fig-e1-test-vx-1.png){#fig-e1-test-vx fig-align='center' width=100%}
:::
:::

::: {#tbl-e1-test-nf-vx .cell layout-ncol="1" layout-align="center" tbl-cap='Testing vx - Empirical Summary'}
::: {.cell-output-display}

`````{=html}
<table>
 <thead>
  <tr>
   <th style="text-align:left;"> Band </th>
   <th style="text-align:left;"> Band Type </th>
   <th style="text-align:right;"> Mean </th>
   <th style="text-align:right;"> Median </th>
   <th style="text-align:right;"> Sd </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> 100-300 </td>
   <td style="text-align:left;"> Extrapolation </td>
   <td style="text-align:right;"> 524 </td>
   <td style="text-align:right;"> 448 </td>
   <td style="text-align:right;"> 327 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 350-550 </td>
   <td style="text-align:left;"> Extrapolation </td>
   <td style="text-align:right;"> 659 </td>
   <td style="text-align:right;"> 624 </td>
   <td style="text-align:right;"> 303 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 600-800 </td>
   <td style="text-align:left;"> Extrapolation </td>
   <td style="text-align:right;"> 770 </td>
   <td style="text-align:right;"> 724 </td>
   <td style="text-align:right;"> 300 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 800-1000 </td>
   <td style="text-align:left;"> Trained </td>
   <td style="text-align:right;"> 1001 </td>
   <td style="text-align:right;"> 940 </td>
   <td style="text-align:right;"> 357 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 1000-1200 </td>
   <td style="text-align:left;"> Extrapolation </td>
   <td style="text-align:right;"> 1167 </td>
   <td style="text-align:right;"> 1104 </td>
   <td style="text-align:right;"> 430 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 1200-1400 </td>
   <td style="text-align:left;"> Extrapolation </td>
   <td style="text-align:right;"> 1283 </td>
   <td style="text-align:right;"> 1225 </td>
   <td style="text-align:right;"> 483 </td>
  </tr>
</tbody>
</table>

`````

:::

::: {.cell-output-display}

`````{=html}
<table>
 <thead>
  <tr>
   <th style="text-align:left;"> Band </th>
   <th style="text-align:left;"> Band Type </th>
   <th style="text-align:right;"> Mean </th>
   <th style="text-align:right;"> Median </th>
   <th style="text-align:right;"> Sd </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> 100-300 </td>
   <td style="text-align:left;"> Extrapolation </td>
   <td style="text-align:right;"> 664 </td>
   <td style="text-align:right;"> 533 </td>
   <td style="text-align:right;"> 448 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 350-550 </td>
   <td style="text-align:left;"> Extrapolation </td>
   <td style="text-align:right;"> 768 </td>
   <td style="text-align:right;"> 677 </td>
   <td style="text-align:right;"> 402 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 600-800 </td>
   <td style="text-align:left;"> Extrapolation </td>
   <td style="text-align:right;"> 876 </td>
   <td style="text-align:right;"> 813 </td>
   <td style="text-align:right;"> 390 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 800-1000 </td>
   <td style="text-align:left;"> Trained </td>
   <td style="text-align:right;"> 1064 </td>
   <td style="text-align:right;"> 1029 </td>
   <td style="text-align:right;"> 370 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 1000-1200 </td>
   <td style="text-align:left;"> Trained </td>
   <td style="text-align:right;"> 1180 </td>
   <td style="text-align:right;"> 1179 </td>
   <td style="text-align:right;"> 372 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 1200-1400 </td>
   <td style="text-align:left;"> Trained </td>
   <td style="text-align:right;"> 1265 </td>
   <td style="text-align:right;"> 1249 </td>
   <td style="text-align:right;"> 412 </td>
  </tr>
</tbody>
</table>

`````

:::
:::

::: {#tbl-e1-bmm-vx .cell layout-align="center" tbl-cap='Experiment 1. Bayesian Mixed Model Predicting Vx as a function of condition (Constant vs. Varied) and Velocity Band' tbl-subcap='["Model fit to all 6 bands","Model fit to 3 extrapolation bands"]'}
::: {.cell-output-display}

`````{=html}
<table>
<caption>Fit to all 6 bands</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> Term </th>
   <th style="text-align:right;"> Estimate </th>
   <th style="text-align:right;"> 95% CrI Lower </th>
   <th style="text-align:right;"> 95% CrI Upper </th>
   <th style="text-align:right;"> pd </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Intercept </td>
   <td style="text-align:right;"> 408.55 </td>
   <td style="text-align:right;"> 327.00 </td>
   <td style="text-align:right;"> 490.61 </td>
   <td style="text-align:right;"> 1.00 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> conditVaried </td>
   <td style="text-align:right;"> 164.05 </td>
   <td style="text-align:right;"> 45.50 </td>
   <td style="text-align:right;"> 278.85 </td>
   <td style="text-align:right;"> 1.00 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Band </td>
   <td style="text-align:right;"> 0.71 </td>
   <td style="text-align:right;"> 0.62 </td>
   <td style="text-align:right;"> 0.80 </td>
   <td style="text-align:right;"> 1.00 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> condit*Band </td>
   <td style="text-align:right;"> -0.14 </td>
   <td style="text-align:right;"> -0.26 </td>
   <td style="text-align:right;"> -0.01 </td>
   <td style="text-align:right;"> 0.98 </td>
  </tr>
</tbody>
</table>

`````

:::

::: {.cell-output-display}

`````{=html}
<table>
<caption>Fit to 3 extrapolation bands</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> Term </th>
   <th style="text-align:right;"> Estimate </th>
   <th style="text-align:right;"> 95% CrI Lower </th>
   <th style="text-align:right;"> 95% CrI Upper </th>
   <th style="text-align:right;"> pd </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Intercept </td>
   <td style="text-align:right;"> 478.47 </td>
   <td style="text-align:right;"> 404.00 </td>
   <td style="text-align:right;"> 551.45 </td>
   <td style="text-align:right;"> 1.00 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> conditVaried </td>
   <td style="text-align:right;"> 142.04 </td>
   <td style="text-align:right;"> 37.17 </td>
   <td style="text-align:right;"> 247.59 </td>
   <td style="text-align:right;"> 1.00 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Band </td>
   <td style="text-align:right;"> 0.50 </td>
   <td style="text-align:right;"> 0.42 </td>
   <td style="text-align:right;"> 0.57 </td>
   <td style="text-align:right;"> 1.00 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> condit*Band </td>
   <td style="text-align:right;"> -0.07 </td>
   <td style="text-align:right;"> -0.17 </td>
   <td style="text-align:right;"> 0.04 </td>
   <td style="text-align:right;"> 0.89 </td>
  </tr>
</tbody>
</table>

`````

:::
:::



See @tbl-e1-bmm-vx for the full model results. The estimated coefficient for training condition ($B$ = 164.05, 95% CrI \[45.5, 278.85\]) suggests that the varied group tends to produce harder throws than the constant group, but is not in and of itself useful for assessing discrimination. Most relevant to the issue of discrimination is the slope on Velocity Band ($B$ = 0.71, 95% CrI \[0.62, 0.8\]). Although the median slope does fall underneath the ideal of value of 1, the fact that the 95% credible interval does not contain 0 provides strong evidence that participants exhibited some discrimination between bands. The estimate for the interaction between slope and condition ($B$ = -0.14, 95% CrI \[-0.26, -0.01\]), suggests that the discrimination was somewhat modulated by training condition, with the varied participants showing less senitivity between vands than the constant condition. This difference is depicted visually in @fig-e1-bmm-vx.@tbl-e1-slope-quartile shows the average slope coefficients for varied and constant participants separately for each quartile. The constant participant participants appear to have larger slopes across quartiles, but the difference between conditions may be less pronounced for the top quartiles of subjects who show the strongest discrimination. Figure @fig-e1-bmm-bx2 shows the distributions of slope values for each participant, and the compares the probability density of slope coefficients between training conditions. @fig-e1-indv-slopes 

The second model, which focused solely on extrapolation bands, revealed similar patterns. The Velocity Band term ($B$ = 0.5, 95% CrI \[0.42, 0.57\]) still demonstrates a high degree of discrimination ability. However, the posterior distribution for interaction term ($B$ = -0.07, 95% CrI \[-0.17, 0.04\] ) does across over 0, suggesting that the evidence for decreased discrimination ability for the varied participants is not as strong when considering only the three extrapolation bands.





::: {#fig-e1-bmm-vx .cell layout-ncol="2" layout-align="center"}
::: {.cell-output-display}
![Model fit to all 6 bands](manuscript_files/figure-odt/fig-e1-bmm-vx-1.png){#fig-e1-bmm-vx-1 fig-align='center' width=100%}
:::

::: {.cell-output-display}
![Model fit to only 3 extrapolation bands](manuscript_files/figure-odt/fig-e1-bmm-vx-2.png){#fig-e1-bmm-vx-2 fig-align='center' width=100%}
:::

Conditional effect of training condition and Band. Ribbons indicate 95% HDI. The steepness of the lines serves as an indicator of how well participants discriminated between velocity bands.
:::

::: {#tbl-e1-slope-quartile .cell layout-align="center" tbl-cap='Slope coefficients by quartile, per condition'}
::: {.cell-output-display}

`````{=html}
<table>
 <thead>
  <tr>
   <th style="text-align:left;"> Condition </th>
   <th style="text-align:right;"> Q_0%_mean </th>
   <th style="text-align:right;"> Q_25%_mean </th>
   <th style="text-align:right;"> Q_50%_mean </th>
   <th style="text-align:right;"> Q_75%_mean </th>
   <th style="text-align:right;"> Q_100%_mean </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Constant </td>
   <td style="text-align:right;"> -0.105 </td>
   <td style="text-align:right;"> 0.482 </td>
   <td style="text-align:right;"> 0.689 </td>
   <td style="text-align:right;"> 0.932 </td>
   <td style="text-align:right;"> 1.39 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Varied </td>
   <td style="text-align:right;"> -0.207 </td>
   <td style="text-align:right;"> 0.265 </td>
   <td style="text-align:right;"> 0.590 </td>
   <td style="text-align:right;"> 0.900 </td>
   <td style="text-align:right;"> 1.30 </td>
  </tr>
</tbody>
</table>

`````

:::
:::

::: {#fig-e1-bmm-bx2 .cell layout-ncol="2" layout-align="center"}
::: {.cell-output-display}
![Slope estimates by participant - ordered from lowest to highest within each condition. ](manuscript_files/figure-odt/fig-e1-bmm-bx2-1.png){#fig-e1-bmm-bx2-1 fig-align='center' width=100%}
:::

::: {.cell-output-display}
![Destiny of slope coefficients by training group](manuscript_files/figure-odt/fig-e1-bmm-bx2-2.png){#fig-e1-bmm-bx2-2 fig-align='center' width=100%}
:::

Slope distributions between condition
:::

::: {#fig-e1-indv-slopes .cell layout-align="center"}
::: {.cell-output-display}
![subset with largest slopes](manuscript_files/figure-odt/fig-e1-indv-slopes-1.png){#fig-e1-indv-slopes-1 fig-align='center' width=100%}
:::

::: {.cell-output-display}
![subset with smallest slopes](manuscript_files/figure-odt/fig-e1-indv-slopes-2.png){#fig-e1-indv-slopes-2 fig-align='center' width=100%}
:::

Subset of Varied and Constant Participants with the smallest and largest estimated slope values. Red lines represent the best fitting line for each participant, gray lines are 200 random samples from the posterior distribution. Colored points and intervals at each band represent the empirical median and 95% HDI.
:::




# Experiment 2







@fig-design-e2 illustrates the design of Experiment 2. The stages of the experiment (i.e. training, testing no-feedback, test with feedback), are identical to that of Experiment 1. The only change is that Experiment 2 participants train, and then test, on bands in the reverse order of Experiment 1 (i.e. training on the softer bands; and testing on the harder bands). 




:::{.cell .column-screen-inset-right fig-width="6" fig-height="2.5" fig-responsive=false}

:::{.cell-output-display}

:::{#fig-design-e2}

:::{}

![](manuscript_files/figure-odt/dot-figure-4.png){width="6in" height="2.5in"}
:::


Experiment 2 Design. Constant and Varied participants complete different training conditions. The training and testing bands are the reverse of Experiment 1. 
:::
:::
:::










## E2 Results

### Testing Phase - No feedback. 

In the first part of the testing phase, participants are tested from each of the velocity bands, and receive no feedback after each throw. 


#### Deviation From Target Band

Descriptive summaries testing deviation data are provided in @tbl-e2-test-nf-deviation and @fig-e2-test-dev. 
To model differences in accuracy between groups, we used Bayesian mixed effects regression models to the trial level data from the testing phase. The primary model predicted the absolute deviation from the target velocity band (dist) as a function of training condition (condit), target velocity band (band), and their interaction, with random intercepts and slopes for each participant (id). 

\begin{equation}
dist_{ij} = \beta_0 + \beta_1 \cdot condit_{ij} + \beta_2 \cdot band_{ij} + \beta_3 \cdot condit_{ij} \cdot band_{ij} + b_{0i} + b_{1i} \cdot band_{ij} + \epsilon_{ij}
\end{equation}



::: {#tbl-e2-test-nf-deviation .cell layout-align="center" tbl-cap='Testing Deviation - Empirical Summary' tbl-subcap='["Constant Testing - Deviation","Varied Testing - Deviation"]'}
::: {.cell-output-display}

`````{=html}
<table>
 <thead>
  <tr>
   <th style="text-align:left;"> Band </th>
   <th style="text-align:left;"> Band Type </th>
   <th style="text-align:right;"> Mean </th>
   <th style="text-align:right;"> Median </th>
   <th style="text-align:right;"> Sd </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> 100-300 </td>
   <td style="text-align:left;"> Extrapolation </td>
   <td style="text-align:right;"> 206 </td>
   <td style="text-align:right;"> 48 </td>
   <td style="text-align:right;"> 317 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 350-550 </td>
   <td style="text-align:left;"> Extrapolation </td>
   <td style="text-align:right;"> 194 </td>
   <td style="text-align:right;"> 86 </td>
   <td style="text-align:right;"> 268 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 600-800 </td>
   <td style="text-align:left;"> Trained </td>
   <td style="text-align:right;"> 182 </td>
   <td style="text-align:right;"> 112 </td>
   <td style="text-align:right;"> 240 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 800-1000 </td>
   <td style="text-align:left;"> Extrapolation </td>
   <td style="text-align:right;"> 200 </td>
   <td style="text-align:right;"> 129 </td>
   <td style="text-align:right;"> 233 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 1000-1200 </td>
   <td style="text-align:left;"> Extrapolation </td>
   <td style="text-align:right;"> 238 </td>
   <td style="text-align:right;"> 190 </td>
   <td style="text-align:right;"> 234 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 1200-1400 </td>
   <td style="text-align:left;"> Extrapolation </td>
   <td style="text-align:right;"> 311 </td>
   <td style="text-align:right;"> 254 </td>
   <td style="text-align:right;"> 288 </td>
  </tr>
</tbody>
</table>

`````

:::

::: {.cell-output-display}

`````{=html}
<table>
 <thead>
  <tr>
   <th style="text-align:left;"> Band </th>
   <th style="text-align:left;"> Band Type </th>
   <th style="text-align:right;"> Mean </th>
   <th style="text-align:right;"> Median </th>
   <th style="text-align:right;"> Sd </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> 100-300 </td>
   <td style="text-align:left;"> Trained </td>
   <td style="text-align:right;"> 153 </td>
   <td style="text-align:right;"> 25 </td>
   <td style="text-align:right;"> 266 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 350-550 </td>
   <td style="text-align:left;"> Trained </td>
   <td style="text-align:right;"> 138 </td>
   <td style="text-align:right;"> 53 </td>
   <td style="text-align:right;"> 233 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 600-800 </td>
   <td style="text-align:left;"> Trained </td>
   <td style="text-align:right;"> 160 </td>
   <td style="text-align:right;"> 120 </td>
   <td style="text-align:right;"> 183 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 800-1000 </td>
   <td style="text-align:left;"> Extrapolation </td>
   <td style="text-align:right;"> 261 </td>
   <td style="text-align:right;"> 207 </td>
   <td style="text-align:right;"> 257 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 1000-1200 </td>
   <td style="text-align:left;"> Extrapolation </td>
   <td style="text-align:right;"> 305 </td>
   <td style="text-align:right;"> 258 </td>
   <td style="text-align:right;"> 273 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 1200-1400 </td>
   <td style="text-align:left;"> Extrapolation </td>
   <td style="text-align:right;"> 363 </td>
   <td style="text-align:right;"> 314 </td>
   <td style="text-align:right;"> 297 </td>
  </tr>
</tbody>
</table>

`````

:::
:::

::: {.cell layout-align="center"}
::: {.cell-output-display}
![E2. Deviations from target band during testing without feedback stage.](manuscript_files/figure-odt/fig-e2-test-dev-1.png){#fig-e2-test-dev fig-align='center' width=100%}
:::
:::

::: {#tbl-e2-bmm-dist .cell layout-align="center" tbl-cap='Experiment 2. Bayesian Mixed Model predicting absolute deviation as a function of condition (Constant vs. Varied) and Velocity Band'}
::: {.cell-output-display}

`````{=html}
<table>
 <thead>
  <tr>
   <th style="text-align:left;"> Term </th>
   <th style="text-align:right;"> Estimate </th>
   <th style="text-align:right;"> 95% CrI Lower </th>
   <th style="text-align:right;"> 95% CrI Upper </th>
   <th style="text-align:right;"> pd </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Intercept </td>
   <td style="text-align:right;"> 151.71 </td>
   <td style="text-align:right;"> 90.51 </td>
   <td style="text-align:right;"> 215.86 </td>
   <td style="text-align:right;"> 1.00 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> conditVaried </td>
   <td style="text-align:right;"> -70.33 </td>
   <td style="text-align:right;"> -156.87 </td>
   <td style="text-align:right;"> 16.66 </td>
   <td style="text-align:right;"> 0.94 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Band </td>
   <td style="text-align:right;"> 0.10 </td>
   <td style="text-align:right;"> 0.02 </td>
   <td style="text-align:right;"> 0.18 </td>
   <td style="text-align:right;"> 1.00 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> condit*Band </td>
   <td style="text-align:right;"> 0.12 </td>
   <td style="text-align:right;"> 0.02 </td>
   <td style="text-align:right;"> 0.23 </td>
   <td style="text-align:right;"> 0.99 </td>
  </tr>
</tbody>
</table>

`````

:::

::: {.cell-output-display}

`````{=html}
<table>
<caption>Contrasts</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> contrast </th>
   <th style="text-align:right;"> Band </th>
   <th style="text-align:right;"> value </th>
   <th style="text-align:right;"> lower </th>
   <th style="text-align:right;"> upper </th>
   <th style="text-align:right;"> pd </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Constant - Varied </td>
   <td style="text-align:right;"> 100 </td>
   <td style="text-align:right;"> 57.6 </td>
   <td style="text-align:right;"> -20.5 </td>
   <td style="text-align:right;"> 135.32 </td>
   <td style="text-align:right;"> 0.93 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Constant - Varied </td>
   <td style="text-align:right;"> 350 </td>
   <td style="text-align:right;"> 26.6 </td>
   <td style="text-align:right;"> -30.9 </td>
   <td style="text-align:right;"> 83.84 </td>
   <td style="text-align:right;"> 0.83 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Constant - Varied </td>
   <td style="text-align:right;"> 600 </td>
   <td style="text-align:right;"> -4.3 </td>
   <td style="text-align:right;"> -46.7 </td>
   <td style="text-align:right;"> 38.52 </td>
   <td style="text-align:right;"> 0.58 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Constant - Varied </td>
   <td style="text-align:right;"> 800 </td>
   <td style="text-align:right;"> -29.3 </td>
   <td style="text-align:right;"> -69.4 </td>
   <td style="text-align:right;"> 11.29 </td>
   <td style="text-align:right;"> 0.92 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Constant - Varied </td>
   <td style="text-align:right;"> 1000 </td>
   <td style="text-align:right;"> -54.6 </td>
   <td style="text-align:right;"> -101.1 </td>
   <td style="text-align:right;"> -5.32 </td>
   <td style="text-align:right;"> 0.98 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Constant - Varied </td>
   <td style="text-align:right;"> 1200 </td>
   <td style="text-align:right;"> -79.6 </td>
   <td style="text-align:right;"> -139.5 </td>
   <td style="text-align:right;"> -15.45 </td>
   <td style="text-align:right;"> 0.99 </td>
  </tr>
</tbody>
</table>

`````

:::
:::




The model predicting absolute deviation showed a modest tendency for the varied training group to have lower deviation compared to the constant training group (β = -70.33, 95% CI \[-156.87, 16.66\]),with 94% of the posterior distribution being less than 0. This suggests a potential benefit of training with variation, though the evidence is not definitive.



# Experiment 3







The major manipulation adjustment of experiment 3 is for participants to receive ordinal feedback during training, in contrast to the continuous feedback of the earlier experiments. Ordinal feedback informs participants whether a throw was too soft, too hard, or fell within the target velocity range. Experiment 3 participants were randomly assigned to both a training condition (Constant vs. Varied) and a Band Order condition (original order used in Experiment 1, or the Reverse order of Experiment 2). 


## Results

### Testing Phase - No feedback. 

In the first part of the testing phase, participants are tested from each of the velocity bands, and receive no feedback after each throw. Note that these no-feedback testing trials are identical to those of Experiment 1 and 2, as the ordinal feedback only occurs during the training phase, and final testing phase, of Experiment 3. 


#### Deviation From Target Band

Descriptive summaries testing deviation data are provided in @tbl-e3-test-nf-deviation and @fig-e3-test-dev. 
To model differences in accuracy between groups, we fit Bayesian mixed effects regression models to the trial level data from the testing phase. The primary model predicted the absolute deviation from the target velocity band (dist) as a function of training condition (condit), target velocity band (band), and their interaction, with random intercepts and slopes for each participant (id). 




::: {#tbl-e3-test-nf-deviation .cell layout-align="center" tbl-cap='Testing Deviation - Empirical Summary' tbl-subcap='["Constant Testing - Deviation","Varied Testing - Deviation"]'}
::: {.cell-output-display}

`````{=html}
<table>
 <thead>
  <tr>
   <th style="text-align:left;"> Band </th>
   <th style="text-align:left;"> Band Type </th>
   <th style="text-align:right;"> Mean </th>
   <th style="text-align:right;"> Median </th>
   <th style="text-align:right;"> Sd </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> 100-300 </td>
   <td style="text-align:left;"> Extrapolation </td>
   <td style="text-align:right;"> 396 </td>
   <td style="text-align:right;"> 325 </td>
   <td style="text-align:right;"> 350 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 350-550 </td>
   <td style="text-align:left;"> Extrapolation </td>
   <td style="text-align:right;"> 278 </td>
   <td style="text-align:right;"> 176 </td>
   <td style="text-align:right;"> 299 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 600-800 </td>
   <td style="text-align:left;"> Extrapolation </td>
   <td style="text-align:right;"> 173 </td>
   <td style="text-align:right;"> 102 </td>
   <td style="text-align:right;"> 215 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 800-1000 </td>
   <td style="text-align:left;"> Trained </td>
   <td style="text-align:right;"> 225 </td>
   <td style="text-align:right;"> 126 </td>
   <td style="text-align:right;"> 284 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 1000-1200 </td>
   <td style="text-align:left;"> Extrapolation </td>
   <td style="text-align:right;"> 253 </td>
   <td style="text-align:right;"> 192 </td>
   <td style="text-align:right;"> 271 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 1200-1400 </td>
   <td style="text-align:left;"> Extrapolation </td>
   <td style="text-align:right;"> 277 </td>
   <td style="text-align:right;"> 210 </td>
   <td style="text-align:right;"> 262 </td>
  </tr>
</tbody>
</table>

`````

:::

::: {.cell-output-display}

`````{=html}
<table>
 <thead>
  <tr>
   <th style="text-align:left;"> Band </th>
   <th style="text-align:left;"> Band Type </th>
   <th style="text-align:right;"> Mean </th>
   <th style="text-align:right;"> Median </th>
   <th style="text-align:right;"> Sd </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> 100-300 </td>
   <td style="text-align:left;"> Extrapolation </td>
   <td style="text-align:right;"> 383 </td>
   <td style="text-align:right;"> 254 </td>
   <td style="text-align:right;"> 385 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 350-550 </td>
   <td style="text-align:left;"> Extrapolation </td>
   <td style="text-align:right;"> 287 </td>
   <td style="text-align:right;"> 154 </td>
   <td style="text-align:right;"> 318 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 600-800 </td>
   <td style="text-align:left;"> Extrapolation </td>
   <td style="text-align:right;"> 213 </td>
   <td style="text-align:right;"> 140 </td>
   <td style="text-align:right;"> 244 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 800-1000 </td>
   <td style="text-align:left;"> Trained </td>
   <td style="text-align:right;"> 199 </td>
   <td style="text-align:right;"> 142 </td>
   <td style="text-align:right;"> 209 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 1000-1200 </td>
   <td style="text-align:left;"> Trained </td>
   <td style="text-align:right;"> 222 </td>
   <td style="text-align:right;"> 163 </td>
   <td style="text-align:right;"> 221 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 1200-1400 </td>
   <td style="text-align:left;"> Trained </td>
   <td style="text-align:right;"> 281 </td>
   <td style="text-align:right;"> 227 </td>
   <td style="text-align:right;"> 246 </td>
  </tr>
</tbody>
</table>

`````

:::

::: {.cell-output-display}

`````{=html}
<table>
 <thead>
  <tr>
   <th style="text-align:left;"> Band </th>
   <th style="text-align:left;"> Band Type </th>
   <th style="text-align:right;"> Mean </th>
   <th style="text-align:right;"> Median </th>
   <th style="text-align:right;"> Sd </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> 100-300 </td>
   <td style="text-align:left;"> Extrapolation </td>
   <td style="text-align:right;"> 403 </td>
   <td style="text-align:right;"> 334 </td>
   <td style="text-align:right;"> 383 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 350-550 </td>
   <td style="text-align:left;"> Extrapolation </td>
   <td style="text-align:right;"> 246 </td>
   <td style="text-align:right;"> 149 </td>
   <td style="text-align:right;"> 287 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 600-800 </td>
   <td style="text-align:left;"> Trained </td>
   <td style="text-align:right;"> 155 </td>
   <td style="text-align:right;"> 82 </td>
   <td style="text-align:right;"> 209 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 800-1000 </td>
   <td style="text-align:left;"> Extrapolation </td>
   <td style="text-align:right;"> 207 </td>
   <td style="text-align:right;"> 151 </td>
   <td style="text-align:right;"> 241 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 1000-1200 </td>
   <td style="text-align:left;"> Extrapolation </td>
   <td style="text-align:right;"> 248 </td>
   <td style="text-align:right;"> 220 </td>
   <td style="text-align:right;"> 222 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 1200-1400 </td>
   <td style="text-align:left;"> Extrapolation </td>
   <td style="text-align:right;"> 322 </td>
   <td style="text-align:right;"> 281 </td>
   <td style="text-align:right;"> 264 </td>
  </tr>
</tbody>
</table>

`````

:::

::: {.cell-output-display}

`````{=html}
<table>
 <thead>
  <tr>
   <th style="text-align:left;"> Band </th>
   <th style="text-align:left;"> Band Type </th>
   <th style="text-align:right;"> Mean </th>
   <th style="text-align:right;"> Median </th>
   <th style="text-align:right;"> Sd </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> 100-300 </td>
   <td style="text-align:left;"> Trained </td>
   <td style="text-align:right;"> 153 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 307 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 350-550 </td>
   <td style="text-align:left;"> Trained </td>
   <td style="text-align:right;"> 147 </td>
   <td style="text-align:right;"> 55 </td>
   <td style="text-align:right;"> 258 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 600-800 </td>
   <td style="text-align:left;"> Trained </td>
   <td style="text-align:right;"> 159 </td>
   <td style="text-align:right;"> 107 </td>
   <td style="text-align:right;"> 192 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 800-1000 </td>
   <td style="text-align:left;"> Extrapolation </td>
   <td style="text-align:right;"> 221 </td>
   <td style="text-align:right;"> 160 </td>
   <td style="text-align:right;"> 235 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 1000-1200 </td>
   <td style="text-align:left;"> Extrapolation </td>
   <td style="text-align:right;"> 244 </td>
   <td style="text-align:right;"> 185 </td>
   <td style="text-align:right;"> 235 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 1200-1400 </td>
   <td style="text-align:left;"> Extrapolation </td>
   <td style="text-align:right;"> 324 </td>
   <td style="text-align:right;"> 264 </td>
   <td style="text-align:right;"> 291 </td>
  </tr>
</tbody>
</table>

`````

:::
:::

::: {.cell layout-align="center"}
::: {.cell-output-display}
![e3. Deviations from target band during testing without feedback stage.](manuscript_files/figure-odt/fig-e3-test-dev-1.png){#fig-e3-test-dev fig-align='center' width=100%}
:::
:::

::: {#tbl-e3-bmm-dist .cell layout-align="center" tbl-cap='Experiment 3. Bayesian Mixed Model predicting absolute deviation as a function of condition (Constant vs. Varied) and Velocity Band'}
::: {.cell-output-display}

`````{=html}
<table>
 <thead>
  <tr>
   <th style="text-align:left;"> Term </th>
   <th style="text-align:right;"> Estimate </th>
   <th style="text-align:right;"> 95% CrI Lower </th>
   <th style="text-align:right;"> 95% CrI Upper </th>
   <th style="text-align:right;"> pd </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Intercept </td>
   <td style="text-align:right;"> 306.47 </td>
   <td style="text-align:right;"> 243.89 </td>
   <td style="text-align:right;"> 368.75 </td>
   <td style="text-align:right;"> 1.00 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> conditVaried </td>
   <td style="text-align:right;"> -90.65 </td>
   <td style="text-align:right;"> -182.79 </td>
   <td style="text-align:right;"> 3.75 </td>
   <td style="text-align:right;"> 0.97 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Band </td>
   <td style="text-align:right;"> -0.07 </td>
   <td style="text-align:right;"> -0.13 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0.97 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> condit*Band </td>
   <td style="text-align:right;"> 0.09 </td>
   <td style="text-align:right;"> -0.01 </td>
   <td style="text-align:right;"> 0.19 </td>
   <td style="text-align:right;"> 0.96 </td>
  </tr>
</tbody>
</table>

`````

:::
:::




The effect of training condition in Experiment 3 showed a similar pattern to Experiment 2, with the varied group tending to have lower deviation than the constant group (β = -90.65, 95% CrI \[-182.79, 3.75\]), with 97% of the posterior distribution falling under 0. 





::: {.cell layout-align="center"}
::: {.cell-output-display}
![e3. Conditioinal Effect of Training Condition and Band. Ribbon indicated 95% Credible Intervals.](manuscript_files/figure-odt/fig-e3-bmm-dist-1.png){#fig-e3-bmm-dist fig-align='center' width=100%}
:::
:::





#### Discrimination between Velocity Bands

In addition to accuracy/deviation. We also assessed the ability of participants to reliably discriminate between the velocity bands (i.e. responding differently when prompted for band 600-800 than when prompted for band 150-350). @tbl-e3-test-nf-vx shows descriptive statistics of this measure, and Figure 1 visualizes the full distributions of throws for each combination of condition and velocity band. To quantify discrimination, we again fit Bayesian Mixed Models as above, but this time the dependent variable was the raw x velocity generated by participants. 

\begin{equation}
vx_{ij} = \beta_0 + \beta_1 \cdot condit_{ij} + \beta_2 \cdot bandInt_{ij} + \beta_3 \cdot condit_{ij} \cdot bandInt_{ij} + b_{0i} + b_{1i} \cdot bandInt_{ij} + \epsilon_{ij}
\end{equation}





::: {.cell .column-screen-inset-right layout-align="center"}
::: {.cell-output-display}
![e3 testing x velocities. Translucent bands with dash lines indicate the correct range for each velocity band.](manuscript_files/figure-odt/fig-e3-test-vx-1.png){#fig-e3-test-vx fig-align='center' width=100%}
:::
:::

::: {#tbl-e3-test-nf-vx .cell layout-align="center" tbl-cap='Testing vx - Empirical Summary' tbl-subcap='["Constant Testing - vx","Varied Testing - vx"]'}
::: {.cell-output-display}

`````{=html}
<table>
 <thead>
  <tr>
   <th style="text-align:left;"> Band </th>
   <th style="text-align:left;"> Band Type </th>
   <th style="text-align:right;"> Mean </th>
   <th style="text-align:right;"> Median </th>
   <th style="text-align:right;"> Sd </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> 100-300 </td>
   <td style="text-align:left;"> Extrapolation </td>
   <td style="text-align:right;"> 680 </td>
   <td style="text-align:right;"> 625 </td>
   <td style="text-align:right;"> 370 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 350-550 </td>
   <td style="text-align:left;"> Extrapolation </td>
   <td style="text-align:right;"> 771 </td>
   <td style="text-align:right;"> 716 </td>
   <td style="text-align:right;"> 357 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 600-800 </td>
   <td style="text-align:left;"> Extrapolation </td>
   <td style="text-align:right;"> 832 </td>
   <td style="text-align:right;"> 786 </td>
   <td style="text-align:right;"> 318 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 800-1000 </td>
   <td style="text-align:left;"> Trained </td>
   <td style="text-align:right;"> 1006 </td>
   <td style="text-align:right;"> 916 </td>
   <td style="text-align:right;"> 417 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 1000-1200 </td>
   <td style="text-align:left;"> Extrapolation </td>
   <td style="text-align:right;"> 1149 </td>
   <td style="text-align:right;"> 1105 </td>
   <td style="text-align:right;"> 441 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 1200-1400 </td>
   <td style="text-align:left;"> Extrapolation </td>
   <td style="text-align:right;"> 1180 </td>
   <td style="text-align:right;"> 1112 </td>
   <td style="text-align:right;"> 443 </td>
  </tr>
</tbody>
</table>

`````

:::

::: {.cell-output-display}

`````{=html}
<table>
 <thead>
  <tr>
   <th style="text-align:left;"> Band </th>
   <th style="text-align:left;"> Band Type </th>
   <th style="text-align:right;"> Mean </th>
   <th style="text-align:right;"> Median </th>
   <th style="text-align:right;"> Sd </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> 100-300 </td>
   <td style="text-align:left;"> Extrapolation </td>
   <td style="text-align:right;"> 667 </td>
   <td style="text-align:right;"> 554 </td>
   <td style="text-align:right;"> 403 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 350-550 </td>
   <td style="text-align:left;"> Extrapolation </td>
   <td style="text-align:right;"> 770 </td>
   <td style="text-align:right;"> 688 </td>
   <td style="text-align:right;"> 383 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 600-800 </td>
   <td style="text-align:left;"> Extrapolation </td>
   <td style="text-align:right;"> 869 </td>
   <td style="text-align:right;"> 814 </td>
   <td style="text-align:right;"> 358 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 800-1000 </td>
   <td style="text-align:left;"> Trained </td>
   <td style="text-align:right;"> 953 </td>
   <td style="text-align:right;"> 928 </td>
   <td style="text-align:right;"> 359 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 1000-1200 </td>
   <td style="text-align:left;"> Trained </td>
   <td style="text-align:right;"> 1072 </td>
   <td style="text-align:right;"> 1066 </td>
   <td style="text-align:right;"> 388 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 1200-1400 </td>
   <td style="text-align:left;"> Trained </td>
   <td style="text-align:right;"> 1144 </td>
   <td style="text-align:right;"> 1093 </td>
   <td style="text-align:right;"> 426 </td>
  </tr>
</tbody>
</table>

`````

:::

::: {.cell-output-display}

`````{=html}
<table>
 <thead>
  <tr>
   <th style="text-align:left;"> Band </th>
   <th style="text-align:left;"> Band Type </th>
   <th style="text-align:right;"> Mean </th>
   <th style="text-align:right;"> Median </th>
   <th style="text-align:right;"> Sd </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> 100-300 </td>
   <td style="text-align:left;"> Extrapolation </td>
   <td style="text-align:right;"> 684 </td>
   <td style="text-align:right;"> 634 </td>
   <td style="text-align:right;"> 406 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 350-550 </td>
   <td style="text-align:left;"> Extrapolation </td>
   <td style="text-align:right;"> 729 </td>
   <td style="text-align:right;"> 679 </td>
   <td style="text-align:right;"> 350 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 600-800 </td>
   <td style="text-align:left;"> Trained </td>
   <td style="text-align:right;"> 776 </td>
   <td style="text-align:right;"> 721 </td>
   <td style="text-align:right;"> 318 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 800-1000 </td>
   <td style="text-align:left;"> Extrapolation </td>
   <td style="text-align:right;"> 941 </td>
   <td style="text-align:right;"> 883 </td>
   <td style="text-align:right;"> 387 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 1000-1200 </td>
   <td style="text-align:left;"> Extrapolation </td>
   <td style="text-align:right;"> 1014 </td>
   <td style="text-align:right;"> 956 </td>
   <td style="text-align:right;"> 403 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 1200-1400 </td>
   <td style="text-align:left;"> Extrapolation </td>
   <td style="text-align:right;"> 1072 </td>
   <td style="text-align:right;"> 1014 </td>
   <td style="text-align:right;"> 442 </td>
  </tr>
</tbody>
</table>

`````

:::

::: {.cell-output-display}

`````{=html}
<table>
 <thead>
  <tr>
   <th style="text-align:left;"> Band </th>
   <th style="text-align:left;"> Band Type </th>
   <th style="text-align:right;"> Mean </th>
   <th style="text-align:right;"> Median </th>
   <th style="text-align:right;"> Sd </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> 100-300 </td>
   <td style="text-align:left;"> Trained </td>
   <td style="text-align:right;"> 392 </td>
   <td style="text-align:right;"> 270 </td>
   <td style="text-align:right;"> 343 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 350-550 </td>
   <td style="text-align:left;"> Trained </td>
   <td style="text-align:right;"> 540 </td>
   <td style="text-align:right;"> 442 </td>
   <td style="text-align:right;"> 343 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 600-800 </td>
   <td style="text-align:left;"> Trained </td>
   <td style="text-align:right;"> 642 </td>
   <td style="text-align:right;"> 588 </td>
   <td style="text-align:right;"> 315 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 800-1000 </td>
   <td style="text-align:left;"> Extrapolation </td>
   <td style="text-align:right;"> 943 </td>
   <td style="text-align:right;"> 899 </td>
   <td style="text-align:right;"> 394 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 1000-1200 </td>
   <td style="text-align:left;"> Extrapolation </td>
   <td style="text-align:right;"> 1081 </td>
   <td style="text-align:right;"> 1048 </td>
   <td style="text-align:right;"> 415 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 1200-1400 </td>
   <td style="text-align:left;"> Extrapolation </td>
   <td style="text-align:right;"> 1185 </td>
   <td style="text-align:right;"> 1129 </td>
   <td style="text-align:right;"> 500 </td>
  </tr>
</tbody>
</table>

`````

:::
:::

::: {#tbl-e3-bmm-vx .cell layout-align="center" tbl-cap='Experiment 3. Bayesian Mixed Model Predicting Vx as a function of condition (Constant vs. Varied) and Velocity Band'}
::: {.cell-output-display}

`````{=html}
<table>
 <thead>
  <tr>
   <th style="text-align:left;"> Term </th>
   <th style="text-align:right;"> Estimate </th>
   <th style="text-align:right;"> 95% CrI Lower </th>
   <th style="text-align:right;"> 95% CrI Upper </th>
   <th style="text-align:right;"> pd </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Intercept </td>
   <td style="text-align:right;"> 607.67 </td>
   <td style="text-align:right;"> 536.02 </td>
   <td style="text-align:right;"> 679.87 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> conditVaried </td>
   <td style="text-align:right;"> -167.76 </td>
   <td style="text-align:right;"> -277.14 </td>
   <td style="text-align:right;"> -64.08 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Band </td>
   <td style="text-align:right;"> 0.44 </td>
   <td style="text-align:right;"> 0.35 </td>
   <td style="text-align:right;"> 0.52 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> condit*Band </td>
   <td style="text-align:right;"> 0.18 </td>
   <td style="text-align:right;"> 0.06 </td>
   <td style="text-align:right;"> 0.31 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
</tbody>
</table>

`````

:::
:::




See @tbl-e3-bmm-vx for the full model results. 

Slope estimates for experiment 3 suggest that participants were capable of distinguishing between velocity bands even when provided only ordinal feedback during training (β = 0.44, 95% CrI \[0.35, 0.52\]). Unlike the previous two experiments, the posterior distribution for the interaction between condition and band was consistently positive, suggestive of superior discrimination for the varied participants 
β = 0.18, 95% CrI \[0.06, 0.31\]. 


# Modeling

In project 1, we applied model-based techniques to quantify and control for the similarity between training and testing experience, which in turn enabled us to account for the difference between varied and constant training via an extended version of a similarity based generalization model. In project 2, we will go a step further, implementing a full process model capable of both 1) producing novel responses and 2) modeling behavior in both the learning and testing stages of the experiment. For this purpose, we will apply the associative learning model (ALM) and the EXAM model of function learning (DeLosh 1997). ALM is a simple connectionist learning model which closely resembles Kruschke's ALCOVE model (Kruscke 1992), with modifications to allow for the generation of continuous responses.

## ALM & Exam Description

@deloshExtrapolationSineQua1997 introduced the associative learning model (ALM), a connectionist model within the popular class of radial-basis networks. ALM was inspired by, and closely resembles Kruschke's influential ALCOVE model of categorization [@kruschkeALCOVEExemplarbasedConnectionist1992]. 

ALM is a localist neural network model, with each input node corresponding to a particular stimulus, and each output node corresponding to a particular response value. The units in the input layer activate as a function of their Gaussian similarity to the input stimulus. So, for example, an input stimulus of value 55 would induce maximal activation of the input unit tuned to 55. Depending on thevalue of the generalization parameter, the nearby units (e.g. 54 and 56; 53 and 57) may also activate to some degree. ALM is structured with input and output nodes that correspond to regions of the stimulus space, and response space, respectively. The units in the input layer activate as a function of their similarity to a presented stimulus. As was the case with the exemplar-based models, similarity in ALM is exponentially decaying function of distance. The input layer is fully connected to the output layer, and the activation for any particular output node is simply the weighted sum of the connection weights between that node and the input activations. The network then produces a response by taking the weighted average of the output units (recall that each output unit has a value corresponding to a particular response). During training, the network receives feedback which activates each output unit as a function of its distance from the ideal level of activation necessary to produce the correct response. The connection weights between input and output units are then updated via the standard delta learning rule, where the magnitude of weight changes are controlled by a learning rate parameter.

See  for a full specification of the equations that define ALM and EXAM.





{{< pagebreak >}}







## Model Table
### ALM Activation & Response
| Step | Equation | Description |
|------|----------|-------------|
| **ALM Activation & Response**  |           |
| Input Activation | $a_i(X) = \frac{e^{-c(X-X_i)^2}}{\sum_{k=1}^M e^{-c(X-X_k)^2}}$ | Activation of each input node $X_i$, is a function of the Gaussian similarity between the node value and stimulus X.  |
| Output Activation                 | $O_j(X) = \sum_{k=1}^M w_{ji} \cdot a_i(X)$ | Activation of each Output unit $O_j$ is the weighted sum of the input activations and association weights.   |
| Output Probability                | $P[Y_j|X] = \frac{O_j(X)}{\sum_{k=1}^M O_k(X)}$ | Each output node has associated response, $Y_j$. The probability of response $Y_j$ is determined by the ratio of output activations.  |
| Mean Output                       | $m(x) = \sum_{j=1}^L Y_j \cdot \frac{O_j(x)}{\sum_{k=1}^M O_k(X)}$   | The response to stimulus x is the weighted average of the response probabilities. |
| **ALM Learning**                  |                                                                                                                                                          |
| Feedback Activation               | $f_j(Z) = e^{-c(Z-Y_j)^2}$                  | After responding, feedback signal Z is presented, activating each output node via the Gaussian similarity to the ideal response. |
| Update Weights                    | $w_{ji}(t + 1) = w_{ji}(t) + \alpha \cdot (f_j(Z(t)) - O_j(X(t)) \cdot a_i(X(t))$    | Delta rule to update weights. Magnitude of weight changes controlled by learning rate parameter alpha.  |  
| **EXAM**                          |                                                                                                                                                          |
| Extrapolation                     | $P[X_i|X] = \frac{a_i(X)}{\sum_{k=1}^M a_k(X)}$ | Novel test stimulus X activates input nodes associated with trained stimuli.     |
|                                   | $E[Y|X_i] = m(X_i) + \frac{m(X_{i+1})-m(X_{i-1})}{X_{i+1}-X_{i-1}} \cdot [X - X_i]$ | Slope value computed from nearest training instances and then added to the response associated with the nearest training instance,m(x)    |               





{{< pagebreak >}}









## Model Fitting and Comparison

Following the procedure used by @mcdanielPredictingTransferPerformance2009, we will assess the ability of both ALM and EXAM to account for the empirical data when fitting the models to 1) only the training data, and 2) both training and testing data. Models will be fit directly to the trial by trial data of each individual participants, both by minimizing the root-mean squared deviation (RMSE), and by maximizing log likelihood. Because ALM has been shown to do poorly at accounting for human patterns extrapolation [@deloshExtrapolationSineQua1997], we will also fit the extended EXAM version of the model, which operates identically to ALM during training, but includes a linear extrapolation mechanism for generating novel responses during testing.

# References