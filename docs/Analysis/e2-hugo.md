---
title: Experiment 2
categories:
  - Analyses
  - R
  - Bayesian
toc: false
format-links: false
page-layout: full
cache: true
code-fold: true
cold-tools: true
lightbox: true
keep-md: true
format:
  html: default
  hugo-md:
    echo: false
    html-math-method: mathjax
    output-file: e2-hugo.md
  gfm:
    echo: true
    output-file: e2-gfm.md
---


<a href="#fig-design-e2" class="quarto-xref">Figure 2</a> illustrates the design of Experiment 2. The stages of the experiment (i.e. training, testing no-feedback, test with feedback), are identical to that of Experiment 1. The only change is that Experiment 2 participants train, and then test, on bands in the reverse order of Experiment 1 (i.e. training on the softer bands; and testing on the harder bands).

<div id="fig-design-e2">

Figure 1: Experiment 2 Design. Constant and Varied participants complete different training conditions. The training and testing bands are the reverse of Experiment 1.
</div>
<img src="../Assets/figs/e2_design.png" id="fig-design-e2"
alt="Figure 2: Experiment 2 Design. Constant and Varied participants complete different training conditions. The training and testing bands are the reverse of Experiment 1." />

<img src="../Assets/figs/e2_train_deviation.png" id="fig-e2-train-dev"
alt="Figure 3: E2. Deviations from target band across training blocks." />
<div id="tbl-e2-train-dist">

| Term         | Estimate | 95% CrI Lower | 95% CrI Upper |  pd |
|:-------------|---------:|--------------:|--------------:|----:|
| Intercept    |    91.01 |         80.67 |        101.26 |   1 |
| conditVaried |    36.15 |         16.35 |         55.67 |   1 |

Table 1: **Experiment 2 - End of training performance**. The Intercept represents the average of the baseline (constant condition), and the conditVaried coefficient reflects the difference between the constant and varied groups. A larger positive coefficient indicates a greater deviation (lower accuracy) for the varied group.
</div>

  

*Training*. <a href="#fig-e2-train-dev" class="quarto-xref">Figure 3</a> presents the deviations across training blocks for both constant and varied training groups. We again compared training performance on the band common to both groups (600-800). The full model results are shown in Table 1. The varied group had a significantly greater deviation than the constant group in the final training block, ( $\beta$ = 36.15, 95% CrI \[16.35, 55.67\]; pd = 99.95%).

<div id="tbl-e2-bmm-dist">

| Term                               | Estimate | 95% CrI Lower | 95% CrI Upper |   pd |
|:-----------------------------|---------:|-------------:|-------------:|-----:|
| Intercept                          |   190.91 |        125.03 |        259.31 | 1.00 |
| conditVaried                       |   -20.58 |        -72.94 |         33.08 | 0.78 |
| bandTypeExtrapolation              |    38.09 |         -6.94 |         83.63 | 0.95 |
| conditVaried:bandTypeExtrapolation |    82.00 |         41.89 |        121.31 | 1.00 |

Table 2: **Experiment 2 testing accuracy**. Main effects of condition and band type (training vs. extrapolation), and the interaction between the two factors. Larger coefficient estimates indicate larger deviations from the baselines (constant & trained bands) - and a positive interaction coefficient indicates disproporionate deviation for the varied condition on the extrapolation bands
</div>

*Testing.* The analysis of testing accuracy examined deviations from the target band as influenced by training condition (Varied vs. Constant) and band type (training vs. extrapolation bands). The results, summarized in <a href="#tbl-e2-bmm-dist" class="quarto-xref">Table 2</a>, reveal no significant main effect of training condition ($\beta$ = -20.58, 95% CrI \[-72.94, 33.08\]; pd = 77.81%). However, the interaction between training condition and band type was significant ($\beta$ = 82, 95% CrI \[41.89, 121.31\]; pd = 100%), with the varied group showing disproportionately larger deviations compared to the constant group on the extrapolation bands (see <a href="#fig-e2-test-dev" class="quarto-xref">Figure 4</a>).

<img src="../Assets/figs/e2_test-dev.png" id="fig-e2-test-dev"
alt="Figure 4: E2. A) Deviations from target band during testing without feedback stage. B) Estimated marginal means for the interaction between training condition and band type. Error bars represent 95% confidence intervals." />

### Discimination

<div id="tbl-e2-bmm-vx">

| Term         | Estimate | 95% CrI Lower | 95% CrI Upper |   pd |
|:-------------|---------:|--------------:|--------------:|-----:|
| Intercept    |   362.64 |        274.85 |        450.02 | 1.00 |
| conditVaried |    -8.56 |       -133.97 |        113.98 | 0.55 |
| Band         |     0.71 |          0.58 |          0.84 | 1.00 |
| condit\*Band |    -0.06 |         -0.24 |          0.13 | 0.73 |

Table 3: Experiment 2. Bayesian Mixed Model Predicting Vx as a function of condition (Constant vs. Varied) and Velocity Band
</div>

Finally, to assess the ability of both conditions to discriminate between velocity bands, we fit a model predicting velocity as a function of training condition and velocity band, with random intercepts and random slopes for each participant. The full model results are shown in <a href="#tbl-e2-bmm-vx" class="quarto-xref">Table 3</a>. The overall slope on target velocity band predictor was significantly positive, ($\beta$ = 0.71, 95% CrI \[0.58, 0.84\]; pd= 100%), indicating that participants exhibited discrimination between bands. The interaction between slope and condition was not significant, ($\beta$ = -0.06, 95% CrI \[-0.24, 0.13\]; pd= `rintCoef2$pd`), suggesting that the two conditions did not differ in their ability to discriminate between bands (see <a href="#fig-e2-test-vx" class="quarto-xref">Figure 5</a>).



<img src="../Assets/figs/e2_test-vx.png" id="fig-e2-test-vx"
alt="Figure 5: Experiment 2. Conditional effect of training condition and Band. Ribbons indicate 95% HDI. The steepness of the lines serves as an indicator of how well participants discriminated between velocity bands." />

## E2 Discussion
