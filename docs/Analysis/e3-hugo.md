---
title: Experiment 3
categories:
  - Analyses
  - R
  - Bayesian
toc: false
format-links: false
page-layout: full
lightbox: true
code-fold: true
cold-tools: true
cache: true
keep-md: true
format:
  html: default
  hugo-md:
    echo: false
    html-math-method: mathjax
    output-file: e3-hugo.md
  gfm:
    echo: true
    output-file: e3-gfm.md
---


The major manipulation adjustment of experiment 3 is for participants to receive ordinal feedback during training, in contrast to the continuous feedback of the earlier experiments. Ordinal feedback informs participants whether a throw was too soft, too hard, or fell within the target velocity range. Experiment 3 participants were randomly assigned to both a training condition (Constant vs. Varied) and a Band Order condition (original order used in Experiment 1, or the Reverse order of Experiment 2).

<div id="tbl-e3-train-dist">

| Term                          | Estimate | 95% CrI Lower | 95% CrI Upper |   pd |
|:---------------------------|---------:|-------------:|-------------:|------:|
| Intercept                     |   121.86 |        109.24 |        134.60 | 1.00 |
| conditVaried                  |    64.93 |         36.99 |         90.80 | 1.00 |
| bandOrderReverse              |     1.11 |        -16.02 |         18.16 | 0.55 |
| conditVaried:bandOrderReverse |   -77.02 |       -114.16 |        -39.61 | 1.00 |

Table 1: **Experiment 3 - End of training performance**. The Intercept represents the average of the baseline (constant condition), and the conditVaried coefficient reflects the difference between the constant and varied groups. A larger positive coefficient indicates a greater deviation (lower accuracy) for the varied group.
</div>
<img src="../Assets/figs/e3_train_deviation.png" id="fig-e3-train-dev"
alt="Figure 1: E3. Deviations from target band across training blocks." />
<div id="tbl-e3-bmm-dist">

| Term                                                | Estimate | 95% CrI Lower | 95% CrI Upper |   pd |
|:------------------------------------|-------:|-----------:|-----------:|-----:|
| Intercept                                           |   288.65 |        199.45 |        374.07 | 1.00 |
| conditVaried                                        |   -40.19 |       -104.68 |         23.13 | 0.89 |
| bandTypeExtrapolation                               |   -23.35 |        -57.28 |         10.35 | 0.92 |
| bandOrderReverse                                    |   -73.72 |       -136.69 |        -11.07 | 0.99 |
| conditVaried:bandTypeExtrapolation                  |    52.66 |         14.16 |         90.23 | 1.00 |
| conditVaried:bandOrderReverse                       |   -37.48 |       -123.28 |         49.37 | 0.80 |
| bandTypeExtrapolation:bandOrderReverse              |    80.69 |         30.01 |        130.93 | 1.00 |
| conditVaried:bandTypeExtrapolation:bandOrderReverse |    30.42 |        -21.00 |         81.65 | 0.87 |

Table 2: **Experiment 3 testing accuracy**. Main effects of condition and band type (training vs. extrapolation), and the interaction between the two factors. Larger coefficient estimates indicate larger deviations from the baselines (constant & trained bands) - and a positive interaction coefficient indicates disproporionate deviation for the varied condition on the extrapolation bands
</div>
<img src="../Assets/figs/e3_cond_effects_dist.png"
id="fig-e3-test-condEffect"
alt="Figure 2: E3. A) Deviations from target band during testing without feedback stage. B) Estimated marginal means for the interaction between training condition and band type. Error bars represent 95% confidence intervals." />

<img src="../Assets/figs/e3_test-dev.png" id="fig-e3-test-dev"
alt="Figure 3: E3. A) Deviations from target band during testing without feedback stage. B) Estimated marginal means for the interaction between training condition and band type. Error bars represent 95% confidence intervals." />

### Discimination

<div id="tbl-e3-bmm-vx">

Table 3: Experiment 3. Bayesian Mixed Model Predicting Vx as a function of condition (Constant vs. Varied) and Velocity Band
</div>
<div id="tbl-e3-bmm-vx">

| Term                                    | Estimate | 95% CrI Lower | 95% CrI Upper |   pd |
|:--------------------------------|--------:|------------:|------------:|-----:|
| b_Intercept                             |   601.83 |        504.75 |        699.42 | 1.00 |
| b_conditVaried                          |    12.18 |       -134.94 |        162.78 | 0.56 |
| b_bandOrderReverse                      |    13.03 |       -123.89 |        144.67 | 0.58 |
| Band                                    |     0.49 |          0.36 |          0.62 | 1.00 |
| b_conditVaried:bandOrderReverse         |  -338.15 |       -541.44 |       -132.58 | 1.00 |
| b_conditVaried:bandInt                  |    -0.04 |         -0.23 |          0.15 | 0.67 |
| b_bandOrderReverse:bandInt              |    -0.10 |         -0.27 |          0.08 | 0.86 |
| b_conditVaried:bandOrderReverse:bandInt |     0.42 |          0.17 |          0.70 | 1.00 |

Table 4: Experiment 3. Bayesian Mixed Model Predicting Vx as a function of condition (Constant vs. Varied) and Velocity Band
</div>



<img src="../Assets/figs/e3_test-vx.png" id="fig-e3-test-vx"
alt="Figure 4: Experiment 3. Conditional effect of training condition and Band. Ribbons indicate 95% HDI. The steepness of the lines serves as an indicator of how well participants discriminated between velocity bands." />
