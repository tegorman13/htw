---
title: Experiment 1 Analysis
date: last-modified
categories: [Analysis, R]
page-layout: article
fig-width: 15
fig-height: 8
code-fold: true
code-tools: true
execute: 
  warning: false
  eval: true
---

| Variable Name | Variable Levels                 | Description                                                                               |
|---------------|---------------------------------|-------------------------------------------------------------------------------------------|
| condit        | Constant, Varied                | Condition of the experiment: constant or varied                                             |
| tOrder        | Test First, Train First         | Order of testing and training stages: test first or train first                            |
| expMode       | train, train-Nf, test-Nf, etc.  | Mode of the experiment: train, train-Nf, test-Nf, etc.                                    |
| trainStage    | Beginning, Middle, End, Test    | Stage of the training: beginning, middle, end, or test                                    |
| expStage      | TrainStart, intTest1, etc.      | Stage of the experiment: TrainStart, intTest1, TrainMid1, etc.                            |
| band          | 1, 2, 3, 4, 5, 6               | Band number                                                                                 |
| vb            | 100-300, 350-550, etc.         | Velocity band range                                                                         |
| lowBound      | 100, 350, 600, etc.            | Lower bound of the velocity band range                                                      |
| feedback      | 0, 1                           | Feedback type: 0 (no feedback), 1 (feedback)                                               |
| stage         | 1, 2, 3, etc.                  | Stage number of the experiment                                                              |

---


```{r}
pacman::p_load(tidyverse,lme4,emmeans,here,knitr,kableExtra,gt)
e1 <- readRDS(here("data/e1_07_04_23.rds"))
source(here("Functions/Display_Functions.R"))
```



### Methods

Participants
A total of `r length(unique(e1$id))` participants (XXX% female, XXX% male) were recruited from the Indiana University Introductory Psychology Course. The average age of participants was XXX years (SD = XXX). Participants were randomly assigned to one of two training conditions: varied training or constant training. 

Design
The experiment employed a 2 (Training Condition: varied vs. constant) x 2 (Order Manipulation: orig vs. rev) x 2 (Feedback Type: continuous vs. ordinal) between-subjects design. 

Procedure
Upon arrival at the laboratory, participants were provided with a description of the experiment and signed informed consent forms. They were then seated in front of a computer equipped with a mouse and were given instructions on how to perform the "Hit The Wall" (HTW) visuomotor extrapolation task.

The HTW task involved launching projectiles to hit a target displayed on the computer screen. Participants completed a total of 90 trials during the training stage. In the varied training condition, participants encountered three velocity bands (800-1000, 1000-1200, and 1200-1400). In contrast, participants in the constant training condition encountered only one velocity band (800-1000).

During the training stage, participants in both conditions also completed "no feedback" trials, where they received no information about their performance. These trials were randomly interleaved with the regular training trials.

Following the training stage, participants proceeded to the testing stage, which consisted of three phases. In the first phase, participants completed "no-feedback" testing from three novel extrapolation bands (100-300, 350-550, and 600-800), with each band consisting of 15 trials.

In the second phase of testing, participants completed "no-feedback" testing from the three velocity bands used during the training stage (800-1000, 1000-1200, and 1200-1400). In the constant training condition, two of these bands were novel, while in the varied training condition, all three bands were encountered during training.

The third and final phase of testing involved "feedback" testing for each of the three extrapolation bands (100-300, 350-550, and 600-800), with each band consisting of 10 trials. Participants received feedback on their performance during this phase.

Throughout the experiment, participants' performance was measured by calculating the distance between the produced x-velocity of the projectiles and the closest edge of the current velocity band. Lower distances indicated better performance.

After completing the experiment, participants were debriefed and provided with an opportunity to ask questions about the study.



```{r design, echo=FALSE, fig.cap="Experiment Procedure"}
DiagrammeR::grViz("digraph {
graph [layout = dot, rankdir = LR]
# define the global styles of the nodes. We can override these in box if we wish
node [shape = rectangle, style = filled, fillcolor = Linen]
data1 [label = <<b>Varied Training</b>  
<br ALIGN = 'LEFT' /> 800-1000
<br/> 1000-1200
<br/>1200-1400>]
data2 [label = <<b>Constant Training</b>  
<br/> 800-1000 >]
Test1 [label = <<b>Testing - No Feedback</b>  
<br ALIGN = 'LEFT' /> 100-300
<br/> 350-550
<br/>600-800>]
Test2 [label = <<b>Test From Train</b>  
<br ALIGN = 'LEFT' /> 800-1000
<br/> 1000-1200
<br/>1200-1400>]
Test3[label = <<b>Testing -Feedback</b>  
<br ALIGN = 'LEFT' /> 100-300
<br/> 350-550
<br/>600-800>]

# edge definitions with the node IDs
{data1 data2}  -> Test1 -> Test2 -> Test3
}
",width=700,height=250)

```

```{r}



```




```{r}

```




```{r}
#| column: page-inset-right
#| fig-cap: "Mean Vx over training blocks"
nb=5
vp1=e1 |> filter(expMode=="train") |> learn_curve_plot(gt.train,vx,condit,facet_var=vb,groupVec=c(gt.train,condit,tOrder,id,vb),nbins=nb) 
vp1
# vpt1=plotWithTable(vp1,vt1,arrange="V") 
```



```{r}
#| column: page-inset-right
#| tbl-cap: "Mean Vx over blocks. Mean (Standard Error)"
#| tbl-cap-location: top

vt1=e1 |> filter(expMode=="train") |> 
  learn_curve_table(gt.train,vx,gw=Trial_Bin,groupVec=c(vb,condit),nbins=nb,prefix="Block_") |>
  rename("Band"=vb,"Group"=condit) 
vt1 %>% gt()

```

```{r}
#| column: page-inset-right
#| fig-cap: "Absolute Deviation from Band over training"
vp2=e1 |> filter(expMode=="train") |> learn_curve_plot(gt.train,dist,condit,facet_var=vb,groupVec=c(gt.train,condit,tOrder,id,vb),nbins=nb) 
vp2

```

```{r}
#| column: page-inset-right
#| tbl-cap: "Mean Vx over blocks. Mean (Standard Error)"
#| tbl-cap-location: bottom

vt2=e1 |> filter(expMode=="train") |> 
  learn_curve_table(gt.train,dist,gw=Trial_Bin,groupVec=c(vb,condit),nbins=nb,prefix="Block_") |>
  rename("Band"=vb,"Group"=condit) 
vt2 %>% gt()

```



## Testing
```{r}



```

