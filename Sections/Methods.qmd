---
title: Methods
date: last-modified
categories: []
code-fold: true
code-tools: true
execute: 
  warning: false
  eval: true
---


```{r}
source(here::here("Functions", "packages.R"))
e1 <- readRDS(here("data/e1_08-04-23.rds")) 
e2 <- readRDS(here("data/e2_08-04-23.rds")) 
e3 <- readRDS(here("data/e3_08-04-23.rds")) 

```



## Methods

Participants
A total of `r length(unique(e1$id))` participants (XXX% female, XXX% male) were recruited from the Indiana University Introductory Psychology Course. The average age of participants was XXX years (SD = XXX). Participants were randomly assigned to one of two training conditions: varied training or constant training. 

Design
The experiment employed a 2 (Training Condition: varied vs. constant). 

Procedure
Upon arrival at the laboratory, participants were provided with a description of the experiment and signed informed consent forms. They were then seated in front of a computer equipped with a mouse and were given instructions on how to perform the "Hit The Wall" (HTW) visuomotor extrapolation task.

The HTW task involved launching projectiles to hit a target displayed on the computer screen. Participants completed a total of 90 trials during the training stage. In the varied training condition, participants encountered three velocity bands (800-1000, 1000-1200, and 1200-1400). In contrast, participants in the constant training condition encountered only one velocity band (800-1000).

During the training stage, participants in both conditions also completed "no feedback" trials, where they received no information about their performance. These trials were randomly interleaved with the regular training trials.

Following the training stage, participants proceeded to the testing stage, which consisted of three phases. In the first phase, participants completed "no-feedback" testing from three novel extrapolation bands (100-300, 350-550, and 600-800), with each band consisting of 15 trials.

In the second phase of testing, participants completed "no-feedback" testing from the three velocity bands used during the training stage (800-1000, 1000-1200, and 1200-1400). In the constant training condition, two of these bands were novel, while in the varied training condition, all three bands were encountered during training.

The third and final phase of testing involved "feedback" testing for each of the three extrapolation bands (100-300, 350-550, and 600-800), with each band consisting of 10 trials. Participants received feedback on their performance during this phase.

Throughout the experiment, participants' performance was measured by calculating the distance between the produced x-velocity of the projectiles and the closest edge of the current velocity band. Lower distances indicated better performance.

After completing the experiment, participants were debriefed and provided with an opportunity to ask questions about the study.


```{dot}
//| label: fig-design-e1
//| fig-cap: "Experiment 1 Design. Constant and Varied participants complete different training conditions."
//| fig-width: 8.0
//| fig-height: 2.5
//| fig-responsive: false
//| column: screen-inset-right
digraph {
  graph [layout = dot, rankdir = LR]

  // define the global styles of the nodes
  node [shape = rectangle, style = filled]

  data1 [label = " Varied Training \n800-1000\n1000-1200\n1200-1400", fillcolor = "#FF0000"]
  data2 [label = " Constant Training \n800-1000", fillcolor = "#00A08A"]
  Test3 [label = "    Final Test \n  Novel With Feedback  \n100-300\n350-550\n600-800", fillcolor = "#ECCBAE"]

  // edge definitions with the node IDs
  data1 -> Test1
  data2 -> Test1

  subgraph cluster {
    label = "Test Phase \n(Counterbalanced Order)"
    Test1 [label = "Test  \nNovel Bands \n100-300\n350-550\n600-800", fillcolor = "#ECCBAE"]
    Test2 [label = "  Test \n  Varied Training Bands  \n800-1000\n1000-1200\n1200-1400", fillcolor = "#ECCBAE"]
    Test1 -> Test2
  }

  Test2 -> Test3
}
#https://renenyffenegger.ch/notes/tools/Graphviz/examples/index
#https://edotor.net/
```







## HTW Task
* need to create a demo version without consent form. And maybe separate windows for the different versions. 

Experimental Task for the HTW Project. Programmed in Javascript, and making use of phaser.js.


```{r out.width = "70%", echo=FALSE, message=FALSE, warning=FALSE, fig.align='center'}
knitr::include_graphics(here::here("images/2023-04-19-16-20-55.png"))
```

<!-- ![](images/2023-04-19-16-20-55.png) -->
<!-- `r knitr::include_graphics(here::here("images/2023-04-19-16-20-55.png"))` -->


## Live Task Demo

Check the box at the bottom of the consent form, and then click Start Experiment. 

![HTW_Task](https://pcl.sitehost.iu.edu/tg/HTW/HTW_Index_Of.html){width="900" height="800"}


