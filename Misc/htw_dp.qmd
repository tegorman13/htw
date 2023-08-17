---
title: "Project 2"
subtitle: " Variability and Function Learning "
bibliography: "../Assets/Bib/HTW.bib"
date: last-modified
keep-tex: true
code-fold: true
---

```{r setup, include = FALSE}
#library(figpatch)
library(kableExtra)
library(knitr)
library(here)

```


# Introduction

In project 1, we applied model-based techniques to quantify and control for the similarity between training and testing experience, which in turn enabled us to account for the difference between varied and constant training via an extended version of a similarity based generalization model. In project 2, we will go a step further, implementing a full process model capable of both 1) producing novel responses and 2) modeling behavior in both the learning and testing stages of the experiment. Project 2 also places a greater emphasis on extrapolation performance following training - as varied training has often been purported to be particularly beneficial in such situations. Extrapolation has long been a focus of the literature on function learning [@brehmerHypothesesRelationsScaled1974; @carrollFunctionalLearningLearning1963]. Central questions of the function learning literature have included the relative difficulties of learning various functional forms (e.g. linear vs.bilinear vs. quadratic), and the relative effectiveness of rule-based vs. association-based exemplar models vs. various hybrid models [@bottNonmonotonicExtrapolationFunction2004; @deloshExtrapolationSineQua1997; @jonesActiveFunctionLearning2018; @kalishPopulationLinearExperts2004; @mcdanielConceptualBasisFunction2005; @mcdanielPredictingTransferPerformance2009]. However the issue of training variation has received surprisingly little attention in this area.

# Methods

## Participants

Data was collected from 647 participants (after exclusions). The results shown below consider data from subjects in our initial experiment, which consisted of 196 participants (106 constant, 90 varied). The follow-up experiments entailed minor manipulations: 1) reversing the velocity bands that were trained on vs. novel during testing; 2) providing ordinal rather than numerical feedback during training (e.g. correct, too low, too high). The data from these subsequent experiments are largely consistently with our initial results shown below.

## Task

We developed a novel visuomotor extrapolation task, termed the ["Hit The Wall" (HTW]{.underline}) task, wherein participants learned to launch a projectile such that it hit a rectangle at the far end of the screen with an appropriate amount of force. Although the projectile had both x and y velocity components, only the x-dimension was relevant for the task.  [Link to task demo](https://pcl.sitehost.iu.edu/tg/HTW/HTW_Index.html?sonaid=){target="_blank"}

## Design

1)  90 training trials split evenly divided between velocity bands. Varied training with 3 velocity bands and Constant training with 1 band.

2)  No-feedback testing from 3 novel extrapolation bands. 15 trials each.  

3)  No-feedbacd testing from the 3 bands used during the training phase (2 of which were novel for the constant group). 9 trials each.

4)  Feedback testing for each of the 3 extrapolation bands. 10 trials each.

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

<!-- {.tabset} -->

# Results

## Training

Training performance is shown in Results Figure 2A. All groups show improvement from each of their training velocity-bands (i.e. decreasing average distance from target). In the velocity band trained at by both groups (800-1000), the constant group maintains a superior level of performance from the early through the final stages of training. This difference is unsurprising given that the constant group had 3x more practice trials from that band.

```{r Training, out.height="100%", out.width="100%"}
#| column: body-outset

#fig.cap="\\label{fig:figs}training performance"

#title=paste0("HTW Training")
#figpatch::fig(here("Assets/Training-1.png"))
knitr::include_graphics(here("Assets/Training-1.png"))
  #fig_lab(.,title,pos="top",fontface='bold',size=12,hjust=.01)

```

## Testing

For evaluating testing performance, we consider 3 separate metrics. 1) The average absolute deviation from the correct velocity, 2) The % of throws in which the wall was hit with the correct velocity and 3) The average x velocity produced.

Results Figure 2B shows the average velocity produced for all 6 bands that were tested. At least at the aggregate level, both conditions were able to differentiate all 6 bands in the correct order, despite only having received training feedback for 1/6 (constant) or 3/6 (varied) bands during training. Participants in both groups also had a bias towards greatly overestimating the correct velocity for band 100-300, for which both groups had an average of greater than 500.

```{r Testing Vx, out.height="100%", out.width="100%" }
#| column: page-inset
# fig2aCap=str_wrap("Figure 2B: Bands 100-300, 350-550 and 600-800 are novel extrapolations for both groups. Band 800-1000 was a training band for both groups. Bands 1000-1200, and 1200-1400 were trained for the varied group, and novel for the constant group.  Top figure displays mean deviation from correct velocity. Bottom figure displays the average % of trials where participants hit the wall with the correct velocity. Error bars indicate standard error of the mean. " ,width=170)


#figpatch::fig(here("Assets/Testing_Vx-1.png"))
knitr::include_graphics(here("Assets/Testing_Vx-1.png"))

```

As is reflected in Results Figure 2C, the constant group performed significantly better than the varied group at the 3 testing bands of greatest interest. Both groups tended to perform worse for testing bands further away from their training conditions. The varied group had a slight advantage for bands 1000-1200 and 1200-1400, which were repeats from training for the varied participants, but novel for the constant participants.

```{r Test Performance, out.height="100%", out.width="100%"}
#| column: page-inset

#figpatch::fig(here("Assets/Test_Performance-1.png"))
knitr::include_graphics(here("Assets/Test_Performance-1.png"))

# 
# gtitle="2C. Testing Performance"
# title = ggdraw()+draw_label(gtitle,fontface = 'bold',x=0,hjust=0)+theme(plot.margin = margin(0, 0, 0, 7))
# captionText=str_wrap("Figure 2C: Bands 100-300, 350-550 and 600-800 are novel extrapolations for both groups. Band 800-1000 was a training band for both groups. Bands 1000-1200, and 1200-1400 were trained for the varied group, and novel for the constant group.  Right side figure displays mean deviation from correct velocity band (lower values correspond to better performance). Bottom Left displays the average % of trials where participants hit the wall with the correct velocity (higher values correspond got better performance). Error bars indicate standard error of the mean. ",150)
# capt=ggdraw()+draw_label(captionText,fontface = 'italic',x=0,hjust=0,size=11)+theme(plot.margin = margin(0, 0, 0, 1))
# 
# plot_grid(title,NULL,leg,NULL,gbDev,gbHit,capt,NULL,ncol=2,rel_heights=c(.1,.1,1,.1),rel_widths=c(1,1))



```

# Modeling

In project 1, we applied model-based techniques to quantify and control for the similarity between training and testing experience, which in turn enabled us to account for the difference between varied and constant training via an extended version of a similarity based generalization model. In project 2, we will go a step further, implementing a full process model capable of both 1) producing novel responses and 2) modeling behavior in both the learning and testing stages of the experiment. For this purpose, we will apply the associative learning model (ALM) and the EXAM model of function learning (DeLosh 1997). ALM is a simple connectionist learning model which closely resembles Kruschke's ALCOVE model (Kruscke 1992), with modifications to allow for the generation of continuous responses.

## ALM & Exam Description

Delosh et al. (1997) introduced the associative learning model (ALM), a connectionist model within the popular class of radial-basis networks. ALM was inspired by, and closely resembles Kruschke's influential ALCOVE model of categorization (Kruscke 1992).

ALM is a localist neural network model, with each input node corresponding to a particular stimulus, and each output node corresponding to a particular response value. The units in the input layer activate as a function of their Gaussian similarity to the input stimulus. So, for example, an input stimulus of value 55 would induce maximal activation of the input unit tuned to 55. Depending on thevalue of the generalization parameter, the nearby units (e.g. 54 and 56; 53 and 57) may also activate to some degree. ALM is structured with input and output nodes that correspond to regions of the stimulus space, and response space, respectively. The units in the input layer activate as a function of their similarity to a presented stimulus. As was the case with the exemplar-based models, similarity in ALM is exponentially decaying function of distance. The input layer is fully connected to the output layer, and the activation for any particular output node is simply the weighted sum of the connection weights between that node and the input activations. The network then produces a response by taking the weighted average of the output units (recall that each output unit has a value corresponding to a particular response). During training, the network receives feedback which activates each output unit as a function of its distance from the ideal level of activation necessary to produce the correct response. The connection weights between input and output units are then updated via the standard delta learning rule, where the magnitude of weight changes are controlled by a learning rate parameter.

See Table 2A for a full specification of the equations that define ALM and EXAM.

## Model Equations

::: column-page-inset-right
```{r ALM_Table, results='asis'}

text_tbl <- data.frame(
    'Step'=c("Input Activation","Output Activation","Output Probability","Mean Output","Feedback Activation","Update Weights","Extrapolation",""),
    'Equation' = c("$a_i(X) = \\frac{e^{-c \\cdot (X-X_i)^2}}{ \\sum_{k=1}^Me^{-c \\cdot (X-X_i)^2}}$", 
                   '$O_j(X) = \\sum_{k=1}^Mw_{ji} \\cdot a_i(X)$',
                   '$P[Y_j | X] = \\frac{O_i(X)}{\\sum_{k=1}^Mo_k(X)}$',
                   "$m(x) = \\sum_{j=1}^LY_j \\cdot \\bigg[\\frac{O_j(X)}{\\sum_{k=1}^Lo_k(X)}\\bigg]$",
                   "$f_j(Z)=e^{-c\\cdot(Z-Y_j)^2}$",
                   "$w_{ji}(t+1)=w_{ji}(t)+\\alpha \\cdot {f_i(Z(t))-O_j(X(t))} \\cdot a_i(X(t))$",
                   "$P[X_i|X] = \\frac{a_i(X)}{\\sum_{k=1}^Ma_k(X)}$",
                   "$E[Y|X_i]=m(X_i) + \\bigg[\\frac{m(X_{i+1})-m(X_{i-1})}{X_{i+1} - X_{i-1}} \\bigg] \\cdot[X-X_i]$"),
    
    'Description'= c(
            "Activation of each input node, $X_i$, is a function of the Gaussian similarity between the node value and stimulus X. ",
            "Activation of each Output unit $O_j$ is the weighted sum of the input activations and association weights",
            "Each output node has associated response, $Y_j$. The probability of response $Y_j$ is determined by the ratio of output activations",
            "The response to stimulus x is the weighted average of the response probabilities",
            "After responding, feedback signal Z is presented, activating each output node via the Gaussian similarity to the ideal response  ",
            "Delta rule to update weights. Magnitude of weight changes controlled by learning rate parameter alpha.",
            "Novel test stimulus X activates input nodes associated with trained stimuli",
            "Slope value computed from nearest training instances and then added to the response associated with the nearest training instance,m(x)")
)
text_tbl$Step=cell_spec(text_tbl$Step,font_size=12)
text_tbl$Equation=cell_spec(text_tbl$Equation,font_size=18)
almTable=kable(text_tbl, 'html', 
  booktabs=T, escape = F, align='l',
  caption = '<span style = "color:black;"><center><strong>Table 1: ALM & EXAM Equations</strong></center></span>',
  col.names=c("","Equation","Description")) %>%
  kable_styling(position="left",bootstrap_options = c("hover")) %>%
  column_spec(1, bold = F,border_right=T) %>%
  column_spec(2, width = '10cm')%>%
  column_spec(3, width = '15cm') %>%
  pack_rows("ALM Activation & Response",1,4,bold=FALSE,italic=TRUE) %>%
  pack_rows("ALM Learning",5,6,bold=FALSE,italic=TRUE) %>%
  pack_rows("EXAM",7,8,bold=FALSE,italic=TRUE)
  #save_kable(file="almTable.html",self_contained=T)
#almTable


cat(almTable)
```
:::

## Model Fitting and Comparison

Following the procedure used by McDaniel & Busemeyer (2009), we will assess the ability of both ALM and EXAM to account for the empirical data when fitting the models to 1) only the training data, and 2) both training and testing data. Models will be fit directly to the trial by trial data of each individual participants, both by minimizing the root-mean squared deviation (RMSE), and by maximizing log likelihood. Because ALM has been shown to do poorly at accounting for human patterns extrapolation [@deloshExtrapolationSineQua1997], we will also fit the extended EXAM version of the model, which operates identically to ALM during training, but includes a linear extrapolation mechanism for generating novel responses during testing.

\
\
\
\
\
\
\

# References

::: {#refs}
:::

```{=tex}
\setlength{\parindent}{-0.5in}
\setlength{\leftskip}{0.5in}
```
<br><br>

<br><br>