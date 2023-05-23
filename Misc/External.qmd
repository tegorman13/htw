---
title: "External"
date: last-modified
---

## Quarto





::: {.callout-note .cBook}
## [External Links](#extSite)


[Regression in Keras](https://www.datatechnotes.com/2019/01/regression-example-with-keras-in-r.html)



[uncertainty in deep learning](http://www.cs.ox.ac.uk/people/yarin.gal/website/blog_2248.html)
<!-- ![uncertainty in deep learning](http://www.cs.ox.ac.uk/people/yarin.gal/website/blog_2248.html){width="1000" height="600"} -->


[Robinson Codebook](https://earobinson95.github.io/log-perception-prolific/analyses/03-estimation/estimation-analysis.html)

[painbrow](https://github.com/steveharoz/painbow)

[NN Flex](https://eranraviv.com/flexible-neural-networks-really/)

[b-splines](https://www.rdatagen.net/post/generating-non-linear-data-using-b-splines/)

::: {.callout-note .cBook}
## Mahr Posts
[Polypoly](https://www.tjmahr.com/polypoly-package-released/){width="1100" height="600"}

[splines](https://www.tjmahr.com/random-effects-penalized-splines-same-thing/){width="1100" height="600"}
:::


[Hock Function Learn with NN](https://github.com/tdhock/2020-yiqi-summer-school)

[Autoencoder with Tensorflow](https://collinerickson.github.io/2018/08/22/creating-an-autoencoder-with-tensorflow-in-r/)

[Connected Scatterplot](http://steveharoz.com/research/connected_scatterplot/)

:::

::: {.callout-note .cTool collapse="true"}
## Quarto Tools to add
[ggpage](https://emilhvitfeldt.github.io/ggpage/) \
[quartoDocSetup](https://quarto.org/docs/websites/website-listings.html) \
[quartoLayout](https://quarto.org/docs/authoring/article-layout.html) \
:::






::: {.callout-note collapse="true"}
## Nutshell items


- [:link to senseless paragraph](#test)  
- [:link to wikipedia article](https://en.wikipedia.org/wiki/Nutshell)
- [:link to invisible sections](#x-invisible)


## [:HtwNN](#htw)


## [:OneNote External](#OneNoteExt)











## [:Website_Add](#webAdd)

[:ggpage](https://emilhvitfeldt.github.io/ggpage/)

[:quartoDocSetup](https://quarto.org/docs/websites/website-listings.html)

[:quartoLayout](https://quarto.org/docs/authoring/article-layout.html)

:::




## another


::: {.callout-note collapse="true"}
## PDF Views



### Main

**Information Sampling Explains Bayesian Learners' biases in correlation judgement **\
<iframe width="560" height="315" src="https://www.youtube.com/embed/quAwHWXwk30 " title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>


### Alt
:::{.column-screen-right}
::: {layout-ncol="2"}
**Neural Network Model of Continual Learning with Cognitive Control**\
<iframe width="560" height="315" src="https://www.youtube.com/embed/lPVQe_kTm28" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

**Geometry of Map-Like Representations Under Dynamic-Cognitive Control**\
<iframe width="560" height="315" src="https://www.youtube.com/embed/vMBS8ckoFl0" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
:::
:::
### Page Embed

[Video](https://www.youtube.com/embed/quAwHWXwk30 )



:::


## Callouts



::: {.callout-note appearance="simple"}

## Pay Attention

Using callouts is an effective way to highlight content that your reader give special consideration or attention.

:::


::: {.callout-note icon=false}

## Pay Attention

Using callouts is an effective way to highlight content that your reader give special consideration or attention.

:::

::: {.callout-note icon=false}
## Video Page


[Video](https://www.youtube.com/embed/quAwHWXwk30 )
:::




## Test Grid

::: {.grid}

::: {.g-col-4}
This column takes 1/3 of the page
:::

::: {.g-col-8}

::: {.callout-important}
## Model To-do

- ALM
- EXAM
- Approximate Bayes?
:::

::: {.callout-important}
## Analysis To-do
- Discrimination
- Mixed Models?
:::



:::

:::


## Alt column approach

::: {layout-ncol="2"}

::: {.callout-important collapse='true'}
## Model To-do
- Train Predict Transfer vs. Full Fit
- Separate ALM and EXAM Fits
- ALM + Prior Knowledge
- Empirical Learning Model
- Noise Parameter
- Model Recovery
- Approximate Bayes?
:::

::: {.callout-important collapse='false'}
## Analysis To-do
- Discrimination
- Mixed Models?

:::


:::






## Test

This is a senseless paragraph

## Testing Links

- [:link to senseless paragraph](#test)  
- [:link to wikipedia article](https://en.wikipedia.org/wiki/Nutshell)
- [:link to invisible sections](#x-invisible)

## end

```{r}
library(flashCard)
df1 <- data.frame(
  front = c("Title front","contentfront", "content second line"),
  back =c("Title back","content back", "second line")
)
flashCard(df1, elementId = "card", front_text_color = "white")
```

