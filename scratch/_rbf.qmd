---
title: "rbf.qmd"
jupyter: julia-1.9
---

```{r}
pacman::p_load(tidyverse,RSNNS,tikzDevice,knitr)
```






```{r fig.width=11, fig.height=9}
#| eval: false 
n_iter=500
n_nodes <- c(c(1,3,6,9),seq(10,300,10))
train <- tibble(input=as.matrix(seq(0,100,1)),output=as.matrix(round(2.2*input + 30,0)))

fits <- map(n_nodes,~rbf(train$input,train$output,size=.x,maxit=n_iter))



crossing(n_nodes,train) %>% 
  cbind(.,fit=unlist(map(fits,fitted))) %T>% 
  {print(ggplot(d,aes(x=input,y=output))+geom_line()+
  geom_line(aes(x=input,y=fit),col="green") + facet_wrap(~n_nodes,scales="free")) } %>% 
  {. ->> d}
  
tibble(n_nodes=as.factor(rep(n_nodes,each=n_iter))) %>% 
  cbind(.,error=unlist(purrr::map(fits, ~.x$IterativeFitError))) %>%
  mutate(it=seq(1,n()),.by=n_nodes) %T>% 
  {print(ggplot(.,aes(it,error,col=n_nodes))+geom_line()+facet_wrap(~n_nodes,scales="free")) } %>% 
  {. ->> errVec}


preds <- map(fits,~predict(.x,as.matrix(seq(0,150,1))))

tibble(n_nodes=n_nodes,fit=map(fits,~predict(.x,as.matrix(seq(0,150,1)))))

tibble(n_nodes=rep(n_nodes,each=150),input=as.matrix(seq(0,150,1)),output=as.matrix(round(2.2*input + 30,0)))

```





```{r}
#| echo: false
#| eval: false 

# d <- crossing(n_nodes,train) %>% 
#   cbind(.,fit=unlist(map(fits,fitted)))

# purrr::map(fits, ~.x$IterativeFitError)
# unlist(purrr::map(fits, ~.x$IterativeFitError))

inputs <- as.matrix(seq(0,100,1))
outputs <- as.matrix(round(2.2*inputs + 30,0))

model <- rbf(inputs, outputs, size=10, maxit=5000,
                     initFuncParams=c(0, 1, 0, 0.01, 0.01),
                     learnFuncParams=c(1e-8, 0, 1e-8, 0.1, 0.8), linOut=TRUE)

ggplot(d,aes(x=input,y=output))+geom_line()+
  geom_line(aes(x=input,y=fit),col="green") + facet_wrap(~n_nodes,scales="free")


inputs <- as.matrix(seq(0,10,0.1))
outputs <- as.matrix(sin(inputs) + runif(inputs*0.2))
outputs <- normalizeData(outputs, "0_1")

model <- rbf(inputs, outputs, size=40, maxit=1000,
                     initFuncParams=c(0, 1, 0, 0.01, 0.01),
                     learnFuncParams=c(1e-8, 0, 1e-8, 0.1, 0.8), linOut=TRUE)

par(mfrow=c(2,1))
plotIterativeError(model)
plot(inputs, outputs)
lines(inputs, fitted(model), col="green")
```

```{r}




```

