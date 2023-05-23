---
title: auto htw
---

```{r}
pacman::p_load(tidyverse,patchwork)
purrr::walk(c("Functions/alm_functions.R","Functions/Display_Functions.R","Functions/Noisy_Functions.R"),source)
```

```{r}
sigmoid_activation <- function(net_input) {
  1 / (1 + exp(-net_input))
}

stimVec=c(100,600,800,1000)

nh=6; 
wi=matrix(runif(nh,-2.0,2.0),nrow=nh)
wh=matrix(runif(nh,-2.0,2.0),ncol=nh)
lr=.010
trainNet <- function(stimVec){
  for (s in stimVec){
    ha <- sigmoid_activation(wi %*% s)
    resp <- sigmoid_activation(t(ha) %*% wi) 
    
    oe=(resp-s)
    he= (t(wh) %*% oe) * (ha * (1-ha))
    
    hd <- lr*(ha %*% t(oe))
    id <- lr*(s %*% t(he))
    
    new_h <- wh - t(hd)
    new_i <- wi - t(id)
    
    wh=new_h; wi=new_i; 
    print(paste0("s=",s,". r=",resp,". oe=",oe))
  }
}

trainNet(rep(stimVec,35))

trainNet(rep(scaled_stimVec,35))

```



```{r}

# Linear activation function
linear_activation <- function(net_input) {
  net_input
}

trainNet_linear <- function(stimVec) {
  wi=matrix(runif(nh*length(stimVec),-.20,.20),nrow=nh)
  wh=matrix(runif(nh*num_bits,-.20,.20),ncol=nh)
  for (s in stimVec) {
    ina <- exp(-c*(stimVec - s)^2)
    ha <- linear_activation(wi %*% ina)
    resp <- linear_activation(t(ha) %*% wi)
    
    oe = (resp - ina)
    he = t(wh) %*% oe
    
    hd <- lr * (ha %*% (oe))
    id <- lr * ((he) %*% (ina))
    
    new_h <- wh - t(hd)
    new_i <- wi - t(id)
    
    wh = new_h
    wi = new_i
    
    print(paste0("s=", s, ". r=", resp, ". oe=", oe))
  }
}
trainNet_linear(rep(stimVec, 35))
```



```{r}

# Gaussian activation function
gaussian_activation <- function(net_input) {
  exp(-net_input^2)
}

trainNet_gaussian <- function(stimVec) {
  wi=matrix(runif(nh*length(stimVec),-.20,.20),nrow=nh)
  wh=matrix(runif(nh*num_bits,-.20,.20),ncol=nh)
  for (s in stimVec) {
    
    ha <- gaussian_activation(wi %*% s)
    resp <- gaussian_activation(t(ha) %*% wi)
    
    oe = (resp - s)
    he = (t(wh) %*% oe) * (-2 * ha * (wi %*% s))
    
    hd <- lr * (ha %*% t(oe))
    id <- lr * (s %*% t(he))
    
    new_h <- wh - t(hd)
    new_i <- wi - t(id)
    
    wh = new_h
    wi = new_i
    
    print(paste0("s=", s, ". r=", resp, ". oe=", oe))
  }
}

trainNet_gaussian(rep(stimVec, 35))

```













```{r}
# Convert the input to binary representation
input_to_binary <- function(input, num_bits = 10) {
  intToBits(input)[1:num_bits]
}

# Train the network with binary input representation
trainNet_binary <- function(stimVec,num_bits = 10) {
  wi=matrix(runif(nh*num_bits,-.20,.20),nrow=nh)
  wh=matrix(runif(nh*num_bits,-.20,.20),ncol=nh)
  for (s in stimVec) {
    binary_input <- as.matrix(as.integer(input_to_binary(s,num_bits)), nrow = 1)

    ha <- sigmoid_activation(wi %*% binary_input)
    resp <- (t(ha) %*% wi)
    
    oe = (resp - t(binary_input))
    he = (t(wh) %*% t(oe)) * (ha * (1 - ha))
    
    hd <- lr * (ha %*% t(oe))
    id <- lr * (binary_input %*% t(he))
    
    new_h <- wh - t(hd)
    new_i <- wi - t(id)
    
    wh = new_h
    wi = new_i
    
    print(paste0("s=", s, ". r=", resp, ". oe=", oe))
  }
}

# Train the network with the binary input representation
trainNet_binary(rep(stimVec, 35))







```

