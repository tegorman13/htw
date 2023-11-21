
```{r setup1, include=FALSE}
pacman::p_load(knitr,kableExtra)
```


```{r ALM_Table1}
#| results: asis

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
