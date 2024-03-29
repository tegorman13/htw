


```{r setup2, include=FALSE}
pacman::p_load(dplyr, knitr,kableExtra)
```




```{r}
#| label: tbl-bmm-structure
#| tbl-cap: Mixed model structure and coefficient descriptions

# Create the data frame for the table
table_data <- data.frame(
  Type = c(
    rep("Population-Level Effects", 4),
    rep("Group-Level Effects", 2),
    "Family Specific Parameters"
  ),
  Parameter = c(
    "\\(\\beta_0\\)", "\\(\\beta_1\\)", "\\(\\beta_2\\)", "\\(\\beta_3\\)",
    "\\(\\sigma_{\\text{Intercept}}\\)", "\\(\\sigma_{\\text{bandInt}}\\)", "\\(\\sigma_{\\text{Observation}}\\)"
  ),
  Term = c(
    "(Intercept)", "conditVaried", "bandInt", "conditVaried:bandInt",
    "sd__(Intercept)", "sd__bandInt", "sd__Observation"
  ),
  Description = c(
    "Intercept representing the baseline deviation", "Effect of condition (Varied vs. Constant) on deviation", "Effect of target velocity band (bandInt) on deviation", "Interaction effect between training condition and target velocity band on deviation",
    "Standard deviation for (Intercept)", "Standard deviation for bandInt", "Standard deviation for Gaussian Family"
  )
) |>   mutate(
    Term = glue::glue("<code>{Term}</code>")
  ) 

# Create the table
kable_out <- table_data %>%
  kbl(format = 'html', escape = FALSE, booktabs = TRUE, 
      #caption = '<span style = "color:black;"><center><strong>Table 1: General Model Structure Information</strong></center></span>',
      col.names = c("Type", "Parameter", "Term", "Description")) %>%
  kable_styling(position="left", bootstrap_options = c("hover"), full_width = FALSE) %>%
  column_spec(1, bold = FALSE, border_right = TRUE) %>%
  column_spec(2, width = '4cm') %>%
  column_spec(3, width = '4cm') %>%
  row_spec(c(4, 7), extra_css = "border-bottom: 2px solid black;") %>%
  pack_rows("", 1, 4, bold = FALSE, italic = TRUE) %>%
  pack_rows("", 5, 6, bold = FALSE, italic = TRUE) %>%
  pack_rows("", 7, 7, bold = FALSE, italic = TRUE)

kable_out

```
