
pacman::p_load(dplyr,purrr,tidyr,ggplot2, data.table, here, patchwork, conflicted, stringr)
conflict_prefer_all("dplyr", quiet = TRUE)

ft <- list.files(here::here("data/abc_tabs"),full.names = TRUE)

pt <- map(ft,readRDS)
pt2 <- imap_dfr(pt, ~.x[["et_sum"]]) |> arrange(Avg_EXAM_error) |> mutate(full_name = factor(full_name, levels = unique(full_name)))
it <- imap_dfr(pt, ~.x[["et_sum_x_indv"]])

pt2 |> filter(Fit_Method=="Test & Train") |> arrange(Avg_EXAM_error)
pt2 |> filter(Fit_Method=="Test & Train", condit=="Varied") |> arrange(Avg_EXAM_error) |> select(-Best_Model,-Fit_Method)
pt2 |> filter(Fit_Method=="Test & Train", condit=="Constant") |> arrange(Avg_EXAM_error) |> select(-Best_Model,-Fit_Method)
pt2 |> filter(Fit_Method=="Test & Train", condit=="Varied") |> arrange(Avg_ALM_error) |> select(-Best_Model,-Fit_Method)
pt2 |> filter(Fit_Method=="Test & Train", condit=="Constant") |> arrange(Avg_ALM_error) |> select(-Best_Model,-Fit_Method)

pt2 |> arrange(Avg_EXAM_error) |> ggplot(aes(x=full_name,y=Avg_EXAM_error,fill=condit)) + geom_col() + coord_flip() 

# use fct_reorder to reorder the x axis
ind_rank <- it |> group_by(id,Fit_Method,condit,full_name) |> summarise(ALM=mean(ALM_error), EXAM=mean(EXAM_error)) |> arrange(id,EXAM) |>
  group_by(id, Fit_Method, condit) |> mutate(rank=rank(EXAM))


it |> filter(id==1, Fit_Method=="Test & Train",x==1200) |> arrange(EXAM_error) 

it |> filter(id==1, Fit_Method=="Test") |> group_by(id,Fit_Method,condit, full_name) |> 
  summarise(ALM=mean(ALM_error), EXAM=mean(EXAM_error)) |> arrange(EXAM)

ind_rank |> filter(rank==1)


library(DT)


dt1 <- pt2 |> 
  mutate(Fit_Method = factor(case_when(
    Fit_Method == "Test Only" ~ "Test",
    Fit_Method == "Train Only" ~ "Train",
    Fit_Method == "Test & Train" ~ "Test_Train", 
    TRUE ~ Fit_Method
  ))) |>
  rename("ALM_Error" = Avg_ALM_error, "EXAM_Error" = Avg_EXAM_error,
         "Best_ALM" = N_Best_ALM, "Best_EXAM"=N_Best_EXAM,
         "ng"=ng_value,"ns"=n_samp,"buf"=buf_value,"type"=run_type) |>
datatable(
  extensions = 'Buttons',
  options = list(
    dom = 'Blfrtip',
    buttons = list('copy', 'csv', 'excel', 'pdf', 'print'),
    pageLength = 100,
    autoWidth = FALSE
  ),
  filter = 'top' # Add filtering at the top of each column
)

htmlwidgets::saveWidget(dt1, "abc_tab_compare.html")

