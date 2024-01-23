
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



ind_rank <- it |> group_by(id,Fit_Method,condit,full_name) |> summarise(ALM=mean(ALM_error), EXAM=mean(EXAM_error)) |> arrange(id,EXAM) |>
  group_by(id, Fit_Method, condit) |> mutate(rank=rank(EXAM))


it |> filter(id==1, Fit_Method=="Test & Train",x==100) |> arrange(EXAM_error) 

ind_rank |> filter(rank==1)
