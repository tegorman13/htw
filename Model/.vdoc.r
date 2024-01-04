#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
pacman::p_load(tidyverse, data.table,dtplyr, here, patchwork, conflicted)
conflict_prefer_all("dplyr", quiet = TRUE)

group_prior=abc_2M_p001 = readRDS(here::here("data/abc_2M_rmse_p001.rds"))
names(group_prior)
str(group_prior$abc_ev[[1]])

#tetr_v <- group_prior$abc_ev$teter_results


teter <- abc_2M_p001 |> map_dfr(~tibble(pluck(.x$teter_results))) 
te <- abc_2M_p001 |> map_dfr(~tibble(pluck(.x$te_results)))
tr <- abc_2M_p001 |> map_dfr(~tibble(pluck(.x$tr_results)))


teter |> select(sim_dat) |> unnest(sim_dat)  %>%  filter(expMode2=="Test") |>
  mutate(facet_label = paste0("rank: ", rank, "\n", "c: ", round(c, 4), "\n", "lr: ", round(lr, 4), "\n", "distance: ", round(distance, 1))) |>
  ggplot(aes(x = x, y = pred, fill=condit)) + 
  stat_summary(fun=mean, geom="bar", position=position_dodge()) +
  stat_summary(fun.data=mean_se, geom="errorbar", position=position_dodge()) +
  facet_wrap(~Model)

teter |> select(sim_dat) |> unnest(sim_dat)  %>%  filter(expMode2=="Test",rank<=10) |>
  ggplot(aes(x = x, y = pred, fill=condit)) + 
  stat_summary(fun=mean, geom="bar", position=position_dodge()) +
  stat_summary(fun.data=mean_se, geom="errorbar", position=position_dodge()) +
  facet_wrap(~Model)



teter |> select(sim_dat) |> unnest(sim_dat)  %>%  filter(expMode2=="Test",rank<=10) |>
  ggplot(aes(x = x, y = pred, fill=condit)) + 
  stat_summary(fun=mean, geom="bar", position=position_dodge()) +
  stat_summary(fun.data=mean_se, geom="errorbar", position=position_dodge()) +
  facet_wrap(rank~Model)


teter |> select(sim_dat) |> unnest(sim_dat)  %>%  filter(condit=="Constant",expMode2=="Test",rank<=20) |>
  mutate(facet_label = paste0("rank: ", rank, "\n", "c: ", round(c, 4), "\n", "lr: ", round(lr, 4), "\n", "distance: ", round(distance, 1))) |>
  ggplot(aes(x = x, y = pred, fill=condit)) + 
  stat_summary(fun=mean, geom="bar", position=position_dodge()) +
  stat_summary(fun.data=mean_se, geom="errorbar", position=position_dodge()) +
  facet_wrap(~facet_label)


tr |> select(sim_dat) |> unnest(sim_dat)  %>%  filter(condit=="Constant",expMode2=="Test",rank<=20) |>
  mutate(facet_label = paste0("rank: ", rank, "\n", "c: ", round(c, 4), "\n", "lr: ", round(lr, 4), "\n", "distance: ", round(distance, 1))) |>
  ggplot(aes(x = x, y = pred, fill=condit)) + 
  stat_summary(fun=mean, geom="bar", position=position_dodge()) +
  stat_summary(fun.data=mean_se, geom="errorbar", position=position_dodge()) +
  facet_wrap(~facet_label)



teter |> select(sim_dat) |> unnest(sim_dat)  %>%  filter(condit=="Varied",expMode2=="Test",rank<=20) |>
  ggplot(aes(x = x, y = pred, fill=condit)) + 
  stat_summary(fun=mean, geom="bar", position=position_dodge()) +
  stat_summary(fun.data=mean_se, geom="errorbar", position=position_dodge()) +
  facet_wrap(rank~Model)





post_data_v %>%  filter(expMode2=="Test") |>
  ggplot(aes(x = x, y = pred, fill=condit)) + 
  stat_summary(fun=mean, geom="bar", position=position_dodge()) +
  stat_summary(fun.data=mean_se, geom="errorbar", position=position_dodge()) 


post_data %>% filter(rank>1970, expMode2=="Test") |> 
  mutate(facet_label = paste0("rank: ", rank, "\n", "c: ", round(c, 4), "\n", "lr: ", round(lr, 4), "\n", "distance: ", round(distance, 1))) |>
  ggplot(aes(x = x, y = pred, fill=condit)) + 
  stat_summary(fun=mean, geom="bar", position=position_dodge()) +
  stat_summary(fun.data=mean_se, geom="errorbar", position=position_dodge()) +
  facet_wrap(~facet_label) 
  







#
#
#
#
#
#

# ggplot(tetr, aes(x = c, y = lr, color = distance)) +
#   geom_point()

#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#

#
#
#
