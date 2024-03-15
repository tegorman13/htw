pacman::p_load(dplyr,purrr,tidyr,tibble,ggplot2,
  stringr, here,conflicted, patchwork, knitr)
walk(c("dplyr"), conflict_prefer_all, quiet = TRUE)
walk(c("Display_Functions"), ~ source(here::here(paste0("Functions/", .x, ".R"))))



e1 <- readRDS(here("data/e1_08-21-23.rds")) 
e1Sbjs <- e1 |> group_by(id,condit) |> summarise(n=n())
testE1 <- e1 |> filter(expMode2 == "Test")

e1 <- readRDS(here("data/e1_08-21-23.rds")) 
e2 <- readRDS(here("data/e2_08-04-23.rds")) 
e3 <- readRDS(here("data/e3_08-04-23.rds")) 
d <- rbind(e1,e2,e3)


testE1 |> group_by(condit,bandType,vb) |> summarize(mean_dist=mean(dist),mean_vx=mean(vx))

wide_data = testE1 %>%
    select(id, condit, vb, vx) %>%
    group_by(id, condit, vb) %>%  # Group by 'id', 'condit', and 'vb' to handle non-unique entries
    summarise(vx = mean(vx), .groups = "drop") %>%  # Summarize 'vx' by taking the mean (or choose another appropriate summary function)
    pivot_wider(names_from = "vb", values_from = "vx")

wide_data_dist = testE1 %>%
    select(id, condit, vb, dist) %>%
    group_by(id, condit, vb) %>%  # Group by 'id', 'condit', and 'vb' to handle non-unique entries
    summarise(dist = mean(dist), .groups = "drop") %>%  # Summarize 'vx' by taking the mean (or choose another appropriate summary function)
    pivot_wider(names_from = "vb", values_from = "dist")    

wide_data_sdist = testE1 %>%
    select(id, condit, vb, sdist) %>%
    group_by(id, condit, vb) %>%  # Group by 'id', 'condit', and 'vb' to handle non-unique entries
    summarise(sdist = mean(sdist), .groups = "drop") %>%  # Summarize 'vx' by taking the mean (or choose another appropriate summary function)
    pivot_wider(names_from = "vb", values_from = "sdist")  


wide_data_dist %>% 
  filter(condit == "Varied") %>% 
  select(-condit,-id) %>% 
  cor()



correlation_matrix <- wide_data %>% 
    select(-id, -condit) %>%
    na.omit() %>%  
    cor(use = "complete.obs")

wide_data_dist %>% 
    select(-id, -condit) %>%
    na.omit() %>%  
    cor(use = "complete.obs")

levels(testE1$vb)
#[1] "100-300"   "350-550"   "600-800"   "800-1000"  "1000-1200" "1200-1400"

correlation_matrix %>%
  as_tibble(rownames = "vb_level") %>%
  pivot_longer(-vb_level, names_to = "vb_level2", values_to = "correlation") %>%
  mutate(vb_level = factor(vb_level, levels = levels(testE1$vb)),
         vb_level2 = factor(vb_level2, levels = levels(testE1$vb))) %>%
  ggplot(aes(x = vb_level, y = vb_level2, fill = correlation)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
  geom_text(aes(label = sprintf("%.2f", correlation)), color = "black", size = 3) +
  scale_x_discrete(limits = levels(testE1$vb)) +
  scale_y_discrete(limits = levels(testE1$vb)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_blank())





#####

correlation_matrices <- wide_data %>%
  group_by(condit) %>%
  summarise(correlation_matrix = list(cor(select(cur_data(), -id), use = "complete.obs"))) %>%
  ungroup()

plot_list <- correlation_matrices %>%
  mutate(plot = map2(condit, correlation_matrix, ~ {
    .y %>%
      as_tibble(rownames = "vb_level") %>%
      pivot_longer(-vb_level, names_to = "vb_level2", values_to = "correlation") %>%
      ggplot(aes(x = vb_level, y = vb_level2, fill = correlation)) +
      geom_tile() +
      scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
    scale_x_discrete(limits = levels(testE1$vb)) +
    scale_y_discrete(limits = levels(testE1$vb)) +
      geom_text(aes(label = sprintf("%.2f", correlation)), color = "black", size = 3) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            axis.title = element_blank()) +
      ggtitle(paste("Correlation Matrix for", .x))
  }))

wrap_plots(plot_list$plot)

##

create_corr_matrices <- function(data, group_var, exclude_cols = NULL) {
  correlation_matrices <- data %>%
    group_by({{ group_var }}) %>%
    summarise(correlation_matrix = list(cor(select(cur_data(), -all_of(exclude_cols)), use = "complete.obs"))) %>%
    ungroup()
  
  plot_list <- correlation_matrices %>%
    mutate(plot = map2({{ group_var }}, correlation_matrix, ~ {
      .y %>%
        as_tibble(rownames = "vb_level") %>%
        pivot_longer(-vb_level, names_to = "vb_level2", values_to = "correlation") %>%
        ggplot(aes(x = vb_level, y = vb_level2, fill = correlation)) +
        geom_tile() +
        scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
        scale_x_discrete(limits = levels(testE1$vb)) +
        scale_y_discrete(limits = levels(testE1$vb)) +
        geom_text(aes(label = sprintf("%.2f", correlation)), color = "black", size = 3) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              axis.title = element_blank()) +
        ggtitle(paste("Correlation Matrix for", .x))
    }))
  
  wrap_plots(plot_list$plot)
}

 #wide_data %>% filter(id %in% unique(testE1$id)) 

# cor mats for indv ids doesn't make sense
 #plot_list <- create_corr_matrices(wide_data |> filter(id %in% levels(testE1$id[1:3])), id, exclude_cols = "condit")

create_corr_matrices(wide_data, condit, exclude_cols = "id")
create_corr_matrices(wide_data_dist, condit, exclude_cols = "id")



## Approach 2

calculate_group_correlations <- function(data, group_var, cor_vars) {
  # Directly use the character vectors for grouping and selecting columns
  data %>%
    select(all_of(group_var), all_of(cor_vars)) %>%
    group_by_at(vars(group_var)) %>%
    group_map(~ {
      # Drop the grouping variable(s) for the correlation calculation
      .x %>% select(-all_of(group_var)) %>% 
        na.omit() %>%
        cor(use = "complete.obs")
    }, .keep = TRUE) %>%
    set_names(paste("Correlation_Matrix", data[[group_var]] %>% unique(), sep = "_"))
}

# Example usage: Calculate correlation matrices for each 'condit'
cor_matrices_condit <- calculate_group_correlations(wide_data, "condit", names(wide_data)[-(1:2)])


 cm_df <- map2(cor_matrices_condit,
    names(cor_matrices_condit),~{as_tibble(.x) |> 
        mutate(condit=stringr::str_remove(.y, "condit_"), 
        vb=c("100-300",   "350-550",   "600-800",   "800-1000",  "1000-1200", "1200-1400"))}) |> 
    data.table::rbindlist() |>
  pivot_longer(cols = -c(condit, vb), names_to = "vb_pair", values_to = "correlation") |>
    mutate(vb_pair = factor(vb_pair, levels = c("100-300", "350-550", "600-800", "800-1000", "1000-1200", "1200-1400")),
           vb = factor(vb, levels = c("100-300", "350-550", "600-800", "800-1000", "1000-1200", "1200-1400")))

 cm_df |> ggplot(aes(x=vb, y = correlation, fill=vb_pair)) + geom_col(stat="identity",position=position_dodge()) + facet_wrap(~condit)





plot_correlation_matrices <- function(cor_matrices_list) {
  plots <- lapply(names(cor_matrices_list), function(matrix_name) {
    # Convert the correlation matrix to a long format for plotting
    cor_data <- cor_matrices_list[[matrix_name]] %>%
      as.data.frame() %>%
      rownames_to_column("Variable1") %>%
      pivot_longer(cols = -Variable1, names_to = "Variable2", values_to = "correlation")
    
    # Plot
    ggplot(cor_data, aes(x = Variable1, y = Variable2, fill = correlation)) +
      geom_tile() +
      geom_text(aes(label = sprintf("%.2f", correlation)), size = 3) +
      scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
        scale_x_discrete(limits = levels(testE1$vb)) +
        scale_y_discrete(limits = levels(testE1$vb)) +
      labs(title = matrix_name, x = NULL, y = NULL) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  })
  
  return(plots)
}


cor_matrices_condit <- calculate_group_correlations(wide_data, "condit", names(wide_data)[-(1:2)])

plots_list <- plot_correlation_matrices(cor_matrices_condit)
#print(plots_list[[1]])
#invisible(lapply(plots_list, print))
wrap_plots(plots_list,ncol=2)


cor_matrices_condit_dist <- calculate_group_correlations(wide_data_dist, "condit", names(wide_data)[-(1:2)])
plots_list <- plot_correlation_matrices(cor_matrices_condit_dist)
wrap_plots(plots_list)


cor_matrices_condit <- calculate_group_correlations(wide_data_sdist, "condit", names(wide_data)[-(1:2)])
plots_list <- plot_correlation_matrices(cor_matrices_condit)
wrap_plots(plots_list)


 cm_df <- map2(cor_matrices_condit_dist,
    names(cor_matrices_condit),~{as_tibble(.x) |> 
        mutate(condit=stringr::str_remove(.y, "condit_"), 
        vb=c("100-300",   "350-550",   "600-800",   "800-1000",  "1000-1200", "1200-1400"))}) |> 
    data.table::rbindlist() |>
  pivot_longer(cols = -c(condit, vb), names_to = "vb_pair", values_to = "correlation") |>
    mutate(vb_pair = factor(vb_pair, levels = c("100-300", "350-550", "600-800", "800-1000", "1000-1200", "1200-1400")),
           vb = factor(vb, levels = c("100-300", "350-550", "600-800", "800-1000", "1000-1200", "1200-1400")))

 cm_df |> ggplot(aes(x=vb, y = correlation, fill=vb_pair)) + geom_col(stat="identity",position=position_dodge()) + facet_wrap(~condit)





cor_matrices_condit <- calculate_group_correlations(wide_data, "condit", names(wide_data)[-(1:2)])

# Function to create an MDS visualization for a correlation matrix
create_mds_plot <- function(cor_matrix, title) {
  mds_result <- cmdscale(1 - cor_matrix, eig = TRUE, k = 2)
  mds_coords <- data.frame(mds_result$points)
  colnames(mds_coords) <- c("Dim1", "Dim2")
  mds_coords$label <- rownames(cor_matrix)
  
  ggplot(mds_coords, aes(x = Dim1, y = Dim2, label = label)) +
    geom_point() +
    geom_text(vjust = -0.5) +
    ggtitle(title) +
    theme_minimal()
}

# Create MDS visualizations for each correlation matrix
mds_plots <- map2(cor_matrices_condit, names(cor_matrices_condit), ~ create_mds_plot(.x, .y))
wrap_plots(mds_plots, ncol = 2)

mds_plots <- map2(cor_matrices_condit_dist, names(cor_matrices_condit_dist), ~ create_mds_plot(.x, .y))
wrap_plots(mds_plots, ncol = 2)

######

cor_diff <- cor_matrices_condit$Correlation_Matrix_Varied - cor_matrices_condit$Correlation_Matrix_Constant
print(cor_diff)

plot_mds <- function(cor_matrix, title) {
  # Convert the correlation matrix to a dissimilarity matrix
  dissimilarity_matrix <- as.dist(1 - abs(cor_matrix))
  # Perform classical MDS
  mds <- cmdscale(dissimilarity_matrix)
  mds_df <- data.frame(Dim1 = mds[,1], Dim2 = mds[,2], Variables = rownames(cor_matrix))
  
  ggplot(mds_df, aes(x = Dim1, y = Dim2, label = Variables)) +
    geom_text() +
    ggtitle(title) +
    theme_minimal() +
    xlab("Dimension 1") +
    ylab("Dimension 2")
}

# Plot MDS for each condition
plot_mds(cor_matrices_condit$Correlation_Matrix_Varied, "MDS for Varied Condition")
plot_mds(cor_matrices_condit$Correlation_Matrix_Constant, "MDS for Constant Condition")



cm_e1 <- wide_data %>% 
    select(-id, -condit) %>%
    cor(use = "complete.obs")
plot_mds(cm_e1, "MDS for All Sbjs.")






cm_all <- wd_all %>% 
    select(-id, -condit) %>%
    cor(use = "complete.obs")
plot_mds(cm_all, "MDS for All Sbjs.")



calculate_group_correlations <- function(data, group_vars, cor_vars) {
  # Directly use the character vectors for grouping and selecting columns
  data %>%
    select(all_of(group_vars), all_of(cor_vars)) %>%
    group_by(across(all_of(group_vars))) %>%
    group_map(~ {
      # Drop the grouping variables for the correlation calculation
      corr_matrix <- .x %>%
        select(-all_of(group_vars)) %>%
        na.omit() %>%
        cor(use = "complete.obs")
      
      # Create a named list with the grouping variable values as names and the correlation matrix as the value
      setNames(list(corr_matrix), paste(group_vars, unlist(.y), sep = "_", collapse = "_"))
    }, .keep = TRUE) %>%
    flatten() %>%  # Flatten the list of lists into a single list
    setNames(names(.))  # Set the names of the list elements
}


wd_all = d %>%
    filter(expMode2 == "Test") %>%
    select(id, condit, vb,bandOrder,fb, vx) %>%
    group_by(id, condit, vb,bandOrder,fb) %>%  # Group by 'id', 'condit', and 'vb' to handle non-unique entries
    summarise(vx = mean(vx), .groups = "drop") %>%  # Summarize 'vx' by taking the mean (or choose another appropriate summary function)
    pivot_wider(names_from = "vb", values_from = "vx")

cm3 <- calculate_group_correlations(
        data=wd_all, 
        group_vars=c("condit","bandOrder","fb"), 
        cor_vars=names(wd_all)[-(1:4)] )

mds_plots <- map2(cm3, names(cm3), ~ create_mds_plot(.x, .y))
wrap_plots(mds_plots, ncol = 2)

plots_list <- plot_correlation_matrices(cm3)
wrap_plots(plots_list,ncol=2)


cm_fb <- calculate_group_correlations(
        data=wd_all, 
        group_vars=c("bandOrder"), 
        cor_vars=names(wd_all)[-(1:4)] )

mds_plots <- map2(cm_fb, names(cm_fb), ~ create_mds_plot(.x, .y))
wrap_plots(mds_plots, ncol = 2)

 
#data %>% select(all_of(group_vars), all_of(cor_vars)) %>% with(split(.,f =c(condit,bandOrder), drop=TRUE))




nbins=4
wd_all_trial =  d |> 
 filter(expMode2 == "Test") |>
group_by(id, condit, vb,bandOrder,fb) |> 
    mutate(Trial_Bin = cut( gt.bandStage, breaks = seq(1, max(gt.bandStage),length.out=nbins+1),
    include.lowest = TRUE, labels=FALSE)) %>%
    group_by(id,Trial_Bin, condit, vb,bandOrder,fb) %>%  # Group by 'id', 'condit', and 'vb' to handle non-unique entries
    summarise(vx = mean(vx), .groups = "drop") %>%  # Summarize 'vx' by taking the mean (or choose another appropriate summary function)
    pivot_wider(names_from = "vb", values_from = "vx")

calculate_group_correlations(
        data=wd_all_trial, 
        group_vars=c("condit","bandOrder","fb"), 
        cor_vars=names(wd_all)[-(1:4)] )





cid <- calculate_group_correlations(
        data=wd_all_trial |> filter(id %in% levels(wd_all_trial$id)[1:3]) , 
        group_vars=c("id"), 
        cor_vars=names(wd_all)[-(1:4)] )


mds_plots <- map2(cid, names(cid), ~ create_mds_plot(.x, .y))
wrap_plots(mds_plots, ncol = 2)

plots_list <- plot_correlation_matrices(cid)
wrap_plots(plots_list,ncol=2)


library(cluster)
individual_dissimilarities <- wd_all_trial %>%
  group_by(id, condit) %>%
  summarise(across(starts_with("dist"), scale)) %>% # Assuming data needs to be scaled
  dist() %>%
  as.matrix()
agnes_result <- agnes(individual_dissimilarities, diss = TRUE)

plot(agnes_result, which.plots = 2)

# Alternatively, perform MDS
mds_result <- cmdscale(individual_dissimilarities)

# Create a data frame for the MDS results
mds_df <- data.frame(Dimension_1 = mds_result[,1], Dimension_2 = mds_result[,2], ID = rownames(mds_result))

# Plot MDS
# ggplot(mds_df, aes(x = Dimension_1, y = Dimension_2, color = ID)) +
#   geom_point() +
#   theme_minimal() +
#   labs(title = "MDS Plot of Individual Differences", x = "Dimension 1", y = "Dimension 2")








####### 

calculate_individual_cor <- function(df,dv=vx) {
  wide_df <- df %>% 
    pivot_wider(names_from = vb, values_from = {{dv}}) %>%
    select(-condit,-id,-gt.bandStage)|>
    #select("100-300",   "350-550",   "600-800",   "800-1000",  "1000-1200", "1200-1400") |>
    relocate("100-300",   "350-550",   "600-800",   "800-1000",  "1000-1200", "1200-1400")

  cor_matrix <- cor(wide_df, use = "pairwise.complete.obs")

  return(cor_matrix)
}

individual_cor_matrices <- list()

# Split the data frame by 'id', then calculate the correlation matrix for each subject
cid <- testE1 %>%
 select(id,gt.bandStage,condit,vb,vx) |>
  group_by(id) %>%
  group_split(.keep = TRUE) %>%  # .keep = FALSE removes the grouping variable from the resulting data frames
  map(~ {
    individual_cor_matrices[[as.character(.x$id[1])]] <- calculate_individual_cor(.x) |> as_tibble() |>
    mutate(id = .x$id[1], condition = .x$condit[1], vb=c("100-300",   "350-550",   "600-800",   "800-1000",  "1000-1200", "1200-1400"),.before=1) 
  })

map(c(1,2,3), ~{cid |> pluck(.x)})


dc <- cid |> data.table::rbindlist() 
head(dc)

plot_data <- dc |>
  pivot_longer(cols = -c(id, condition, vb), names_to = "vb_pair", values_to = "correlation") |>
    mutate(vb_pair = factor(vb_pair, levels = c("100-300", "350-550", "600-800", "800-1000", "1000-1200", "1200-1400")),
           vb = factor(vb, levels = c("100-300", "350-550", "600-800", "800-1000", "1000-1200", "1200-1400")))

plot_data |> 
filter(correlation<.95) |>
ggplot(aes(x = vb, y = correlation, fill = vb_pair)) +
    stat_bar + facet_wrap(~condition)




cid_dist <- testE1 %>%
 select(id,gt.bandStage,condit,vb,dist) |>
  group_by(id) %>%
  group_split(.keep = TRUE) %>%  # .keep = FALSE removes the grouping variable from the resulting data frames
  map(~ {
    individual_cor_matrices[[as.character(.x$id[1])]] <- calculate_individual_cor(.x,dist) |> as_tibble() |>
    mutate(id = .x$id[1], condition = .x$condit[1], vb=c("100-300",   "350-550",   "600-800",   "800-1000",  "1000-1200", "1200-1400"),.before=1) 
  })

map(c(1,2,3), ~{cid_dist |> pluck(.x)})
dc_dist <- cid_dist |> data.table::rbindlist() 

plot_data <- dc_dist |>
  pivot_longer(cols = -c(id, condition, vb), names_to = "vb_pair", values_to = "correlation") |>
    mutate(vb_pair = factor(vb_pair, levels = c("100-300", "350-550", "600-800", "800-1000", "1000-1200", "1200-1400")),
           vb = factor(vb, levels = c("100-300", "350-550", "600-800", "800-1000", "1000-1200", "1200-1400")))

plot_data |> 
filter(correlation<.95) |>
ggplot(aes(x = vb, y = correlation, fill = vb_pair)) +
    stat_bar + facet_wrap(~condition)



#######    


testE1 |>  select(id,gt.bandStage,condit,vb,vx) |> mutate(k=1, .before=names(cur_data()[1]))
testE1 |>  select(id,gt.bandStage,condit,vb,vx) |> mutate(k=1, .before=1)

correlations_of_interest <- map_df(cid, ~ {
  training_item <- "800-1000"
  extrapolation_items <- c("100-300", "350-550", "600-800")
  
  correlations <- .x[extrapolation_items, training_item]
  


  tibble(
    id = names(.),
    `800-1000_100-300` = correlations["100-300"],
    `800-1000_350-550` = correlations["350-550"],
    `800-1000_600-800` = correlations["600-800"]
  )
})


cdf <- cid |> data.table::rbindlist()

names(cid)

#ind_fits_df <- ind_fits_df |> map(~rbindlist(.x$dat) |> mutate(Model = .x$Model, Fit_Method = .x$Fit_Method)) |> rbindlist() 