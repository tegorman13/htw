
# 
# create_summary_table <- function(data, DV, mfun = list(mean = mean, sd = sd)) {
#   # Grouping by required variables and computing statistics specified in mfun
#   summarization <- data %>%
#     group_by(vb, bandType, condit) %>%
#     summarise(across(.data[[DV]], mfun, .names = "{.fn}_{.col}"), .groups = 'drop')
#   
#   # Correctly renaming the columns to have clear stat names
#   cols_to_rename <- grep(paste0("^.*_", DV), names(summarization), value = TRUE)
#   new_colnames <- sapply(strsplit(cols_to_rename, "_"), `[`, 1)
#   names(summarization)[names(summarization) %in% cols_to_rename] <- new_colnames
#   
#   return(summarization)
# }
# 
# create_summary_table(test, "vx", mfun = list(mean = mean, median = median, sd = sd))


create_summary_table <- function(data, DV, mfun = list(mean = mean, sd = sd)) {
  # Grouping by required variables and computing statistics specified in mfun
  summarization <- data %>%
    group_by(vb, bandType, condit) %>%
    summarise(across(all_of(DV), mfun, .names = "{.fn}_{.col}"), .groups = 'drop') %>%
    mutate(across(matches(DV), ~ round(.x, 0)))
  
  # Correctly renaming the columns to have clear stat names
  cols_to_rename <- grep(paste0("^.*_", DV), names(summarization), value = TRUE)
  new_colnames <- sapply(strsplit(cols_to_rename, "_"), `[`, 1)
  new_colnames <- sapply(new_colnames, function(name) paste0(toupper(substr(name, 1, 1)), substr(name, 2, nchar(name))))
  names(summarization)[names(summarization) %in% cols_to_rename] <- new_colnames
  summarization <- summarization |> rename("Band"=vb, "Band Type"=bandType)
  
  
  # Splitting the table into constant and varied groups
  constant_table <- summarization %>%
    filter(condit == "Constant") %>%
    select(-condit) %>%
    kable(caption = paste("Summary of", DV, "- Constant")) #|>
    #kable_minimal(full_width = FALSE) |>
    add_header_above(c("Constant Testing " = ncol(summarization) - 1))
  
  varied_table <- summarization %>%
    filter(condit == "Varied") %>%
    select(-condit) %>%
    kable(caption = paste("Summary of", DV, "- Varied")) #|>
    #kable_minimal(full_width = FALSE) |>
    #add_header_above(c("Varied Testing " = ncol(summarization) - 1))
  
  return(list(constant = constant_table, varied = varied_table))
}

#result <- create_summary_table(test, "vx", mfun = list(mean = mean, median = median, sd = sd))
#result$constant



create_table <- function(data, DV) {
  compute_summary_label <- function(data, DV) {
    if (DV == "Percent_Hit") {
      data %>%
        group_by(id, condit, vb, bandType) %>%
        summarise(nHits = sum(dist == 0), Percent_Hit = nHits / n(), .groups = 'drop') %>%
        group_by(vb, condit, bandType) %>%
        summarise(Percent_HitMean = mean(Percent_Hit), Percent_HitSd = sd(Percent_Hit), .groups = 'drop') %>%
        mutate(
          meanLab = paste0("Mean=", round(Percent_HitMean, 3)),
          sdLab = paste0("SD=", round(Percent_HitSd, 2)),
          sumStatLab = paste0(meanLab, "\n", sdLab)
        ) %>%
        select(vb, condit, bandType, sumStatLab)
    } else {
      summarization <-
        data %>%
        group_by(id, vb, bandType, condit) %>%
        summarise(value = mean(.data[[DV]]), .groups = 'drop') %>%
        group_by(vb, bandType, condit) %>%
        summarise(meanValue = mean(value), sdValue = sd(value), .groups = 'drop') %>%
        mutate(
          meanLab = paste0("Mean=", round(meanValue, 0)),
          sdLab = paste0("SD=", round(sdValue, 0)),
          sumStatLab = paste0(meanLab, "\n", sdLab)
        ) %>%
        select(vb, condit, bandType, sumStatLab)
    }
  }
  
  summary_label <- compute_summary_label(data, DV)
  constant_table <- summary_label %>%
    filter(condit == "Constant") %>%
    rename(Constant = sumStatLab) %>%
    select(-condit)
  
  varied_table <- summary_label %>%
    filter(condit == "Varied") %>%
    rename(Varied = sumStatLab) %>%
    select(-condit)
  
  final_table <- full_join(constant_table, varied_table, by = "vb") |>
    rename("Band Type"=bandType.x, "Band Type "=bandType.y,"Band"=vb)
  
  final_table %>%
    kbl(digits = c(0, 0, 0, 0, 0),
        caption = paste("Summary of", DV),escape=F,booktabs = T) %>%
    # kable_minimal(full_width = FALSE,
    #               position = "left") %>%
    add_header_above(c(" " = 1, " " = 2, " " = 2))
  #add_header_above(c(" " = 1, "Constant" = 2, "Varied" = 2))
  
}









# 
# 
# create_table <- function(data, DV, IdAgg) {
#   compute_summary_label <- function(data, DV, IdAgg) {
#     if (DV == "Percent_Hit") {
#       data %>%
#         group_by(id, condit, vb, bandType) %>%
#         summarise(nHits = sum(dist == 0), Percent_Hit = nHits / n(), .groups = 'drop') %>%
#         group_by(vb, condit, bandType) %>%
#         summarise(Percent_HitMean = mean(Percent_Hit), Percent_HitSd = sd(Percent_Hit), .groups = 'drop') %>%
#         mutate(
#           meanLab = paste0("Mean=", round(Percent_HitMean, 3)),
#           sdLab = paste0("SD=", round(Percent_HitSd, 2)),
#           sumStatLab = paste0(meanLab, "\n", sdLab)
#         ) %>%
#         select(vb, condit, bandType, sumStatLab)
#     } else {
#       summarization <- if (IdAgg) {
#         data %>%
#           group_by(id, vb, bandType, condit) %>%
#           summarise(value = mean(.data[[DV]]), .groups = 'drop') %>%
#           group_by(vb, bandType, condit) %>%
#           summarise(meanValue = mean(value), sdValue = sd(value), .groups = 'drop')
#       } else {
#         data %>%
#           group_by(vb, bandType, condit) %>%
#           summarise(meanValue = mean(.data[[DV]]), sdValue = sd(.data[[DV]]), .groups = 'drop')
#       }
#       summarization %>%
#         mutate(
#           meanLab = paste0("Mean=", round(meanValue, 0)),
#           sdLab = paste0("SD=", round(sdValue, 0)),
#           sumStatLab = paste0(meanLab, "\n", sdLab)
#         ) %>%
#         select(vb, condit, bandType, sumStatLab)
#     }
#   }
#   
#   summary_label <- compute_summary_label(data, DV, IdAgg)
#   constant_table <- summary_label %>%
#     filter(condit == "Constant") %>%
#     rename(Constant = sumStatLab) %>%
#     select(-condit)
#   
#   varied_table <- summary_label %>%
#     filter(condit == "Varied") %>%
#     rename(Varied = sumStatLab) %>%
#     select(-condit)
#   
#   final_table <- full_join(constant_table, varied_table, by = "vb") |>
#     rename("Band Type"=bandType.x, "Band Type "=bandType.y,"Band"=vb)
#   
#   final_table %>%
#     kbl(digits = c(0, 0, 0, 0, 0),
#         caption = paste("Summary of", DV),escape=F,booktabs = T) %>%
#     kable_minimal(full_width = FALSE,
#                   position = "left") %>%
#     add_header_above(c(" " = 1, " " = 2, " " = 2))
#   #add_header_above(c(" " = 1, "Constant" = 2, "Varied" = 2))
# }
# 
# 
# create_table(test, "vx", IdAgg=FALSE)
# create_table(test, "Percent_Hit", IdAgg=TRUE)