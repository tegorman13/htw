#
#
#
#
#
#
#


process_dataset <- function(dataset) {
  # Process the dataset and return the combined results
  teter_combined <- map_dfr(dataset, "teter_results", .id = "run_name") 
  te_combined <- map_dfr(dataset, "te_results", .id = "run_name")
  tr_combined <- map_dfr(dataset, "tr_results", .id = "run_name") 
  
  return(list(teter_combined = teter_combined, te_combined = te_combined, tr_combined = tr_combined))
}

create_density_plots <- function(combined_data, plot_title_suffix) {
  # Create and return density plot objects with dynamic titles
  plots <- map(combined_data, function(data, name) {
    ggplot(data, aes(x = c)) + geom_density(aes(color = Model)) +
      facet_wrap(~Group, scales = "free") + ggtitle(sprintf("c posterior - %s - %s", name, plot_title_suffix))
  }, .id = "Group Type")
  
  # Combine density plots
  density_combined <- wrap_plots(plots, ncol = 2)
  
  return(density_combined)
}

create_distance_plots <- function(combined_data, plot_title_suffix) {
  # Create and return distance plot objects with dynamic titles
  dist_plots <- map(combined_data, function(data, name) {
    ggplot(data, aes(x = Group, y = distance, fill = Model)) + 
      stat_summary(fun = mean, geom = "bar", position = position_dodge()) +
      stat_summary(fun.data = mean_se, geom = "errorbar", position = position_dodge()) +
      ggtitle(sprintf("%s - %s", name, plot_title_suffix))
  }, .id = "Group Type")
  
  # Combine distance plots
  distance_combined <- wrap_plots(dist_plots, ncol = 2)
  
  return(distance_combined)
}

save_plots <- function(plot, filename) {
  # Save the plot to the specified filename
  ggsave(filename = here::here(filename), plot = plot, width = 10, height = 6)
}

# Example usage
name <- "example_dataset_name"
sample_size <- str_extract(name, "\\d+p?\\d*?M")
posterior_cutoff <- str_extract(name, "p\\d+$")
plot_title_suffix <- sprintf("%s - %s", sample_size, posterior_cutoff)

# Process the dataset
processed_data <- process_dataset(dataset)

# Create density plots
density_combined <- create_density_plots(processed_data, plot_title_suffix)
density_filename <- sprintf("assets/tmp_plots/density_plots_combo_abc_%s_rmse_%s.png", sample_size, posterior_cutoff)
save_plots(density_combined, density_filename)

# Create distance plots
distance_combined <- create_distance_plots(processed_data, plot_title_suffix)
distance_filename <- sprintf("assets/tmp_plots/distance_plots_combo_abc_%s_rmse_%s.png", sample_size, posterior_cutoff)
save_plots(distance_combined, distance_filename)


#
#
#
