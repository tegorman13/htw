pacman::p_load(dplyr,purrr,tidyr,ggplot2, here, patchwork, conflicted, stringr, quarto)
conflict_prefer_all("dplyr", quiet = TRUE)

# currently works running in vs code when in the tmp_plot folder. Doesn't work from terminal

# find Data/indv_sim -size +18M | cat
#list.files("../../data/indv_sim/", pattern = "ind_abc_ss*")

#ind_abc_ss_1500_ng100_buf10.rds
#ind_abc_ss_1500_ng75_buf20.rds
#ind_abc_ss_1600_g150_bufp25.rds


abc_files <- c("ind_abc_ss_1500_ng100_buf10.rds", "ind_abc_ss_1500_ng75_buf20.rds", 
  "ind_abc_ss_1600_g150_bufp25.rds", "ind_abc_ss_1500_ng75_buf20.rds", "ind_abc_ss_1600_g150_bufp25.rds",
  "ind_abc_ss_2000_ng100_buf25.rds","ind_abc_ss_2500_ng100_buf10.rds","ind_abc_ss_2000_ng100_buf20.rds")

#abc_files <- list.files("../../data/indv_sim/", pattern = "ind_abc_ss_*")
abc_files <- list.files("../../data/indv_sim",pattern="*.rds")

abc_files <- c("ind_abc_1000_ng150_buf20.rds","ind_abc_ss_1500_ng100_buf10.rds" )
#abc_files <- c("ind_abc_1000_ngrid150_buf20.rds")
# abc_files <- c("ind_abc_ss_1500_ng100_buf10.rds")

# print messaging listing files in abc_files that report will be rendered for
print(paste0("Rendering reports for ", length(abc_files), " files"))
print(abc_files)


# check if .qmd files already exist for items in abc_files





render_report <- function(input) {
  tmp_template <- paste0(tools::file_path_sans_ext(input),".qmd")
  title <- tools::file_path_sans_ext(input)
  tmp_template_path <- paste0("hier_reports/",tmp_template)
  file.copy("_temp_hier.qmd",tmp_template_path)
  output_file <- paste0("hier_reports/", paste0(tools::file_path_sans_ext(input),".html"))
  output_file_name <- basename(output_file)
  quarto::quarto_render(input = tmp_template_path,
                #output_file=output_file_name,
                execute_params = list("file_name" = paste0("../../../data/indv_sim/",input),
                                      "title"=title))
  #file.copy(from = output_file_name, to = output_file)
  #file.remove(output_file_name)
}





# print("removing existing .qmd files in hier_reports")
# file.remove(list.files("hier_reports", pattern = ".qmd", full.names = TRUE))




qmd_files <- list.files("hier_reports", pattern = "\\.qmd$")
abc_files_modified <- sub("\\.rds$", ".qmd", abc_files)
files_to_remove <- abc_files_modified %in% qmd_files
abc_files <- abc_files[!files_to_remove]






purrr::map(abc_files, render_report)




# delete all .qmd files in hier_reports - won't work in R project





#render_report <- render_report("ind_abc_ss_1000_ng100_bufp55.rds", "ind_abc_ss_1000_ng100_bufp55_2.html")


# render_report <- function(input, output_file) {
#   output_file <- paste0("Assets/tmp_plots/", output_file)
#   output_file_name <- basename(output_file)
#   quarto_render(input = here::here("Assets/tmp_plots/_temp_hier.qmd"),
#                 output_file=output_file_name,
#                         execute_params = list("file_name" = paste0("data/indv_sim/",input),
#                                               "title"=paste0(input,"Plots")))
#   file.copy(from = output_file_name, to = output_file)
#   file.remove(output_file_name)
# }







# purrr::walk(all_data, ~quarto::quarto_render(
#   input = here::here("Model/_temp_hier.qmd"),
#   execute_params = list("year" = .x,
#                         "state" = .y),
#   output_file = glue::glue("{.y}_{.x}.html")
# ))


#paste0(tools::file_path_sans_ext(basename("data/indv_sim/ind_abc_ss_1000_ng100_bufp55.rds")),".html")
#raw_name <- tools::file_path_sans_ext(basename(params$file_name))
#fn="data/indv_sim/ind_abc_ss_1000_ng100_bufp55.rds"
# get basename
#basename(fn)
#"ind_abc_ss_1000_ng100_bufp55.rds"
# get filename without extension



