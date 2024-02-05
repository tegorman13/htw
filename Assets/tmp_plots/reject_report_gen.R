pacman::p_load(dplyr,purrr,tidyr,ggplot2, here, patchwork, conflicted, stringr, quarto)
conflict_prefer_all("dplyr", quiet = TRUE)




render_report <- function(input) {
 title <- tools::file_path_sans_ext(basename(input))
  tmp_template <- paste0(tools::file_path_sans_ext(title),".qmd")
  tmp_template_path <- paste0("reject_reports/",tmp_template)
  file.copy("_temp_reject.qmd",tmp_template_path)
  output_file <- paste0("reject_reports/", paste0(tools::file_path_sans_ext(input),".html"))
  output_file_name <- basename(output_file)
  quarto::quarto_render(input = tmp_template_path,
                execute_params = list("file_name" = paste0(input),
                                      "title"=title))
}



folders <- tibble(path=list.files('../../data/abc_reject',full.names=TRUE)) |> 
    mutate(mtime = file.info(path)$mtime) |> 
    arrange(desc(mtime)) 

folder_2_gen <- folders |> slice(1:35) |> pull(path)


complete_reports <- tools::file_path_sans_ext(list.files("reject_reports", pattern = "\\.qmd$"))
# abc_files_modified <- sub("\\.rds$", ".qmd", basename(complete_reports))
files_to_remove <- which(basename(folder_2_gen) %in% complete_reports)
abc_files <- folder_2_gen[-files_to_remove] # files that still need to be processed


# file.remove(list.files("reject_reports", pattern = ".qmd", full.names = TRUE))

purrr::map(abc_files, render_report)


source("ind_fit_inspect.R")
# run "python gen_index.py with system"

#system("python gen_index.py")














# list.files(here('data/abc_reject'))
# (grep("Train",list.files(here(paste0('data/abc_reject/',file_name)),
#                                            pattern="EXAM_Test",full.names = TRUE),
#                                 invert=TRUE, value=TRUE))
# files <- cbind(files, m= file.info(files$name)$mtime) %>%
#  arrange(desc(m)) 