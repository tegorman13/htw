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


folder_2_gen <- folders |> slice(1:5) |> pull(path)



purrr::map(folder_2_gen, render_report)















# list.files(here('data/abc_reject'))
# (grep("Train",list.files(here(paste0('data/abc_reject/',file_name)),
#                                            pattern="EXAM_Test",full.names = TRUE),
#                                 invert=TRUE, value=TRUE))
# files <- cbind(files, m= file.info(files$name)$mtime) %>%
#  arrange(desc(m)) 