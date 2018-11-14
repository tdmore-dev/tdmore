setwd(devtools::package_file("bookdown/"))
output <- devtools::package_file("docs/book/")
bookdown::render_book(input="index.Rmd",
                      output_format="bookdown::gitbook",
                      output_dir = output,
                      new_session=TRUE)
unlink("_bookdown_files", recursive=TRUE) #clean up
unlink("_main.rds")
setwd(devtools::package_file("."))
