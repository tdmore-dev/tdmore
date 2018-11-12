library(bookdown)

if(requireNamespace("rstudioapi")) setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
bookdown::render_book("index.Rmd", "bookdown::gitbook", output_dir = "./../")
