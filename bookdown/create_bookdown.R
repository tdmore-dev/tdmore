setwd(devtools::package_file("bookdown/"))
output <- devtools::package_file("docs/book/")
bookdown::render_book(input="index.Rmd", output_format="bookdown::gitbook", output_dir = output)
