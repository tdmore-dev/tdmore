setwd(devtools::package_file("bookdown/"))
bookoutputdir <- devtools::package_file("docs/book/")
bookoutputdir <- force(bookoutputdir)
message("Creating bookdown in ", bookoutputdir)

## Because bookdown::render_book is executed in a new session, it does not have access to
## variables defined in the current R script. This is why we need to use do.call(),
## and cannot call bookdown::render_book directly (it would otherwise complain that 'bookoutputdir'
## is an unknown variable )
arguments <- list(input="index.Rmd",
                      output_format="bookdown::gitbook",
                      output_dir = bookoutputdir,
                      new_session=TRUE)
do.call(bookdown::render_book, arguments)
unlink("_bookdown_files", recursive=TRUE) #clean up
unlink("_main.rds")
setwd(devtools::package_file("."))
