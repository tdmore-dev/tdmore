if(packageVersion('bookdown') != "0.12") stop("Expecting bookdown version 0.12")


build <- function() {
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
}

buildWithPackage <- function() {
  pkg <- devtools::as.package(".")
  tmp_lib <- tempfile("R_LIBS")
  dir.create(tmp_lib)
  on.exit({
    unlink(tmp_lib, recursive = TRUE)
  }, add=T)
  utils::install.packages(repos = NULL,
                          lib = tmp_lib,
                          pkg$path,
                          type = "source",
                          INSTALL_opts = c("--example",
                                           "--install-tests",
                                           "--with-keep.source",
                                           "--with-keep.parse.data",
                                           "--no-multiarch"),
                          quiet = FALSE)
  withr::with_libpaths(tmp_lib, {
    build()
  }, action="prefix")
}

buildWithPackage()
