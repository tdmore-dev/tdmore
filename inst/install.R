if(R.version$major != "3" || as.numeric(R.version$minor) < 5.2) {
  stop("tdmore requires at least R version 3.5.2
Please install a more recent version from https://cran.r-project.org/bin/windows/base/
You may have to adapt the settings of RStudio in Tools -> Global Options as well.")
}
options(install.packages.check.source = "no")
options("install.packages.compile.from.source"="no")
ensureInstalled <- function(x) {tryCatch( find.package(x), error=function(e) install.packages(x) )}
ensureInstalled("devtools")
ensureInstalled("tidyverse")

Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")
devtools::install_github("nlmixrdevelopment/RxODE@937597db213208ccef737f2bb532a88416ea139e")

PAT_TOKEN <- "eab9e908f5073f3ff82a7cb77e5fced10ed01d5d"
devtools::install_github("tdmore-dev/tdmore", user="tdmore-training", auth_token=PAT_TOKEN, upgrade="never")

tryCatch(RxODE::RxODE( "d/dt(A0) = -0.5 * A0;" ), error=function(e) message("
= Warning: could not compile model
We failed to compile an RxODE model.
This does not stop tdmore from working, nor will it stop you from completing the training.

However, you will not be able to use RxODE or nlmixr models.
To use nlmixr-style models, you will have to install Rtools.
You can find this on https://cran.r-project.org/bin/windows/Rtools/Rtools35.exe
"
))

if(library(tdmore, logical.return=T)) {
    message("===========================================")
    message("==== TDMore was installed succesfully! ====")
    message("===========================================")
} else {
    stop("Failed to install tdmore...")
}
