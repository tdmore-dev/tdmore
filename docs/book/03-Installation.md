---
output: html_document
editor_options: 
  chunk_output_type: console
---
# Installation  {#install}

TDMore depends on package `RxODE`, which requires a working C and fortran compiler to work. Installation procedure can be found [here](https://github.com/nlmixrdevelopment/RxODE). Once `RxODE` installed, simply execute the following command in the R console:


```r
devtools::install_github("tdmore-dev/tdmore")
```

If you would like to install the package with all its vignettes, please run the following command instead:

```r
devtools::install_github("tdmore-dev/tdmore", build_vignettes = TRUE)
```

# Installation script

The below script can be used to automatically install TDMore. It uses an API key that still needs to be filled in.


```r
if(R.version$major != "3" || as.numeric(R.version$minor) < 5.2) {
  stop("tdmore requires at least R version 3.5.2
Please install a more recent version from https://cran.r-project.org/bin/windows/base/
You may have to adapt the settings of RStudio in Tools -> Global Options as well.")
}
if(!require(devtools)) install.packages('devtools')
if(!require(tidyverse)) install.packages('tidyverse')

PAT_TOKEN <- "SECRET" #request the PAT_TOKEN from us to access the repository
devtools::install_github("tdmore-dev/tdmore", user="tdmore-training", auth_token=PAT_TOKEN,
                         upgrade="never")

tryCatch(
  {
    RxODE::RxODE( "d/dt(A0) = -0.5 * A0;" )
  }, error=function(e) {
    stop(
"RxODE could not compile a model.
This does not stop tdmore from working, nor will it stop you from completing the training.

However, you will not be able to use RxODE or nlmixr models.
To use nlmixr-style models, you will have to install Rtools.
You can find this on https://cran.r-project.org/bin/windows/Rtools/Rtools35.exe
"
)
  })
```