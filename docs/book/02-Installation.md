# Installation  {#install}

TDMore depends on package `RxODE`, which requires a working C and fortran compiler to work. Installation procedure can be found [here](https://github.com/nlmixrdevelopment/RxODE). Once `RxODE` installed, simply execute the following command in the R console:


```r
devtools::install_github("tdmore-dev/tdmore")
```

If you would like to install the package with all its vignettes, please run the following command instead:

```r
devtools::install_github("tdmore-dev/tdmore", build_vignettes = TRUE)
```