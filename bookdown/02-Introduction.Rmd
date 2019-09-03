# Introduction {#intro}
## Prerequisites {#prereq}
The `TDMore` framework is based on the [`R`](https://cran.r-project.org/) programming language, a free sofware environment for statistical computing.

`TDMore` relies on a series of R open-source packages. In particular, it depends on:

- [`RxODE`](https://cran.r-project.org/package=RxODE)(An R package for simulating PK/PD models. This package is intensively used by `TDMore` to predict data)
- [`nlmixr`](https://cran.r-project.org/package=nlmixr)(An R package for estimating PK/PD model parameters. Although `TDMore` also estimates model parameters, it does not rely on the `nlmixr` for the estimation part. However, this package is  used for its nice way of defining PK/PD models.)

More information on the installation of packages and software can be found in Chapter \@ref(install).
