# Preface {-}
`TDMore` provides an easy interface to execute post-hoc bayesian estimation of individual profiles, and to find the best dose to put the patient on target.

The package is intended to make it easy to define your own models, and easily create a dose decision support tool for physicians.

## About this book {-}
This book provides guidance on the use of `TDMore` and serves as a basic reference manual.

We assume that users have some notions in PK/PD modeling and simulation and are familiar with the R language. The installation of `TDMore` and related software and packages is described in Chapter \@ref(install).

## Acknowledgements {-}
We would like to thank all the people who are contributing to this project.

# Introduction {#intro}
## Prerequisites {#prereq}
The `TDMore` framework is based on the [`R`](https://cran.r-project.org/) programming language, a free sofware environment for statistical computing.

`TDMore` relies on a series of R open-source packages. In particular, it depends on:

- [`RxODE`](https://cran.r-project.org/package=RxODE)(An R package for simulating PK/PD models. This package is intensively used by `TDMore` to predict data)
- [`nlmixr`](https://cran.r-project.org/package=nlmixr)(An R package for estimating PK/PD model parameters. Although `TDMore` also estimates model parameters, it does not rely on the `nlmixr` for the estimation part. However, this package is  used for its nice way of defining PK/PD models.)

More information on the installation of packages and software can be found in Chapter \@ref(install).
