% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/algebraicLibrary.R
\name{linCmt}
\alias{linCmt}
\title{This function guesses which model should be calculated, based on the available variables
in the caller environment}
\usage{
linCmt(env = parent.frame())
}
\arguments{
\item{env}{environment in which to search for parameters}
}
\description{
Please note this function also modifies the caller environment with the
appropriate values for the call to the algebraic function.
As an example, it will define 'K'.
}
