% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/modelMixture.R
\name{tdmore_mixture}
\alias{tdmore_mixture}
\title{Create a TDMore mixture model.}
\usage{
tdmore_mixture(..., probs)
}
\arguments{
\item{...}{2 or more tdmore models that describe different subpopulations}

\item{probs}{'a priori' probabilities for belonging to the different subpopulations, numeric vector}
}
\value{
a tdmore_mixture object
}
\description{
Create a TDMore mixture model.
}
