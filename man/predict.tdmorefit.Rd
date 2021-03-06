% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ebe.R
\name{predict.tdmorefit}
\alias{predict.tdmorefit}
\title{Predict new data using a model fit}
\usage{
\method{predict}{tdmorefit}(
  object,
  newdata = NULL,
  regimen = NULL,
  parameters = NULL,
  covariates = NULL,
  se.fit = FALSE,
  level = 0.95,
  mc.maxpts = 100,
  ...,
  .progress = interactive()
)
}
\arguments{
\item{object}{A tdmorefit object}

\item{newdata}{A data.frame with new data and the columns to predict,
or a numeric vector to specify times, and predict all model output
or NULL to interpolate between 0 and the maximum known times}

\item{regimen}{Treatment regimen}

\item{parameters}{named numeric vector of fixed parameters}

\item{covariates}{the model covariates, named vector, or data.frame with column 'TIME', and at least TIME 0}

\item{se.fit}{TRUE to provide a confidence interval on the prediction, adding columns xxx.median, xxx.upper and xxx.lower
FALSE to show the model prediction (IPRED)}

\item{level}{The confidence interval, or NA to return all mc.maxpts results}

\item{mc.maxpts}{Maximum number of points to sample in Monte Carlo simulation}

\item{...}{ignored}

\item{.progress}{Allows to specify a plyr-like progress object
A plyr progress object is a list with 3 function definitions: `init(N)`, `step()` and `term()`.
This can also be specified as a boolean. TRUE uses the default dplyr progress_estimated.}
}
\value{
A data.frame
}
\description{
Predict new data using a model fit
}
