% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ggplot.autolayer.R
\name{parameterPlot.tdmorefit}
\alias{parameterPlot.tdmorefit}
\title{Generate a plot of the parameters}
\usage{
parameterPlot.tdmorefit(x, newdata, vars, ...)
}
\arguments{
\item{x}{tdmorefit object}

\item{newdata}{A data.frame with new data and the columns to predict,
or a numeric vector to specify times, and predict all model output
or NULL to interpolate between 0 and the maximum known times}

\item{vars}{the variables to plot, as a character vector. If missing, this is taken from the observed_variables metadata}

\item{...}{ignored}
}
\description{
Generate a plot of the parameters
}
