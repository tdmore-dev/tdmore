#' Get an overview of the log-likelihood for varying parameter values.
#' If arguments `fix` and `limits` are not set for a specific parameter, individual eta's will be generated to cover 95 percent of the population distribution.
#'
#' @param fitted A tdmorefit object
#' @param fix Which parameters to fix? Named vector of the parameters that are fixed and should not be profiled
#' Specify a FIX of NA to estimate the optimal value at that specific grid-point.
#' @param maxpts Maximum number of points per parameter
#' @param limits limits to explore (numeric vector of form c(min, max) or specific limits per parameter in the form list(ETA1=c(min,max), ETA2=c(min,max), etc))
#' @param type log-lokelihood function type, 3 possible values: 'pop', 'pred' or 'll' (= pop + pred)
#' @param .progress Allows to specify a plyr-like progress object
#' A plyr progress object is a list with 3 function definitions: `init(N)`, `step()` and `term()`.
#' This can also be specified as a boolean. TRUE uses the default dplyr progress_estimated.
#'
#' @param ... ignored
#'
#' @return a tdmore profile object. It namely contains a data.frame with each parameter value tested, and an additional `logLik` column with the log-likelihood for each parameter combination.
#' @engine
profile.tdmorefit <- function(fitted, fix=NULL, maxpts = 50, limits=NULL, type=c('ll', 'pop', 'pred'), .progress=TRUE,
                              ...) {
  tdmorefit <- fitted
  model <- tdmorefit$tdmore
  omegas <- diag(tdmorefit$tdmore$omega)
  profiledParameters <- model$parameters
  if(!is.null(fix)) profiledParameters <- profiledParameters[ !(profiledParameters %in% names(fix)) ]
  if(length(profiledParameters)==0) stop("No parameters to profile")

  limitsAsList <- !is.null(limits) && is.list(limits)
  if(limitsAsList) stopifnot(all(names(limits) %in% profiledParameters))

  limitsAsNumeric <- !is.null(limits) && is.numeric(limits)
  if(limitsAsNumeric) stopifnot(length(limits)==2)

  fun <- getLikelihoodFun(type)

  list <- lapply(
    profiledParameters,
    FUN = function(profiledParameter) {
      if (limitsAsNumeric) {
        return(seq(limits[1], limits[2], length.out=maxpts))

      } else if (limitsAsList){
        parameterRange <- unlist(limits[profiledParameter])
        if(!is.null(parameterRange)) {
          cat()
          return(seq(parameterRange[1], parameterRange[2], length.out=maxpts))
        }
      }
      # Default case, take omega value to generate a 95% range of possible eta's
      omega <- omegas[profiledParameter]
      return(seq(-1.96 * sqrt(omega), 1.96 * sqrt(omega), length.out = maxpts))
    }
  )
  grid <- expand.grid(list)
  colnames(grid) <- profiledParameters
  for(i in names(fix)) grid[,i] <- fix[i]
  grid <- grid[, model$parameters, drop=FALSE] # Reorder the columns

  cModel <- model
  cModel$cache <- model_prepare(model=model$model,
                                times=tdmorefit$observed$TIME,
                                regimen=tdmorefit$regimen,
                                parameters=coef(tdmorefit),
                                covariates=tdmorefit$covariates,
                                iov=model$iov,
                                extraArguments=model$extraArguments)
  tdmorefit$model <- cModel

  omega <- expandOmega(model, getMaxOccasion(tdmorefit$regimen))
  omega <- chol(omega) #performance improvement

  p <- to_progress(.progress)
  p$initialize(total=nrow(grid))

  profile <- apply(grid, 1, function(estimate) {
    p$tick()
    eta <- as.numeric(estimate)
    names(eta) <- model$parameters
    if(anyNA(eta)) {
      i <- is.na(eta)
      eta[i] <- 0
      eta <- estimate(tdmorefit, par=eta, fix=eta[!i])$res
    }
    c(eta,
      logLik=fun(par=eta,
                 omega=omega,
                 fix=NULL,
                 tdmore=cModel,
                 observed=tdmorefit$observed,
                 regimen=tdmorefit$regimen,
                 covariates=tdmorefit$covariates,
                 isChol=TRUE))
  })

  return(structure(
    list(
      profile = as.data.frame(t(profile)),
      profiledParameters = profiledParameters,
      tdmorefit = tdmorefit
    ),
    class = c("tdmoreprofile")
  ))
}



#' Plot a tdmoreprofile object.
#'
#' @param x the tdmoreprofile object
#' @param parameters the parameters to plot
#' @param contour use the 'geom_contour' method to draw the profile, logical value
#' @param raster use the 'raster method' to draw the profile, logical value
#' @param ... unused arguments
#' @importFrom ggplot2 autoplot ggplot aes_string geom_raster geom_contour geom_point geom_line labs
#'
#' @engine
autoplot.tdmoreprofile <- function(x, parameters=NULL, contour=T, raster=T, ...) {
  profile <- x
  tdmorefit <- profile$tdmorefit

  if(is.null(parameters)) {
    selectedParameters <- profile$profiledParameters
  } else {
    selectedParameters <- parameters
  }

  stopifnot(all(selectedParameters %in% colnames(profile$profile)))
  length <- length(selectedParameters)

  if(length > 2) {
    stop("Cannot produce plot in more than 2 dimensions... Use tdmore::shinyProfile() to explore the log-likelihood space instead")
  } else if(length == 2) {
    if(length(profile$profiledParameters) > 2) {
      warning(paste0("More than 2 profiled parameters, leading to multiple logLik-scores for a given parameter combination.\n",
           "Execute the profile() again with just 2 parameters, or use facets to ensure unique logLik scores per point."))
    }
    plot <- ggplot(data = profile$profile, aes_string(x=selectedParameters[1], y=selectedParameters[2], z="exp(logLik)"))
    if(raster) {plot <- plot + geom_raster(aes_string(fill = "exp(logLik)"))}
    if(contour){plot <- plot + geom_contour(colour="grey")}
    res <- coef(profile$tdmorefit)
    plot + geom_point(x=c(res[[selectedParameters[1]]]), y=c(res[[selectedParameters[2]]]), colour="white", inherit.aes=FALSE)

  } else if(length == 1) {
    notSelected <- setdiff(profile$profiledParameters, selectedParameters)
    group <- paste0("interaction(", paste(notSelected, collapse=","), ")")
    if(length(notSelected)==0) group <- NULL
    ggplot(data = profile$profile, aes_string(x=selectedParameters[1], y="exp(logLik)", group=group)) +
      geom_line()
  } else {
    stop("No plot can be produced. There should be at least 1 or 2 parameters to draw.")
  }
}

#' @export
plot.tdmoreprofile <- autoplot.tdmoreprofile
