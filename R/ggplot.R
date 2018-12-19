## Strategy:
## We use fortify to convert the tdmorefit to a special data.frame
## It remembers all arguments passed to fortify.tdmorefit
##
## Any stat_predict or other functions will generate a layer,
## but the layer will have a special extra GGTdmoreLayer class.
## When it is added to a plot using ggplot_add.GGTdmoreLayer,
## we use the original plot data to fill in missing arguments for the layer stat methods.
##
## There is one potential bug:
## If %+% is used, the new tdmorefit is not applied.
## To avoid this, we mark the resulting plot from 'ggplot_add'
## with a 'ggtdmore'
## The %+%.ggtdmore function will always generate an error.

#' @export
#' @importFrom ggplot2 fortify
fortify.tdmorefit <- function(model, newdata=NULL, regimen=NULL, covariates=NULL, parameters=NULL, ...) {
  observed <- model.frame(model, data=newdata)
  observed <- ggplot2::fortify(observed)
  attr(observed, 'tdmoreArgs') <- list(
    tdmorefit=model,
    regimen=regimen,
    covariates=covariates,
    parameters=parameters,
    extraArgs=list(...)
    )
  observed
}
has_tdmoreArgs <- function(x) {"tdmoreArgs" %in% names(attributes(x)) }

`%||%` <- function(x,y) {if(!is.null(x)) x else y}

#'@export
ggplot_add.GGTdmoreLayer <- function(object, plot, object_name) {
  class(object) <- setdiff(class(object), "GGTdmoreLayer") #unset class

  data <- plot$data
  layer <- object

  #params in layer: geom_params, stat_params, aes_params,
  list_union <- function(foo, bar){ c(foo, bar)[ union(names(foo), names(bar)) ] }
  if(has_tdmoreArgs(data)) {
    tdmoreArgs <- attr(data, 'tdmoreArgs')
    for(i in intersect(names(tdmoreArgs), names(layer$stat_params))) {
      value <- object$stat_params[[i]] %||% tdmoreArgs[[i]]
      #beware, assigning NULL to a list element removes that element!
      if(!is.null(value)) object$stat_params[[i]] <- value
    }
    object$stat_params$extraArgs <- tdmoreArgs$extraArgs
  }

  NextMethod()
}


#' Calculate predictions using the provided tdmore model. This is very similar to stat_function,
#' but the calculated values are calculated using `predict.tdmorefit`.
#'
#' It can be compared to `stat_function(fun=predict, object, newdata, regimen, parameters, covariates, ...)`.
#'
#' Not all arguments have to be provided. Missing arguments can be inherited from the original `ggplot()` call,
#' in a similar way as aesthetics.
#'
#' @param data
#' The data to be displayed in this layer. There are three options:
#' If NULL, data will be generated on a grid of evenly spaced values along the x axis.
#' A data.frame, or other object, will override the plot data.
#' A numeric vector will be coered to a data.frame with column TIME.
#' All objects will be fortified to produce a data frame.
#' See fortify() for which variables will be created.
#' A function will be called with a single argument, the plot data.
#' The return value must be a data.frame, and will be used as the layer data.
#'
#' This data is used in the `predict()` call.
#'
#' @inheritParams predict.tdmorefit
#' @inheritParams ggplot2::stat_function
#' @param tdmorefit if specified, overrides the fit object
#' @export
#'
#' @details
#' This is implemented as a `stat` function.
#'
#' If the model predicts an `y` or `x` output, those are passed to the geom instead.
#' This may be outside of the range of the scales.
#'
#' @importFrom ggplot2 layer
#' @examples
#' m1 <- tdmore(theopp_nlmixr)
#' pred <- estimate(m1, regimen = data.frame(TIME=0, AMT=5))
#' ggplot2::ggplot(pred, ggplot2::aes(TIME)) + ggplot2::geom_point() + stat_predict(data=seq(0, 24))
#'
#' ggplot2::ggplot(data.frame(TIME=c(0, 4), CONC=c(12,14)), mapping=ggplot2::aes(x=TIME, y=CONC)) +
#'     ggplot2::geom_point() +
#'     stat_predict(tdmorefit=pred, xlim=c(NA, NA))
stat_predict <- function(mapping=NULL, data=NULL, geom="line",
                         position="identity",
                         tdmorefit=NULL, regimen=NULL, parameters=NULL, covariates=NULL,
                         xlim=NULL, n=101,
                         se.fit=FALSE, level=0.95,
                         ...,
                         show.legend=NA, inherit.aes=TRUE) {
  # ggplot complains if we provide `data=seq(0, 10)`, as that is not a data.frame
  if(is.numeric(data)) {
    data <- data.frame(TIME=data)
  }
  z <- ggplot2::layer(
    data = data, mapping = mapping, stat = StatPredict,
    geom = geom, position = position, show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      xlim=xlim,
      n=n,
      tdmorefit=tdmorefit,
      regimen=regimen,
      parameters=parameters,
      covariates=covariates,
      ...
      )
  )
  #allow override of layer parameters by ggtdmore object in ggplot_add.GGTdmoreLayer
  class(z) <- union("GGTdmoreLayer", class(z))
  z
}

geom_fit <- function(mapping=NULL, data=NULL,
                         tdmorefit=NULL, regimen=NULL, parameters=NULL, covariates=NULL,
                         xlim=NULL, n=101,
                         ...,
                     color="red",
                         show.legend=NA, inherit.aes=TRUE) {
  stat_predict(mapping=mapping,
               data=data,
               geom="line",
               color=color,
               xlim=xlim,
               n=n,
               show.legend=show.legend,
               inherit.aes=inherit.aes)
}

## Ggplot_build does the following with the 'data' argument before we get to compute_statistic:
## layers <- plot$layers
## layer_data <- lapply(layers, function(y) y$layer_data(plot$data))
##    allow the layer to override `data` from the original plot, or filter using a function
## layout <- create_layout(plot$facet, plot$coordinates)
## data <- layout$setup(layer_data, plot$data, plot$plot_env)
##    data <- c(list(plot_data), data)
##    data <- self$facet$setup_data(data, self$facet_params)
##      Split the data.frame for every panel, and add a PANEL column.
##    data <- self$coord$setup_data(data, self$coord_params)
##      Never actually used by any default coord- in ggplot
## data <- by_layer(function(l, d) l$compute_aesthetics(d, plot))
##    Create a new data.frame with the computed values from the aesthetics.
##    If this is completely empty, or completely determined by stat(), then create a data.frame with length 1
## data <- lapply(data, scales_transform_df, scales = scales)
## layout$train_position(data, scale_x(), scale_y())
## data <- layout$map_position(data)
##    Execute the 'map' function on all x-y values
## data <- by_layer(function(l, d) l$compute_statistic(d, layout))
##    Finally, we compute the statistic
##
## We need to ensure that the data.frame is not empty at this point, otherwise
## setup_data never gets called by compute_statistic!
StatPredict <- ggplot2::ggproto("StatPredict", ggplot2::Stat,
  #default_aes = ggplot2::aes(x=TIME), #Does not work! This is applied only _after_ stats are calculated
  compute_group = function(data, scales,
                           tdmorefit, regimen, parameters, covariates, extraArgs = list(),
                           xlim = NULL, n = 101) {
    if(!is.null(xlim)) {
      # Evenly spaced grid, determined by xlim
      # NA is replaced by the axis limits
      range <- xlim

      if(anyNA(range)) {
        #Can we actually replace NA by the axis limits?
        if(is.null(scales$x)) stop("NA was used in xlim, but the x axis does not contain any real data...")
        if(scales$x$is_discrete()) stop("stat_predict requires a continuous (time) x axis")
        xrange <- scales$x$dimension()  # dimension() returns transformed values
        xrange <- scales$x$trans$inverse(xrange)
        range[is.na(range)] <- xrange[is.na(range)]
      }
      xseq <- seq(range[1], range[2], length.out = n)

      # re-initialize `data` with evenly spaced grid
      data <- data.frame(
        x=if(is.null(scales$x)) xseq else scales$x$trans$inverse(xseq)
      )
    } else {
      if("x" %in% names(data)) {
        xseq <- scales$x$trans$inverse( data$x ) #data$x is already transformed
      }
      else stop("No `x` value available for stat_predict. Either specify `xlim`, or provide an aesthetic that maps `x`.")
    }

    result <- do.call(
      predict,
      c(list(object=tdmorefit,
           newdata=xseq,
           regimen=regimen,
           covariates=covariates),
        extraArgs)
    )

    # Push new columns into `data`
    for(i in setdiff(names(result), names(data))) data[,i] <- result[,i]

    # Assign 'y' as first variable with residual error
    vars <- lapply(tdmorefit$tdmore$res_var, function(x) x$var)
    if(length(vars) > 0 && !"y" %in% colnames(result))
      data$y <- result[, vars[[1]] ] #first output column

    data
  }
)
