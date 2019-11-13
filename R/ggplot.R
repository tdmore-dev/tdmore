## Strategy:
## We use fortify to convert the tdmorefit to a special data.frame
## It remembers all arguments passed to fortify.tdmorefit
##
## Any stat_predict or other functions will generate a layer with a special layer_class
## This layer_class can change the data.frame to the predicted values instead

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
#' @inheritParams ggplot2::stat_function
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
                         xlim=NULL, n=101,
                         se.fit=FALSE, level=0.95,
                         ...,
                         show.legend=NA, inherit.aes=TRUE) {
  params <- list(
    xlim=xlim,
    n=n,
    se.fit=se.fit, level=level,
    ...
  )
  z <- ggplot2::layer(
    data = data, mapping = mapping, stat = StatPredict,
    geom = geom, position = position, show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = params,
    layer_class=TDMoreInheritLayer
  )
  z
}

StatPredict <- ggplot2::ggproto("StatPredict", ggplot2::Stat,
  default_aes = ggplot2::aes(x=stat(TIME), y=stat(CONC)), #Does not work! This is applied only _after_ stats are calculated
  compute_group = function(data, scales,
                           tdmorefit, regimen=NULL, parameters=NULL, covariates=NULL,
                           se.fit=FALSE, level=0.95, mc.maxpts=100,
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
           parameters=parameters,
           covariates=covariates,
           se.fit=se.fit,
           level=level,
           mc.maxpts=mc.maxpts))
    )

    # Push new columns into `data`
    for(i in setdiff(names(result), names(data))) data[,i] <- result[,i]

    # FIXME: if the aesthetics contained 'aes(y=CONC)', then it is a pity this is not maintained...


    # FIXME
    # Assign 'y' as first variable with residual error
    #vars <- lapply(tdmorefit$tdmore$res_var, function(x) x$var)
    #vars <- unname(vars)
    #if(length(vars) > 0 && !"y" %in% colnames(result)) {
    #  data$y <- result[, vars[[1]] ] #first output column
    #}

    data
  }
)

## This layer takes care of inheriting the provided parameters of ggplot() into the arguments of geom and stat
TDMoreInheritLayer <- ggplot2::ggproto(
  "TDMoreInheritLayer",
  ggplot2:::Layer,
  setup_layer = function(self, data, plot) {
    plotData <- plot$data
    if(has_tdmoreArgs(plotData)) { ## if it is a tdmore data object, it can be used to override missing parameters
      tdmoreArgs <- attr(plotData, 'tdmoreArgs')

      # a layer has geom_params and stat_params
      # we can nicely fill missing values using the tdmoreArgs
      i <- intersect(names(tdmoreArgs), self$geom$parameters(extra=TRUE) )
      self$geom_params[i] <- tdmoreArgs[i]

      i <- intersect(names(tdmoreArgs), self$stat$parameters(extra=TRUE) )
      self$stat_params[i] <- tdmoreArgs[i]
    }

    data
  }
)

##### ----- The below code is experimental, and relies heavily on ggplot2 internals! -----
## This layer replaces the data with a predicted data.frame
## In layer_data, we set up the tdmoreArgs frame
TDMorePredictLayer <- ggplot2::ggproto(
  "TDMorePredictLayer",
  ggplot2:::Layer,
  tdmoreArgs=NULL, ## cache parameter
  layer_data = function(self, plot_data) {
    super <- ggproto_parent(ggplot2:::Layer, self)
    layerData <- super$layer_data(plot_data)

    ## First try self$data. Maybe it has some 'tdmoreArgs'?
    params <- list()
    if(has_tdmoreArgs(self$data)) {
      tdmoreArgs <- attr(self$data, 'tdmoreArgs')
      ## What is not yet specified comes from tdmoreArgs
      i <- setdiff(names(tdmoreArgs), names(params))
      params[i] <- tdmoreArgs[i]
    } else if (!is.waive(self$data)) {
      ### It is not a waiver; someone really specified something!
      ### Use it as newdata argument
      if(!empty(layerData) && is.null(params$newdata)) params$newdata <- layerData
    }

    ## Also check the plot_data!
    if(has_tdmoreArgs(plot_data)) {
      tdmoreArgs <- attr(plot_data, 'tdmoreArgs')
      ## What is not yet specified comes from tdmoreArgs
      i <- setdiff(names(tdmoreArgs), names(params))
      params[i] <- tdmoreArgs[i]
    }

    if(is.null(params$object)) params$object <- params$tdmorefit #ensure the 'object' parameter is defined
    self$tdmoreArgs <- params

    do.call(what=predict, args=params) ## Already pass a prediction, even though it may be a prediction without observed data
  },
  setup_layer = function(self, data, plot) {
    self$plot <- plot ## Store the plot for later use
    data
  },
  compute_statistic = function(self, data, layout) {
    ### Now calculate the data.frame
    if(is.null(self$tdmoreArgs)) stop()
    params <- self$tdmoreArgs
    browser()

    data <- plyr::ddply(data, "PANEL", function(data) {
      if(is.null(params$newdata)) {
        scales <- layout$get_scales(data$PANEL[1])
        range <- range( c(scales$x$dimension(), layout$coord$limits$x) )
        params$newdata <- seq(range[1], range[2], length.out=101)
      }
      result <- do.call(what=predict, args=params)
      result$PANEL <- data$PANEL[1]

      browser()
      data <- self$compute_aesthetics(result, self$plot) ## re-compute previous aes
      data <- scales_transform_df(scales, data)
      data <- layout$map_position(data)

      for(i in setdiff(names(result), names(data))) data[,i] <- result[,i] ## add missing columns

      data
    })

    super <- ggproto_parent(ggplot2:::Layer, self)
    data <- super$compute_statistic(data, layout)
    data
  }
)

#' A geom to display a fit object as a line
#' @export
geom_fit <- function(mapping=NULL, data=NULL,
                     stat="identity",
                     position="identity",
                     ...,
                     se.fit=FALSE, level=0.95,
                     color="tomato1",
                     show.legend=NA, inherit.aes=TRUE) {
  ## Add the tdmoreArgs to the list
  tdmoreArgs <- list(
    se.fit=se.fit, level=level
  )
  data <- fortify(data)
  attr(data, 'tdmoreArgs') <- c( tdmoreArgs, attr(data, 'tdmoreArgs'))

  params <- list(
    color=color,
    ...
  )

  z <- ggplot2::layer(
    data = data, mapping = mapping, stat = stat,
    geom = "line", position = position, show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = params,
    layer_class=TDMorePredictLayer
  )
  z
}

