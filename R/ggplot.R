#' Convert a tdmorefit object to a data.frame using `model.frame`.
#' The data.frame has a special 'tdmoreArgs' attribute that is used to capture all extra arguments to the original ggplot call.
#' The attribute contains a list with element 'tdmorefit' and all extra named arguments to the fortify call.
#'
#' @param model the tdmorefit object
#' @param ... extra arguments passed to `predict`
#'
#' @importFrom ggplot2 fortify
#' @export
#' @noRd
fortify.tdmorefit <- function(model, ...) {
  observed <- model.frame(model) %>% ggplot2::fortify()
  ## We ignore 'newdata' as an argument
  attr(observed, 'tdmoreArgs') <- list(tdmorefit=model, ...)
  observed
}

#' @importFrom ggplot2 fortify
#' @export
#' @noRd
fortify.tdmorefit_mixture <- fortify.tdmorefit

has_tdmoreArgs <- function(x) {"tdmoreArgs" %in% names(attributes(x)) }
`%||%` <- function(x,y) {if(!is.null(x)) x else y}


#' This creates a special type of ggplot2 Layer object.
#' The layer searches for a tdmorefit object in either the plot, or its own arguments.
#' It then calls the `predict()` function and provides this data.frame instead of the original data
#'
#' You can use this just like a regular `layer` call.
#' Beware that this function does not warn you about spurious or misspelled arguments, as
#' these will be passed along to `predict()` instead.
#'
#' @inheritParams ggplot2::layer
#' @param n Specifies whether the prediction must be performed on interpolated points along the x-axis.
#' @param ... extra arguments, which go partly into the geom/stat parameters, and partly into extra_params (which go into the `predict` call)
#'
#' @export
#' @examples
#' m1 <- tdmore(theopp_nlmixr)
#' db <- dataTibble(
#'      observed=data.frame(ID=c(1,2), TIME=5, CONC=c(10,15)),
#'      regimen=data.frame(TIME=0, AMT=5)
#' )
#' pred <- posthoc(object=m1, db)
#' ggplot2::ggplot(pred$fit[[1]], ggplot2::aes(x=TIME, y=CONC)) +
#'     ggplot2::geom_point() + predictionLayer(newdata=seq(0, 10, 0.1)) +
#'     predictionLayer(data=pred$fit[[2]], newdata=seq(0, 10, 0.1), color=2)
predictionLayer <- function(mapping=NULL, data=NULL, geom="line",
                            stat="identity",
                            position="identity",
                            ...,
                            n=101,
                            show.legend=NA, inherit.aes=TRUE) {
  params <- list(
    ...
  )

  geom <- ggplot2:::check_subclass(geom, "Geom", env = parent.frame())
  stat <- ggplot2:::check_subclass(stat, "Stat", env = parent.frame())
  params <- ggplot2:::rename_aes(params) #standardisze color and colour
  aes_params  <- params[intersect(names(params), geom$aesthetics())]
  geom_params <- params[intersect(names(params), geom$parameters(TRUE))]
  stat_params <- params[intersect(names(params), stat$parameters(TRUE))]
  all <- c("key_glyph", geom$parameters(TRUE), stat$parameters(TRUE), geom$aesthetics())

  isLayerParam <- names(params) %in% all
  layer_params <- params[isLayerParam]
  extra_param <- params[!isLayerParam]

  z <- ggplot2::layer(
    data = data, mapping = mapping,
    stat = stat,
    geom = geom, position = position, show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = layer_params,
    layer_class=PredictionLayer
  )
  z$extra_param <- extra_param
  if(is.null(n)) {
    #OK
  } else {
    #interpolation
    if(!is.numeric(n) || length(n) < 1) stop("the `n` argument should either be a single number, or a numeric vector")
  }
  z$interpolation_n <- n
  z
}

PredictionLayer <- ggplot2::ggproto(
  "PredictionLayer",
  ggplot2:::Layer,
  extra_param = list(), #extra parameters passed
  interpolation_n = NULL, #interpolate?
  layer_data = function(self, plot_data) {
    super <- ggplot2::ggproto_parent(ggplot2:::Layer, self)
    layerData <- super$layer_data(plot_data)

    params <- self$extra_param

    if(has_tdmoreArgs(layerData)) {
      extra <- attr(layerData, 'tdmoreArgs')
      i <- setdiff(names(extra), names(params))
      params[i] <- extra[i]
    }
    if(has_tdmoreArgs(plot_data)) {
      extra <- attr(plot_data, 'tdmoreArgs')
      i <- setdiff(names(extra), names(params))
      params[i] <- extra[i]
    }

    if(is.null(params$object)) params$object <- params$tdmorefit #ensure the 'object' parameter is defined
    if(is.null(params$object)) stop("Could not find a tdmorefit object as 'tdmorefit' or 'object' argument")

    self$params <- params
    self$predict <- function(...) {
      params <- c(self$params, list(...) ) #give original parameters preference
      params <- params[ unique(names(params)) ]
      do.call(what=stats::predict, args=params)
    }
    data <- self$predict(newdata=layerData)  ## Already pass a prediction, even though it may be a prediction no data
    if(nrow(data)==0) {
      #Ggplot handles empty data poorly. We add dummy NA data instead. This will get replaced by actual data later
      dummy <- as.list( rep(NA_real_, length.out=ncol(data)) )
      names(dummy) <- colnames(data)
      data <- rbind(data, dummy)
    }
    as.data.frame(data)
  },
  setup_layer = function(self, data, plot) {
    self$plot <- plot ## Store the plot for later use
    data
  },
  compute_statistic = function(self, data, layout) {
    ## First interpolate the data-frame if required
    n <- self$interpolation_n
    calcPanel <- function(panel) {
      # per panel
      data <- data[data$PANEL == panel, ]
      if(length(n)==1) {
        scales <- layout$get_scales(data$PANEL[1])
        xLimits <- c(scales$x$dimension(), layout$coord$limits$x)
        xLimits[is.infinite(xLimits)] <- NA
        range <- range(xLimits, na.rm=TRUE )
        if(any(!is.finite(range))) stop("Prediction layer does not know where to interpolate.\nPlease specify an x-coordinate range for this plot using e.g. coord_cartesian(), or add some actual data.")
        xseq <- seq(range[1], range[2], length.out=n)
        events <- getEvents(self$params$tdmorefit, self$params$object, self$params$regimen, self$params$covariates)
        xseq <- includeEvents(xseq, events)
      } else if (is.numeric(n)) {
        xseq <- n
      } else if (is.null(n)) {
        xseq <- data$x
        if(all(is.na(xseq))) stop("`n`=NULL specified, but there is no x data to use!")
      } else {
        stop("Argument `n` is unsupported: ", n)
      }
      result <- self$predict(newdata=xseq)
      result$PANEL <- data$PANEL[1]

      data <- self$compute_aesthetics(result, self$plot)

      ## For aesthetics mapped to NULL, we will auto-guess their mapping here
      vars <- c("x", "y", "ymin", "ymax")
      missing <- setdiff(vars, colnames(data)) %>%
        getAes(self$params$tdmorefit, vars=., as.aes=FALSE, data=result)
      missing <- purrr::compact(missing)
      if(length(missing)>0) data[, names(missing)] <- result[, unlist(missing) ]

      data <- ggplot2:::scales_transform_df(self$plot$scales, data)

      for(i in setdiff(names(result), names(data))) data[,i] <- result[,i] ## add missing columns

      data
    }
    result <- lapply(unique(data$PANEL), calcPanel)
    data <- do.call(dplyr::bind_rows, args=result)

    super <- ggplot2::ggproto_parent(ggplot2:::Layer, self)
    data <- super$compute_statistic(data, layout)
    as.data.frame(data)
  }
)

# Includes the regimen, observed and covariate times of the tdmorefit into
# the xseq sequence, as well as an event just before and just after.
includeEvents <- function(xseq, events) {
  range <- range(xseq)
  events <- c(events, modifyMantissa(events, 2^-52), modifyMantissa(events, -2^-52))
  events <- events[ events > range[1] & events < range[2] ]
  xseq <- unique( sort( c(xseq, events) ) )
}

getEvents <- function(...) {
  result <- lapply(list(...), function(i) {
    if(is.data.frame(i) & "TIME" %in% colnames(i)) return(i$TIME)
    if(is.tdmorefit(i)) return(getEvents(i$regimen, i$observed, i$covariates))
    NULL
  })
  unlist(result)
}

#' This function tries to construct an aes
#' by guessing the required variables
#'
#' If no 'x' is defined, we put a NULL in the aesthetic
#' It will then get filled in by the predictionLayer
#' @param x tdmorefit object
#' @param vars the aes vars to guess
#' @param as.aes TRUE to return an aes object, FALSE to return a named list
#' @param data data.frame to check whether the values actually exist
#' @noRd
getAes <- function(x, vars=c("x", "y", "ymin", "ymax"), as.aes=TRUE, data = NULL) {
  if(!is.tdmorefit(x)) stop()
  varNames <- lapply(x$tdmore$res_var, function(x){x$var}) %>% unlist

  if(length(varNames)==0) {
    var <- NULL
  } else if (is.null(data)) {
    var <- varNames[1]
  } else {
    ## Consult data to only suggest available column names
    varNames <- intersect(varNames, colnames(data) )
    if(length(varNames)==0) var <- NULL #nothing left; too bad!
    else var <- varNames[1]
  }

  mapping <- lapply(vars, function(i) {
    switch (i,
            x = "TIME",
            y = var,
            ymin = if(is.null(var)) NULL else paste0(var, ".lower"),
            ymax = if(is.null(var)) NULL else paste0(var, ".upper"),
            stop("Cannot guess aes for aesthetic ", i)
    )
  })
  names(mapping) <- vars
  mapping <- purrr::keep(mapping, function(i) {
    i %in% colnames(data) #remove columns that are not available in source dataset
  })

  if(as.aes) mapping <- do.call( ggplot2::aes_string, args=mapping)

  return(mapping)
}
