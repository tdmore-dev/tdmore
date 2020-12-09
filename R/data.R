detectDataType <- function(data) {
  if(all(c("time") %in% colnames(data)) ) "monolix"
  else if (all(c("TIME") %in% colnames(data)) ) "NONMEM"
  else stop("Could not detect data type based on column names; cannot automatically split data into observed/regimen/covariates")
}

loadObserved <- function(tdmore, data, type=detectDataType(data)) {
  if(type == "monolix") {
    data <- data %>%
      dplyr::select(.data$time, .data$observation)
    colnames(data) <- c("TIME", tdmore$res_var[[1]]$var )
    return(data)
  } else {
    stop("Data file type ", type, " not supported")
  }
}

loadRegimen <- function(tdmore, data, type=detectDataType(data)) {
  if(type == "monolix") {
    data %>%
      dplyr::select(.data$time, .data$amount) %>%
      dplyr::rename(TIME=.data$time, AMT=.data$amount)
  } else {
    stop("Data file type ", type, " not supported")
  }
}

loadCovariates <- function(tdmore, data, type=detectDataType(data)) {
  if(type == "monolix") {
    data %>% dplyr::select(.data$time, !! tdmore$covariates) %>% dplyr::rename(TIME=.data$time)
  } else {
    stop("Data file type ", type, " not supported")
  }
}

#' Utility function to create a data tibble.
#'
#' A data tibble is a tibble with an ID column, and separate columns
#' with the arguments for `estimate` function calls.
#' @param ... each of these arguments will be included in the resulting tibble.
#'
#' @return All arguments passed to the dataTibble function will be present in the resulting tibble as columns.
#' @export
dataTibble <- function(...) {
  qargs <- rlang::enquos(..., .named=TRUE)
  args <- list(...)
  names(args) <- names(qargs)

  hasId <- function(.x) !is.null(.x) && is.data.frame(.x) && "ID" %in% colnames(.x)
  argsWithId <- purrr::map(args, hasId) %>% unlist

  accId <- function(acc, .x){
    if(!hasId(.x)) return(acc)
    if(length(acc) == 0) return(unique(.x$ID))
    if( !setequal(acc, .x$ID) ) stop("ID vector not the same across all arguments: data-frame with columns ",
                                     paste(names(.x), collapse=","),
                                     " has new IDs ",
                                     paste(setdiff(acc, .x$ID), collapse=","))
    return(acc)
  }
  ID <- purrr::reduce( args, accId, .init=c())

  argsForTibble <- c(
    purrr::map(args[ !argsWithId ], ~list(.x)),
    purrr::map(args[argsWithId], ~split(.x[, colnames(.x)!="ID"], .x$ID))
  )
  myTibble <- c(
    if(is.null(ID)) NULL else list(ID=ID),
    argsForTibble[ names(args) ] #use same order as specified in arguments
  )
  tibble::as_tibble(myTibble)
}
