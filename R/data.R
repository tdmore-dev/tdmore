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
