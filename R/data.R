detectDataType <- function(data) {
  if(all(c("time") %in% colnames(data)) ) "monolix"
  else if (all(c("TIME") %in% colnames(data)) ) "NONMEM"
  else stop("Could not detect data type based on column names; cannot automatically split data into observed/regimen/covariates")
}

loadObserved <- function(tdmore, data, type=detectDataType(data)) {
  if(type == "monolix") {
    data <- data %>% select(time, observation)
    colnames(data) <- c("TIME", tdmore$res_var[[1]]$var )
    return(data)
  } else {
    stop("Data file type ", type, " not supported")
  }
}

loadRegimen <- function(tdmore, data, type=detectDataType(data)) {
  if(type == "monolix") {
    data %>% select(time, amount) %>% rename(TIME=time, AMT=amount)
  } else {
    stop("Data file type ", type, " not supported")
  }
}

loadCovariates <- function(tdmore, data, type=detectDataType(data)) {
  if(type == "monolix") {
    data %>% select(time, !! tdmore$covariates) %>% rename(TIME=time)
  } else {
    stop("Data file type ", type, " not supported")
  }
}
