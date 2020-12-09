#' Transform a specific treatment regimen using a modifier
#'
#' @param regimen a treatment regimen
#' @param x A row of a possibilitygrid
#'
#' @return The modified treatment regimen
#' @noRd
transformRegimen <- function(regimen, x) {
  for(j in names(x)) {
    rIndex <- as.numeric(strsplit(j, "\\.")[[1]][1]) ## Which row of the regimen
    rCol <- sub("^\\d+\\.", "", j) ## Which column
    regimen[rIndex, rCol] <- x[j]
  }
  regimen
}

#' Flatten a treatment regimen
#'
#' @param regimen data.frame with treatment regimen
#'
#' @return The equivalent treatment regimen with all ADDL objects moved to 0
#' @noRd
#'
#' @examples flatten(data.frame(TIME=10, AMT=10, ADDL=8, II=24))
flatten <- function(regimen) {
  if(! "ADDL" %in% colnames(regimen)) return(regimen)
  if(! "II" %in% colnames(regimen)) stop("ADDL column detected, but II not present")
  myList <- apply(regimen, 1, function(row) {
    row <- as.list(row)
    row$TIME <- row$TIME + seq(0, row$ADDL)*row$II
    row$ADDL <- NULL
    row$II <- NULL
    row
  })
  dplyr::bind_rows( lapply(myList, as.data.frame) )
}

#' Get the times related to each occasion in a regimen.
#'
#' @param regimen the specified regimen
#' @return a numeric vector with the times
#' @noRd
getOccasionTimes <- function(regimen) {
  if(hasOccasion(regimen)) {
    temp <- regimen[!duplicated(regimen$OCC),]
    retValue <- temp$TIME
  } else {
    retValue <- NULL
  }
  return(retValue)
}

#' Tell if a regimen has an occasion column.
#'
#' @param regimen the specified regimen
#' @return TRUE if the regimen has a 'OCC' column, false otherwise
#' @noRd
hasOccasion <- function(regimen) {
  return("OCC" %in% colnames(regimen))
}

#' Get the highest occasion number in a regimen.
#' If the regimen does not contain an 'OCC' column, the number of rows in the regimen is returned (each regimen a separate occasion)
#'
#' @param regimen the specified regimen
#' @return the highest occasion number
#' @noRd
getMaxOccasion <- function(regimen) {
  if(hasOccasion(regimen)) {
    retValue <- max(regimen$OCC)
  } else {
    retValue <- nrow(regimen)
  }
  return(retValue)
}
