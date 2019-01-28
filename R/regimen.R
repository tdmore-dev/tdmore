#' Transform a specific treatment regimen using a modifier
#'
#' @param regimen a treatment regimen
#' @param x A row of a possibilitygrid
#'
#' @return The modified treatment regimen
#' @export
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
#' @export
#'
#' @examples flatten(data.frame(TIME=10, AMT=10, ADDL=8, II=24))
flatten <- function(regimen) {
  newRegimen <- data.frame()
  for(i in seq_len(nrow(regimen))) {
    tmt <- data.frame( regimen[i, ] )
    addl <- tmt$ADDL
    if(is.null(addl)) addl <- 0
    tmt$ADDL <- 0
    newRegimen <- rbind(newRegimen, tmt)
    for(i in seq_len(addl)) {
      tmt$TIME <- tmt$TIME + tmt$II
      newRegimen <- rbind(newRegimen, tmt)
    }
  }
  newRegimen[, ! colnames(newRegimen) %in% c("ADDL", "II")] # Without columns ADDL and II
}

#' Get the times related to each occasion in a regimen.
#'
#' @param regimen the specified regimen
#' @return a numeric vector with the times
#'
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
#'
hasOccasion <- function(regimen) {
  return("OCC" %in% colnames(regimen))
}

#' Get the highest occasion number in a regimen.
#' If the regimen does not contain an 'OCC' column, 1 is returned.
#'
#' @param regimen the specified regimen
#' @return the highest occasion number
#'
getMaxOccasion <- function(regimen) {
  if(hasOccasion(regimen)) {
    retValue <- max(regimen$OCC)
  } else {
    retValue <- 1
  }
  return(retValue)
}
