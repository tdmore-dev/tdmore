# Defines the search space for dose finding
# This search space is added as a column to the regimen table, identifying which regimen attributes can be modified

#' Define a blank search space
#'
#' @return a searchSpace object
#' @export
searchspace <- function() {
  structure( list(
    options=data.frame()
  ), class="searchspace")
}

#' Add discrete options to a search space
#'
#' @param searchspace Defined search space
#' @param options data.frame with discrete options
#'
#' @return A search space with the discrete options in `options` added
#' @export
discrete.searchspace <- function(searchspace, options=NULL) {
  assert_that(is.data.frame(options))
  assert_that(nrow(options) > 0)
  assert_that(ncol(options) > 0)

  if(nrow(searchspace$options) == 0) {
    searchspace$options <- options
    return(searchspace)
  }

  result <- data.frame()
  grid <- expand.grid(seq_len(nrow(searchspace$options)), seq_len(nrow(options)))
  for(i in seq_len(nrow(grid))) {
    row1 <- searchspace$options[ grid[i, 1], , drop=FALSE]
    row2 <- options[ grid[i, 2], ,drop=FALSE]
    result <- rbind(result, c(row1, row2))
  }
  searchspace$options <- result
  searchspace
}

#' Calculate all possible regimen combinations, as rows in a data.frame
#'
#' @param searchspace searchspace object
#'
#' @return A data.frame with e.g. columns 1.Foo, 1.Bar 2.Foo,
#' with the rows containing all discrete combinations for the Foo value in the first row,
#' the Bar value in the first row and the Foo value in the second row.
#' @export
possibilitygrid <- function(searchspace) {
  assert_that(inherits(searchspace, "list"))
  sapply(searchspace, function(x) {assert_that(is.null(x) || inherits(x, "searchspace"))})

  grid <- NULL
  for(i in seq_along(searchspace)) {
    ss <- searchspace[[i]]
    if(is.null(ss)) next
    thisGrid <- ss$options
    names(thisGrid) <- paste0(i, ".", names(thisGrid))
    grid <- if(is.null(grid)) { thisGrid } else { expand.grid(grid, thisGrid) }
  }
  grid
}
