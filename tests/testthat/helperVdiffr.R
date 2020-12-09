## Goal: fix bug https://github.com/r-lib/vdiffr/issues/96

expect_doppelganger <- function(title, fig, path = NULL, ...) {
  vdiffr::expect_doppelganger(title, fig, path="", ...)
}
