## Given a shinyProfile interface
app <- ShinyDriver$new("../../", loadTimeout = 30000, seed=1234)
app$snapshotInit("test")

env <- new.env()
source("../../app.R", local = env)

expectUpdate <- function(...) {
  app$expectUpdate(...)
  # The progress bar has a fade-out transition of 250ms
  # It may also have a width transition of 600ms
  # To be sure, we use a sleep of 600ms
  #
  #https://github.com/rstudio/shiny/blob/44843a7768f899cf623297a3d4ae8cac5e0ac6f8/inst/www/shared/shiny.css#L179
  #https://github.com/rstudio/shiny/blob/c332c051f33fe325f6c2e75426daaabb6366d50a/inst/www/shared/bootstrap/css/bootstrap.css#L5146
  #See also: https://github.com/rstudio/shinytest/commit/29f3a23bcf97a102f805a369df1621e0aef26aaa
  Sys.sleep(0.6)
}

## Then by default the first parameter is profiled, and all others are fixed
Sys.sleep(5) ## first wait 5 seconds for rendering of plot to be complete...
app$snapshot(filename="EV1-fixed-ECL-profiled")
expect_equal(app$getValue("var"), list("ECL"))
expect_equal(app$getValue("EV1"), unname(coef(env$ipred)["EV1"]) )

## When I remove the fix for EV1
expectUpdate(timeout=30*1000, output="plot", EV1="")
## Then EV1 is optimized for every profiled ECL value
app$snapshot(filename="EV1-estimated-ECL-profiled")

## When I profile both ECL and EV1
expectUpdate(timeout=60*1000, output="plot", var=c("ECL", "EV1")) #long wait!
## Then a 2d-profile for EV1 and ECL is shown
app$snapshot(filename="EV1,ECL-profiled")

## When I profile only EV1
expectUpdate(timeout=10*1000, output="plot", var="EV1")
## Then EV1 is profiled in 1D, and ECL is fixed
app$snapshot(filename="EV1-profiled-ECL-fixed")
expect_equal(app$getValue("ECL"), unname(coef(env$ipred)["ECL"]))

## When I remove the fix for ECL
expectUpdate(timeout=30*1000, output="plot", ECL="")
## Then ECL is optimized for every profiled EV1 value
app$snapshot(filename="EV1-profile-ECL-estimated")

## When ECL is fixed to -5
expectUpdate(timeout=10*1000, output="plot", ECL=-5)
## Then EV1 is profiled for that ECL value
app$snapshot(filename="EV1-profile-ECL-fixedDiffValue")

### TODO: We do not test the plot zoom
p <- app$.__enclos_env__$private$shinyProcess
p$interrupt()
p$wait()

