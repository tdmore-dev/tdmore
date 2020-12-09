## ShinyProfile runs a shiny gadget for a given profile

#' Open a Shiny Gadget to explore the log-likelihood profile
#' of a specific fit
#' @inheritParams shinyProfileApp
#' @engine
shinyProfile <- function(fitted, fix=NULL, ...) {
  app <- shinyProfileApp(fitted, fix, ...)
  shiny::runGadget(app)
}

#' Build a Shiny app to explore the log-likelihood profile
#' of a specific fit
#' @inheritParams profile.tdmorefit
#' @param ... Extra arguments passed to `profile`
#' @noRd
shinyProfileApp <- function(fitted, fix=NULL, ...) {
  tdmorefit <- fitted
  if (!requireNamespace("shiny", quietly = TRUE)) {
    warning("The shiny package must be installed to use this functionality")
    return(NULL)
  }
  if (!requireNamespace("miniUI", quietly = TRUE)) {
    warning("The miniUI package must be installed to use this functionality")
    return(NULL)
  }
  names <- names(coef(tdmorefit))
  omega <- diag(tdmorefit$tdmore$omega)
  init <- coef(tdmorefit)
  initNames <- names[1]
  if(!is.null(fix)) {
    init[names(fix)] <- fix
    initNames <- setdiff(names, names(fix))
  }
  allButtons <- lapply(names, function(n){
    shiny::conditionalPanel( condition = paste0("input.var.indexOf('",n,"') == -1"),
      shiny::numericInput(inputId=n, label=n, value=init[n],
                        step=sqrt(omega[n])/10)
    )
  })
  allButtons <- c(allButtons, list(
    shiny::selectInput("var", label=NULL, choices=names, multiple=TRUE,
                       selected=initNames)
  ))
  ui <- miniUI::miniPage(
    do.call(miniUI::miniButtonBlock, allButtons),
    miniUI::miniContentPanel(
      shiny::plotOutput("plot",
                        height = "100%",
                        brush=shiny::brushOpts(
                          id="brush",
                          resetOnNew = TRUE
                          ), click = "click", dblclick="dblclick")
    )
  )

  server <- function(input, output, session) {
    prefs <- shiny::reactiveValues(
      zoom=list(xmin=-1, xmax=1, ymin=-1, ymax=1)
    )
    selected <- shiny::reactive({
      res <- sapply(names, function(n){ input[[n]] })
      names(res) <- names
      res
    })

    myProfile <- shiny::reactive({
      plotVars <- input$var
      limits=list(
        c(prefs$zoom$xmin, prefs$zoom$xmax),
        c(prefs$zoom$ymin, prefs$zoom$ymax)
      )
      if(length(plotVars) == 1) limits=limits[1]
      if(length(plotVars) == 3) stop("Cannot explore more than 2 variables at the same time")
      names(limits) <- plotVars

      progress <- shiny::Progress$new()
      on.exit(progress$close())
      progress$set(message="Calculating")
      profile <- stats::profile(fitted=tdmorefit,
                                  fix=selected()[! names %in% plotVars],
                                  limits=limits,
                                  .progress=progress,
                                  ...)
      profile
    })
    shiny::exportTestValues( profile = {
      foo <- myProfile()
      foo$tdmorefit <- NULL
      foo$profile <- signif( foo$profile, digits=6 ) #round, so tests always get same result
      foo
      })
    output$plot <- shiny::snapshotExclude(shiny::renderPlot({
      autoplot(myProfile())
    }))

    shiny::observeEvent(input$dblclick, {
      ## zoom in, or zoom out
      if(is.null(input$brush)) {
      width <- prefs$zoom$xmax - prefs$zoom$xmin
      height <- prefs$zoom$ymax - prefs$zoom$ymin
      width <- width * 2
      height <- height * 2
      prefs$zoom <- list(
        xmin= prefs$zoom$xmin -width/2,
        xmax= prefs$zoom$xmax + width/2,
        ymin= prefs$zoom$ymin -height/2,
        ymax= prefs$zoom$ymax + height/2
      )
      } else {
        prefs$zoom <- list(
          xmin=input$brush$xmin,
          xmax=input$brush$xmax,
          ymin=input$brush$ymin,
          ymax=input$brush$ymax
        )
      }
    })
    shiny::observeEvent(input$done, {
      shiny::stopApp(NULL)
    })
  }
  shiny::shinyApp(ui=ui, server=server)
}
