## ShinyProfile runs a shiny gadget for a given profile

to_plyr <- function(x) {
  ## x is a shiny::Progress object
  N <- 1
  list(
    init=function(N) {N <<- N},
    step=function() {x$inc(1/N, detail="Rendering...")},
    term=function() {}
  )
}

#' Open a Shiny Gadget to explore the log-likelihood profile
#' of a specific fit
#' @inheritParams shinyProfileApp
#' @export
shinyProfile <- function(fitted, fix=NULL, ...) {
  app <- shinyProfileApp(fitted, fix, ...)
  shiny::runGadget(app)
}

#' Build a Shiny app to explore the log-likelihood profile
#' of a specific fit
#' @inheritParams profile.tdmorefit
#' @param ... Extra arguments passed to `profile`
#' @export
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
    shiny::conditionalPanel( condition = paste0("! input.var.includes('",n,"')"),
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

    output$plot <- shiny::renderPlot({
      plotVars <- input$var
      limits=list(
        c(prefs$zoom$xmin, prefs$zoom$xmax),
        c(prefs$zoom$ymin, prefs$zoom$ymax)
      )
      if(length(plotVars) == 1) limits=limits[1]
      names(limits) <- plotVars

      progress <- shiny::Progress$new()
      on.exit(progress$close())
      myProfile <- stats::profile(fitted=tdmorefit,
                           fix=selected()[! names %in% plotVars],
                           limits=limits,
                           .progress=to_plyr(progress),
                           ...)
      plot(myProfile)
    })

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
