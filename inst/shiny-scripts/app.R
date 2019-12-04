# A shiny app for hotgenes

library(shiny)
library(shinyWidgets)
library(hotgenes)

#DNA strand to be rendered start bp
dnaStart <- 1
#DNA strand to be rendered end bp
dnaEnd <- 195471971
#1 is most flexible, might speed up if use 1000; the scaling options must be divisible by this
baseScale <- 1000
dnaChr <- "chr1"

ui <- fluidPage(
  titlePanel("Hotgenes"),
  fluidRow(
    column(12,
      sliderInput("range", "Start-End bases", min = dnaStart, max = dnaEnd,
                  value = c(dnaStart,dnaEnd), width="99%")
    )
  ),
  fluidRow(
    column(4,
           sliderTextInput("scaling","Scaling",
                           #  choices=c(1, 100, 1000, 10000, 100000, 1000000, 1000000, 100000000),
                           choices=c(1000, 10000, 100000, 1000000, 1000000, 100000000),
                           selected=100000, grid = T)
           ),
    column(1,
           radioButtons("strand", "Strand", c("+", "-"))
           ),
    column(4,
           numericInput("iterations", "Number of iterations to simulate", value=10)
           ),
    column(3,
           numericInput("conductivity", "Heat 'conductivity'", min= 0, max = 1, value = 0.1, step = 0.05)
           )
  ),
  tags$b("Heated Maps"),
  tags$br(),
  tags$br(),
  plotOutput("heatedmap")
)

server <- function(input, output) {

  smpos <- reactive(
  {
    #use shiny Progress class as per https://shiny.rstudio.com/articles/progress.html
    progress <- shiny::Progress$new(min=0, max=length(musCh1fc[["counts"]][,1]))
    on.exit(progress$close())
    progress$set(value=0,  message="Loading positive strand...")
    #generate model and cache
    hotgenes::generateStrandModel(startBase=dnaStart, endBase=dnaEnd,
                                  fc=musCh1fc, chr=dnaChr, strand="+", scaling=baseScale,
                                  updateProgressBar=function(){progress$inc(1000)})
  })
  smneg <- reactive(
  {
    progress <- shiny::Progress$new(min=0, max=length(musCh1fc[["counts"]][,1]))
    on.exit(progress$close())
    progress$set(value=0, message="Loading negative strand...")
    #generate model and cache
    hotgenes::generateStrandModel(startBase=dnaStart, endBase=dnaEnd,
                                  fc=musCh1fc, chr=dnaChr, strand="-", scaling=baseScale,
                                  updateProgressBar=function(){progress$inc(1000)})

  })
  smr <- reactive(
    {
      if (input$strand == "+")
      {
        hotgenes::rebuildStrandModel(strandModel=smpos(), newStartBase=input$range[1], newEndBase=input$range[2], newScaling=input$scaling,
                                                              startBase=dnaStart, endBase=dnaEnd, scaling=baseScale)
      }
      else if (input$strand == "-")
      {
        hotgenes::rebuildStrandModel(strandModel=smneg(), newStartBase=input$range[1], newEndBase=input$range[2], newScaling=input$scaling,
                                     startBase=dnaStart, endBase=dnaEnd, scaling=baseScale)
      }
    })
  smt <- reactive(
    {
      hotgenes::simulateHeatSpread(smr(), input$conductivity, input$iterations)
    })
  locations <- reactive({
     hotgenes::generateLocationModel(startBase=input$range[1], endBase=input$range[2], scaling=input$scaling)
  })
  output$heatedmap <- renderPlot({
    hotgenes::plotHeatedMap(smt(), locations())


  })
}

shinyApp(ui = ui, server = server)
