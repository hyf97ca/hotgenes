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

ui <- fluidPage(
  sliderInput("range", "Start-End bases",
              min = dnaStart, max = dnaEnd,
              value = c(dnaStart,dnaEnd)),
  sliderTextInput("scaling","Scaling",
#  choices=c(1, 100, 1000, 10000, 100000, 1000000, 1000000, 100000000),
  choices=c(1000, 10000, 100000, 1000000, 1000000, 100000000),
  selected=100000, grid = T),
  radioButtons("strand", "Strand (+/-)", c("+", "-")),
  numericInput("iterations", "Number of iterations to simulate", value=5),
  numericInput("conductivity", "Heat 'conductivity'", min= 0, max = 1, value = 0.001, step = 0.001),
  plotOutput("heatedmap")
)

server <- function(input, output) {

  sm <- reactive(
    {
      hotgenes::generateStrandModel(startBase=dnaStart, endBase=dnaEnd,
                                    fc=musCh1fc, chr="chr1", strand=input$strand, scaling=baseScale)
    })
  smr <- reactive(
    {
      hotgenes::rebuildStrandModel(strandModel=sm(), newStartBase=input$range[1], newEndBase=input$range[2], newScaling=input$scaling,
                                   startBase=dnaStart, endBase=dnaEnd, scaling=baseScale)
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
