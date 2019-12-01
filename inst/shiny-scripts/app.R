# A shiny app for hotgenes

library(shiny)
library(shinyWidgets)
library(hotgenes)

ui <- fluidPage(
  sliderInput("range", "Start-End bases",
              min = 1, max = 195471971,
              value = c(1,195471971)),
  sliderTextInput("scaling","Scaling",
  choices=c(1, 100, 1000, 10000, 100000, 1000000, 1000000, 100000000),
  selected=100000, grid = T),
  radioButtons("strand", "Strand (+/-)", c("+", "-")),
  numericInput("iterations", "Number of iterations to simulate", value=5),
  numericInput("conductivity", "Heat 'conductivity'", min= 0, max = 1, value = 0.001),
  plotOutput("heatedmap")
)

server <- function(input, output) {

  sm <- reactive(
    {
      hotgenes::generateStrandModel(startBase=input$range[1], endBase=input$range[2],
                                    fc=musCh1fc, chr="chr1", strand=input$strand, scaling=input$scaling)
    })
  sm1 <- reactive(
    {
      hotgenes::simulateHeatSpread(sm(), input$conductivity, input$iterations)
    })
  locations <- reactive({
     hotgenes::generateLocationModel(startBase=input$range[1], endBase=input$range[2], scaling=input$scaling)
  })
  output$heatedmap <- renderPlot({
    hotgenes::plotHeatedMap(sm1(), locations())


  })
}

shinyApp(ui = ui, server = server)
