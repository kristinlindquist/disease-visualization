library(shiny);
source('./disease.R');

ui <- fluidPage(
   
   titlePanel("Disease Model"),
   
   sidebarLayout(
      sidebarPanel(
         numericInput("R0", "R0:", min = 0.1, max = 50, value = 2.0),
         numericInput("CFR", "Case fatality rate:", min = 0.0, max = 1.0, value = 0.30),
         sliderInput("Size", "Size:", min = 5, max = 100, value = 30),
         sliderInput("S0", "Population initially susceptible:", min = 0.0, max = 1.0, value = 0.9),
         sliderInput("initialInfected", "Infected @ t = 0:", min = 1, max = 25, value = 1)
      ),
      
      mainPanel(
         plotOutput("distPlot")
      )
   )
)

server <- function(input, output) {
   
   output$distPlot <- renderPlot({
      x    <- faithful[, 2] 
      bins <- seq(min(x), max(x), length.out = input$bins + 1)
      hist(x, breaks = bins, col = 'darkgray', border = 'white')
   })
}

shinyApp(ui = ui, server = server)

