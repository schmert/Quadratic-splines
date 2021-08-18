library(shiny)

ui <- fluidPage( 
  shiny::numericInput(inputId = 'n',label = 'n',value = 2),
  shiny::textOutput('nthValue'),
  shiny::textOutput('nthValueInv')   
)

fib <- function(n) {return(list(x=n, xx=n^2))}

server<-shinyServer(function(input, output, session) {   
  currentFib         <- reactive({ fib(as.numeric(input$n)) })  
  output$nthValue    <- renderText({ currentFib()$x })
  output$nthValueInv <- renderText({ currentFib()$xx })   
})

shinyApp(ui = ui, server = server)