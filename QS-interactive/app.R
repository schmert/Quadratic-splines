library(shiny)
library(tidyverse)

# Define UI for application that draws a histogram ----
ui <- fluidPage(

    # Application title
    titlePanel("Quadratic Spline Fertility Model: Illustration"),
    titlePanel(h6("This app illustrates the Quadratic Spline (QS) model proposed in Schmertmann (2003)")),
    titlePanel(h6("Parameters \u03b1, P, H, and R control the shape and level.")),
    titlePanel(h6("Colored sections illustrate the piecewise quadratic curves that comprise the spline schedule")),
    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(width=2,
            h4('SHAPE'),         
            numericInput(inputId = "alpha",  
                         HTML("&alpha;<br>Initial Age"), value = 14,
                         width='75%'),
            numericInput(inputId = "P",  
                         HTML("P<br>Peak Age"), value = 29,
                         width='75%'),
            numericInput(inputId = "H",  
                         HTML("H<br>Half-Peak Age"), value = 36,
                         width='75%'),
            h4('LEVEL'),
            numericInput(inputId = "R",  
                         HTML("R<br>Peak Level"), value = 200,
                         width='75%',
                         min=50, max=350, step=10)
        ),

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("fx_plot")
        )
    )
)

# Define server logic required to draw a histogram ----
server <- function(input, output) {

    QS = function(x, R,alpha,P,H) {
        
        w = min(.75, .25+.025*(P-alpha))
        beta = max (  min(50, 4*H-3*P) , (4*H-P)/3 )
        
        D = P-20          # delay index
        C = (P+50)/2 - H  # control index
        
        knot    = vector("numeric",5)
        knot[1] = alpha
        knot[2] = alpha + w*(P-alpha)
        knot[3] = P
        knot[4] = (P+H)/2
        knot[5] = (H+beta)/2
        
        A      = matrix(NA,5,5)
        target = vector("numeric",5)
        
        # target value at P=1
        A[1,]     = (pmax(P-knot, 0))^2
        target[1] = 1
        
        # target value at H=1/2
        A[2,]     = (pmax(H-knot, 0))^2
        target[2] = 1/2
        
        # target value at beta=0
        A[3,]     = (pmax(beta-knot, 0))^2
        target[3] = 0
        
        # target slope at P=0
        A[4,]     = 2*pmax(P-knot, 0) 
        target[4] = 0
        
        # target slope at beta=0
        A[5,]     = 2*pmax(beta-knot, 0) 
        target[5] = 0
        
        # calculate thetas
        theta = solve(A,target)
        
        tmp = (pmax( outer(x,knot,"-") , 0))^2
        
        y = (x >= alpha) * (x <= beta) * R * (tmp %*% theta)  # ASFRs
        
        return(list(R=R,alpha=alpha,P=P,H=H,
                    theta=theta,knot=knot,w=w,beta=beta,delay=D,control=C,
                    x=as.vector(x) ,y=as.vector(y) ))
        
    } #QS
  
  
    output$fx_plot <- renderPlot({

        x       = seq(10,55,.25)
        this_QS = QS(x, input$R, input$alpha, input$P, input$H)
        fx      = this_QS$y
      
        dividers = c(-Inf,this_QS$knot,this_QS$beta,Inf)
        model = tibble(x,fx, section=cut(x,breaks=dividers))
        
        G = ggplot(data=model) + 
            aes(x=x,y=fx, color=section) +
            geom_line(size=5, alpha=.60) +
            geom_line(aes(x=x,y=fx),size=1.5, color='black') +
            labs(x='Age',y='Fertility Rate per 1000') +
            theme_bw() +
            theme(axis.text  = element_text(size=14,face='bold'),
                  axis.title = element_text(size=14,face='bold')) +
            geom_point(x=input$P, y=0, shape=16, size=4,color='black') +
            geom_text(x=input$P+1,y=10, label='P', size=6,color='black') +
            geom_segment(x=input$P, y=0, xend=input$P, yend = input$R, size=1,color='black' ) +
            geom_point(x=input$H, y=0, shape=16, size=4,color='black') +
            geom_text(x=input$H+1,y=10, label='H', size=6,color='black') +
            geom_segment(x=input$H, y=0, xend=input$H, yend = input$R/2, size=1,color='black' ) +
            geom_point(x=input$alpha, y=0, shape=16, size=4,color='black') +
            geom_hline(yintercept = c(0,input$R,input$R/2), 
                       lty='dashed', size=0.8,color='black') +
            geom_text(x=10,y=0.95*input$R, label='R', size=6,color='black') +
            geom_text(x=10,y=0.45*input$R, label='R/2', size=6,color='black') +
            geom_text(x=input$alpha,y=0.05*input$R, label="\u03b1", size=6,color='black') +
            scale_y_continuous(expand=c(0.02,0.02)) +
            scale_x_continuous(breaks=seq(10,55,5),minor_breaks = NULL) +
            guides(color=FALSE)
            
        
        
        print(G)
        
    })
}

# Run the application ----
shinyApp(ui = ui, server = server)
