library(shiny)
library(tidyverse)
library(bslib)

# Define UI for application that draws a histogram ----
ui <- fluidPage(
  
    theme = bslib::bs_theme(
              bootswatch = "simplex",
              base_font = font_google('Fira Sans')
              ),
  
    # Application title
    titlePanel("Quadratic Spline Fertility Model: What do the parameters do?"),
    titlePanel(HTML("<h6>This app illustrates the Quadratic Spline (QS) model proposed in <a href='https://www.demographic-research.org/volumes/vol9/5/' target=_blank'>Schmertmann (2003)</a>. 
                    Colored sections illustrate the piecewise quadratic curves that comprise the spline schedule.</br>
                    Parameters (\u03b1, P, H, R) control the shape and level. Change them to see their effects on the schedule.</h6>")),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(width=3,
            h4('SHAPE'),     
            sliderInput(inputId='alpha', label=HTML('&alpha; (Initial Age)'), 
                        min=12, max=20, value=15, 
                        step = 0.25 ),
            sliderInput(inputId='P', label='P (Peak Age)', 
                        min=20, max=48, value=26, 
                        step = 0.25 ),
            sliderInput(inputId='H', label='H (Half-peak Age)', 
                        min=20, max=48, value=34, 
                        step = 0.25 ),

            h4('LEVEL'),
            numericInput(inputId = "R",  
                         HTML("R per 1000 (Peak Level)"), value = 200,
                         width='75%',
                         min=50, max=350, step=10)
        ),

        # Show a plot of the generated distribution
        mainPanel(
           fluidRow(
             checkboxInput("show_coefs", 
                           label='Show Spline Knots and Coefs'),
             checkboxInput("show_quadratics", 
                           label='Plot All Quadratic Functions')
           ),
           
           tableOutput("coef_table"),
           
           plotOutput("fx_plot"),
        )
    )
)


QS_calc = function(x, alpha, P, H, R) {
  
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

  # calculate thetas
  theta = solve(A,target)
  
  # fitted spline
  tmp = (pmax( outer(x,knot,"-") , 0))^2
  
  y = (x >= alpha) * (x <= beta) * R * (tmp %*% theta)  # ASFRs
  
  fitted_spline = tibble(x  = x,
                         fx = as.vector(y))

  tmp   = (pmax( outer(knot,knot,"-") , 0))^2
  fknot = (knot >= alpha) * (knot <= beta) * R * (tmp %*% theta)  # ASFRs

  knot_info = tibble(x  = c(knot,beta),                    
                     fx = c(fknot,0),
                     quadratic = c(paste0('f',1:5), NA))

  # calculate the a,b,c for each section in  [a x^2 + b x + c] form
  aa =    R * cumsum(theta)
  bb = -2*R*cumsum(theta*knot)
  cc =    R*cumsum(theta*knot^2)
  
  # make an |x| by 5 data frame with the 5 complete
  # quadratic functions (each function over the whole
  # range of x, not just the interval where it *is*
  # the spline value)
  
  fvals = sapply(x, function(z) aa*z^2 + bb*z + cc) %>% 
    t() %>% 
    as_tibble() %>% 
    add_column(x, .before=1)
  
  names(fvals) = c('x',paste0('f',1:5))
  
  dividers = c(-Inf,knot,beta,Inf)
  
  sel = function(x,q) {
    ok = (q=='f1') & (x >= knot[1]) & (x < knot[2])
    ok = ok | (q=='f2') & (x >= knot[2]) & (x < knot[3])
    ok = ok | (q=='f3') & (x >= knot[3]) & (x < knot[4])
    ok = ok | (q=='f4') & (x >= knot[4]) & (x < knot[5])
    ok = ok | (q=='f5') & (x >= knot[5]) & (x < beta)
    return(ok)
  }
  
  fvals = fvals %>% 
    pivot_longer(cols=starts_with('f'),
                 names_to='quadratic',
                 values_to='fx') %>% 
    mutate(section = cut(x, breaks=dividers),
           include = sel(x,quadratic))

  
  return(list( knot          = knot, 
               theta         = theta, 
               beta          = beta,
               deltax        = diff(x)[1],
               fvals         = fvals,
               fitted_spline = fitted_spline,
               knot_info     = knot_info))
} #QS


# Define server logic ----
server <- function(input, output) {

  deltax = .25
  x      = seq(from=10, to=55, by=deltax)
  
  QS = reactive( QS_calc(x, input$alpha, 
                            input$P, 
                            input$H, 
                            input$R) )
  
  output$coef_table = renderTable({ 
    
    if (input$show_coefs) {
      knot_text = sprintf('% .3f',c(QS()$knot, QS()$beta))
      coef_text = c(sprintf('%+.3f', QS()$theta) , '--')
      
      tab        = rbind('Knot'=knot_text, 
                         'Coef'=coef_text)  %>% 
                    as_tibble()
      names(tab) = paste(0:5)
      vnames = c('Knot (t)','Coef (&#x3B8;)')
      tab        = add_column(tab, 'Variable'=vnames, .before=1)
      
      tab    
    } # if tshow_coefs
  }, align='r',sanitize.text.function = function(x) x   
  )
  
    output$fx_plot <- renderPlot({

      shiny::validate(
        need(input$alpha < input$P, '\n\n\n\n\u03b1 must be < P'),
        need(input$P     < input$H, HTML('\n\n\n\nP must be < H'))
      )

      #computation

        TFR  = 0.001 * sum( head(QS()$fitted_spline$fx,-1)+
                            tail(QS()$fitted_spline$fx,-1))/2 * QS()$deltax
        
        macb = weighted.mean( QS()$fitted_spline$x, 
                              QS()$fitted_spline$fx) + QS()$deltax

        info = paste0('TFR = '       , sprintf(fmt='%4.2f',TFR),
                      '\nMean Age = ', sprintf(fmt='%4.2f',macb))

        selected_data = QS()$fvals
        
        if (!(input$show_quadratics)) {
           selected_data = filter(selected_data, include)  
        }
        
        print(dim(data))
        
        selected_size = ifelse(input$show_quadratics,3.5, 6)
        
      #plotting
        G = ggplot(data=selected_data) +
            aes(x=x,y=fx, color=quadratic) +
            geom_line(size=selected_size, alpha=.60) +
            geom_line(data=QS()$fitted_spline,
                      aes(x=x,y=fx),size=1.5, 
                      color='black', inherit.aes = FALSE) +
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
            scale_y_continuous(expand=c(0.02,0.02), limits=c(0,1.20*input$R)) +
            scale_x_continuous(breaks=seq(10,55,5),minor_breaks = NULL) +
            guides(color=FALSE) +
            geom_text(x=45, y=.90*input$R, label=info,
                      hjust=0,color='black', size=6)

       G = G + geom_point(data=QS()$knot_info,
                          aes(x=x,y=fx),
                          shape=16, size=8, alpha=.50)
       
       print(G)

    })


    

  
}

# Run the application ----
shinyApp(ui = ui, server = server)
