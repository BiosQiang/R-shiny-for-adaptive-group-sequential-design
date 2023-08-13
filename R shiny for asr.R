
library(gsDesign)

# install.packages('ggdark')
library(ggdark)

library(shiny)

# install.packages("shinythemes")


library(shinythemes)


ui <- fluidPage(themeSelector(),
  
  # App title ----
  titlePanel("Adaptive group sequential with sample size re-estimation"),
  
  # Sidebar layout with input and output definitions ----
  fluidRow(
   
     column(2,
      
      sliderInput(inputId = "CP",
                  label = "adjusting zone for CP:",
                  min = 0,
                  max = 1,
                  value =c(0.1,0.9)),
      br(),

      numericInput(inputId = 'alpha',
                   label='familiwise type I error',
                   min = 0,
                   max = 1,
                   value =0.025)
     ),

    
    column(2,
      numericInput(inputId = 'beta',
                   label='type II error',
                   min = 0,
                   max = 1,
                   value =0.1),
      br(),
      
      numericInput(inputId = 'delta0',
                   label='absulate difference under H0',
                   value =0)
    ),
  
    
    column(2,

      numericInput(inputId = 'delta1',
                   label='absulate difference under H1',
                   value =1),
      br(),
      
      numericInput(inputId = 'SD',
                   label='pooled standard devarience',
                   value =1)
    ),
  
    column(2,
      numericInput(inputId = 'maxinf',
                   label='maximum inflation ratio',
                   value=2,
                   min=1,
                   max=Inf),
      br(),
      
      numericInput(inputId = 'z1',
                   label='statistics for stage 1',
                   value=1,
                   min=0,
                   max=4,
                   step=0.1)
  ),

column(2,
      numericInput(inputId = 'timing',
                   label='timing for SSR',
                   value=0.5,
                   min=0.000001,
                   max=0.9999999),
      br(),
      
      numericInput(inputId = 'upar',
                   label='parameter for upper boundary',
                   value =-4,
                   min=-40,
                   max=40)
      ),

column(2,
      numericInput(inputId = 'lpar',
                   label='parameter for lower boundary',
                   value =-4,
                   min=-40,
                   max=40),
      br(),
      
      selectInput(inputId = 'end',
                    label='type of endpoint',
                    c('normal'='normal','binary'='binary','time-to-event'='time-to-event'))
     )),
    
    # Main panel for displaying outputs ----
    

     column(6, 
      # Output: Histogram ----
      textOutput(outputId = "zone"),
      h4('observed effect size vs.stage II sample size (using observed effect size at IA)'),
      plotOutput(outputId = "adaptive_sample"),
      h4('critical value for stage II under different combination methods'),
      plotOutput(outputId = "combin_z"),
      h4('observed effect size vs.stage II sample size under different combination methods'),
      plotOutput(outputId = "combin_n2")
      ),
     column(6,
      h4('total sample size under different assmpution on future effect size'),
      plotOutput(outputId = "hypoth"),
      h4('expected sample size under different assmpution on future effect size'),
      plotOutput(outputId = "expect"),
      h4('power changing under different assmpution on future effect size'),
      plotOutput(outputId = "power"))
      
    )
  



server <- function(input, output) {
  
  library(gsDesign)
  library(ggplot2)
  
n.fix=reactive({
  if (input$end=='normal') {n.fix=nNormal(delta1=input$delta1,alpha = input$alpha,beta=input$beta,sd=input$SD)}
  else if (input$end=='binary') {n.fix=nBinomial(p1=input$delta1,p2=0.0001,alpha = input$alpha,beta=input$beta)}
  else if (input$end=='time-to-event') {n.fix=nEvents(hr=input$delta1,alpha = input$alpha,beta=input$beta)}
})

design=reactive({
  gsDesign(k=2,test.type=3,alpha =input$alpha,beta=input$beta,n.fix=n.fix(),timing = c(input$timing,1),delta0 = input$delta0,delta1 = input$delta1,sfu=sfHSD,sfupar = input$upar,sfl=sfHSD,sflpar = input$lpar)
})

ssrNC=reactive(
  ssrCP(
    x = design(), z1 = seq(0,4,0.05), overrun = 0, beta = input$beta, cpadj = input$CP,
    maxinc = input$maxinf, z2 = z2NC
  )
  
)

ssrFisher=reactive(
  ssrCP(
    x = design(), z1 = seq(0,4,0.05), overrun = 0, beta = input$beta, cpadj = input$CP,
    maxinc = input$maxinf, z2 = z2Fisher)
)
  
ssrSuf=reactive(
    ssrCP(
      x = design(), z1 = seq(0,4,0.05), overrun = 0, beta = input$beta, cpadj = input$CP,
      maxinc = input$maxinf, z2 = z2Z)
)

combin=reactive(
  rbind(
    data.frame(cbind(ssrNC()$dat, Test = "Normal combination")),
    data.frame(cbind(ssrSuf()$dat, Test = "Sufficient statistic")),
    data.frame(cbind(ssrFisher()$dat, Test = "Fisher combination"))
  )
)

altdelta=reactive({
  
 ssrCP(
    x = design(), z1 = seq(0,4,0.05), overrun = 0, beta = input$beta, cpadj = input$CP,
    maxinc = input$maxinf, z2 = z2NC, theta = design()$delta)
  # combine data frames for the 2 designs
  
})

combin_del=reactive({
  rbind(
    data.frame(cbind(ssrNC()$dat, "CP effect size" = "Obs. at IA")),
    data.frame(cbind(altdelta()$dat, "CP effect size" = "Alt. hypothesis"))
  )
})

power=reactive({
  y1 <- Power.ssrCP(x = ssrNC())
  y2 <- Power.ssrCP(x = altdelta())
  # combine data frames for the 2 designs
  y3 <- rbind(
    data.frame(cbind(y1, "CP effect size" = "Obs. at IA")),
    data.frame(cbind(y2, "CP effect size" = "Alt. hypothesis"))
  )
})

  output$adaptive_sample <- renderPlot(
    plot(ssrNC())
  )
  
  output$combin_z<-renderPlot(
    ggplot(data = combin(), aes(x = z1, y = z2, col = Test)) +
    geom_line()+theme_classic()


  )
  output$combin_n2<-renderPlot(
    ggplot(data = combin(), aes(x = z1, y = n2, col = Test)) +
      geom_line()+theme_classic()

  )
  output$hypoth<-renderPlot(
    ggplot(data = combin_del(), aes(x = z1, y = n2, col = CP.effect.size)) +
      geom_line()+theme_classic()

  )
  
  output$expect<-renderPlot(
    ggplot(data = power(), aes(x = delta, y = en, col = CP.effect.size)) +
      geom_line()+
      xlab(expression(delta)) +
      ylab("Expected sample size")+theme_classic()
  )
  
  output$power<-renderPlot(
    ggplot(data = power(), aes(x = delta, y = Power, col = CP.effect.size)) +
      geom_line()+theme_classic()
    
  )
  
}


shinyApp(ui = ui, server = server)

