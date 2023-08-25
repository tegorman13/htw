
pacman::p_load(tidyverse,shiny,reactable)


input.activation<-function(x.target, association.parameter){
  #return(exp(-1*association.parameter*(100*x.target-100*x.plotting)^2))
  return(exp(-1*association.parameter*(x.target-x.plotting)^2))
}

output.activation<-function(x.target, weights, association.parameter){
  return(weights%*%input.activation(x.target, association.parameter))
}

mean.prediction<-function(x.target, weights, association.parameter){
  probability<-output.activation(x.target, weights, association.parameter)/sum(output.activation(x.target, weights, association.parameter))
  return(y.plotting%*%probability)
}
# function to generate exam predictions
exam.prediction<-function(x.target, weights, association.parameter){
  trainVec = sort(unique(x.learning))
  nearestTrain = trainVec[which.min(abs(trainVec-x.target))]
  aresp = mean.prediction(nearestTrain, weights, association.parameter)
  xUnder = ifelse(min(trainVec) == nearestTrain, nearestTrain, trainVec[which(trainVec == nearestTrain) - 1])
  xOver = ifelse(max(trainVec) == nearestTrain, nearestTrain, trainVec[which(trainVec == nearestTrain) + 1])
  mUnder = mean.prediction(xUnder, weights, association.parameter)
  mOver = mean.prediction(xOver, weights, association.parameter)
  exam.output = round(aresp + ((mOver - mUnder) / (xOver - xUnder)) * (x.target - nearestTrain), 3)
  exam.output
}

update.weights<-function(x.new, y.new, weights, association.parameter, update.parameter){
  y.feedback.activation<-exp(-1*association.parameter*(y.new-y.plotting)^2)
  x.feedback.activation<-output.activation(x.new, weights, association.parameter)
  return(weights+update.parameter*(y.feedback.activation-x.feedback.activation)%*%t(input.activation(x.new, association.parameter)))
}

learn.alm<-function(y.learning, association.parameter=0.05, update.parameter=0.5){
  weights<-matrix(rep(0.00, length(y.plotting)*length(x.plotting)), nrow=length(y.plotting), ncol=length(x.plotting))
  for (i in 1:length(y.learning)){
    weights<-update.weights(x.learning[i], y.learning[i], weights, association.parameter, update.parameter)
    weights[weights<0]=0
  }
  alm.predictions<-sapply(x.plotting, mean.prediction, weights=weights, association.parameter=association.parameter)
  exam.predictions <- sapply(x.plotting, exam.prediction, weights=weights, association.parameter=association.parameter)
  return(list(alm.predictions=alm.predictions, exam.predictions=exam.predictions))
}


x.plotting<<-seq(0,90, .5)
y.plotting<<-seq(0, 210, by=2)
trainOptions=round(seq(1,length(x.plotting),length.out=50),0)
trainItems=trainOptions[c(10,11,12)]


ui <- fluidPage(
  
  titlePanel("ALM & EXAM Predictions"),
  
  sidebarLayout(
    sidebarPanel(
      sliderInput("assoc", "Association Parameter (c):",
                  min = 0, max = 1, value = 0.03, step = 0.01),
      sliderInput("update", "Update Parameter (lr):",
                  min = 0, max = 1, value = 0.2, step = 0.05),
      sliderInput("trainRep", "Training Repetitions:",
                  min = 1, max = 10, value = 1, step = 1),
      sliderInput("Noise","Noise Level:",
                  min = 0, max = 2, value = 0.00, step = 0.1),
      checkboxGroupInput("trainItems", "Training Items:", choices = trainOptions, selected = trainOptions[c(10,15,35)],inline=TRUE),
      # radio buttons for selecting function form
      radioButtons("functionForm", "Function Form:", 
                   choices = c("Linear", "Quadratic", "Exponential"), 
                   selected = "Linear"),
      numericInput("nRep", "Number of Replications:", value = 1, min = 1, max = 100),
      actionButton("run", "Run Simulation")
    ),
    mainPanel(
      plotOutput("plot"),
      h4("Average Model Performance"),
      reactableOutput("table"),
      h4("Model Performance by Item Type"),
      reactableOutput("table2")
    )
  )
)

server <- function(input, output) {
  
  observeEvent(input$run, {
    assoc = input$assoc
    update = input$update
    nRep = input$nRep
    Noise=input$Noise
    functionForm=input$functionForm
    trainRep=input$trainRep
    trainItems <- input$trainItems
    trainItems <- as.numeric(trainItems)

    if (functionForm=="Linear"){
      f.plotting<<-as.numeric(x.plotting*2.2+30)
    } else if (functionForm=="Quadratic"){
      f.plotting<<-as.numeric(210 - ((x.plotting-50)^2)/12)
    } else if (functionForm=="Exponential"){
      #f.plotting<<-as.numeric(scale(200*(1-exp(-x.plotting/25))))
      f.plotting<<-as.numeric(200*(1-exp(-x.plotting/25)))
    }
    
    y.plotting<<-seq(0, max(f.plotting), by=1)
    x.learning<<-rep(x.plotting[trainItems],times=trainRep)
    f.learning<<-rep(f.plotting[trainItems],times=trainRep)
    
    output_list <- replicate(nRep, list(learn.alm(f.learning + rnorm(length(f.learning), sd=Noise), 
                                                  association.parameter=assoc, update.parameter=update)))
    
    output_df <- lapply(output_list, function(x) as.data.frame(x))
    output_df <- Reduce(rbind, output_df) %>% mutate(x=x.plotting,y=f.plotting)
    output_df <- output_df %>% pivot_longer(names_to = "Model", values_to = "Prediction", cols = c(alm.predictions, exam.predictions))
    
    output$plot <- renderPlot({
      ggplot(data = output_df, aes(x = x, y = Prediction,color=Model),alpha=.2) + 
        geom_line() + 
        geom_point(data = data.frame(x.learning, f.learning), 
                   aes(x = x.learning,y = f.learning),color="black",size=4,shape=4) +
        geom_line(data = data.frame(x.plotting, f.plotting), 
                  aes(x = x.plotting, y = f.plotting), color = "black",alpha=.3) + 
        ggtitle(paste("Association Parameter:", assoc, " Update Parameter:", update, " Train Reps:", trainRep, " Noise:", Noise))
    }) 
    # table 1 reports the summary stats for all items. Table uses GT library to make gt table
    output$table <- renderReactable({
      output_df %>% group_by(Model) %>%
        summarise(MeanDeviation = mean(abs(Prediction - y)), 
                  RMSD = sqrt(mean((Prediction -y)^2)),Correlation = cor(Prediction, y)) %>%
        reactable::reactable(compact=TRUE)
    })
    # table 2 reports the summary stats separately for training items, interpolation items, and extrapolation items
    output$table2 <- renderReactable({
      output_df %>% # add a column to indicate whether the item is a training item, interpolation item, or extrapolation 
        mutate(ItemType = ifelse(x %in% x.learning, "Training", ifelse(x > min(x.learning) & x < max(x.learning), "Interpolation", "Extrapolation"))) %>%
        group_by(ItemType,Model) %>%
        summarise(MeanDeviation = mean(abs(Prediction - y)), 
                  RMSD = sqrt(mean((Prediction -y)^2)),Correlation = cor(Prediction, y)) %>% 
        reactable::reactable() 
    })},ignoreNULL = FALSE)
}

shinyApp(ui, server)