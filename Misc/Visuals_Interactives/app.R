#library(plotly)
#library(gridlayout)

pacman::p_load(tidyverse,shiny,reactable,shinydashboard,shinydashboardPlus)

#library(shinydashboard)
#library(shinyWidgets)
#library(shinyjs)
#library(shinythemes)
#library(shinydashboardPlus)
#library(bslib)


## create navbarPage for ALM app ui 



input.activation<-function(x.target, c){
  return(exp(-1*c*(x.target-x.plotting)^2))
}

output.activation<-function(x.target, weights, c){
  return(weights%*%input.activation(x.target, c))
}

mean.prediction<-function(x.target, weights, c){
  probability<-output.activation(x.target, weights, c)/sum(output.activation(x.target, weights, c))
  return(y.plotting%*%probability)
}
# function to generate exam predictions
exam.prediction<-function(x.target, weights, c){
  trainVec = sort(unique(x.learning))
  nearestTrain = trainVec[which.min(abs(trainVec-x.target))]
  aresp = mean.prediction(nearestTrain, weights, c)
  xUnder = ifelse(min(trainVec) == nearestTrain, nearestTrain, trainVec[which(trainVec == nearestTrain) - 1])
  xOver = ifelse(max(trainVec) == nearestTrain, nearestTrain, trainVec[which(trainVec == nearestTrain) + 1])
  mUnder = mean.prediction(xUnder, weights, c)
  mOver = mean.prediction(xOver, weights, c)
  exam.output = round(aresp + ((mOver - mUnder) / (xOver - xUnder)) * (x.target - nearestTrain), 3)
  exam.output
}

update.weights<-function(x.new, y.new, weights, c, lr){
  y.feedback.activation<-exp(-1*c*(y.new-y.plotting)^2)
  x.feedback.activation<-output.activation(x.new, weights, c)
  return(weights+lr*(y.feedback.activation-x.feedback.activation)%*%t(input.activation(x.new, c)))
}

update.weights.with_noise <- function(x.new, y.new,weights, c, lr, noise_sd){
  y.feedback.activation <- exp(-1 * c * (y.new - y.plotting)^2)
  x.feedback.activation <- output.activation(x.new, weights, c)
  weight_updates <- lr * (y.feedback.activation - x.feedback.activation) %*% t(input.activation(x.new, c))
  # Add random noise to the weight updates
  noiseMat <- matrix(rnorm(nrow(weight_updates) * ncol(weight_updates), mean = noise_sd), 
                  nrow = nrow(weight_updates), ncol = ncol(weight_updates))
  updated_weights <- weights + weight_updates + noiseMat
  return(updated_weights)
}



learn.alm<-function(y.learning, c=0.05, lr=0.5,noise_sd=0){
  weights<-matrix(rep(0.00, length(y.plotting)*length(x.plotting)), nrow=length(y.plotting), ncol=length(x.plotting))
  alm.train <- rep(NA, length(y.learning))
  for (i in 1:length(y.learning)){
    weights<-update.weights.with_noise(x.learning[i], y.learning[i], weights, c, lr,noise_sd)
    weights[weights<0]=0
    resp = mean.prediction(x.learning[i], weights, c)
    alm.train[i] = resp
  }
  alm.predictions<-sapply(x.plotting, mean.prediction, weights=weights, c=c)
  exam.predictions <- sapply(x.plotting, exam.prediction, weights=weights, c=c)
  return(list(alm.train=alm.train, alm.predictions=alm.predictions, exam.predictions=exam.predictions))
  #return(list(alm.predictions=alm.predictions, exam.predictions=exam.predictions,wmFinal=weights))
}



x.plotting<<-seq(0,90, .5)
y.plotting<<-seq(0, 210, by=2)
#trainOptions=round(seq(1,length(x.plotting),length.out=21),0)
trainOptions=x.plotting[seq(1,181,by=4)]
trainItems=trainOptions[c(10,11,12)]



### An alternate, shiny dashboard version of the app. Which has the main app in one tab, and a separate tab showing the functions used to generate the data and the model, and briefly describing how
## the code works to implement the models. 


# Define UI for application
# 
ui <- dashboardPage(

  skin = "black",
  dashboardHeader(title = "ALM Simulation App"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Home", tabName = "home", icon = icon("home")),
      menuItem("Code", tabName = "code", icon = icon("code"))
    )
  ),
  dashboardBody(
    tabItems(
      tabItem(tabName = "home",
              fluidRow(
                column(4,
                       box(
                         title = "Simulation Parameters",
                         status = "primary",
                         solidHeader = TRUE,
                         collapsible = TRUE,
                         collapsed = FALSE,
                         width = 12,
                         sliderInput("assoc", "Generalization Parameter (c):",
                                     min = .0001, max = 1, value = 0.5, step = 0.01),
                         sliderInput("update", "Learning Rate:",
                                     min = 0.00001, max = 1, value = 0.5, step = 0.01),
                         sliderInput("trainRep", "Training Repetitions Per Item:",
                                     min = 1, max = 200, value = 3, step = 1),
                         sliderInput("Noise","Output Noise:",
                                     min = 0, max = 100, value = 0.00, step = 1),
                         sliderInput("Learn_Noise","Learning Noise:",
                                     min = 0, max = 100, value = 0.00, step = 1),
                         checkboxGroupInput("trainItems", "Training Items:", choices = trainOptions, selected = trainOptions[c(10,15,35)],inline=TRUE),
                         # radio buttons for selecting function form
                         radioButtons("functionForm", "Function Form:",
                                      choices = c("Linear", "Quadratic", "Exponential"),
                                      selected = "Quadratic"),
                        # numericInput("nRep", "Number of Replications:", value = 1, min = 1, max = 100),
                         actionButton("run", "Run Simulation")
                       )
                ),
                column(8,
                       box(
                         title = "Model Performance",
                         status = "primary",
                         solidHeader = TRUE,
                         collapsible = TRUE,
                         collapsed = FALSE,
                          width = 12,
                         plotOutput("trainPlot"),
                         plotOutput("plot"),
                         h5("*Dashed line shows true function. Red shows ALM, and blue depicts EXAM predictions*"),
                         h4("Average Model Performance"),
                         reactableOutput("table"),
                         h4("Model Performance by Item Type"),
                         reactableOutput("table2")
                       )
                )
              )
      ),
      tabItem(tabName = "code",
              fluidRow(
                column(12,
                       box(
                         title = "Code",
                         status = "primary",
                         solidHeader = TRUE,
                         collapsible = TRUE,
                         collapsed = FALSE,
                         width = 12,
                         verbatimTextOutput("code")
                       )
                )
                )
        )
    )
    )
)

# Define server 



server <- function(input, output, session) {
  
  nRep=1
  user_choice <- eventReactive(input$run, {
    return(list(assoc = input$assoc, update = input$update, Noise=input$Noise,
                learn_noise=input$Learn_Noise,
                functionForm=input$functionForm,trainRep = as.numeric(input$trainRep),
                trainItems = input$trainItems))
    
  }, ignoreNULL = FALSE)
  

    output_df <- eventReactive(input$run, {
      uc <- reactive({user_choice()})
    if (uc()$functionForm == "Linear") {
      f.plotting <<- as.numeric(x.plotting * 2.2 + 30)
    } else if (uc()$functionForm == "Quadratic") {
      f.plotting <<- as.numeric(210 - ((x.plotting - 50)^2) / 12)
    } else if (uc()$functionForm == "Exponential") {
      # f.plotting<<-as.numeric(scale(200*(1-exp(-x.plotting/25))))
      f.plotting <<- as.numeric(200 * (1 - exp(-x.plotting / 25)))
    }
    trainItems <- as.numeric(uc()$trainItems)
    y.plotting <<- seq(0, max(f.plotting), by = 1)
    x.learning <<- rep(trainItems, times = uc()$trainRep)
    f.learning <<- rep(f.plotting[which(x.plotting %in% trainItems)], times = uc()$trainRep)
    # print(x.learning)
    # print(f.learning)
    # print(uc()$trainRep)
    # print(trainItems)
    # print(uc()$functionForm)
    
    # output_list <- replicate(nRep, list(learn.alm(f.learning + rnorm(length(f.learning), sd = Noise),
    #                                                   c = assoc, lr = update)))

    output_list <- replicate(nRep, list(learn.alm(f.learning + rnorm(length(f.learning), sd = uc()$Noise),
                                                  c = uc()$assoc, lr = uc()$update,noise_sd=uc()$learn_nois)))
    
    almTrain <- tibble(trial=seq(1,length(f.learning)),Input=as.factor(x.learning),y=f.learning,pred=map(output_list, "alm.train") %>% unlist(),
    error=abs(y-pred))  

   output_list <- lapply(output_list, "[", 2:3)
    output_df <- lapply(output_list, function(x) as.data.frame(x))
    #output_df <- lapply(output_list, function(x) lapply(x, as.data.frame)) # 10 dfs x 9 lists
    output_df <- Reduce(rbind, output_df) %>% mutate(x = x.plotting, y = f.plotting)
    #output_df <- lappl y(output_df, function(x) Reduce(rbind,x))# 1 df x 9 lists
    output_df <- output_df %>%
      pivot_longer(names_to = "Model", values_to = "Prediction", cols = c(alm.predictions, exam.predictions)) %>%
      rbind(data.frame(data.frame(x = x.plotting, y = f.plotting, Model = "True Function", Prediction = f.plotting)), .)
    #str(output_df)

    output_df <- list(output_df=output_df,almTrain=almTrain)
    return(output_df)
    
    }, ignoreNULL = FALSE)
    

  output$trainPlot <- renderPlot({
      almTrain <- reactive({output_df() %>% pluck("almTrain")})

       ggplot(almTrain(),aes(x = trial, y = error, color = Input)) + geom_line()+geom_point()
  })


    output$plot <- renderPlot({
       
      output_df2 <- reactive({output_df() %>% pluck("output_df")})
      ggplot(data = output_df2(), aes(x = x, y = Prediction,color=Model),alpha=.2) + 
        geom_line(aes(linetype=Model,alpha=Model)) + 
        geom_point(data = data.frame(x.learning, f.learning), 
                   aes(x = x.learning,y = f.learning),color="black",size=4,shape=4) +
        # geom_line(data = data.frame(x.plotting, f.plotting), 
        #           aes(x = x.plotting, y = f.plotting),linetype=2, color = "black",alpha=.3) + 
        scale_color_manual(values = c("red", "blue", "black"))+
        scale_alpha_manual(values=c(.8,.8,.4))+
        scale_linetype_manual(values=c(1,1,2))+
        ylim(c(0,250))#+
        # ggtitle(paste("Association Parameter:", user_choice()$assoc, " Update Parameter:", 
        #               uc$update, "Train Reps:", 
        #               uc$trainRep, " Noise:", uc$Noise))
    })  
    # table 1 reports the summary stats for all items. Table uses GT library to make gt table
    output$table <- renderReactable({
      output_df <- output_df() %>% pluck("output_df") 
      output_df %>% group_by(Model) %>% filter(Model !="True Function") %>%
        summarise(MeanDeviation = mean(abs(Prediction - y)), 
                  RMSD = sqrt(mean((Prediction -y)^2)),Correlation = cor(Prediction, y)) %>%
        mutate(across(where(is.numeric), round, 1)) %>%
        reactable::reactable(compact=TRUE,bordered = TRUE, highlight = TRUE, resizable=TRUE)
    })
    # table 2 reports the summary stats separately for training items, interpolation items, and extrapolation items # nolint
    output$table2 <- renderReactable({
      uc <- reactive({user_choice()})
      output_df() %>% pluck("output_df") %>% filter(Model !="True Function") %>% 
        mutate(ItemType = ifelse(x %in% x.learning, "Training", ifelse(x > min(x.learning) & x < max(x.learning), "Interpolation", "Extrapolation"))) %>%
        group_by(ItemType,Model) %>%
        summarise(MeanDeviation = mean(abs(Prediction - y)), 
                  RMSD = sqrt(mean((Prediction -y)^2)),Correlation = cor(Prediction, y), 
                  .groups="keep") %>% 
        mutate(across(where(is.numeric), round, 1)) %>%
        reactable::reactable(compact=TRUE,bordered = TRUE, highlight = TRUE, resizable=TRUE) 
    })
    
    
    output$code <- renderPrint({
      # code to implement the ALM and EXAM models
      # code to generate data
      # code to run models
      # code to format output
      cat(" input.activation<-function(x.target, c){
  return(exp(-1*c*(x.target-x.plotting)^2))
}

output.activation<-function(x.target, weights, c){
  return(weights%*%input.activation(x.target, c))
}

mean.prediction<-function(x.target, weights, c){
  probability<-output.activation(x.target, weights, c)/sum(output.activation(x.target, weights, c))
  return(y.plotting%*%probability)
}
# function to generate exam predictions
exam.prediction<-function(x.target, weights, c){
  trainVec = sort(unique(x.learning))
  nearestTrain = trainVec[which.min(abs(trainVec-x.target))]
  aresp = mean.prediction(nearestTrain, weights, c)
  xUnder = ifelse(min(trainVec) == nearestTrain, nearestTrain, trainVec[which(trainVec == nearestTrain) - 1])
  xOver = ifelse(max(trainVec) == nearestTrain, nearestTrain, trainVec[which(trainVec == nearestTrain) + 1])
  mUnder = mean.prediction(xUnder, weights, c)
  mOver = mean.prediction(xOver, weights, c)
  exam.output = round(aresp + ((mOver - mUnder) / (xOver - xUnder)) * (x.target - nearestTrain), 3)
  exam.output
}

update.weights<-function(x.new, y.new, weights, c, lr){
  y.feedback.activation<-exp(-1*c*(y.new-y.plotting)^2)
  x.feedback.activation<-output.activation(x.new, weights, c)
  return(weights+lr*(y.feedback.activation-x.feedback.activation)%*%t(input.activation(x.new, c)))
}

learn.alm<-function(y.learning, c=0.05, lr=0.5){
  weights<-matrix(rep(0.00, length(y.plotting)*length(x.plotting)), nrow=length(y.plotting), ncol=length(x.plotting))
  for (i in 1:length(y.learning)){
    weights<-update.weights(x.learning[i], y.learning[i], weights, c, lr)
    weights[weights<0]=0
  }
  alm.predictions<-sapply(x.plotting, mean.prediction, weights=weights, c=c)
  exam.predictions <- sapply(x.plotting, exam.prediction, weights=weights, c=c)
  return(list(alm.predictions=alm.predictions, exam.predictions=exam.predictions))
  #return(list(alm.predictions=alm.predictions, exam.predictions=exam.predictions,wmFinal=weights))
}


    ")
    })
    
}



# Run the application


shinyApp(ui, server)






# 
# ui <- navbarPage(
#   title = "ALM EXAM",
#   theme = bs_theme(bootswatch = "morph"),
#   collapsible = FALSE,
#   tabPanel(
#     title = "ALM & EXAM",
#     sidebarLayout(
#       sidebarPanel(
#         sliderInput("assoc",
#           "Association Parameter",
#           min = 0,
#           max = 1,
#           value = 0.5,
#           step = 0.1
#         ),
#         sliderInput("update",
#           "Update Parameter",
#           min = 0,
#           max = 1,
#           value = 0.5,
#           step = 0.1
#         ),
#         sliderInput("trainRep",
#           "Training Repetitions",
#           min = 1,
#           max = 100,
#           value = 10,
#           step = 1
#         ),
#         sliderInput("Noise",
#           "Noise",
#           min = 0,
#           max = 1,
#           value = 0.1,
#           step = 0.1
#         ),
#         checkboxGroupInput("trainItems",
#           "Training Items",
#           choices = trainOptions,
#           selected = trainOptions[c(
#             10,
#             15,
#             35
#           )],
#           inline = TRUE
#         ),
#         radioButtons("functionForm",
#           "Function Form",
#           choices = c(
#             "Linear",
#             "Quadratic",
#             "Exponential"
#           ),
#           selected = "Linear"
#         ),
#         numericInput("nRep",
#           "Number of Repetitions",
#           value = 1,
#           min = 1,
#           max = 100
#         ),
#         actionButton(
#           "run",
#           "Run Simulation"
#         )
#       ),
#       mainPanel(
#         plotOutput("plot"),
#         h3("Average Model Performance"),
#         reactableOutput("table"),
#         h3("Model Performance by Item Type"),
#         reactableOutput("table2")
#       )
#     )
#   ),
#   tabPanel(
#     title = "Model Definition",
#     grid_container(
#       layout = c(
#         "area0 .",
#         ".     ."
#       ),
#       row_sizes = c(
#         "1.73fr",
#         "0.27fr"
#       ),
#       col_sizes = c(
#         "1.73fr",
#         "0.27fr"
#       ),
#       gap_size = "10px",
#       grid_card_text(
#         content = htmltools::includeMarkdown("modelCode.md"),
#         alignment = "start",
#         area = "area0"
#       )
#     )
#   )
# 
# )




# server <- function(input, output, session) {
#   
# observeEvent(input$run, {
#      assoc = input$assoc
#      update = input$update
#      #nRep = input$nRep
#      nRep=1
#      Noise = input$Noise
#      functionForm = input$functionForm
#      trainRep = input$trainRep
#      trainItems <- input$trainItems
#      trainItems <- as.numeric(trainItems)
# 
#      if (functionForm == "Linear") {
#          f.plotting <<- as.numeric(x.plotting * 2.2 + 30)
#      } else if (functionForm == "Quadratic") {
#          f.plotting <<- as.numeric(210 - ((x.plotting - 50)^2) / 12)
#      } else if (functionForm == "Exponential") {
#          # f.plotting<<-as.numeric(scale(200*(1-exp(-x.plotting/25))))
#          f.plotting <<- as.numeric(200 * (1 - exp(-x.plotting / 25)))
#      }
# 
#      y.plotting <<- seq(0, max(f.plotting), by = 1)
#      x.learning <<- rep(x.plotting[trainItems], times = trainRep)
#      f.learning <<- rep(f.plotting[trainItems], times = trainRep)
# 
#      output_list <- replicate(nRep, list(learn.alm(f.learning + rnorm(length(f.learning), sd = Noise),
#          c = assoc, lr = update)))
# 
# 
#      
#      # weight.mat = output_list[[1]]$wmFinal
#      # output_list = list(output_list[[1]][1], output_list[[1]][2])
# 
#      output_df <- lapply(output_list, function(x) as.data.frame(x))
#      #output_df <- lapply(output_list, function(x) lapply(x, as.data.frame)) # 10 dfs x 9 lists
#      
#      output_df <- Reduce(rbind, output_df) %>% mutate(x = x.plotting, y = f.plotting)
#      #output_df <- lapply(output_df, function(x) Reduce(rbind,x))# 1 df x 9 lists
#      
#      output_df <- output_df %>%
#          pivot_longer(names_to = "Model", values_to = "Prediction", cols = c(alm.predictions, exam.predictions)) %>%
#          rbind(data.frame(data.frame(x = x.plotting, y = f.plotting, Model = "True Function", Prediction = f.plotting)), .)
#     
#    output$plot <- renderPlot({
#       ggplot(data = output_df, aes(x = x, y = Prediction,color=Model),alpha=.2) + 
#         geom_line(aes(linetype=Model,alpha=Model)) + 
#         geom_point(data = data.frame(x.learning, f.learning), 
#                    aes(x = x.learning,y = f.learning),color="black",size=4,shape=4) +
#         # geom_line(data = data.frame(x.plotting, f.plotting), 
#         #           aes(x = x.plotting, y = f.plotting),linetype=2, color = "black",alpha=.3) + 
#         scale_color_manual(values = c("red", "blue", "black"))+
#         scale_alpha_manual(values=c(.8,.8,.4))+
#         scale_linetype_manual(values=c(1,1,2))+
#         ggtitle(paste("Association Parameter:", assoc, " Update Parameter:", update, " Train Reps:", trainRep, " Noise:", Noise))
#     }) 
#     # table 1 reports the summary stats for all items. Table uses GT library to make gt table
#     output$table <- renderReactable({
#       output_df %>% group_by(Model) %>% filter(Model !="True Function") %>%
#         summarise(MeanDeviation = mean(abs(Prediction - y)), 
#                   RMSD = sqrt(mean((Prediction -y)^2)),Correlation = cor(Prediction, y)) %>%
#         mutate(across(where(is.numeric), round, 1)) %>%
#         reactable::reactable(compact=TRUE,bordered = TRUE, highlight = TRUE, resizable=TRUE)
#     })
#     # table 2 reports the summary stats separately for training items, interpolation items, and extrapolation items
#     output$table2 <- renderReactable({
#       output_df %>% filter(Model !="True Function") %>% 
#         mutate(ItemType = ifelse(x %in% x.learning, "Training", ifelse(x > min(x.learning) & x < max(x.learning), "Interpolation", "Extrapolation"))) %>%
#         group_by(ItemType,Model) %>%
#         summarise(MeanDeviation = mean(abs(Prediction - y)), 
#                   RMSD = sqrt(mean((Prediction -y)^2)),Correlation = cor(Prediction, y)) %>% 
#         mutate(across(where(is.numeric), round, 1)) %>%
#         reactable::reactable(compact=TRUE,bordered = TRUE, highlight = TRUE, resizable=TRUE) 
#     })
#     
#     
#     output$code <- renderPrint({
#         # code to implement the ALM and EXAM models
#         # code to generate data
#         # code to run models
#         # code to format output
#         cat(" input.activation <- function(x,assoc,update,trainRep){
#   # activation function for ALM
#   # x is the input vector
#   # assoc is the association parameter
#   # update is the update parameter
#   # trainRep is the number of training repetitions
#   # returns the activation vector
#   # initialize activation vector
#   activation <- rep(0,length(x))
#   # loop through training repetitions
#   for (i in 1:trainRep){
#     # loop through items
#     for (j in 1:length(x)){
#       # update activation vector
#       activation[j] <- activation[j] + assoc*x[j] + update*(1-activation[j])
#     }
#   }
#   return(activation)
#   }
#   output.activation <- function(x,assoc,update,trainRep){
#     # activation function for EXAM
#     # x is the input vector
#     # assoc is the association parameter
#     ")
#     })
#     
#     
#     
#     },ignoreNULL = FALSE)
# 
# 
# 
# }




#runApp(list(ui=ui,server=server),launch.browser=TRUE)


#library(markdown)


#runApp(list(server=server),launch.browser=TRUE)


