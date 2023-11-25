


packages <-  c("tidyverse","foreach","doParallel","doSNOW","data.table")

have = packages %in% rownames(installed.packages())
if ( any(!have) ) { install.packages(packages[!have]) }
invisible(lapply(packages, require, character.only = TRUE))
source("Exam_Functions3.R")
d <- readRDS("data/htw-04-07-22.rds")
d$trainStage=factor(d$trainStage,ordered = TRUE)
d <- d %>% filter(catOrder=="orig",feedbackType=="continuous",goodThrow=="TRUE") %>% droplevels()
dt <- d %>% filter(expMode=="train") %>% select(sbjCode,condit,expMode,trial,trialType,band,throwCategory,vx,dist,bandInt,highBound,feedback)
rm(d)


runInfoBase=list(functions=list(inputActivation=inputActivation,outputActivation=outputActivation,
                                almResponse=almResponse,computeError=computeError,updateWeights=updateWeights,
                                almLikelihood=almLikelihood,trainAlm=trainAlm),
                 runScripts=list(baseScript=readLines("fit_alm_ll.R"),
                                 exam_functions=readLines("Exam_Functions3.R"),
                                 environmentVars=ls(),
                                 path=getwd()),
                 systemInfo=Sys.info(),
                 sessionInfo=sessionInfo(),timeInit=Sys.time())


note1 = "HTW E1 - optimizing training fits with LL - ALM Only, no test trials of any kind - EXP choice"
runInfo1 <-list(almTrainFitsLL=almTrainFitsLL,
                inputs=list(indvCodes=indvCodes,parmStart=parStart.c.lr,e1train=dfit),
                runInfo=c(runInfoBase,saveTime=format(Sys.time(), "%a %b %d %Y %X "),fitDuration=t1),
                notes=note1,
                minutes.to.fit=(t1[3]/60))

saveRDS(runInfo1,paste0("model_fits/opt.indv.E1_train_ll_EXP_Choice",format(Sys.time(), "%H_%M_%OS"),".rds"))

in.opts <- list(objFun="LL",parmLowLimit=list(cLow=.0000000001,lrLow=.00000001,returnExtreme=100000))

coreCheck=detectCores(all.tests = FALSE, logical = TRUE)
ncores   <- coreCheck
n_reps <- 1
n_d= length(indvCodes)
cl       <- makeSOCKcluster(ncores)
pb       <- txtProgressBar(max = n_d*n_reps, style = 3)
progress <- function(n){setTxtProgressBar(pb, n)}
sopts     <- list(progress = progress)
registerDoSNOW(cl)

t1 <- system.time({ 
  almTrainFitsLL <- foreach(i=1:length(indvCodes),.packages = packages, .combine= "rbind",.options.snow = sopts,.verbose=T) %dopar% {
    source("Exam_Functions3.R")
    sbj=indvCodes[i]
    
    dfit <- dt %>% filter() %>% 
      mutate(bandF=band,band=as.numeric(band),vxFit=ifelse(vx>max(outputNodes),max(outputNodes)+10,vx)) %>% setDT
    
    sbjResults <- fitAlmIndv_ll(sbj,parStart.c.lr,fd=dfit,in.opts)
    sbj2 <- sbjResults %>% as_tibble()
    sbj2
  }
})
stopCluster(cl)


mFiles = list.files("model_fits",".rds",full.names=TRUE)
allM <- lapply(mFiles,readRDS)
#res <- lapply(allM, summary)
names(allM) <- substr(mFiles , 12, 40)
#format(object.size(allM),units="Mb")

mNotes <- c(sapply(allM,"[[","note"),sapply(allM,"[[","notes"))
mNotes
mInclude=c("fitSbjTrain","noW_trainTrain")
mInclude=c("fitSbjTrain","noW_trainTrain","noW_trainTrainNf")

mpl= list()
for (i in 1:length(allM))
{
  fname=names(allM)[i]
  name=stringr::str_sub(fname,1,nchar(fname)-4)
  mp.i=allM %>% extract2(fname) %>% extract2("modelPredictions") %>%
    mutate(model=factor(name),almError=almPred.vx-vxb,examError=examPred.vx-vxb,
           almErrorAbs=abs(almError),examErrorAbs=abs(examError),
           bestModel=ifelse(abs(almError)<=abs(examError),"ALM","EXAM"),
           bestPred=ifelse(bestModel=="ALM",almPred.vx,examPred.vx),
           bestError=ifelse(bestModel=="ALM",almError,examError)) %>% 
    relocate(model,expMode,nGoodTrial,bestModel,bestPred,bestError,almError,examError,almPred.vx,examPred.vx,.after=condit)
  
  
  mpl[[i]]=mp.i
}


mp <- ldply(mpl,data.frame) %>% filter(model %in% mInclude) %>%
  mutate(model=recode_factor(model,"fitSbjTrain"="lr_c_w","noW_trainTrain"="lr_c","noW_trainTrainNf"="tnf_lr_c"))


predict_alm_exam_weighted_hybrid <- function(input, c, weight_exam, trainData, input.layer, output.layer, weight.mat) {
  
  alm_pred <- predict_alm_only(input, c, trainData, input.layer, output.layer, weight.mat)
  exam_pred <- predict_alm_exam(input, c, trainData, input.layer, output.layer, weight.mat)
  
  
  # Weighted combination of the predictions from ALM and EXAM
  (1 - weight_exam) * alm_pred + weight_exam * exam_pred
}

optimize_params_weighted <- function(dsAvg, c_values, lr_values, weight_exam_values, input.layer, output.layer) {
  best_error <- Inf
  best_params <- list(c = NA, lr = NA, weight_exam = NA)
  error_data <- data.frame() # Empty data frame to store errors
  
  # Split the data into training and testing sets
  train_data <- dsAvg[dsAvg$expMode2 == "Train", ]
  test_data <- dsAvg[dsAvg$expMode2 == "Test", ] |> group_by(x) |> summarise(y=mean(y))
  
  # Grid search over c, lr, and weight_exam values
  for (c in c_values) {
    for (lr in lr_values) {
      for (weight_exam in weight_exam_values) {
        # Train ALM with current parameters
        sim_results <- alm.sim(train_data, c, lr, input.layer, output.layer)
        weight.mat <- sim_results$wm
        
        # Test on the testing data using the weighted hybrid prediction function
        predictions <- sapply(test_data$x, function(input) {
          predict_alm_exam_weighted_hybrid(input, c, weight_exam, train_data, input.layer, output.layer, weight.mat)
        })
        
        # Compute mean squared error
        mse <- mean((predictions - test_data$y)^2)
        error_data <- rbind(error_data, data.frame(c=c, lr=lr, weight_exam=weight_exam, mse=mse))
        
        # Update best parameters if current error is lower
        if (!is.na(mse) && mse < best_error) {
          best_error <- mse
          best_params$c <- c
          best_params$lr <- lr
          best_params$weight_exam <- weight_exam
        }
      }
    }
  }
  
  list(best_params=best_params, best_error=best_error, error_data=error_data)
}