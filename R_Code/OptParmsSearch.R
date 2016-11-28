setwd("/home/nmoniz/TResampling/Code 2")

library(lubridate)
library(xts)
library(performanceEstimation)
library(uba)
library(UBL) # only for distance functions
library(DMwR)

# data sest generated with target variable always on the last column.
# This is mandatory for smote versions to work correctly.
# Function to create the embed
create.data<-function(ts,embed){
  t<-index(ts)[-(1:(embed-1))]
  e<-embed(ts,embed)[,embed:1]
  colnames(e)<-paste('V',1:embed,sep='')
  d <-xts(e,t)
  as.data.frame(d)
}

#EVALUATION STATISTICS (Utility-Based Regression Framework)
eval.stats <- function(form,train,test,preds,ph,ls) {
  trues <- resp(form,test)
  trainY <- resp(form,train)
  
  prec <- util(preds,trues,ph,ls,util.control(umetric="P",event.thr=0.9))
  rec  <- util(preds,trues,ph,ls,util.control(umetric="R",event.thr=0.9))
  F05  <- util(preds,trues,ph,ls,util.control(umetric="Fm",beta=0.5,event.thr=0.9))
  F1   <- util(preds,trues,ph,ls,util.control(umetric="Fm",beta=1,event.thr=0.9))
  F2   <- util(preds,trues,ph,ls,util.control(umetric="Fm",beta=2,event.thr=0.9))
  
  mad=mean(abs(trues-preds))
  mse=mean((trues-preds)^2)
  mape= mean((abs(trues-preds)/trues))*100
  rmse= sqrt(mean((trues-preds)^2))
  mae_phi= mean(phi(trues,control.parms=ph)*(abs(trues-preds)))
  mape_phi= mean(phi(trues,control.parms=ph)*(abs(trues-preds)/trues))*100
  mse_phi= mean(phi(trues,control.parms=ph)*(trues-preds)^2)
  rmse_phi= sqrt(mean(phi(trues,control.parms=ph)*(trues-preds)^2))
  prec=prec
  rec=rec
  F05=F05
  F1=F1
  F2=F2
  
  c(
    prec=prec,rec=rec,F1=F1
  )
  
}

mc.lm <- function(form,train,test,...) {
  ph <- uba::phi.control(train[,ncol(train)], method="extremes")
  ls <- uba::loss.control(train[,ncol(train)])
  m <- lm(form,train,na.action = NULL,...)
  p <- predict(m,test)
  eval <- eval.stats(form,train,test,p,ph,ls)
  res <- list(trues=responseValues(form,test),preds=p,evaluation=eval)
  res
}

mc.svm <- function(form,train,test,cost,gamma,...) {
  require(e1071)
  ph <- uba::phi.control(train[,ncol(train)], method="extremes")
  ls <- uba::loss.control(train[,ncol(train)])
  m <- svm(form,train,cost=cost,gamma=gamma,...)
  p <- predict(m,test)
  eval <- eval.stats(form,train,test,p,ph,ls)
  res <- list(trues=responseValues(form,test),preds=p,evaluation=eval)
  res
}

mc.mars <- function(form,train,test,nk,degree,thresh,...) {
  require(earth)
  ph <- uba::phi.control(train[,ncol(train)], method="extremes")
  ls <- uba::loss.control(train[,ncol(train)])
  m <- earth(form,train,nk=nk,degree=degree,thresh=thresh)
  p <- predict(m,test)
  p <- as.vector(p)
  names(p) <- rownames(test)
  eval <- eval.stats(form,train,test,p,ph,ls)
  res <- list(trues=responseValues(form,test),preds=p,evaluation=eval)
  res
}

mc.rf <- function(form,train,test,mtry,ntree,...) {
  require(randomForest)
  ph <- uba::phi.control(train[,ncol(train)], method="extremes")
  ls <- uba::loss.control(train[,ncol(train)])
  m <- randomForest(form,train,mtry=mtry,ntree=ntree,...)
  p <- predict(m,test)
  eval <- eval.stats(form,train,test,p,ph,ls)
  res <- list(trues=responseValues(form,test),preds=p,evaluation=eval)
  res
}

mc.rpart <- function(form,train,test,minsplit,cp,...) {
  require(rpart)
  ph <- uba::phi.control(train[,ncol(train)], method="extremes")
  ls <- uba::loss.control(train[,ncol(train)])
  m <- rpart(form,train,control=rpart.control(minsplit=minsplit,cp=cp),...)
  p <- predict(m,test)
  eval <- eval.stats(form,train,test,p,ph,ls)
  res <- list(trues=responseValues(form,test),preds=p,evaluation=eval)
  res
}

#####################

mc.svm_UNDERB <- function(form,train,test,cost,gamma,un,...) {
  require(e1071)
  ph <- uba::phi.control(train[,ncol(train)], method="extremes")
  ls <- uba::loss.control(train[,ncol(train)])
  train <- randUnderRegressB(form, train, rel="auto", thr.rel=0.9, C.perc=list(un), repl=FALSE)
  m <- svm(form,train,cost=cost,gamma=gamma,...)
  p <- predict(m,test)
  eval <- eval.stats(form,train,test,p,ph,ls)
  res <- list(evaluation=eval)
  res
}

mc.svm_UNDERT <- function(form,train,test,cost,gamma,un,...) {
  require(e1071)
  ph <- uba::phi.control(train[,ncol(train)], method="extremes")
  ls <- uba::loss.control(train[,ncol(train)])
  train <- randUnderRegressT(form, train, rel="auto", thr.rel=0.9, C.perc=list(un), repl=FALSE)
  m <- svm(form,train,cost=cost,gamma=gamma,...)
  p <- predict(m,test)
  eval <- eval.stats(form,train,test,p,ph,ls)
  res <- list(evaluation=eval)
  res
}

mc.svm_UNDERTPhi <- function(form,train,test,cost,gamma,un,...) {
  require(e1071)
  ph <- uba::phi.control(train[,ncol(train)], method="extremes")
  ls <- uba::loss.control(train[,ncol(train)])
  train <- randUnderRegressTPhi(form, train, rel="auto", thr.rel=0.9, C.perc=list(un), repl=FALSE)
  m <- svm(form,train,cost=cost,gamma=gamma,...)
  p <- predict(m,test)
  eval <- eval.stats(form,train,test,p,ph,ls)
  res <- list(evaluation=eval)
  res
}

mc.svm_OVERB <- function(form,train,test,cost,gamma,ov,...) {
  require(e1071)
  ph <- uba::phi.control(train[,ncol(train)], method="extremes")
  ls <- uba::loss.control(train[,ncol(train)])
  train <- randOverRegressB(form, train, rel="auto", thr.rel=0.9, C.perc=list(ov), repl=TRUE)
  m <- svm(form,train,cost=cost,gamma=gamma,...)
  p <- predict(m,test)
  eval <- eval.stats(form,train,test,p,ph,ls)
  res <- list(evaluation=eval)
  res
}

mc.svm_OVERT <- function(form,train,test,cost,gamma,ov,...) {
  require(e1071)
  ph <- uba::phi.control(train[,ncol(train)], method="extremes")
  ls <- uba::loss.control(train[,ncol(train)])
  train <- randOverRegressT(form, train, rel="auto", thr.rel=0.9, C.perc=list(ov), repl=TRUE)
  m <- svm(form,train,cost=cost,gamma=gamma,...)
  p <- predict(m,test)
  eval <- eval.stats(form,train,test,p,ph,ls)
  res <- list(evaluation=eval)
  res
}

mc.svm_OVERTPhi <- function(form,train,test,cost,gamma,ov,...) {
  require(e1071)
  ph <- uba::phi.control(train[,ncol(train)], method="extremes")
  ls <- uba::loss.control(train[,ncol(train)])
  train <- randOverRegressTPhi(form, train, rel="auto", thr.rel=0.9, C.perc=list(ov), repl=TRUE)
  m <- svm(form,train,cost=cost,gamma=gamma,...)
  p <- predict(m,test)
  eval <- eval.stats(form,train,test,p,ph,ls)
  res <- list(evaluation=eval)
  res
}

mc.svm_SMOTEB <- function(form,train,test,cost,gamma,un,ov,...) {
  require(e1071)
  ph <- uba::phi.control(train[,ncol(train)], method="extremes")
  ls <- uba::loss.control(train[,ncol(train)])
  train <- smoteRegressB(form, train, rel="auto", thr.rel=0.9, C.perc=list(un,ov), k=5, repl=TRUE, dist="Euclidean", p=2)
  m <- svm(form,train,cost=cost,gamma=gamma,...)
  p <- predict(m,test)
  eval <- eval.stats(form,train,test,p,ph,ls)
  res <- list(evaluation=eval)
  res
}

mc.svm_SMOTET <- function(form,train,test,cost,gamma,un,ov,...) {
  require(e1071)
  ph <- uba::phi.control(train[,ncol(train)], method="extremes")
  ls <- uba::loss.control(train[,ncol(train)])
  train <- smoteRegressT(form, train, rel="auto", thr.rel=0.9, C.perc=list(un,ov), k=5, repl=TRUE, dist="Euclidean", p=2)
  m <- svm(form,train,cost=cost,gamma=gamma,...)
  p <- predict(m,test)
  eval <- eval.stats(form,train,test,p,ph,ls)
  res <- list(evaluation=eval)
  res
}

mc.svm_SMOTETPhi <- function(form,train,test,cost,gamma,un,ov,...) {
  require(e1071)
  ph <- uba::phi.control(train[,ncol(train)], method="extremes")
  ls <- uba::loss.control(train[,ncol(train)])
  train <- smoteRegressTPhi(form, train, rel="auto", thr.rel=0.9, C.perc=list(un,ov), k=5, repl=TRUE, dist="Euclidean", p=2)
  m <- svm(form,train,cost=cost,gamma=gamma,...)
  p <- predict(m,test)
  eval <- eval.stats(form,train,test,p,ph,ls)
  res <- list(evaluation=eval)
  res
}

exp.best <- performanceEstimation(PredTask(form,ds),
                             c(Workflow("mc.svm",cost=150,gamma=0.01),
                               workflowVariants("mc.svm_UNDERB",cost=c(10,150,300),gamma=c(0.01,0.001),un=c(.1,.2,.4,.6,.8)),
                               workflowVariants("mc.svm_UNDERT",cost=c(10,150,300),gamma=c(0.01,0.001),un=c(.1,.2,.4,.6,.8)),
                               workflowVariants("mc.svm_UNDERTPhi",cost=c(10,150,300),gamma=c(0.01,0.001),un=c(.1,.2,.4,.6,.8)),
                               workflowVariants("mc.svm_OVERB",cost=c(10,150,300),gamma=c(0.01,0.001),ov=c(2,3,5,10)),
                               workflowVariants("mc.svm_OVERT",cost=c(10,150,300),gamma=c(0.01,0.001),ov=c(2,3,5,10)),
                               workflowVariants("mc.svm_OVERTPhi",cost=c(10,150,300),gamma=c(0.01,0.001),ov=c(2,3,5,10)),
                               workflowVariants("mc.svm_SMOTEB",cost=c(10,150,300),gamma=c(0.01,0.001),un=c(.05,.1,.2,.4,.6,.8),ov=c(2,3,5,10)),
                               workflowVariants("mc.svm_SMOTET",cost=c(10,150,300),gamma=c(0.01,0.001),un=c(.05,.1,.2,.4,.6,.8),ov=c(2,3,5,10)),
                               workflowVariants("mc.svm_SMOTETPhi",cost=c(10,150,300),gamma=c(0.01,0.001),un=c(.05,.1,.2,.4,.6,.8),ov=c(2,3,5,10))),
                             EstimationTask("totTime",method=MonteCarlo(nReps=10,szTrain=.5,szTest=.25))
)

globalres <- data.frame(model=character(0),prec=numeric(0),rec=numeric(0),F1=numeric(0))

models <- c("mc.svm",paste0("mc.svm_UNDERB",seq(1,30,1)),paste0("mc.svm_UNDERT",seq(1,30,1)),paste0("mc.svm_UNDERTPhi",seq(1,30,1)),
            paste0("mc.svm_OVERB",seq(1,24,1)),paste0("mc.svm_OVERT",seq(1,24,1)),paste0("mc.svm_OVERTPhi",seq(1,24,1)),
            paste0("mc.svm_SMOTEB",seq(1,144,1)),paste0("mc.svm_SMOTET",seq(1,144,1)),paste0("mc.svm_SMOTETPhi",seq(1,144,1)))

for(wf in 1:length(models)) {
  res <- c()
  for(i in 1:10) {
    res <- rbind(res,getIterationsInfo(exp.best,workflow=wf,task=1,it=i)$evaluation)
  }
  
  exp.best[[1]][[wf]]@iterationsScores <- res
  exp.best[[1]][[wf]]@estTask@metrics <- c("prec","rec","F1")
  
  globalres.row <- data.frame(model=models[wf],prec=mean(res[,1]),rec=mean(res[,2]),F1=median(res[,3]))  
  globalres <- rbind(globalres,globalres.row)
}

globalres
