#WORKFLOWS
mc.lm <- function(form,train,test,...) {
  ph <- uba::phi.control(train[,ncol(train)], method="extremes")
  ls <- uba::loss.control(train[,ncol(train)])
  ptm <- proc.time()
  m <- lm(form,train,...)
  train.time <- proc.time() - ptm
  p <- predict(m,test)
  trainpred.time <- proc.time() - ptm
  eval <- eval.stats(form,train,test,p,ph,ls)
  res <- list(evaluation=eval,traintime=train.time,trainpredtime=trainpred.time)
  res
}

mc.lm_UNDERB <- function(form,train,test,...) {
  ph <- uba::phi.control(train[,ncol(train)], method="extremes")
  ls <- uba::loss.control(train[,ncol(train)])
  train <- randUnderRegressB(form, train, rel="auto", thr.rel=0.9, C.perc="balance", repl=FALSE)
  ptm <- proc.time()
  m <- lm(form,train,...)
  train.time <- proc.time() - ptm
  p <- predict(m,test)
  trainpred.time <- proc.time() - ptm
  eval <- eval.stats(form,train,test,p,ph,ls)
  res <- list(evaluation=eval,traintime=train.time,trainpredtime=trainpred.time)
  res
}

mc.lm_UNDERT <- function(form,train,test,...) {
  ph <- uba::phi.control(train[,ncol(train)], method="extremes")
  ls <- uba::loss.control(train[,ncol(train)])
  train <- randUnderRegressT(form, train, rel="auto", thr.rel=0.9, C.perc="balance", repl=FALSE)
  ptm <- proc.time()
  m <- lm(form,train,...)
  train.time <- proc.time() - ptm
  p <- predict(m,test)
  trainpred.time <- proc.time() - ptm
  eval <- eval.stats(form,train,test,p,ph,ls)
  res <- list(evaluation=eval,traintime=train.time,trainpredtime=trainpred.time)
  res
}

mc.lm_UNDERTPhi <- function(form,train,test,...) {
  ph <- uba::phi.control(train[,ncol(train)], method="extremes")
  ls <- uba::loss.control(train[,ncol(train)])
  train <- randUnderRegressTPhi(form, train, rel="auto", thr.rel=0.9, C.perc="balance", repl=FALSE)
  ptm <- proc.time()
  m <- lm(form,train,...)
  train.time <- proc.time() - ptm
  p <- predict(m,test)
  trainpred.time <- proc.time() - ptm
  eval <- eval.stats(form,train,test,p,ph,ls)
  res <- list(evaluation=eval,traintime=train.time,trainpredtime=trainpred.time)
  res
}

mc.lm_OVERB <- function(form,train,test,...) {
  ph <- uba::phi.control(train[,ncol(train)], method="extremes")
  ls <- uba::loss.control(train[,ncol(train)])
  train <- randOverRegressB(form, train, rel="auto", thr.rel=0.9, C.perc="balance", repl=TRUE)
  ptm <- proc.time()
  m <- lm(form,train,...)
  train.time <- proc.time() - ptm
  p <- predict(m,test)
  trainpred.time <- proc.time() - ptm
  eval <- eval.stats(form,train,test,p,ph,ls)
  res <- list(evaluation=eval,traintime=train.time,trainpredtime=trainpred.time)
  res
}

mc.lm_OVERT <- function(form,train,test,...) {
  ph <- uba::phi.control(train[,ncol(train)], method="extremes")
  ls <- uba::loss.control(train[,ncol(train)])
  train <- randOverRegressT(form, train, rel="auto", thr.rel=0.9, C.perc="balance", repl=TRUE)
  ptm <- proc.time()
  m <- lm(form,train,...)
  train.time <- proc.time() - ptm
  p <- predict(m,test)
  trainpred.time <- proc.time() - ptm
  eval <- eval.stats(form,train,test,p,ph,ls)
  res <- list(evaluation=eval,traintime=train.time,trainpredtime=trainpred.time)
  res
}

mc.lm_OVERTPhi <- function(form,train,test,...) {
  ph <- uba::phi.control(train[,ncol(train)], method="extremes")
  ls <- uba::loss.control(train[,ncol(train)])
  train <- randOverRegressTPhi(form, train, rel="auto", thr.rel=0.9, C.perc="balance", repl=TRUE)
  ptm <- proc.time()
  m <- lm(form,train,...)
  train.time <- proc.time() - ptm
  p <- predict(m,test)
  trainpred.time <- proc.time() - ptm
  eval <- eval.stats(form,train,test,p,ph,ls)
  res <- list(evaluation=eval,traintime=train.time,trainpredtime=trainpred.time)
  res
}

mc.lm_SMOTEB <- function(form,train,test,...) {
  ph <- uba::phi.control(train[,ncol(train)], method="extremes")
  ls <- uba::loss.control(train[,ncol(train)])
  train <- smoteRegressB(form, train, rel="auto", thr.rel=0.9, C.perc="balance", k=5, repl=TRUE, dist="Euclidean", p=2)
  ptm <- proc.time()
  m <- lm(form,train,...)
  train.time <- proc.time() - ptm
  p <- predict(m,test)
  trainpred.time <- proc.time() - ptm
  eval <- eval.stats(form,train,test,p,ph,ls)
  res <- list(evaluation=eval,traintime=train.time,trainpredtime=trainpred.time)
  res
}

mc.lm_SMOTET <- function(form,train,test,...) {
  ph <- uba::phi.control(train[,ncol(train)], method="extremes")
  ls <- uba::loss.control(train[,ncol(train)])
  train <- smoteRegressT(form, train, rel="auto", thr.rel=0.9, C.perc="balance", k=5, repl=TRUE, dist="Euclidean", p=2)
  ptm <- proc.time()
  m <- lm(form,train,...)
  train.time <- proc.time() - ptm
  p <- predict(m,test)
  trainpred.time <- proc.time() - ptm
  eval <- eval.stats(form,train,test,p,ph,ls)
  res <- list(evaluation=eval,traintime=train.time,trainpredtime=trainpred.time)
  res
}

mc.lm_SMOTETPhi <- function(form,train,test,...) {
  ph <- uba::phi.control(train[,ncol(train)], method="extremes")
  ls <- uba::loss.control(train[,ncol(train)])
  train <- smoteRegressTPhi(form, train, rel="auto", thr.rel=0.9, C.perc="balance", k=5, repl=TRUE, dist="Euclidean", p=2)
  ptm <- proc.time()
  m <- lm(form,train,...)
  train.time <- proc.time() - ptm
  p <- predict(m,test)
  trainpred.time <- proc.time() - ptm
  eval <- eval.stats(form,train,test,p,ph,ls)
  res <- list(evaluation=eval,traintime=train.time,trainpredtime=trainpred.time)
  res
}

mc.rpart <- function(form,train,test,minsplit,cp,...) {
  require(rpart)
  ph <- uba::phi.control(train[,ncol(train)], method="extremes")
  ls <- uba::loss.control(train[,ncol(train)])
  ptm <- proc.time()
  m <- rpart(form,train,control=rpart.control(minsplit=minsplit,cp=cp),...)
  train.time <- proc.time() - ptm
  p <- predict(m,test)
  trainpred.time <- proc.time() - ptm
  eval <- eval.stats(form,train,test,p,ph,ls)
  res <- list(evaluation=eval,traintime=train.time,trainpredtime=trainpred.time)
  res
}

mc.rpart_UNDERB <- function(form,train,test,minsplit,cp,...) {
  require(rpart)
  ph <- uba::phi.control(train[,ncol(train)], method="extremes")
  ls <- uba::loss.control(train[,ncol(train)])
  train <- randUnderRegressB(form, train, rel="auto", thr.rel=0.9, C.perc="balance", repl=FALSE)
  ptm <- proc.time()
  m <- rpart(form,train,control=rpart.control(minsplit=minsplit,cp=cp),...)
  train.time <- proc.time() - ptm
  p <- predict(m,test)
  trainpred.time <- proc.time() - ptm
  eval <- eval.stats(form,train,test,p,ph,ls)
  res <- list(evaluation=eval,traintime=train.time,trainpredtime=trainpred.time)
  res
}

mc.rpart_UNDERT <- function(form,train,test,minsplit,cp,...) {
  require(rpart)
  ph <- uba::phi.control(train[,ncol(train)], method="extremes")
  ls <- uba::loss.control(train[,ncol(train)])
  train <- randUnderRegressT(form, train, rel="auto", thr.rel=0.9, C.perc="balance", repl=FALSE)
  ptm <- proc.time()
  m <- rpart(form,train,control=rpart.control(minsplit=minsplit,cp=cp),...)
  train.time <- proc.time() - ptm
  p <- predict(m,test)
  trainpred.time <- proc.time() - ptm
  eval <- eval.stats(form,train,test,p,ph,ls)
  res <- list(evaluation=eval,traintime=train.time,trainpredtime=trainpred.time)
  res
}

mc.rpart_UNDERTPhi <- function(form,train,test,minsplit,cp,...) {
  require(rpart)
  ph <- uba::phi.control(train[,ncol(train)], method="extremes")
  ls <- uba::loss.control(train[,ncol(train)])
  train <- randUnderRegressTPhi(form, train, rel="auto", thr.rel=0.9, C.perc="balance", repl=FALSE)
  ptm <- proc.time()
  m <- rpart(form,train,control=rpart.control(minsplit=minsplit,cp=cp),...)
  train.time <- proc.time() - ptm
  p <- predict(m,test)
  trainpred.time <- proc.time() - ptm
  eval <- eval.stats(form,train,test,p,ph,ls)
  res <- list(evaluation=eval,traintime=train.time,trainpredtime=trainpred.time)
  res
}

mc.rpart_OVERB <- function(form,train,test,minsplit,cp,...) {
  require(rpart)
  ph <- uba::phi.control(train[,ncol(train)], method="extremes")
  ls <- uba::loss.control(train[,ncol(train)])
  train <- randOverRegressB(form, train, rel="auto", thr.rel=0.9, C.perc="balance", repl=TRUE)
  ptm <- proc.time()
  m <- rpart(form,train,control=rpart.control(minsplit=minsplit,cp=cp),...)
  train.time <- proc.time() - ptm
  p <- predict(m,test)
  trainpred.time <- proc.time() - ptm
  eval <- eval.stats(form,train,test,p,ph,ls)
  res <- list(evaluation=eval,traintime=train.time,trainpredtime=trainpred.time)
  res
}

mc.rpart_OVERT <- function(form,train,test,minsplit,cp,...) {
  require(rpart)
  ph <- uba::phi.control(train[,ncol(train)], method="extremes")
  ls <- uba::loss.control(train[,ncol(train)])
  train <- randOverRegressT(form, train, rel="auto", thr.rel=0.9, C.perc="balance", repl=TRUE)
  ptm <- proc.time()
  m <- rpart(form,train,control=rpart.control(minsplit=minsplit,cp=cp),...)
  train.time <- proc.time() - ptm
  p <- predict(m,test)
  trainpred.time <- proc.time() - ptm
  eval <- eval.stats(form,train,test,p,ph,ls)
  res <- list(evaluation=eval,traintime=train.time,trainpredtime=trainpred.time)
  res
}

mc.rpart_OVERTPhi <- function(form,train,test,minsplit,cp,...) {
  require(rpart)
  ph <- uba::phi.control(train[,ncol(train)], method="extremes")
  ls <- uba::loss.control(train[,ncol(train)])
  train <- randOverRegressTPhi(form, train, rel="auto", thr.rel=0.9, C.perc="balance", repl=TRUE)
  ptm <- proc.time()
  m <- rpart(form,train,control=rpart.control(minsplit=minsplit,cp=cp),...)
  train.time <- proc.time() - ptm
  p <- predict(m,test)
  trainpred.time <- proc.time() - ptm
  eval <- eval.stats(form,train,test,p,ph,ls)
  res <- list(evaluation=eval,traintime=train.time,trainpredtime=trainpred.time)
  res
}

mc.rpart_SMOTEB <- function(form,train,test,minsplit,cp,...) {
  require(rpart)
  ph <- uba::phi.control(train[,ncol(train)], method="extremes")
  ls <- uba::loss.control(train[,ncol(train)])
  train <- smoteRegressB(form, train, rel="auto", thr.rel=0.9, C.perc="balance", k=5, repl=TRUE, dist="Euclidean", p=2)
  ptm <- proc.time()
  m <- rpart(form,train,control=rpart.control(minsplit=minsplit,cp=cp),...)
  train.time <- proc.time() - ptm
  p <- predict(m,test)
  trainpred.time <- proc.time() - ptm
  eval <- eval.stats(form,train,test,p,ph,ls)
  res <- list(evaluation=eval,traintime=train.time,trainpredtime=trainpred.time)
  res
}

mc.rpart_SMOTET <- function(form,train,test,minsplit,cp,...) {
  require(rpart)
  ph <- uba::phi.control(train[,ncol(train)], method="extremes")
  ls <- uba::loss.control(train[,ncol(train)])
  train <- smoteRegressT(form, train, rel="auto", thr.rel=0.9, C.perc="balance", k=5, repl=TRUE, dist="Euclidean", p=2)
  ptm <- proc.time()
  m <- rpart(form,train,control=rpart.control(minsplit=minsplit,cp=cp),...)
  train.time <- proc.time() - ptm
  p <- predict(m,test)
  trainpred.time <- proc.time() - ptm
  eval <- eval.stats(form,train,test,p,ph,ls)
  res <- list(evaluation=eval,traintime=train.time,trainpredtime=trainpred.time)
  res
}

mc.rpart_SMOTETPhi <- function(form,train,test,minsplit,cp,...) {
  require(rpart)
  ph <- uba::phi.control(train[,ncol(train)], method="extremes")
  ls <- uba::loss.control(train[,ncol(train)])
  train <- smoteRegressTPhi(form, train, rel="auto", thr.rel=0.9, C.perc="balance", k=5, repl=TRUE, dist="Euclidean", p=2)
  ptm <- proc.time()
  m <- rpart(form,train,control=rpart.control(minsplit=minsplit,cp=cp),...)
  train.time <- proc.time() - ptm
  p <- predict(m,test)
  trainpred.time <- proc.time() - ptm
  eval <- eval.stats(form,train,test,p,ph,ls)
  res <- list(evaluation=eval,traintime=train.time,trainpredtime=trainpred.time)
  res
}

mc.svm <- function(form,train,test,cost,gamma,...) {
  require(e1071)
  ph <- uba::phi.control(train[,ncol(train)], method="extremes")
  ls <- uba::loss.control(train[,ncol(train)])
  ptm <- proc.time()
  m <- svm(form,train,cost=cost,gamma=gamma,...)
  train.time <- proc.time() - ptm
  p <- predict(m,test)
  trainpred.time <- proc.time() - ptm
  eval <- eval.stats(form,train,test,p,ph,ls)
  res <- list(evaluation=eval,traintime=train.time,trainpredtime=trainpred.time)
  res
}

mc.svm_UNDERB <- function(form,train,test,cost,gamma,...) {
  require(e1071)
  ph <- uba::phi.control(train[,ncol(train)], method="extremes")
  ls <- uba::loss.control(train[,ncol(train)])
  train <- randUnderRegressB(form, train, rel="auto", thr.rel=0.9, C.perc="balance", repl=FALSE)
  ptm <- proc.time()
  m <- svm(form,train,cost=cost,gamma=gamma,...)
  train.time <- proc.time() - ptm
  p <- predict(m,test)
  trainpred.time <- proc.time() - ptm
  eval <- eval.stats(form,train,test,p,ph,ls)
  res <- list(evaluation=eval,traintime=train.time,trainpredtime=trainpred.time)
  res
}

mc.svm_UNDERT <- function(form,train,test,cost,gamma,...) {
  require(e1071)
  ph <- uba::phi.control(train[,ncol(train)], method="extremes")
  ls <- uba::loss.control(train[,ncol(train)])
  train <- randUnderRegressT(form, train, rel="auto", thr.rel=0.9, C.perc="balance", repl=FALSE)
  ptm <- proc.time()
  m <- svm(form,train,cost=cost,gamma=gamma,...)
  train.time <- proc.time() - ptm
  p <- predict(m,test)
  trainpred.time <- proc.time() - ptm
  eval <- eval.stats(form,train,test,p,ph,ls)
  res <- list(evaluation=eval,traintime=train.time,trainpredtime=trainpred.time)
  res
}

mc.svm_UNDERTPhi <- function(form,train,test,cost,gamma,...) {
  require(e1071)
  ph <- uba::phi.control(train[,ncol(train)], method="extremes")
  ls <- uba::loss.control(train[,ncol(train)])
  train <- randUnderRegressTPhi(form, train, rel="auto", thr.rel=0.9, C.perc="balance", repl=FALSE)
  ptm <- proc.time()
  m <- svm(form,train,cost=cost,gamma=gamma,...)
  train.time <- proc.time() - ptm
  p <- predict(m,test)
  trainpred.time <- proc.time() - ptm
  eval <- eval.stats(form,train,test,p,ph,ls)
  res <- list(evaluation=eval,traintime=train.time,trainpredtime=trainpred.time)
  res
}

mc.svm_OVERB <- function(form,train,test,cost,gamma,...) {
  require(e1071)
  ph <- uba::phi.control(train[,ncol(train)], method="extremes")
  ls <- uba::loss.control(train[,ncol(train)])
  train <- randOverRegressB(form, train, rel="auto", thr.rel=0.9, C.perc="balance", repl=TRUE)
  ptm <- proc.time()
  m <- svm(form,train,cost=cost,gamma=gamma,...)
  train.time <- proc.time() - ptm
  p <- predict(m,test)
  trainpred.time <- proc.time() - ptm
  eval <- eval.stats(form,train,test,p,ph,ls)
  res <- list(evaluation=eval,traintime=train.time,trainpredtime=trainpred.time)
  res
}

mc.svm_OVERT <- function(form,train,test,cost,gamma,...) {
  require(e1071)
  ph <- uba::phi.control(train[,ncol(train)], method="extremes")
  ls <- uba::loss.control(train[,ncol(train)])
  train <- randOverRegressT(form, train, rel="auto", thr.rel=0.9, C.perc="balance", repl=TRUE)
  ptm <- proc.time()
  m <- svm(form,train,cost=cost,gamma=gamma,...)
  train.time <- proc.time() - ptm
  p <- predict(m,test)
  trainpred.time <- proc.time() - ptm
  eval <- eval.stats(form,train,test,p,ph,ls)
  res <- list(evaluation=eval,traintime=train.time,trainpredtime=trainpred.time)
  res
}

mc.svm_OVERTPhi <- function(form,train,test,cost,gamma,...) {
  require(e1071)
  ph <- uba::phi.control(train[,ncol(train)], method="extremes")
  ls <- uba::loss.control(train[,ncol(train)])
  train <- randOverRegressTPhi(form, train, rel="auto", thr.rel=0.9, C.perc="balance", repl=TRUE)
  ptm <- proc.time()
  m <- svm(form,train,cost=cost,gamma=gamma,...)
  train.time <- proc.time() - ptm
  p <- predict(m,test)
  trainpred.time <- proc.time() - ptm
  eval <- eval.stats(form,train,test,p,ph,ls)
  res <- list(evaluation=eval,traintime=train.time,trainpredtime=trainpred.time)
  res
}

mc.svm_SMOTEB <- function(form,train,test,cost,gamma,...) {
  require(e1071)
  ph <- uba::phi.control(train[,ncol(train)], method="extremes")
  ls <- uba::loss.control(train[,ncol(train)])
  train <- smoteRegressB(form, train, rel="auto", thr.rel=0.9, C.perc="balance", k=5, repl=TRUE, dist="Euclidean", p=2)
  ptm <- proc.time()
  m <- svm(form,train,cost=cost,gamma=gamma,...)
  train.time <- proc.time() - ptm
  p <- predict(m,test)
  trainpred.time <- proc.time() - ptm
  eval <- eval.stats(form,train,test,p,ph,ls)
  res <- list(evaluation=eval,traintime=train.time,trainpredtime=trainpred.time)
  res
}

mc.svm_SMOTET <- function(form,train,test,cost,gamma,...) {
  require(e1071)
  ph <- uba::phi.control(train[,ncol(train)], method="extremes")
  ls <- uba::loss.control(train[,ncol(train)])
  train <- smoteRegressT(form, train, rel="auto", thr.rel=0.9, C.perc="balance", k=5, repl=TRUE, dist="Euclidean", p=2)
  ptm <- proc.time()
  m <- svm(form,train,cost=cost,gamma=gamma,...)
  train.time <- proc.time() - ptm
  p <- predict(m,test)
  trainpred.time <- proc.time() - ptm
  eval <- eval.stats(form,train,test,p,ph,ls)
  res <- list(evaluation=eval,traintime=train.time,trainpredtime=trainpred.time)
  res
}

mc.svm_SMOTETPhi <- function(form,train,test,cost,gamma,...) {
  require(e1071)
  ph <- uba::phi.control(train[,ncol(train)], method="extremes")
  ls <- uba::loss.control(train[,ncol(train)])
  train <- smoteRegressTPhi(form, train, rel="auto", thr.rel=0.9, C.perc="balance", k=5, repl=TRUE, dist="Euclidean", p=2)
  ptm <- proc.time()
  m <- svm(form,train,cost=cost,gamma=gamma,...)
  train.time <- proc.time() - ptm
  p <- predict(m,test)
  trainpred.time <- proc.time() - ptm
  eval <- eval.stats(form,train,test,p,ph,ls)
  res <- list(evaluation=eval,traintime=train.time,trainpredtime=trainpred.time)
  res
}

mc.mars <- function(form,train,test,nk,degree,thresh,...) {
  require(earth)
  ph <- uba::phi.control(train[,ncol(train)], method="extremes")
  ls <- uba::loss.control(train[,ncol(train)])
  ptm <- proc.time()
  m <- earth(form,train,nk=nk,degree=degree,thresh=thresh)
  train.time <- proc.time() - ptm
  p <- predict(m,test)
  trainpred.time <- proc.time() - ptm
  p <- as.vector(p)
  names(p) <- rownames(test)
  eval <- eval.stats(form,train,test,p,ph,ls)
  res <- list(evaluation=eval,traintime=train.time,trainpredtime=trainpred.time)
  res
}

mc.mars_UNDERB <- function(form,train,test,nk,degree,thresh,...) {
  require(earth)
  ph <- uba::phi.control(train[,ncol(train)], method="extremes")
  ls <- uba::loss.control(train[,ncol(train)])
  train <- randUnderRegressB(form, train, rel="auto", thr.rel=0.9, C.perc="balance", repl=FALSE)
  ptm <- proc.time()
  m <- earth(form,train,nk=nk,degree=degree,thresh=thresh)
  train.time <- proc.time() - ptm
  p <- predict(m,test)
  trainpred.time <- proc.time() - ptm
  p <- as.vector(p)
  names(p) <- rownames(test)
  eval <- eval.stats(form,train,test,p,ph,ls)
  res <- list(evaluation=eval,traintime=train.time,trainpredtime=trainpred.time)
  res
}

mc.mars_UNDERT <- function(form,train,test,nk,degree,thresh,...) {
  require(earth)
  ph <- uba::phi.control(train[,ncol(train)], method="extremes")
  ls <- uba::loss.control(train[,ncol(train)])
  train <- randUnderRegressT(form, train, rel="auto", thr.rel=0.9, C.perc="balance", repl=FALSE)
  ptm <- proc.time()
  m <- earth(form,train,nk=nk,degree=degree,thresh=thresh)
  train.time <- proc.time() - ptm
  p <- predict(m,test)
  trainpred.time <- proc.time() - ptm
  p <- as.vector(p)
  names(p) <- rownames(test)
  eval <- eval.stats(form,train,test,p,ph,ls)
  res <- list(evaluation=eval,traintime=train.time,trainpredtime=trainpred.time)
  res
}

mc.mars_UNDERTPhi <- function(form,train,test,nk,degree,thresh,...) {
  require(earth)
  ph <- uba::phi.control(train[,ncol(train)], method="extremes")
  ls <- uba::loss.control(train[,ncol(train)])
  train <- randUnderRegressTPhi(form, train, rel="auto", thr.rel=0.9, C.perc="balance", repl=FALSE)
  ptm <- proc.time()
  m <- earth(form,train,nk=nk,degree=degree,thresh=thresh)
  train.time <- proc.time() - ptm
  p <- predict(m,test)
  trainpred.time <- proc.time() - ptm
  p <- as.vector(p)
  names(p) <- rownames(test)
  eval <- eval.stats(form,train,test,p,ph,ls)
  res <- list(evaluation=eval,traintime=train.time,trainpredtime=trainpred.time)
  res
}

mc.mars_OVERB <- function(form,train,test,nk,degree,thresh,...) {
  require(earth)
  ph <- uba::phi.control(train[,ncol(train)], method="extremes")
  ls <- uba::loss.control(train[,ncol(train)])
  train <- randOverRegressB(form, train, rel="auto", thr.rel=0.9, C.perc="balance", repl=TRUE)
  ptm <- proc.time()
  m <- earth(form,train,nk=nk,degree=degree,thresh=thresh)
  train.time <- proc.time() - ptm
  p <- predict(m,test)
  trainpred.time <- proc.time() - ptm
  p <- as.vector(p)
  names(p) <- rownames(test)
  eval <- eval.stats(form,train,test,p,ph,ls)
  res <- list(evaluation=eval,traintime=train.time,trainpredtime=trainpred.time)
  res
}

mc.mars_OVERT <- function(form,train,test,nk,degree,thresh,...) {
  require(earth)
  ph <- uba::phi.control(train[,ncol(train)], method="extremes")
  ls <- uba::loss.control(train[,ncol(train)])
  train <- randOverRegressT(form, train, rel="auto", thr.rel=0.9, C.perc="balance", repl=TRUE)
  ptm <- proc.time()
  m <- earth(form,train,nk=nk,degree=degree,thresh=thresh)
  train.time <- proc.time() - ptm
  p <- predict(m,test)
  trainpred.time <- proc.time() - ptm
  p <- as.vector(p)
  names(p) <- rownames(test)
  eval <- eval.stats(form,train,test,p,ph,ls)
  res <- list(evaluation=eval,traintime=train.time,trainpredtime=trainpred.time)
  res
}

mc.mars_OVERTPhi <- function(form,train,test,nk,degree,thresh,...) {
  require(earth)
  ph <- uba::phi.control(train[,ncol(train)], method="extremes")
  ls <- uba::loss.control(train[,ncol(train)])
  train <- randOverRegressTPhi(form, train, rel="auto", thr.rel=0.9, C.perc="balance", repl=TRUE)
  ptm <- proc.time()
  m <- earth(form,train,nk=nk,degree=degree,thresh=thresh)
  train.time <- proc.time() - ptm
  p <- predict(m,test)
  trainpred.time <- proc.time() - ptm
  p <- as.vector(p)
  names(p) <- rownames(test)
  eval <- eval.stats(form,train,test,p,ph,ls)
  res <- list(evaluation=eval,traintime=train.time,trainpredtime=trainpred.time)
  res
}

mc.mars_SMOTEB <- function(form,train,test,nk,degree,thresh,...) {
  require(earth)
  ph <- uba::phi.control(train[,ncol(train)], method="extremes")
  ls <- uba::loss.control(train[,ncol(train)])
  train <- smoteRegressB(form, train, rel="auto", thr.rel=0.9, C.perc="balance", k=5, repl=TRUE, dist="Euclidean", p=2)
  ptm <- proc.time()
  m <- earth(form,train,nk=nk,degree=degree,thresh=thresh)
  train.time <- proc.time() - ptm
  p <- predict(m,test)
  trainpred.time <- proc.time() - ptm
  p <- as.vector(p)
  names(p) <- rownames(test)
  eval <- eval.stats(form,train,test,p,ph,ls)
  res <- list(evaluation=eval,traintime=train.time,trainpredtime=trainpred.time)
  res
}

mc.mars_SMOTET <- function(form,train,test,nk,degree,thresh,...) {
  require(earth)
  ph <- uba::phi.control(train[,ncol(train)], method="extremes")
  ls <- uba::loss.control(train[,ncol(train)])
  train <- smoteRegressT(form, train, rel="auto", thr.rel=0.9, C.perc="balance", k=5, repl=TRUE, dist="Euclidean", p=2)
  ptm <- proc.time()
  m <- earth(form,train,nk=nk,degree=degree,thresh=thresh)
  train.time <- proc.time() - ptm
  p <- predict(m,test)
  trainpred.time <- proc.time() - ptm
  p <- as.vector(p)
  names(p) <- rownames(test)
  eval <- eval.stats(form,train,test,p,ph,ls)
  res <- list(evaluation=eval,traintime=train.time,trainpredtime=trainpred.time)
  res
}

mc.mars_SMOTETPhi <- function(form,train,test,nk,degree,thresh,...) {
  require(earth)
  ph <- uba::phi.control(train[,ncol(train)], method="extremes")
  ls <- uba::loss.control(train[,ncol(train)])
  train <- smoteRegressTPhi(form, train, rel="auto", thr.rel=0.9, C.perc="balance", k=5, repl=TRUE, dist="Euclidean", p=2)
  ptm <- proc.time()
  m <- earth(form,train,nk=nk,degree=degree,thresh=thresh)
  train.time <- proc.time() - ptm
  p <- predict(m,test)
  trainpred.time <- proc.time() - ptm
  p <- as.vector(p)
  names(p) <- rownames(test)
  eval <- eval.stats(form,train,test,p,ph,ls)
  res <- list(evaluation=eval,traintime=train.time,trainpredtime=trainpred.time)
  res
}


mc.rf <- function(form,train,test,mtry,ntree,...) {
  require(randomForest)
  ph <- uba::phi.control(train[,ncol(train)], method="extremes")
  ls <- uba::loss.control(train[,ncol(train)])
  ptm <- proc.time()
  m <- randomForest(form,train,mtry=mtry,ntree=ntree,...)
  train.time <- proc.time() - ptm
  p <- predict(m,test)
  trainpred.time <- proc.time() - ptm
  eval <- eval.stats(form,train,test,p,ph,ls)
  res <- list(evaluation=eval,traintime=train.time,trainpredtime=trainpred.time)
  res
}

mc.rf_UNDERB <- function(form,train,test,mtry,ntree,...) {
  require(randomForest)
  ph <- uba::phi.control(train[,ncol(train)], method="extremes")
  ls <- uba::loss.control(train[,ncol(train)])
  train <- randUnderRegressB(form, train, rel="auto", thr.rel=0.9, C.perc="balance", repl=FALSE)
  ptm <- proc.time()
  m <- randomForest(form,train,mtry=mtry,ntree=ntree,...)
  train.time <- proc.time() - ptm
  p <- predict(m,test)
  trainpred.time <- proc.time() - ptm
  eval <- eval.stats(form,train,test,p,ph,ls)
  res <- list(evaluation=eval,traintime=train.time,trainpredtime=trainpred.time)
  res
}

mc.rf_UNDERT <- function(form,train,test,mtry,ntree,...) {
  require(randomForest)
  ph <- uba::phi.control(train[,ncol(train)], method="extremes")
  ls <- uba::loss.control(train[,ncol(train)])
  train <- randUnderRegressT(form, train, rel="auto", thr.rel=0.9, C.perc="balance", repl=FALSE)
  ptm <- proc.time()
  m <- randomForest(form,train,mtry=mtry,ntree=ntree,...)
  train.time <- proc.time() - ptm
  p <- predict(m,test)
  trainpred.time <- proc.time() - ptm
  eval <- eval.stats(form,train,test,p,ph,ls)
  res <- list(evaluation=eval,traintime=train.time,trainpredtime=trainpred.time)
  res
}

mc.rf_UNDERTPhi <- function(form,train,test,mtry,ntree,...) {
  require(randomForest)
  ph <- uba::phi.control(train[,ncol(train)], method="extremes")
  ls <- uba::loss.control(train[,ncol(train)])
  train <- randUnderRegressTPhi(form, train, rel="auto", thr.rel=0.9, C.perc="balance", repl=FALSE)
  ptm <- proc.time()
  m <- randomForest(form,train,mtry=mtry,ntree=ntree,...)
  train.time <- proc.time() - ptm
  p <- predict(m,test)
  trainpred.time <- proc.time() - ptm
  eval <- eval.stats(form,train,test,p,ph,ls)
  res <- list(evaluation=eval,traintime=train.time,trainpredtime=trainpred.time)
  res
}

mc.rf_OVERB <- function(form,train,test,mtry,ntree,...) {
  require(randomForest)
  ph <- uba::phi.control(train[,ncol(train)], method="extremes")
  ls <- uba::loss.control(train[,ncol(train)])
  train <- randOverRegressB(form, train, rel="auto", thr.rel=0.9, C.perc="balance", repl=TRUE)
  ptm <- proc.time()
  m <- randomForest(form,train,mtry=mtry,ntree=ntree,...)
  train.time <- proc.time() - ptm
  p <- predict(m,test)
  trainpred.time <- proc.time() - ptm
  eval <- eval.stats(form,train,test,p,ph,ls)
  res <- list(evaluation=eval,traintime=train.time,trainpredtime=trainpred.time)
  res
}

mc.rf_OVERT <- function(form,train,test,mtry,ntree,...) {
  require(randomForest)
  ph <- uba::phi.control(train[,ncol(train)], method="extremes")
  ls <- uba::loss.control(train[,ncol(train)])
  train <- randOverRegressT(form, train, rel="auto", thr.rel=0.9, C.perc="balance", repl=TRUE)
  ptm <- proc.time()
  m <- randomForest(form,train,mtry=mtry,ntree=ntree,...)
  train.time <- proc.time() - ptm
  p <- predict(m,test)
  trainpred.time <- proc.time() - ptm
  eval <- eval.stats(form,train,test,p,ph,ls)
  res <- list(evaluation=eval,traintime=train.time,trainpredtime=trainpred.time)
  res
}

mc.rf_OVERTPhi <- function(form,train,test,mtry,ntree,...) {
  require(randomForest)
  ph <- uba::phi.control(train[,ncol(train)], method="extremes")
  ls <- uba::loss.control(train[,ncol(train)])
  train <- randOverRegressTPhi(form, train, rel="auto", thr.rel=0.9, C.perc="balance", repl=TRUE)
  ptm <- proc.time()
  m <- randomForest(form,train,mtry=mtry,ntree=ntree,...)
  train.time <- proc.time() - ptm
  p <- predict(m,test)
  trainpred.time <- proc.time() - ptm
  eval <- eval.stats(form,train,test,p,ph,ls)
  res <- list(evaluation=eval,traintime=train.time,trainpredtime=trainpred.time)
  res
}

mc.rf_SMOTEB <- function(form,train,test,mtry,ntree,...) {
  require(randomForest)
  ph <- uba::phi.control(train[,ncol(train)], method="extremes")
  ls <- uba::loss.control(train[,ncol(train)])
  train <- smoteRegressB(form, train, rel="auto", thr.rel=0.9, C.perc="balance", k=5, repl=TRUE, dist="Euclidean", p=2)
  ptm <- proc.time()
  m <- randomForest(form,train,mtry=mtry,ntree=ntree,...)
  train.time <- proc.time() - ptm
  p <- predict(m,test)
  trainpred.time <- proc.time() - ptm
  eval <- eval.stats(form,train,test,p,ph,ls)
  res <- list(evaluation=eval,traintime=train.time,trainpredtime=trainpred.time)
  res
}

mc.rf_SMOTET <- function(form,train,test,mtry,ntree,...) {
  require(randomForest)
  ph <- uba::phi.control(train[,ncol(train)], method="extremes")
  ls <- uba::loss.control(train[,ncol(train)])
  train <- smoteRegressT(form, train, rel="auto", thr.rel=0.9, C.perc="balance", k=5, repl=TRUE, dist="Euclidean", p=2)
  ptm <- proc.time()
  m <- randomForest(form,train,mtry=mtry,ntree=ntree,...)
  train.time <- proc.time() - ptm
  p <- predict(m,test)
  trainpred.time <- proc.time() - ptm
  eval <- eval.stats(form,train,test,p,ph,ls)
  res <- list(evaluation=eval,traintime=train.time,trainpredtime=trainpred.time)
  res
}

mc.rf_SMOTETPhi <- function(form,train,test,mtry,ntree,...) {
  require(randomForest)
  ph <- uba::phi.control(train[,ncol(train)], method="extremes")
  ls <- uba::loss.control(train[,ncol(train)])
  train <- smoteRegressTPhi(form, train, rel="auto", thr.rel=0.9, C.perc="balance", k=5, repl=TRUE, dist="Euclidean", p=2)
  ptm <- proc.time()
  m <- randomForest(form,train,mtry=mtry,ntree=ntree,...)
  train.time <- proc.time() - ptm
  p <- predict(m,test)
  trainpred.time <- proc.time() - ptm
  eval <- eval.stats(form,train,test,p,ph,ls)
  res <- list(evaluation=eval,traintime=train.time,trainpredtime=trainpred.time)
  res
}

mc.arima <- function(form,train,test,...) {
  require(forecast)
  ph <- uba::phi.control(train[,ncol(train)], method="extremes")
  ls <- uba::loss.control(train[,ncol(train)])
  
  trainY <- resp(form,train)
  trues <- resp(form,test)
  
  ptm <- proc.time()
  m <- auto.arima(trainY)
  data <- c(trainY,trues)
  
  train.time <- proc.time() - ptm
  p <- fitted(Arima(data,model=m))[(length(trainY)+1):length(data)]
  trainpred.time <- proc.time() - ptm
  eval <- eval.stats(form,train,test,p,ph,ls)
  res <- list(evaluation=eval,traintime=train.time,trainpredtime=trainpred.time)
  res
}

#############################

embedStats <- function(data, stats = c(mean, var)) {
  if (methods::is(data, "vector")) return(c(data, mean=mean(data), var=var(data)))
  
  dstats <- vapply(stats, function(stat) apply(data, 1, stat), double(NROW(data)))
  colnames(dstats) <- c("mean", "var")
  
  cbind.data.frame(data, dstats)
}

##############################

bootstrap <- function(n) sample(n, n, replace = TRUE)

##############################

baggedtrees <- function(form,
                        data,
                        embedding.dimension,
                        nstats,
                        learner.pars, ...) {
  
  require(rpart)
  
  embed.split.at <- c(embedding.dimension, embedding.dimension / 2, embedding.dimension / 4)
  
  ntrees <- learner.pars[["baggedtrees"]][["ntrees"]]
  n <- NROW(data)
  .i <- 0L
  .n <- floor(ntrees/(2 * length(embed.split.at))); .seqn <- seq_len(.n)
  
  .Models <- lapply(embed.split.at, function(K) {
    predictors <- seq_len(K)
    t1 <- lapply(.seqn, function(i) {
      do.call('rpart', c(list(form, data[bootstrap(n), predictors])))
    })
    predictors <- c(seq_len(K),(NCOL(data) - nstats + 1L):NCOL(data))
    t2 <- lapply(.seqn, function(i) {
      do.call('rpart', c(list(form, data[bootstrap(n), predictors])))
    })
    c(t1, t2)
  })
  names(.Models) <- paste("decisiontree", seq_along(.Models), sep = "_")
  
  unlist(.Models, recursive = FALSE)
}

##############################

mc.BDES <- function(form,train,test,...) {
  require(rpart)
  ph <- uba::phi.control(train[,ncol(train)], method="extremes")
  ls <- uba::loss.control(train[,ncol(train)])
  
  train <- train[,c(10,9,8,7,6,5,4,3,2,1)]
  train <- embedStats(train)
  
  test <- test[,c(10,9,8,7,6,5,4,3,2,1)]
  test <- embedStats(test)
  
  ptm <- proc.time()
  m <- baggedtrees(form,train,embedding.dimension=10,nstats=2,learner.pars = list(baggedtrees = list(ntrees = 500)))
  train.time <- proc.time() - ptm
  nullcheck <- vapply(m, is.null, logical(1))
  if (any(nullcheck)) {
    null.models <- which(nullcheck)
    models[null.models] <- NULL
  }
  predictions <- sapply(m, function(model) predict(model, test))
  
  p <- rowMeans(predictions)
  trainpred.time <- proc.time() - ptm
  names(p) <- rownames(test)
  eval <- eval.stats(form,train,test,p,ph,ls)
  res <- list(evaluation=eval,traintime=train.time,trainpredtime=trainpred.time)
  res
}

exp.time <- performanceEstimation(PredTask(form,ds),
                             c(Workflow("mc.lm"),
                               Workflow("mc.lm_UNDERB"),
                               Workflow("mc.lm_UNDERT"),
                               Workflow("mc.lm_UNDERTPhi"),
                               Workflow("mc.lm_OVERB"),
                               Workflow("mc.lm_OVERT"),
                               Workflow("mc.lm_OVERTPhi"),
                               Workflow("mc.lm_SMOTEB"),
                               Workflow("mc.lm_SMOTET"),
                               Workflow("mc.lm_SMOTETPhi"),
                               Workflow("mc.svm",cost=cost,gamma=gamma),
                               Workflow("mc.svm_UNDERB",cost=cost,gamma=gamma),
                               Workflow("mc.svm_UNDERT",cost=cost,gamma=gamma),
                               Workflow("mc.svm_UNDERTPhi",cost=cost,gamma=gamma),
                               Workflow("mc.svm_OVERB",cost=cost,gamma=gamma),
                               Workflow("mc.svm_OVERT",cost=cost,gamma=gamma),
                               Workflow("mc.svm_OVERTPhi",cost=cost,gamma=gamma),
                               Workflow("mc.svm_SMOTEB",cost=cost,gamma=gamma),
                               Workflow("mc.svm_SMOTET",cost=cost,gamma=gamma),
                               Workflow("mc.svm_SMOTETPhi",cost=cost,gamma=gamma),
                               Workflow("mc.mars",nk=nk,degree=degree,thresh=thresh),
                               Workflow("mc.mars_UNDERB",nk=nk,degree=degree,thresh=thresh),
                               Workflow("mc.mars_UNDERT",nk=nk,degree=degree,thresh=thresh),
                               Workflow("mc.mars_UNDERTPhi",nk=nk,degree=degree,thresh=thresh),
                               Workflow("mc.mars_OVERB",nk=nk,degree=degree,thresh=thresh),
                               Workflow("mc.mars_OVERT",nk=nk,degree=degree,thresh=thresh),
                               Workflow("mc.mars_OVERTPhi",nk=nk,degree=degree,thresh=thresh),
                               Workflow("mc.mars_SMOTEB",nk=nk,degree=degree,thresh=thresh),
                               Workflow("mc.mars_SMOTET",nk=nk,degree=degree,thresh=thresh),
                               Workflow("mc.mars_SMOTETPhi",nk=nk,degree=degree,thresh=thresh),
                               Workflow("mc.rf",mtry=mtry,ntree=ntree),
                               Workflow("mc.rf_UNDERB",mtry=mtry,ntree=ntree),
                               Workflow("mc.rf_UNDERT",mtry=mtry,ntree=ntree),
                               Workflow("mc.rf_UNDERTPhi",mtry=mtry,ntree=ntree),
                               Workflow("mc.rf_OVERB",mtry=mtry,ntree=ntree),
                               Workflow("mc.rf_OVERT",mtry=mtry,ntree=ntree),
                               Workflow("mc.rf_OVERTPhi",mtry=mtry,ntree=ntree),
                               Workflow("mc.rf_SMOTEB",mtry=mtry,ntree=ntree),
                               Workflow("mc.rf_SMOTET",mtry=mtry,ntree=ntree),
                               Workflow("mc.rf_SMOTETPhi",mtry=mtry,ntree=ntree),
                               Workflow("mc.rpart",minsplit=minsplit,cp=cp),
                               Workflow("mc.rpart_UNDERB",minsplit=minsplit,cp=cp),
                               Workflow("mc.rpart_UNDERT",minsplit=minsplit,cp=cp),
                               Workflow("mc.rpart_UNDERTPhi",minsplit=minsplit,cp=cp),
                               Workflow("mc.rpart_OVERB",minsplit=minsplit,cp=cp),
                               Workflow("mc.rpart_OVERT",minsplit=minsplit,cp=cp),
                               Workflow("mc.rpart_OVERTPhi",minsplit=minsplit,cp=cp),
                               Workflow("mc.rpart_SMOTEB",minsplit=minsplit,cp=cp),
                               Workflow("mc.rpart_SMOTET",minsplit=minsplit,cp=cp),
                               Workflow("mc.rpart_SMOTETPhi",minsplit=minsplit,cp=cp),
                               Workflow("mc.arima"),
                               Workflow("mc.BDES")),
                             EstimationTask("totTime",method=MonteCarlo(nReps=50,szTrain=.5,szTest=.25))
)


