library(performanceEstimation)

#EVAL
models <- c("mc.lm","mc.lm_UNDERB","mc.lm_UNDERT","mc.lm_UNDERTPhi","mc.lm_OVERB","mc.lm_OVERT","mc.lm_OVERTPhi","mc.lm_SMOTEB","mc.lm_SMOTET","mc.lm_SMOTETPhi",
            "mc.svm","mc.svm_UNDERB","mc.svm_UNDERT","mc.svm_UNDERTPhi","mc.svm_OVERB","mc.svm_OVERT","mc.svm_OVERTPhi","mc.svm_SMOTEB","mc.svm_SMOTET","mc.svm_SMOTETPhi",
            "mc.mars","mc.mars_UNDERB","mc.mars_UNDERT","mc.mars_UNDERTPhi","mc.mars_OVERB","mc.mars_OVERT","mc.mars_OVERTPhi","mc.mars_SMOTEB","mc.mars_SMOTET","mc.mars_SMOTETPhi",
            "mc.rf","mc.rf_UNDERB","mc.rf_UNDERT","mc.rf_UNDERTPhi","mc.rf_OVERB","mc.rf_OVERT","mc.rf_OVERTPhi","mc.rf_SMOTEB","mc.rf_SMOTET","mc.rf_SMOTETPhi",
            "mc.rpart","mc.rpart_UNDERB","mc.rpart_UNDERT","mc.rpart_UNDERTPhi","mc.rpart_OVERB","mc.rpart_OVERT","mc.rpart_OVERTPhi","mc.rpart_SMOTEB","mc.rpart_SMOTET","mc.rpart_SMOTETPhi",
            "mc.arima","mc.BDES")

globalres.time <- data.frame(model=character(0),tr.time=numeric(0),pr.time=numeric(0),tot.time=numeric(0))

for(wf in 1:length(models)) {
	tr.time <- c()
	pr.time <- c()
	tot.time <- c()
	for(i in 1:10) {
		tr.time_row <- as.numeric(getIterationsInfo(exp.time,workflow=wf,task=1,it=i)$traintime[1])
		tot.time_row <- as.numeric(getIterationsInfo(exp.time,workflow=wf,task=1,it=i)$trainpredtime[1])
    pr.time_row <- tot.time_row - tr.time_row
		tr.time <- c(tr.time,tr.time_row)
		tot.time <- c(tot.time,tot.time_row)
    pr.time <- c(pr.time,pr.time_row)
	}
	
	globalres.row <- data.frame(model=models[wf],tr.time=mean(tr.time),pr.time=mean(pr.time),tot.time=mean(tot.time))
	globalres.time <- rbind(globalres.time,globalres.row)
}

globalres.time